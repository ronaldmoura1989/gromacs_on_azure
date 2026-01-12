#!/bin/bash
# Main resume job script for GROMACS simulations
# This script can be used both for auto-resume on reboot and manual resume after updates

MOUNT_POINT="/data"
DISK_DEVICE="/dev/nvme0n2"
MAX_WAIT=120 # Wait up to 2 minutes for disk
SCRIPTS_DIR="/data/scripts"

# Parse command line arguments
STOP_CURRENT=false
if [ "$1" = "--stop" ]; then
    STOP_CURRENT=true
    echo "Will stop current running jobs before resuming..."
fi

echo "Starting resume job at $(date)"

# If --stop flag is provided, stop all running GROMACS jobs
if [ "$STOP_CURRENT" = true ]; then
    echo "Stopping current GROMACS jobs..."
    pkill -SIGTERM gmx_mpi || true
    pkill -SIGTERM gmx || true

    # Wait for graceful shutdown
    sleep 10

    # Force kill if still running
    pkill -SIGKILL gmx_mpi || true
    pkill -SIGKILL gmx || true

    echo "All GROMACS jobs stopped."
fi

# 1. Wait for disk to auto-mount (from /etc/fstab)
echo "Waiting for data disk to mount..."
count=0
while [ $count -lt $MAX_WAIT ]; do
    if grep -qs "$MOUNT_POINT" /proc/mounts; then
        echo "✓ Disk mounted at $MOUNT_POINT"
        break
    fi
    sleep 1
    ((count++))
done

# 2. If not mounted, try manual mount
if ! grep -qs "$MOUNT_POINT" /proc/mounts; then
    echo "Disk not auto-mounted, attempting manual mount..."

    # Detect partition
    if [ -b "${DISK_DEVICE}p1" ]; then
        PARTITION="${DISK_DEVICE}p1"
    elif [ -b "${DISK_DEVICE}1" ]; then
        PARTITION="${DISK_DEVICE}1"
    else
        echo "ERROR: Cannot find data disk partition!"
        exit 1
    fi

    mkdir -p $MOUNT_POINT
    mount $PARTITION $MOUNT_POINT || {
        echo "ERROR: Failed to mount data disk!"
        exit 1
    }
fi

# 3. Source GROMACS environment
source /usr/local/gromacs/bin/GMXRC

# 4. Ensure resume scripts are executable
chmod +x $SCRIPTS_DIR/*.sh 2>/dev/null || true

# 5. Resume GROMACS simulations
# Find all simulation directories and resume each one
for sim_dir in /data/simulations/*/; do
    if [ -d "$sim_dir" ]; then
        echo "Checking simulation in $sim_dir"
        cd "$sim_dir"

        # Determine which step to resume (MD > NPT > NVT)
        if [ -f "md.tpr" ] && [ -f "md.cpt" ]; then
            echo "Resuming Production MD in $sim_dir"
            $SCRIPTS_DIR/resume_md.sh "$sim_dir"
        elif [ -f "npt.tpr" ] && [ -f "npt.cpt" ]; then
            echo "Resuming NPT Equilibration in $sim_dir"
            $SCRIPTS_DIR/resume_npt.sh "$sim_dir"
        elif [ -f "nvt.tpr" ] && [ -f "nvt.cpt" ]; then
            echo "Resuming NVT Equilibration in $sim_dir"
            $SCRIPTS_DIR/resume_nvt.sh "$sim_dir"
        else
            echo "No resumable checkpoint found in $sim_dir, skipping"
        fi
    fi
done

echo "Resume job completed at $(date)"
