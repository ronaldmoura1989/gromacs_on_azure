#!/bin/bash
set -e

# Configuration
MOUNT_POINT="/data"
USERNAME="${1:-ronaldmoura}"
DISK_DEVICE="/dev/nvme0n2" # Updated to identified NVMe device path

echo "Starting data disk configuration on $(hostname)..."

# 1. Verify Disk Existence
if [ ! -b "$DISK_DEVICE" ]; then
    echo "Error: Device $DISK_DEVICE not found."
    exit 1
fi

# 2. Check if already mounted
if grep -qs "$MOUNT_POINT" /proc/mounts; then
    echo "Disk already mounted at $MOUNT_POINT."
else 
    # 3. Partitioning (if partition 1 does not exist)
    # 3. Partitioning (if partition 1 does not exist)
    if [ ! -b "${DISK_DEVICE}p1" ] && [ ! -b "${DISK_DEVICE}1" ]; then
        echo "Partitioning $DISK_DEVICE..."
        # Use parted to create GPT label and primary partition
        parted $DISK_DEVICE --script mklabel gpt mkpart primary ext4 0% 100%
        partprobe $DISK_DEVICE
        sleep 5 # Wait for kernel to register partition
    fi
    
    # Detect existing partition or partition just created
    if [ -b "${DISK_DEVICE}p1" ]; then
        PARTITION="${DISK_DEVICE}p1"
    elif [ -b "${DISK_DEVICE}1" ]; then
        PARTITION="${DISK_DEVICE}1"
    else
        # Partition likely doesn't exist yet, retry probing
        partprobe $DISK_DEVICE
        sleep 2
        if [ -b "${DISK_DEVICE}p1" ]; then
             PARTITION="${DISK_DEVICE}p1"
        elif [ -b "${DISK_DEVICE}1" ]; then
             PARTITION="${DISK_DEVICE}1"
        else
             echo "Error: Partition detection failed on $DISK_DEVICE"
             exit 1
        fi
    fi
    
    # 4. Formatting - Check for valid filesystem
    # blkid returns 0 if attributes found, 2 if not found.
    # We specifically check if a TYPE (filesystem) exists.
    FS_TYPE=$(blkid -o value -s TYPE $PARTITION || true)
    
    if [ -z "$FS_TYPE" ]; then
        echo "Formatting $PARTITION with ext4..."
        mkfs.ext4 -F $PARTITION
    else
        echo "Partition $PARTITION already has filesystem: $FS_TYPE"
    fi
    
    # 5. Mounting
    echo "Mounting at $MOUNT_POINT..."
    mkdir -p $MOUNT_POINT
    mount $PARTITION $MOUNT_POINT
fi

# 6. Persistence (Add to /etc/fstab)
# Get UUID of the partition
# PARTITION variable is preserved from previous steps
UUID=$(blkid -s UUID -o value $PARTITION)

if [[ -z "$UUID" ]]; then
    echo "Error: Could not obtain UUID for $PARTITION"
    exit 1
fi

if ! grep -qs "$UUID" /etc/fstab; then
    echo "Adding entry to /etc/fstab..."
    echo "UUID=$UUID $MOUNT_POINT ext4 defaults,nofail 0 2" >> /etc/fstab
fi

# 7. Permissions
echo "Setting ownership to $USERNAME..."
chown -R $USERNAME:$USERNAME $MOUNT_POINT
chmod 755 $MOUNT_POINT

echo "✓ Data disk successfully mounted at $MOUNT_POINT"
