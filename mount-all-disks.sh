#!/bin/bash
set -e

# Load environment variables
if [ -f .env ]; then
    source .env
else
    echo "Error: .env file not found."
    exit 1
fi

TARGET_RG="rg-md-eastus"
SCRIPT_FILE="mount_data_disk.sh"

# Verify script exists
if [[ ! -f "$SCRIPT_FILE" ]]; then
    echo "Error: $SCRIPT_FILE not found in current directory."
    exit 1
fi

echo "=== Mounting Data Disks on all VMs in $TARGET_RG ==="
echo "Using script: $SCRIPT_FILE"
echo "Target User: $SSH_USER"

# Get list of VM names
vm_list=$(az vm list --resource-group $TARGET_RG --query "[].name" -o tsv)

# Arrays to track jobs (Parallel arrays for compatibility with Bash 3.2)
vm_names=()
job_pids=()

for vm_name in $vm_list; do
    echo "Triggering disk mount on $vm_name in background..."
    
    # Run the command in background
    az vm run-command invoke \
        --resource-group $TARGET_RG \
        --name $vm_name \
        --command-id RunShellScript \
        --scripts "@$SCRIPT_FILE" \
        --parameters "$SSH_USER" \
        --output none &
        
    # Store PID and Name in separate arrays by index
    job_pids+=($!)
    vm_names+=($vm_name)
done

echo "Waiting for all operations to complete..."

# Wait for all background jobs
for i in "${!job_pids[@]}"; do
    pid=${job_pids[$i]}
    name=${vm_names[$i]}
    
    if wait $pid; then
        echo "✓ Disk mounted successfully on $name"
    else
        echo "✗ Disk mount FAILED on $name"
    fi
done

echo "=== All operations completed ==="
