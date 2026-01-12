# Using Azure Spot Instances for Molecular Dynamics

This guide explains how to use **Azure Spot Instances** to reduce computing costs by up to **90%**.

## 1. What are Spot Instances?
Spare capacity that Azure offers at a deep discount. 
*   **Pros**: Extremely cheap (~$0.30/hr instead of $3.05/hr).
*   **Cons**: Azure can shut them down (evict) with 30 seconds of notice if they need the capacity back.

**Ideal Workflow**: Run simulation -> Checkpoint every 15 mins -> Evicted -> Wait for capacity -> Resume from last checkpoint.

---

## 2. Converting to Spot Instances

You cannot simply "toggle" a button to convert an existing Standard VM to Spot. You must:
1.  **Delete** the existing VM (but keep the 2TB data disk and Public IP/NIC).
2.  **Recreate** the VM as a Spot Instance, attaching the *existing* components.

### Script: `migrate_to_spot.sh`
Save the following script as `migrate_to_spot.sh` in your `azure_md_infrastructure` folder. Ideally, ensure you have backed up any data on the **OS Disk** (the 30GB one) because it will be **reset** (GROMACS needs re-installation). Your 2TB Data Disk is safe.

```bash
#!/bin/bash
# migrate_to_spot.sh
# Deletes existing Standard VMs and recreates them as SPOT instances.
# PRESERVES: 2TB Data Disks, Public IPs, NICs.
# DESTROYS: OS Disks (You must run install-gromacs again).

set -e

# Load credentials
if [ -f .env ]; then
    source .env
else
    echo "Error: .env file not found."
    exit 1
fi

TARGET_RG="rg-md-eastus"
VM_SIZE="Standard_D64alds_v6"
ADMIN_USER="${SSH_USER}"
ADMIN_PASS="${SSH_PASSWORD}"

echo "=== MIGRATING TO SPOT INSTANCES ==="
echo "WARNING: This will DELETE the VM object and OS Disk."
echo "Your 2TB Data Disks will be PRESERVED and re-attached."
echo "Press Ctrl+C to cancel within 5 seconds..."
sleep 5

for i in {1..9}; do
  vm_num=$(printf "%02d" $i)
  vm_name="vm-md-${vm_num}"
  nic_name="${vm_name}-nic"
  data_disk_name="${vm_name}-datadisk"
  os_disk_name="${vm_name}-osdisk"

  echo "Processing $vm_name..."

  # 1. DELETE EXISTING VM (If exists)
  echo "  [1/4] Deleting existing Standard VM..."
  # --yes to confirm non-interactive
  # Only deletes the VM computes, keeps NICs/Disks by default unless specified
  az vm delete --resource-group $TARGET_RG --name $vm_name --yes --output none || true
  
  # Delete the old OS disk to avoid conflict/ensure clean slate
  echo "  [2/4] Cleaning old OS disk..."
  az disk delete --resource-group $TARGET_RG --name $os_disk_name --yes --output none || true

  # 2. CREATE SPOT VM
  echo "  [3/4] Creating SPOT VM (Eviction Policy: Deallocate)..."
  az vm create \
    --resource-group $TARGET_RG \
    --name $vm_name \
    --location eastus \
    --size $VM_SIZE \
    --image Ubuntu2204 \
    --priority Spot \
    --eviction-policy Deallocate \
    --max-price -1 \
    --os-disk-size-gb 30 \
    --os-disk-name $os_disk_name \
    --storage-sku Premium_LRS \
    --admin-username $ADMIN_USER \
    --admin-password $ADMIN_PASS \
    --authentication-type password \
    --nics $nic_name \
    --tags environment=production workload=molecular-dynamics type=spot \
    --output none

  # 3. ATTACH EXISTING DATA DISK
  echo "  [4/4] Re-attaching existing 2TB Data Disk..."
  az vm disk attach \
    --resource-group $TARGET_RG \
    --vm-name $vm_name \
    --name $data_disk_name \
    --output none

  echo "  ✓ $vm_name is now a Spot Instance."
done

echo "=== Migration Complete ==="
echo "Next Steps:"
echo "1. Run ./install-gromacs-on-all-vms.sh"
echo "2. Run ./mount-all-disks.sh"
```

---

## 3. Handling Eviction (The Watchdog)

Since Azure Spot instances can be stopped at any time, you need an automated logical loop to handle this.

### A. Inside the VM (Graceful Shutdown)
GROMACS writes checkpoints (`.cpt`) automatically. To ensure the *very last* bits are saved when Azure warns us of shutdown (30 seconds warning), run this script in the background of every simulation.

**Script**: `watch_preemption.py` (Deploy this to `/data/scripts/` on the VM)
```python
import requests
import time
import os
import smtplib
from email.mime.text import MIMEText
from dotenv import load_dotenv

# Load environment variables
load_dotenv("/data/.env") # Assumes .env is in /data

# Email Configuration
EMAIL_HOST = os.getenv("EMAIL_HOST", "smtp.gmail.com")
EMAIL_PORT = int(os.getenv("EMAIL_PORT", 587))
EMAIL_HOST_USER = os.getenv("EMAIL_HOST_USER")
EMAIL_HOST_PASSWORD = os.getenv("GOOGLE_APP_PASSWORD") 
EMAIL_FROM = os.getenv("DEFAULT_FROM_EMAIL")
EMAIL_TO = os.getenv("EMAIL_TO_SEND")

def get_vm_name():
    try:
        # Get VM Name from Azure Metadata Service
        url = "http://169.254.169.254/metadata/instance/compute/name?api-version=2021-02-01&format=text"
        headers = {"Metadata": "true"}
        return requests.get(url, headers=headers, timeout=2).text
    except:
        return "Unknown-VM"

def send_email_alert(vm_name):
    subject = f"Azure Spot Eviction Alert: {vm_name}"
    body = f"URGENT: The Azure Spot Instance '{vm_name}' has received a preemption notice and will be shut down in 30 seconds.\n\nGROMACS checkpointing has been triggered."
    
    msg = MIMEText(body)
    msg['Subject'] = subject
    msg['From'] = EMAIL_FROM
    msg['To'] = EMAIL_TO

    try:
        server = smtplib.SMTP(EMAIL_HOST, EMAIL_PORT)
        server.starttls()
        server.login(EMAIL_HOST_USER, EMAIL_HOST_PASSWORD)
        server.sendmail(EMAIL_FROM, EMAIL_TO, msg.as_string())
        server.quit()
        print(f"Email alert sent for {vm_name}")
    except Exception as e:
        print(f"Failed to send email: {e}")

def check_preemption():
    # Azure Metadata Service
    url = "http://169.254.169.254/metadata/scheduledevents?api-version=2020-07-01"
    headers = {"Metadata": "true"}
    try:
        resp = requests.get(url, headers=headers, timeout=2).json()
        for event in resp.get('Events', []):
            if event['EventType'] == 'Preempt':
                return True
    except:
        pass
    return False

print("Watcher started. Monitoring for Azure Preemption...")
while True:
    if check_preemption():
        vm_name = get_vm_name()
        print(f"PREEMPTION DETECTED for {vm_name}! Sending alert and killing GROMACS...")
        
        # 1. Send Email (Do this first as it takes <1s)
        send_email_alert(vm_name)
        
        # 2. Cleanly kill GROMACS.
        os.system("pkill -SIGINT gmx_mpi")
        os.system("pkill -SIGINT gmx")
        break
    time.sleep(5)
```

### B. Outside the VM (Auto-Restart)
You need a "Manager" (your laptop or a tiny separate VM) to check if the Spot VMs have been stopped and try to start them again.

**Script**: `spot_manager.sh` (Run continuously on your laptop)
```bash
#!/bin/bash
# Checks for STOPPED VMs every 15 minutes and tries to start them.

RG="rg-md-eastus"

while true; do
    echo "[$(date)] Checking for stopped Spot instances..."
    
    # Get IDs of all STOPPED VMs
    IDS=$(az vm list -g $RG -d --query "[?powerState=='VM stopped'].id" -o tsv)
    
    if [ ! -z "$IDS" ]; then
        echo "Found stopped VMs. Attempting to restart..."
        # Try to start them. If capacity is unavailable, this command will fail/error, which is fine.
        az vm start --ids $IDS || echo "  > Start failed (likely no capacity). Retrying in 15 mins."
    else
        echo "  > All VMs are running."
    fi
    
    # Wait 15 minutes
    sleep 900
done
```

### C. Auto-Resume Logic
When the VM starts up again, it needs to automatically mount the disk and resume the run.
Add this line to your `crontab` on the VM (`crontab -e`):
```bash
@reboot /data/resume_job.sh
```

**`resume_job.sh`**:
```bash
#!/bin/bash
sleep 30 # Wait for network/disks

# 1. Mount Disk
mount /dev/disk/by-uuid/<UUID> /data

# 2. Source GROMACS environment
source /usr/local/gromacs/bin/GMXRC

# 3. Resume GROMACS simulations
# Find all simulation directories and resume each one
for sim_dir in /data/simulations/*/; do
    if [ -d "$sim_dir" ]; then
        echo "Checking simulation in $sim_dir"
        cd "$sim_dir"

        # Check if there's a checkpoint file and TPR file
        if [ -f "state.cpt" ] && [ -f "topol.tpr" ]; then
            echo "Resuming simulation in $sim_dir"
            # -cpi checks for checkpoint. If found, it appends.
            # Run with nohup to keep it running in background
            nohup /usr/local/gromacs/bin/gmx mdrun -s topol.tpr -cpi state.cpt -maxh 24 >> run.log 2>&1 &
            echo "Simulation resumed with PID $!"
        else
            echo "No checkpoint or TPR file found in $sim_dir, skipping"
        fi
    fi
done
```
