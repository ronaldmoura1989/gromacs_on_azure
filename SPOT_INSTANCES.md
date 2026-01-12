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

### B. Outside the VM (Auto-Restart Manager)
You need a "Manager" VM to continuously monitor and restart evicted Spot instances. While you could run this on your laptop, deploying a dedicated small Azure VM ensures 24/7 monitoring without keeping your local machine running.

#### Option 1: Deploy a Dedicated Manager VM (Recommended)
Create a small, always-on VM in the same resource group to manage your Spot instances. This VM will be extremely cheap (~$10-15/month) compared to your compute cluster.

**Script**: `create-manager-vm.sh`
```bash
#!/bin/bash
# Creates a small manager VM to monitor and restart Spot instances

set -e

# Load credentials
if [ -f .env ]; then
    source .env
else
    echo "Error: .env file not found."
    exit 1
fi

TARGET_RG="rg-md-eastus"
TARGET_LOCATION="eastus"
MANAGER_NAME="vm-spot-manager"
VNET_NAME="vnet-md-eastus"
SUBNET_NAME="subnet-compute"
ADMIN_USER="${SSH_USER}"
ADMIN_PASS="${SSH_PASSWORD}"

echo "=== Creating Spot Manager VM ==="

# Create Public IP
echo "Creating Public IP..."
az network public-ip create \
  --resource-group $TARGET_RG \
  --name ${MANAGER_NAME}-pip \
  --sku Standard \
  --location $TARGET_LOCATION \
  --allocation-method Static \
  --output none

# Create NIC
echo "Creating NIC..."
az network nic create \
  --resource-group $TARGET_RG \
  --name ${MANAGER_NAME}-nic \
  --vnet-name $VNET_NAME \
  --subnet $SUBNET_NAME \
  --location $TARGET_LOCATION \
  --public-ip-address ${MANAGER_NAME}-pip \
  --output none

# Create small VM (B2s: 2 vCPUs, 4GB RAM - ~$0.05/hr = $36/month)
echo "Creating Manager VM..."
az vm create \
  --resource-group $TARGET_RG \
  --name $MANAGER_NAME \
  --location $TARGET_LOCATION \
  --size Standard_B2s \
  --image Ubuntu2204 \
  --os-disk-size-gb 30 \
  --os-disk-name "${MANAGER_NAME}-osdisk" \
  --storage-sku StandardSSD_LRS \
  --admin-username $ADMIN_USER \
  --admin-password $ADMIN_PASS \
  --authentication-type password \
  --nics ${MANAGER_NAME}-nic \
  --tags environment=production workload=spot-manager \
  --output none

# Assign VM Contributor role to the Manager VM (so it can start/stop VMs)
echo "Assigning VM Contributor role..."
MANAGER_PRINCIPAL_ID=$(az vm identity assign \
  --resource-group $TARGET_RG \
  --name $MANAGER_NAME \
  --query principalId -o tsv)

az role assignment create \
  --assignee $MANAGER_PRINCIPAL_ID \
  --role "Virtual Machine Contributor" \
  --scope "/subscriptions/${AZURE_SUBSCRIPTION}/resourceGroups/${TARGET_RG}"

echo "✓ Manager VM created successfully"

# Get Public IP
MANAGER_IP=$(az network public-ip show \
  --resource-group $TARGET_RG \
  --name ${MANAGER_NAME}-pip \
  --query ipAddress -o tsv)

echo ""
echo "Manager VM IP: $MANAGER_IP"
echo "SSH: ssh $ADMIN_USER@$MANAGER_IP"
echo ""
echo "Next Steps:"
echo "1. SSH into the manager VM"
echo "2. Install Azure CLI: curl -sL https://aka.ms/InstallAzureCLIDeb | sudo bash"
echo "3. Login with managed identity: az login --identity"
echo "4. Deploy spot_manager.sh script (see below)"
```

#### Deploying the Monitor Script on Manager VM
Once the manager VM is created, SSH into it and set up the monitoring script:

**Script**: `spot_manager.sh` (Deploy to `/home/[user]/spot_manager.sh` on manager VM)
```bash
#!/bin/bash
# Checks for STOPPED VMs every 15 minutes and tries to start them.

RG="rg-md-eastus"

# Log file
LOG_FILE="/home/$(whoami)/spot_manager.log"

# Login using managed identity (no credentials needed!)
az login --identity > /dev/null 2>&1

while true; do
    echo "[$(date)] Checking for stopped Spot instances..." | tee -a $LOG_FILE

    # Get IDs of all STOPPED VMs (excluding the manager itself)
    IDS=$(az vm list -g $RG -d \
      --query "[?powerState=='VM stopped' && name!='vm-spot-manager'].id" -o tsv)

    if [ ! -z "$IDS" ]; then
        echo "Found stopped VMs. Attempting to restart..." | tee -a $LOG_FILE
        # Try to start them. If capacity is unavailable, this command will fail/error, which is fine.
        az vm start --ids $IDS 2>&1 | tee -a $LOG_FILE || \
          echo "  > Start failed (likely no capacity). Retrying in 15 mins." | tee -a $LOG_FILE
    else
        echo "  > All VMs are running." | tee -a $LOG_FILE
    fi

    # Wait 15 minutes
    sleep 900
done
```

**Setup on Manager VM**:
```bash
# 1. SSH into manager VM
ssh ${SSH_USER}@<MANAGER_IP>

# 2. Install Azure CLI
curl -sL https://aka.ms/InstallAzureCLIDeb | sudo bash

# 3. Login with managed identity
az login --identity

# 4. Create the script
nano ~/spot_manager.sh
# (Paste the script above)
chmod +x ~/spot_manager.sh

# 5. Run in background with nohup
nohup ~/spot_manager.sh &

# 6. (Optional) Add to crontab for auto-restart on reboot
crontab -e
# Add: @reboot nohup /home/[username]/spot_manager.sh &
```

#### Option 2: Run on Your Laptop
If you prefer to run the manager on your local machine (must stay on 24/7):

```bash
# Ensure you're logged in to Azure CLI
az login

# Run the spot_manager.sh script (without the 'az login --identity' line)
./spot_manager.sh
```

### C. Auto-Resume Logic
When the VM starts up again, it needs to automatically mount the disk and resume the run.

#### Setup Resume Scripts
All resume scripts are located in the `scripts/` folder of this repository. Deploy them to `/data/scripts/` on each VM:

```bash
# From your local machine, copy scripts to all VMs
scp scripts/*.sh user@vm-ip:/data/scripts/
ssh user@vm-ip "chmod +x /data/scripts/*.sh"
```

Add this line to your `crontab` on each VM (`crontab -e`):
```bash
@reboot /data/scripts/resume_job.sh
```

#### Manual Resume After Updates
If you need to update scripts and resume simulations manually, use the automated update script:

```bash
# Option 1: Use update_scripts.py (RECOMMENDED - updates all VMs from gromacs_runs.tsv)
python3 update_scripts.py --restart

# Option 2: Manual update for a single VM
scp scripts/*.sh user@vm-ip:/data/scripts/
scp gromacs_azure.sh user@vm-ip:/data/simulations/{basename}/
scp gromacs_parameters/*.mdp user@vm-ip:/data/simulations/{basename}/
ssh user@vm-ip "/data/scripts/resume_job.sh --stop"
```

The `update_scripts.py` script:
- Reads VMs from `gromacs_runs.tsv` (same as `run_simulations.py`)
- Uploads `gromacs_azure.sh` and `.mdp` files to simulation directories
- Uploads resume scripts to `/data/scripts/`
- With `--restart` flag: stops jobs and resumes from checkpoints
- Ensures all VMs get updates simultaneously

See README.md Step 6 for full documentation.

#### Resume Scripts Overview

All scripts are available in the `scripts/` folder. See the actual files for full implementation details.

**`scripts/resume_job.sh`** - Main orchestration script:
- Supports `--stop` flag to stop current jobs before resuming
- Waits for disk auto-mount (up to 2 minutes) or manually mounts `/dev/nvme0n2`
- Scans `/data/simulations/*/` for checkpoints
- Determines simulation step (MD > NPT > NVT) and calls appropriate resume script
- Used both for auto-resume on reboot and manual resume after updates

**`scripts/resume_nvt.sh`** - Resume from NVT step:
- Resumes NVT equilibration from checkpoint
- Continues full workflow: NVT → NPT → MD → Analysis
- Uses optimized GMX flags with dynamic load balancing
- Runs in background with nohup

**`scripts/resume_npt.sh`** - Resume from NPT step:
- Resumes NPT equilibration from checkpoint
- Continues workflow: NPT → MD → Analysis
- Uses optimized GMX flags with dynamic load balancing

**`scripts/resume_md.sh`** - Resume from MD step:
- Resumes Production MD from checkpoint
- Continues to Analysis after MD completes
- Uses PME auto-tuning for optimal performance
