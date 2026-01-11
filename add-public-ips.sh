#!/bin/bash
set -e

TARGET_RG="rg-md-eastus"
TARGET_LOCATION="eastus"

echo "=== Adding Public IPs to existing VMs in $TARGET_RG ==="

# Loop through the VM numbers 1 to 9
for i in {1..9}; do
  # Force 2-digit zero padding (compatibility with macOS Bash 3.2)
  vm_num=$(printf "%02d" $i)
  
  vm_name="vm-md-${vm_num}"
  # Naming convention matches the other scripts
  nic_name="${vm_name}-nic"
  public_ip_name="${vm_name}-pip"

  echo "Processing $vm_name..."

  # 1. Check if NIC exists first to avoid errors if VM doesn't exist
  nic_exists=$(az network nic show --resource-group $TARGET_RG --name $nic_name --query id -o tsv 2>/dev/null || true)
  if [[ -z "$nic_exists" ]]; then
    echo "  ⚠️  NIC $nic_name not found. Skipping VM..."
    continue
  fi

  # 2. Create Public IP
  echo "  Creating Public IP resource: $public_ip_name"
  az network public-ip create \
    --resource-group $TARGET_RG \
    --name $public_ip_name \
    --sku Standard \
    --allocation-method Static \
    --location $TARGET_LOCATION \
    --output none

  # 3. Get the IP Configuration name (usually ipconfig1, but safest to query)
  ip_config_name=$(az network nic show \
    --resource-group $TARGET_RG \
    --name $nic_name \
    --query "ipConfigurations[0].name" \
    -o tsv)

  # 4. Update the NIC IP Config to attach the Public IP
  echo "  Attaching Public IP to NIC ($nic_name / $ip_config_name)..."
  az network nic ip-config update \
    --resource-group $TARGET_RG \
    --nic-name $nic_name \
    --name $ip_config_name \
    --public-ip-address $public_ip_name \
    --output none

  echo "  ✓ Public IP attached successfully."
  echo ""
done

echo "=== Operation Complete ==="
echo "Current Public IPs:"
az network public-ip list \
  --resource-group $TARGET_RG \
  --query "[].{VM_Suffix:name, IP:ipAddress}" \
  --output table
