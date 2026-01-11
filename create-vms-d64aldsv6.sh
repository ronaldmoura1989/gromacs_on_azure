#### login on Azure ####

# Azure CLI (versão 2.50+)
az --version

# Login no Azure
az login

# Selecionar subscription
az account set --subscription "${AZURE_SUBSCRIPTION}"

# Instalar extensões necessárias
az extension add --name account
az extension add --name storage-preview

#### verify quotas ####

# 1. Verificar quota em US East para D-series
az vm list-usage \
  --location eastus \
  --output table | grep "standard Daldv6 Family"

#### create resource group ####

az group create \
  --name rg-md-eastus \
  --location eastus \
  --tags project=molecular-dynamics environment=production

#### create vnet ####

# Criar VNet
az network vnet create \
  --resource-group rg-md-eastus \
  --name vnet-md-eastus \
  --address-prefix 10.1.0.0/16 \
  --subnet-name subnet-compute \
  --subnet-prefix 10.1.1.0/24 \
  --location eastus

# Criar subnet adicional para futuro uso
az network vnet subnet create \
  --resource-group rg-md-eastus \
  --vnet-name vnet-md-eastus \
  --name subnet-storage \
  --address-prefix 10.1.2.0/24

#### network security group ####

# Criar NSG
az network nsg create \
  --resource-group rg-md-eastus \
  --name nsg-md-compute \
  --location eastus

# Regra SSH (ajuste source IP para seu range)
az network nsg rule create \
  --resource-group rg-md-eastus \
  --nsg-name nsg-md-compute \
  --name allow-ssh-brazil \
  --priority 1000 \
  --source-address-prefixes 200.0.0.0/8 45.4.61.238 \
  --destination-port-ranges 22 \
  --protocol Tcp \
  --access Allow
 
# Regra outbound para atualizações
az network nsg rule create \
  --resource-group rg-md-eastus \
  --nsg-name nsg-md-compute \
  --name allow-outbound-http \
  --priority 1001 \
  --direction Outbound \
  --source-address-prefixes '*' \
  --destination-port-ranges 80 443 \
  --protocol Tcp \
  --access Allow

# Associar NSG à subnet
az network vnet subnet update \
  --resource-group rg-md-eastus \
  --vnet-name vnet-md-eastus \
  --name subnet-compute \
  --network-security-group nsg-md-compute

#### deploy machines #####

# Load environment variables
if [ -f .env ]; then
    source .env
else
    echo "Error: .env file not found."
    exit 1
fi

TARGET_RG="rg-md-eastus"
TARGET_LOCATION="eastus"
VM_SIZE="Standard_D64alds_v6"
VNET_NAME="vnet-md-eastus"
SUBNET_NAME="subnet-compute"
# CAUTION: Change this password before running in production!
ADMIN_USER="${SSH_USER}"
ADMIN_PASS="${SSH_PASSWORD}" 

echo "=== Creating D64alds_v6 VMs ==="

# Loop using printf for robust zero-padding on macOS Bash 3.2
for i in {1..9}; do
  vm_num=$(printf "%02d" $i)
  vm_name="vm-md-${vm_num}"
  nic_name="${vm_name}-nic"
  data_disk_name="${vm_name}-datadisk"
  public_ip_name="${vm_name}-pip"
  
  echo "Processing VM: $vm_name"

  # Create Data disk (2 TB HDD -> Standard_LRS)
  echo "  Creating Data disk: $data_disk_name (2TB)"
  az disk create \
    --resource-group $TARGET_RG \
    --name $data_disk_name \
    --size-gb 2048 \
    --sku Standard_LRS \
    --location $TARGET_LOCATION \
    --output none

  # Create Public IP
  echo "  Creating Public IP: $public_ip_name"
  az network public-ip create \
    --resource-group $TARGET_RG \
    --name $public_ip_name \
    --sku Standard \
    --location $TARGET_LOCATION \
    --allocation-method Static \
    --output none

  # Create NIC with Public IP
  echo "  Creating NIC: $nic_name"
  az network nic create \
    --resource-group $TARGET_RG \
    --name $nic_name \
    --vnet-name $VNET_NAME \
    --subnet $SUBNET_NAME \
    --location $TARGET_LOCATION \
    --public-ip-address $public_ip_name \
    --output none
  
  # Create VM with Image and Password Auth (Creates 30GB OS Disk Automatically)
  # Note: Standard_D64alds_v6 is an updated AMD generation.
  echo "  Creating VM $vm_name"
  az vm create \
    --resource-group $TARGET_RG \
    --name $vm_name \
    --location $TARGET_LOCATION \
    --size $VM_SIZE \
    --image Ubuntu2204 \
    --os-disk-size-gb 30 \
    --os-disk-name "${vm_name}-osdisk" \
    --storage-sku Premium_LRS \
    --admin-username $ADMIN_USER \
    --admin-password $ADMIN_PASS \
    --authentication-type password \
    --nics $nic_name \
    --tags environment=production workload=molecular-dynamics \
    --output none
  
  # Attach data disk
  echo "  Attaching data disk: $data_disk_name"
  az vm disk attach \
    --resource-group $TARGET_RG \
    --vm-name $vm_name \
    --name $data_disk_name \
    --lun 0 \
    --output none
    
  echo "  ✓ VM $vm_name created successfully"
  echo ""
done

echo "=== All VMs created ==="
echo "Fetching Public IPs..."

# List name and Public IP
az network public-ip list \
  --resource-group $TARGET_RG \
  --query "[].{Name:name, IP:ipAddress}" \
  --output table

echo ""
echo "Access instructions:"
echo "ssh $ADMIN_USER@<PUBLIC-IP>"
echo "Example: ssh $ADMIN_USER@20.x.x.x"

#chmod +x create-vms-d64aldsv6.sh
#./create-vms-d64aldsv6.sh