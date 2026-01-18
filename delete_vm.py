import argparse
import subprocess
import sys
import os

# Configuration
DEFAULT_RG = "rg-md-eastus"

def run_command(cmd, verbose=True):
    if verbose:
        print(f"Running: {' '.join(cmd)}")
    try:
        subprocess.check_call(cmd)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        return False

def check_resource_exists(name, resource_type, resource_group):
    # Check if resource exists using az resource show or specific list commands
    # Simpler to just try delete or check with list
    # az <type> show --name <name> -g <rg>
    cmd = ["az", resource_type, "show", "--name", name, "--resource-group", resource_group, "--output", "none"]
    try:
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False

def delete_vm_and_resources(vm_name, resource_group):
    print(f"\n!!! WARNING: You are about to DELETE {vm_name} and its associated resources in {resource_group} !!!")
    print(f"The following will be DELETED:")
    print(f"  - Virtual Machine: {vm_name}")
    print(f"  - OS Disk:         {vm_name}-osdisk")
    print(f"  - Network Interface: {vm_name}-nic")
    print(f"  - Public IP:       {vm_name}-pip")
    print(f"\nThe following will be PRESERVED (Backup):")
    print(f"  + Data Disk:       {vm_name}-datadisk (2TB HDD)")
    print("-" * 60)
    
    confirm = input("Type 'DELETE' to confirm: ")
    if confirm != "DELETE":
        print("Operation cancelled.")
        return

    print("\n=== Starting Deletion Process ===")

    # 1. Delete VM
    print(f"-> Deleting VM: {vm_name}...")
    # --yes to skip confirmation prompt from CLI since we handled it
    run_command(["az", "vm", "delete", "--resource-group", resource_group, "--name", vm_name, "--yes"])

    # 2. Delete NIC
    # NICs usually can't be deleted if attached. Since VM is deleted, it should be detached.
    nic_name = f"{vm_name}-nic"
    print(f"-> Deleting NIC: {nic_name}...")
    run_command(["az", "network", "nic", "delete", "--resource-group", resource_group, "--name", nic_name])

    # 3. Delete Public IP
    pip_name = f"{vm_name}-pip"
    print(f"-> Deleting Public IP: {pip_name}...")
    run_command(["az", "network", "public-ip", "delete", "--resource-group", resource_group, "--name", pip_name])

    # 4. Delete OS Disk
    os_disk_name = f"{vm_name}-osdisk"
    print(f"-> Deleting OS Disk: {os_disk_name}...")
    run_command(["az", "disk", "delete", "--resource-group", resource_group, "--name", os_disk_name, "--yes"])

    print("\n=== Deletion Complete ===")
    print(f"✅ Resources for {vm_name} deleted.")
    print(f"💾 Preserved Data Disk: {vm_name}-datadisk")

def main():
    parser = argparse.ArgumentParser(description="Delete an Azure VM and its resources, preserving the Data Disk.")
    parser.add_argument("vm_name", help="The name of the VM to delete (e.g., vm-md-09)")
    parser.add_argument("--rg", default=DEFAULT_RG, help=f"Resource Group (default: {DEFAULT_RG})")
    
    args = parser.parse_args()

    # Simple check if az is logged in
    try:
        subprocess.check_call(["az", "account", "show"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print("Error: Please login to Azure CLI using 'az login' first.")
        sys.exit(1)

    delete_vm_and_resources(args.vm_name, args.rg)

if __name__ == "__main__":
    main()
