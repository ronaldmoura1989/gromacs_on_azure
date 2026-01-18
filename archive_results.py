import csv
import os
import subprocess
import sys
import datetime
import json
import argparse

# Configuration
ENV_FILE = ".env"
TSV_FILE = "gromacs_runs.tsv"
CONTAINER_NAME = "simulations"

def load_env(env_path):
    config = {}
    if os.path.exists(env_path):
        with open(env_path) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    # Handle cases with multiple = in value (like connection strings)
                    parts = line.strip().split('=', 1)
                    if len(parts) == 2:
                        config[parts[0]] = parts[1]
    return config

def run_command(cmd, capture_output=False):
    if capture_output:
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return result.stdout.strip()
        except subprocess.CalledProcessError as e:
            print(f"DEBUG: Command failed: {cmd}")
            if e.stderr:
                print(f"DEBUG: Stderr: {e.stderr.strip()}")
            return None
    else:
        try:
            subprocess.check_call(cmd)
            return True
        except subprocess.CalledProcessError:
            return False

def get_storage_info(config):
    account_name = config.get("STORAGE_ACCOUNT_NAME")
    account_key = config.get("STORAGE_ACCOUNT_KEY")
    resource_group = "rg-md-eastus" # Default from infrastructure

    if not account_name:
        print("Storage account name not found in .env, attempting to find one in rg-md-eastus...")
        # Dictionary query to find the account starting with stmdsimulations
        cmd = ["az", "storage", "account", "list", "--resource-group", resource_group, "--query", "[?starts_with(name, 'stmdsimulations')].name | [0]", "-o", "tsv"]
        account_name = run_command(cmd, capture_output=True)
        
        if not account_name:
            print("❌ Could not find a storage account. Please run ./setup_storage.sh first.")
            sys.exit(1)
        print(f"Found Storage Account: {account_name}")

    if not account_key:
        print(f"Fetching key for {account_name}...")
        cmd = ["az", "storage", "account", "keys", "list", "--account-name", account_name, "--resource-group", resource_group, "--query", "[0].value", "-o", "tsv"]
        account_key = run_command(cmd, capture_output=True)
        
        if not account_key:
            print("❌ Failed to fetch storage account key. Please check your Azure login.")
            sys.exit(1)

    return account_name, account_key

def generate_sas_token(account_name, account_key, container_name):
    # Expiry: tomorrow
    expiry = (datetime.datetime.now() + datetime.timedelta(days=1)).strftime('%Y-%m-%dT%H:%M:%SZ')
    
    cmd = [
        "az", "storage", "container", "generate-sas",
        "--account-name", account_name,
        "--account-key", account_key,
        "--name", container_name,
        "--permissions", "racwdl",
        "--expiry", expiry,
        "--output", "tsv"
    ]
    
    sas = run_command(cmd, capture_output=True)
    if not sas:
        print("❌ Failed to generate SAS token.")
        sys.exit(1)
    return sas

def archive_vm(vm_id, basename, ip_address, ssh_user, ssh_pass, has_sshpass, container_url, sas_token, file_pattern):
    print(f"\n--------------------------------------------------")
    print(f"Archiving {vm_id} ({ip_address}) -> {basename}")
    print(f"--------------------------------------------------")

    # SSH Base
    ssh_base = ["ssh", "-o", "StrictHostKeyChecking=no", f"{ssh_user}@{ip_address}"]
    if has_sshpass and ssh_pass:
        ssh_wrapper = ["sshpass", "-p", ssh_pass]
        ssh_base = ssh_wrapper + ssh_base

    # Remote Directory
    # User confirmed directory does NOT have .pdb at the end
    folder_name = basename.replace(".pdb", "").replace(".PDB", "")
    remote_base = "/data/simulations"
    
    # Destination URL
    # https://<account>.blob.../<container>/<folder>?<sas>
    # Note: AzCopy syntax is `azcopy copy "src" "https://.../container/folder?sas"`
    # So we append `/{folder_name}?{sas_token}` to container_url
    
    remote_script = f"""
    # 1. Setup AzCopy
    if [ ! -f ./azcopy ]; then
        echo "Installing AzCopy..."
        wget -q -O azcopy_v10.tar.gz https://aka.ms/downloadazcopy-v10-linux
        tar -xf azcopy_v10.tar.gz --strip-components=1
        chmod +x azcopy
    fi

    # 2. Run Copy
    SOURCE_DIR="{remote_base}/{folder_name}"
    DEST_URL="{container_url}/{folder_name}?{sas_token}"
    INCLUDE_PATTERN="{file_pattern}"
    
    echo "Source: $SOURCE_DIR"
    echo "Dest:   {container_url}/{folder_name}?[SAS_HIDDEN]"
    
    # Check if source exists
    if [ ! -d "$SOURCE_DIR" ]; then
        echo "❌ Source directory not found: $SOURCE_DIR"
        exit 1
    fi

    # Run AzCopy
    # recursive=true is needed even for flat files if we point to a dir?
    # include-pattern filters filenames.
    
    ./azcopy copy "$SOURCE_DIR" "$DEST_URL" \\
        --recursive=true \\
        --include-pattern "$INCLUDE_PATTERN" \\
        --output-type=text
    """
    
    # Run via SSH
    if not run_command(ssh_base + [remote_script], capture_output=False):
        print(f"❌ Archive failed for {basename}")
    else:
        print(f"✅ Archive successful for {basename}")

def main():
    parser = argparse.ArgumentParser(description="Archive simulation results from Azure VMs to Blob Storage.")
    parser.add_argument("--ip", help="Specific VM IP Address to archive. If omitted, all VMs in TSV are processed.")
    parser.add_argument("--all", action="store_true", help="If set, archive ALL files in the simulation folder instead of the selected subset.")
    args = parser.parse_args()

    # 1. Load Env
    if not os.path.exists(ENV_FILE):
        print(f"Error: {ENV_FILE} not found.")
        sys.exit(1)
    config = load_env(ENV_FILE)
    
    ssh_user = config.get("SSH_USER")
    ssh_pass = config.get("SSH_PASSWORD")
    
    has_sshpass = False
    if subprocess.call(["which", "sshpass"], stdout=subprocess.DEVNULL) == 0:
        has_sshpass = True
    else:
        print("Warning: sshpass not found.")

    # 2. Get Storage Info & SAS
    account_name, account_key = get_storage_info(config)
    sas_token = generate_sas_token(account_name, account_key, CONTAINER_NAME)
    
    # Base Container URL (No SAS yet)
    container_url = f"https://{account_name}.blob.core.windows.net/{CONTAINER_NAME}"
    
    # 3. File Pattern
    if args.all:
        print("ℹ️  Mode: Archiving ALL files.")
        file_pattern = "*"
    else:
        # Default subset
        # npt.gro, md.cpt, md.edr, md.gro, md.log, md.mdp, md.tpr, md.trr, md_noPBC.xtc, md_prev.cpt, mdout.mdp, *.xvg, *.log
        file_list = [
            "npt.gro", "md.cpt", "md.edr", "md.gro", "md.log", "md.mdp", 
            "md.tpr", "md.trr", "md_noPBC.xtc", "md_prev.cpt", "mdout.mdp", "*.xvg", "*.log"
        ]
        file_pattern = ";".join(file_list)

    # 4. Iterate VMs
    if not os.path.exists(TSV_FILE):
        print(f"Error: {TSV_FILE} not found.")
        sys.exit(1)

    print(f"Reading from {TSV_FILE}...")
    
    found_vm = False
    with open(TSV_FILE) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        
        for row in reader:
            vm_id = row['vm_id']
            basename = row['basename']
            vm_public_ip = row.get('vm_public_ip')

            if not vm_public_ip:
                print(f"Skipping {vm_id}: No IP.")
                continue

            # Filter by IP if argument provided
            if args.ip and args.ip != vm_public_ip:
                continue

            found_vm = True
            # Pass container URL and SAS token separately
            archive_vm(vm_id, basename, vm_public_ip, ssh_user, ssh_pass, has_sshpass, container_url, sas_token, file_pattern)

    if args.ip and not found_vm:
        print(f"❌ No VM found in TSV with IP: {args.ip}")
    else:
        print("\n=== Archive Process Complete ===")

if __name__ == "__main__":
    main()
