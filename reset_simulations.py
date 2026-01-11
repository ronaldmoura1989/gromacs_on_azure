import csv
import os
import subprocess
import sys
import argparse

# Configuration
ENV_FILE = ".env"
TSV_FILE = "gromacs_runs.tsv"

def load_env(env_path):
    config = {}
    if os.path.exists(env_path):
        with open(env_path) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    key, value = line.strip().split('=', 1)
                    config[key] = value
    return config

def run_command(cmd):
    try:
        subprocess.check_call(cmd)
        return True
    except subprocess.CalledProcessError:
        return False

def main():
    parser = argparse.ArgumentParser(description="Reset simulations on Azure VMs (Kill GMX & Delete Data).")
    parser.add_argument("--confirm", action="store_true", help="Confirm execution without prompt.")
    args = parser.parse_args()

    if not args.confirm:
        print("WARNING: This will KILL all running simulations and DELETE ALL DATA in /data/simulations/ on all VMs.")
        response = input("Are you sure? (type 'yes' to proceed): ")
        if response.lower() != "yes":
            print("Aborted.")
            sys.exit(0)

    if not os.path.exists(ENV_FILE):
        print(f"Error: {ENV_FILE} not found.")
        sys.exit(1)
    
    config = load_env(ENV_FILE)
    ssh_user = config.get("SSH_USER")
    ssh_pass = config.get("SSH_PASSWORD")
    
    has_sshpass = False
    if subprocess.call(["which", "sshpass"], stdout=subprocess.DEVNULL) == 0:
        has_sshpass = True
    
    if not os.path.exists(TSV_FILE):
        print(f"Error: {TSV_FILE} not found.")
        sys.exit(1)

    with open(TSV_FILE) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        
        for row in reader:
            vm_id = row['vm_id']
            ip_address = row.get('vm_public_ip')

            if not ip_address:
                continue
            
            print(f"Resetting {vm_id} ({ip_address})...")
            
            ssh_base = ["ssh", "-o", "StrictHostKeyChecking=no", "-o", "ConnectTimeout=5", f"{ssh_user}@{ip_address}"]
            if has_sshpass and ssh_pass:
                ssh_wrapper = ["sshpass", "-p", ssh_pass]
                ssh_base = ssh_wrapper + ssh_base
            
            # The Reset Command
            # 1. Kill gmx_mpi and gmx
            # 2. Wait a second
            # 3. Delete the simulation folder contents (but keep the top level /data/simulations dir)
            reset_cmd = "pkill gmx_mpi; pkill gmx; echo 'Killed processes'; sleep 2; rm -rf /data/simulations/*; echo 'Cleaned data'"
            
            run_command(ssh_base + [reset_cmd])
            print("✓ Done")

if __name__ == "__main__":
    main()
