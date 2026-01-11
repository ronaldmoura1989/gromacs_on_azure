import csv
import os
import subprocess
import time
import sys

# Configuration
ENV_FILE = ".env"
TSV_FILE = "gromacs_runs.tsv"
SCRIPT_PATH = "gromacs_azure.sh"

def load_env(env_path):
    config = {}
    if os.path.exists(env_path):
        with open(env_path) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    key, value = line.strip().split('=', 1)
                    config[key] = value
    return config

def run_command(cmd, shell=False):
    # Helper to run shell command and print output
    print(f"CMD: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    try:
        subprocess.check_call(cmd, shell=shell)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        return False
    return True

def main():
    # 1. Load Environment
    if not os.path.exists(ENV_FILE):
        print(f"Error: {ENV_FILE} not found.")
        sys.exit(1)
    
    config = load_env(ENV_FILE)
    project_path = config.get("PROJECT_PATH")
    ssh_user = config.get("SSH_USER")
    ssh_pass = config.get("SSH_PASSWORD") 
    
    # Check for sshpass
    has_sshpass = False
    if subprocess.call(["which", "sshpass"], stdout=subprocess.DEVNULL) == 0:
        has_sshpass = True
    else:
        print("Warning: 'sshpass' not found. You will be prompted for passwords for each VM unless SSH keys are set up.")
    
    # 2. Read TSV and Launch Jobs
    with open(TSV_FILE) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        
        for row in reader:
            vm_id = row['vm_id']
            basename = row['basename']
            filepath = row['filepath']
            
            # Use IP directly from TSV
            if 'vm_public_ip' not in row or not row['vm_public_ip']:
                print(f"Skipping {basename}: No IP found in TSV for {vm_id}.")
                continue
                
            ip_address = row['vm_public_ip']
            print(f"\n=== Launching {basename} on {vm_id} ({ip_address}) ===")
            
            # Local Paths
            local_pdb = os.path.join(project_path, filepath)
            # Use gromacs_parameters from the directory where this script resides
            script_dir = os.path.dirname(os.path.abspath(__file__))
            local_params_dir = os.path.join(script_dir, "gromacs_parameters/")
            
            # Remote Paths
            # Strip .pdb from folder name if present in basename
            folder_name = basename.replace(".pdb", "").replace(".PDB", "")
            remote_base = "/data/simulations"
            remote_dir = f"{remote_base}/{folder_name}"
            
            # Construct SSH/SCP base commands
            ssh_base = ["ssh", "-o", "StrictHostKeyChecking=no", f"{ssh_user}@{ip_address}"]
            scp_base = ["scp", "-o", "StrictHostKeyChecking=no"]
            
            if has_sshpass and ssh_pass:
                ssh_wrapper = ["sshpass", "-p", ssh_pass]
                ssh_base = ssh_wrapper + ssh_base
                scp_base = ssh_wrapper + scp_base
            
            # A. Create Remote Directory
            print("Creating remote directory...")
            run_command(ssh_base + [f"mkdir -p {remote_dir}"])
            
            # B. Upload PDB
            print(f"Uploading PDB: {filepath}...")
            if not run_command(scp_base + [local_pdb, f"{ssh_user}@{ip_address}:{remote_dir}/"]):
                continue

            # C. Upload Parameters (.mdp)
            print("Uploading Parameter files...")
            local_mdp_wildcard = os.path.join(local_params_dir, "*.mdp")
            cmd_upload_mdp = f"{' '.join(scp_base)} {local_mdp_wildcard} {ssh_user}@{ip_address}:{remote_dir}/"
            run_command(cmd_upload_mdp, shell=True)

            # D. Upload Script
            print("Uploading Job Script...")
            run_command(scp_base + [SCRIPT_PATH, f"{ssh_user}@{ip_address}:{remote_dir}/"])
            
            # E. Execute Job (nohup)
            print("Starting Simulation...")
            pdb_filename = os.path.basename(local_pdb)
            
            remote_cmd = (
                f"source /usr/local/gromacs/bin/GMXRC; "
                f"cd {remote_dir}; "
                f"chmod +x {os.path.basename(SCRIPT_PATH)}; "
                f"nohup ./{os.path.basename(SCRIPT_PATH)} {basename} {pdb_filename} > run.log 2>&1 & "
                f"echo 'Job submitted with PID $!'"
            )
            
            run_command(ssh_base + [remote_cmd])
            
            print(f"Job launched for {basename} on {vm_id}")

    print("\n=== All Jobs Submitted ===")

if __name__ == "__main__":
    main()
