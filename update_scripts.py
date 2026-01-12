#!/usr/bin/env python3
"""
Update Scripts on Remote VMs

This script updates configuration files and scripts on running VMs:
- Uploads updated gromacs_azure.sh to simulation directories
- Uploads updated .mdp files to simulation directories
- Uploads resume scripts to /data/scripts/
- Optionally stops current jobs and resumes from checkpoints

Usage:
    python3 update_scripts.py                    # Just upload files
    python3 update_scripts.py --restart          # Upload and restart jobs from checkpoints
"""

import csv
import os
import subprocess
import sys
import argparse

# Configuration
ENV_FILE = ".env"
TSV_FILE = "gromacs_runs.tsv"
SCRIPT_PATH = "gromacs_azure.sh"
PARAMS_DIR = "gromacs_parameters"
RESUME_SCRIPTS_DIR = "scripts"

def load_env(env_path):
    """Load environment variables from .env file"""
    config = {}
    if os.path.exists(env_path):
        with open(env_path) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    key, value = line.strip().split('=', 1)
                    config[key] = value
    return config

def run_command(cmd, shell=False, ignore_errors=False):
    """Helper to run shell command and print output"""
    print(f"CMD: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    try:
        subprocess.check_call(cmd, shell=shell)
    except subprocess.CalledProcessError as e:
        if not ignore_errors:
            print(f"Error running command: {e}")
            return False
    return True

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Update scripts on remote VMs')
    parser.add_argument('--restart', action='store_true',
                       help='Stop current jobs and restart from checkpoints')
    args = parser.parse_args()

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

    # Verify local files exist
    script_dir = os.path.dirname(os.path.abspath(__file__))
    local_params_dir = os.path.join(script_dir, PARAMS_DIR)
    local_resume_dir = os.path.join(script_dir, RESUME_SCRIPTS_DIR)

    if not os.path.exists(SCRIPT_PATH):
        print(f"Error: {SCRIPT_PATH} not found in current directory.")
        sys.exit(1)
    if not os.path.exists(local_params_dir):
        print(f"Error: {PARAMS_DIR}/ directory not found.")
        sys.exit(1)
    if not os.path.exists(local_resume_dir):
        print(f"Error: {RESUME_SCRIPTS_DIR}/ directory not found.")
        sys.exit(1)

    # 2. Read TSV and Update VMs
    if not os.path.exists(TSV_FILE):
        print(f"Error: {TSV_FILE} not found.")
        sys.exit(1)

    with open(TSV_FILE) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        for row in reader:
            vm_id = row['vm_id']
            basename = row['basename']

            # Use IP directly from TSV
            if 'vm_public_ip' not in row or not row['vm_public_ip']:
                print(f"Skipping {basename}: No IP found in TSV for {vm_id}.")
                continue

            ip_address = row['vm_public_ip']
            print(f"\n{'='*70}")
            print(f"=== Updating {vm_id} ({ip_address}) ===")
            print(f"{'='*70}")

            # Remote Paths
            folder_name = basename.replace(".pdb", "").replace(".PDB", "")
            remote_sim_dir = f"/data/simulations/{folder_name}"
            remote_scripts_dir = "/data/scripts"

            # Construct SSH/SCP base commands
            ssh_base = ["ssh", "-o", "StrictHostKeyChecking=no", f"{ssh_user}@{ip_address}"]
            scp_base = ["scp", "-o", "StrictHostKeyChecking=no"]

            if has_sshpass and ssh_pass:
                ssh_wrapper = ["sshpass", "-p", ssh_pass]
                ssh_base = ssh_wrapper + ssh_base
                scp_base = ssh_wrapper + scp_base

            # A. Create /data/scripts directory if it doesn't exist
            print("Ensuring /data/scripts directory exists...")
            run_command(ssh_base + [f"mkdir -p {remote_scripts_dir}"], ignore_errors=True)

            # B. Upload Resume Scripts to /data/scripts/
            print("Uploading resume scripts to /data/scripts/...")
            local_resume_wildcard = os.path.join(local_resume_dir, "*.sh")
            cmd_upload_resume = f"{' '.join(scp_base)} {local_resume_wildcard} {ssh_user}@{ip_address}:{remote_scripts_dir}/"
            if not run_command(cmd_upload_resume, shell=True):
                print(f"Warning: Failed to upload resume scripts to {vm_id}")

            # Make resume scripts executable
            run_command(ssh_base + [f"chmod +x {remote_scripts_dir}/*.sh"], ignore_errors=True)

            # C. Upload gromacs_azure.sh to simulation directory
            print(f"Uploading {SCRIPT_PATH} to {remote_sim_dir}/...")
            if not run_command(scp_base + [SCRIPT_PATH, f"{ssh_user}@{ip_address}:{remote_sim_dir}/"]):
                print(f"Warning: Failed to upload {SCRIPT_PATH} to {vm_id}")

            # Make gromacs_azure.sh executable
            run_command(ssh_base + [f"chmod +x {remote_sim_dir}/{SCRIPT_PATH}"], ignore_errors=True)

            # D. Upload Parameter files (.mdp) to simulation directory
            print(f"Uploading .mdp files to {remote_sim_dir}/...")
            local_mdp_wildcard = os.path.join(local_params_dir, "*.mdp")
            cmd_upload_mdp = f"{' '.join(scp_base)} {local_mdp_wildcard} {ssh_user}@{ip_address}:{remote_sim_dir}/"
            if not run_command(cmd_upload_mdp, shell=True):
                print(f"Warning: Failed to upload .mdp files to {vm_id}")

            # E. Optionally restart simulations
            if args.restart:
                print(f"Stopping current jobs and resuming from checkpoint on {vm_id}...")
                remote_cmd = f"{remote_scripts_dir}/resume_job.sh --stop"
                run_command(ssh_base + [remote_cmd], ignore_errors=True)

            print(f"✓ Update complete for {vm_id}")

    print("\n" + "="*70)
    print("=== All VMs Updated ===")
    print("="*70)

    if args.restart:
        print("\nSimulations have been restarted from checkpoints with updated scripts.")
    else:
        print("\nFiles updated. To restart simulations with new settings, run:")
        print("  python3 update_scripts.py --restart")
        print("\nOr manually on each VM:")
        print("  ssh user@vm-ip '/data/scripts/resume_job.sh --stop'")

if __name__ == "__main__":
    main()
