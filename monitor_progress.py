import csv
import os
import subprocess
import sys
import re
from datetime import datetime
import argparse
import json
import time

# Configuration
ENV_FILE = ".env"
TSV_FILE = "gromacs_runs.tsv"
CACHE_FILE = ".progress_cache.json"

def load_env(env_path):
    config = {}
    if os.path.exists(env_path):
        with open(env_path) as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    key, value = line.strip().split('=', 1)
                    config[key] = value
    return config

def load_cache():
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE, 'r') as f:
                return json.load(f)
        except:
            pass
    return {}

def save_cache(data):
    with open(CACHE_FILE, 'w') as f:
        json.dump(data, f)

def run_command(cmd, capture_output=False):
    if capture_output:
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return result.stdout.strip()
        except subprocess.CalledProcessError:
            return None
    else:
        try:
            subprocess.check_call(cmd)
            return True
        except subprocess.CalledProcessError:
            return False

def check_ssh_connection(ssh_base):
    return run_command(ssh_base + ["echo 'Connection OK'"], capture_output=True) is not None

def parse_last_step_time(log_content):
    """
    Parses the GROMACS log tail to find the last occurrence of:
           Step           Time
         290500      581.00000
    Returns ("Display String", FloatValue)
    """
    if not log_content or "File not found" in log_content:
        return ("Not Started", None)
    
    # regex to match the pattern. Since we get a tail, allow searching from end.
    # Pattern explanation:
    # Step\s+Time matches the header
    # \s*(\d+)\s+([\d\.]+) matches the values on the next line
    matches = re.findall(r"Step\s+Time\s*\n\s*(\d+)\s+([\d\.]+)", log_content)
    if matches:
        last_match = matches[-1]
        return (f"{last_match[1]} ps", float(last_match[1])) 
    
    # If process finished, it might look different, check for "Finished"
    if "Finished" in log_content:
        return ("Done", None)
        
    return ("Running (No Data)", None)

def calculate_eta_from_snapshot(vm_id, log_name, current_sim_time, cache):
    """
    Calculates ETA by comparing current sim time with previous cached time.
    """
    target_time_ps = 10000.0 # Default NVT/NPT
    if log_name == "md.log":
        target_time_ps = 100000.0

    if current_sim_time is None:
        return "N/A"

    current_wall_time = time.time()
    
    # Key for cache includes VM ID and Log Name (Phase)
    cache_key = f"{vm_id}_{log_name}"
    
    entry = cache.get(cache_key)
    eta_str = "Calc..."
    
    if entry:
        old_sim_time = entry.get('sim_time', 0.0)
        old_wall_time = entry.get('wall_time', 0.0)
        
        delta_sim = current_sim_time - old_sim_time
        delta_wall = current_wall_time - old_wall_time
        
        # Only calculate if time has moved forward and enough wall time passed (>1s)
        if delta_sim > 0 and delta_wall > 1.0:
            speed_ps_per_sec = delta_sim / delta_wall
            
            remaining_ps = target_time_ps - current_sim_time
            if remaining_ps <= 0:
                eta_str = "Done"
            else:
                remaining_seconds = remaining_ps / speed_ps_per_sec
                hours_left = remaining_seconds / 3600.0
                eta_str = f"{hours_left:.1f}h"
        else:
            # Keep old ETA if we haven't moved or checked too fast
            eta_str = entry.get('last_eta', "Calc...")
            
    # Update Cache (Always define new state)
    cache[cache_key] = {
        'sim_time': current_sim_time,
        'wall_time': current_wall_time,
        'last_eta': eta_str
    }
    
    return eta_str

def get_vm_status(vm_id, basename, ip_address, ssh_user, ssh_pass, has_sshpass, cache):
    # SSH Command Construction
    ssh_base = ["ssh", "-o", "StrictHostKeyChecking=no", "-o", "ConnectTimeout=5", f"{ssh_user}@{ip_address}"]
    if has_sshpass and ssh_pass:
        ssh_wrapper = ["sshpass", "-p", ssh_pass]
        ssh_base = ssh_wrapper + ssh_base

    # Remote Directory
    # Ensure we strip .pdb extension to match run_simulations.py behavior
    folder_name = basename.replace(".pdb", "").replace(".PDB", "")
    remote_base = "/data/simulations"
    remote_dir = f"{remote_base}/{folder_name}"
    
    # We check all logs
    logs_to_check = ["nvt.log", "npt.log", "md.log"]
    
    # Delimiter for separating log outputs
    delimiter = "||||"
    
    # Smaller tail is fine now since we parse the Step header
    remote_script = f"""
    if cd {remote_dir} 2>/dev/null; then
        for logfile in {' '.join(logs_to_check)}; do
            if [ -f "$logfile" ]; then
                tail -n 20 "$logfile"
            else
                echo "File not found"
            fi
            echo "{delimiter}" 
        done
    else
        echo "DIR_NOT_FOUND"
    fi
    """
    
    output = run_command(ssh_base + [remote_script], capture_output=True)
    
    results = {}
    active_eta = "N/A"

    if output and "DIR_NOT_FOUND" not in output:
        parts = output.split(delimiter)
        for i, log_name in enumerate(logs_to_check):
            if i < len(parts):
                status_str, time_val = parse_last_step_time(parts[i])
                results[log_name] = status_str
                
                # Determine ETA based on which log is active (Has time and not Done)
                # Priority: MD > NPT > NVT
                # We overwrite active_eta as we process the list in order (nvt -> npt -> md)
                # If a step is running (has time_val and not Done), we calculate its ETA.
                
                if time_val is not None:
                     eta_val = calculate_eta_from_snapshot(vm_id, log_name, time_val, cache)
                     if eta_val and eta_val != "Done" and eta_val != "N/A":
                         # Overwrite previous phases (NVT->NPT->MD)
                         active_eta = f"{log_name[:3].upper()}:{eta_val}"
            else:
                results[log_name] = "Error"
    else:
        for log_name in logs_to_check:
            results[log_name] = "Missing Dir" if output and "DIR_NOT_FOUND" in output else "Conn Error"

    return results, active_eta

def main():
    parser = argparse.ArgumentParser(description="Monitor GROMACS simulations on Azure VMs.")
    args = parser.parse_args()

    if not os.path.exists(ENV_FILE):
        print(f"Error: {ENV_FILE} not found.")
        sys.exit(1)
    
    config = load_env(ENV_FILE)
    ssh_user = config.get("SSH_USER")
    ssh_pass = config.get("SSH_PASSWORD")
    
    # Load Cache
    progress_cache = load_cache()
    
    has_sshpass = False
    if subprocess.call(["which", "sshpass"], stdout=subprocess.DEVNULL) == 0:
        has_sshpass = True
    else:
        print("Warning: 'sshpass' not found.")

    if not os.path.exists(TSV_FILE):
        print(f"Error: {TSV_FILE} not found.")
        sys.exit(1)

    print(f"\n{'VM ID':<10} | {'IP Address':<15} | {'Basename':<28} | {'NVT':<14} | {'NPT':<14} | {'Produc MD':<14} | {'Estimate Time':<15} | {'Last Check':<20}")
    print("-" * 155)
    
    with open(TSV_FILE) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        
        for row in reader:
            vm_id = row['vm_id']
            basename = row['basename']
            vm_public_ip = row.get('vm_public_ip')

            if not vm_public_ip:
                continue
            
            # Shorten basename for display
            display_name = (basename[:25] + '..') if len(basename) > 27 else basename
            
            metrics, eta = get_vm_status(vm_id, basename, vm_public_ip, ssh_user, ssh_pass, has_sshpass, progress_cache)
            now = datetime.now().strftime("%H:%M:%S")
            
            print(f"{vm_id:<10} | {vm_public_ip:<15} | {display_name:<28} | {metrics['nvt.log']:<14} | {metrics['npt.log']:<14} | {metrics['md.log']:<14} | {eta:<15} | {now:<20}")

    # Save Cache at end
    save_cache(progress_cache)

if __name__ == "__main__":
    main()
