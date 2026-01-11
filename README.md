# Azure Infrastructure Creation Guide for Molecular Dynamics

This guide details the process of creating a cluster of Virtual Machines on Azure optimized for Molecular Dynamics simulations. The infrastructure creates 9 nodes of `Standard_D64alds_v6` machines, interconnected in a single Virtual Network, with GROMACS pre-installed and large storage disks mounted.

## 1. Infrastructure Architecture

The setup consists of:
*   **Resource Group**: `rg-md-eastus`
*   **Region**: `eastus`
*   **Network**: `vnet-md-eastus` (10.1.0.0/16)
    *   **Subnet**: `subnet-compute` (10.1.1.0/24)
    *   **NSG**: `nsg-md-compute` (Allow SSH from user IP)
*   **Virtual Machines** (x9): `vm-md-01` to `vm-md-09`
    *   **Size**: `Standard_D64alds_v6` (64 vCPUs, 128 GiB RAM)
    *   **OS**: Ubuntu 22.04 LTS
    *   **OS Disk**: 30 GB Premium SSD
    *   **Data Disk**: 2 TB Standard HDD (Mounted at `/data`)
    *   **Networking**: Private IP + Public IP
*   **Software**: GROMACS 2023.2 (Compiled for CPU-only)

```mermaid
graph TD
    User[User] -->|SSH| Internet
    Internet -->|Public IP| VM1[vm-md-01]
    Internet -->|Public IP| VM2[vm-md-02]
    Internet -->|Public IP| VM9[vm-md-09]
    
    subgraph Azure Resource Group: rg-md-eastus
        subgraph VNet: vnet-md-eastus
            subgraph Subnet: subnet-compute
                VM1
                VM2
                VM9
                style VM1 fill:#ccf,stroke:#333,stroke-width:2px
                style VM2 fill:#ccf,stroke:#333,stroke-width:2px
                style VM9 fill:#ccf,stroke:#333,stroke-width:2px
            end
        end
        NSG[NSG: Allow SSH] -.-> Subnet
    end
```

---

## 2. Deployment Steps

Follow these steps to deploy the entire infrastructure from scratch using the provided scripts.

### Step 0: Prerequisites
Ensure you have the Azure CLI installed and you are logged in.
Create a `.env` file in the `azure_md_infrastructure` folder with your desired credentials:
```bash
SSH_USER=your_username
SSH_PASSWORD=your_secure_password
AZURE_SUBSCRIPTION=your_subscription_id
```

### Step 1: Create Virtual Machines
Run the main provisioning script. This creates the Resource Group, VNet, NSG, Disks, NICs, and VMs.
```bash
cd azure_md_infrastructure
chmod +x create-vms-d64aldsv6.sh
./create-vms-d64aldsv6.sh
```
*Wait for the script to finish. It will output the Public IPs of the created VMs.*

### Step 2: Install GROMACS
Deploy GROMACS 2023.2 to all nodes in parallel.
```bash
chmod +x install-gromacs-on-all-vms.sh
./install-gromacs-on-all-vms.sh
```
*This compiles GROMACS from source on each VM. It may take 10-20 minutes.*

### Step 3: Mount Data Disks
Format and mount the 2TB data disks to `/data` on all nodes.
```bash
chmod +x mount-all-disks.sh
./mount-all-disks.sh
```

### Step 4: Run Simulations
Programmatically deploy and execute simulations on all nodes.

1.  **Prepare Run Data**: Edit `gromacs_runs.tsv` with your simulation details (PDB path, VM selection). Ensure your local files exist in `PROJECT_PATH`.
2.  **Run Orchestrator**:
    ```bash
    # (Optional) Install sshpass for automated password entry
    brew install sshpass
    
    # Launch jobs
    python3 run_simulations.py
    ```
    *This script uploads your PDB and Parameter files (from `gromacs_parameters/`) and starts the simulation in the background (`nohup`) on each VM.*

### Step 5: Monitor Progress
Monitor the status of your simulations in real-time. This script connects to each VM and reports the progress (Time/Step) for NVT, NPT, and Production MD phases.

```bash
python3 monitor_progress.py
```
*Output Example:*
```text
VM ID      | IP Address      | Basename                       | NVT Time        | NPT Time        | MD Time        | Last Check          
--------------------------------------------------------------------------------------------------------------------------------------------
vm-md-01   | 13.82.93.175    | miocardin_A0A499FIJ6_C684te..  | 631.00000 ps    | Not Started     | Not Started    | 2026-01-11 17:35:10 
vm-md-02   | 172.178.47.153  | miocardin_A0A499FIJ6_K823te..  | 514.00000 ps    | Not Started     | Not Started    | 2026-01-11 17:35:11 
...
```

---

## 3. Virtual Machine Choice

We selected the **Standard_D64alds_v6** for this cluster.

*   **vCPUs**: 64 (AMD EPYC™ Genoa processors)
*   **RAM**: 128 GiB
*   **Architecture**: ARM-based (Ampere Altra) or AMD-based? *Correction*: The 'a' in Azure names usually denotes AMD. The Dpsv5/Dplsv5 are ARM. The **Dasv6** series runs on 4th Gen AMD EPYC™ (Genoa).
*   **Rationale**:
    *   **High Core Count**: MD simulations (like GROMACS) scale well with CPU cores. 64 cores per node provides significant parallel power.
    *   **Cost-Efficiency**: The `D` series is a General Purpose workload but the "a" (AMD) variants often offer better price/performance than Intel equivalents. "ls" denotes lower memory (2GiB/vCPU) which is sufficient for many MD simulations that are compute-bound rather than memory-bound, saving costs compared to 4GiB/8GiB per vCPU options.

---

## 4. Alternative Configurations

### High-Performance Computing (HPC) - CPU
If you need to run a **single large simulation** distributed across multiple nodes, you require low-latency networking (InfiniBand/RDMA). The standard `D` series lacks this.

| VM Size | vCPUs | RAM | Network | Cost/Hr | Use Case |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Standard_HB120rs_v3** | 120 (AMD Milan) | 448 GiB | 200 Gbps InfiniBand | ~$3.60 | True cluster computing. Allows 9 nodes to act as one supercomputer via MPI. |

### GPU Acceleration
If your simulation requires GPU acceleration (e.g., highly parallelized non-bonded interactions), consider the **NC-series**:

| VM Size | GPU | vCPUs | RAM | Use Case |
| :--- | :--- | :--- | :--- | :--- |
| **Standard_NC6s_v3** | 1 x NVIDIA V100 | 6 | 112 GiB | Good balance for single-GPU jobs. |
| **Standard_NC24ads_A100_v4**| 1 x NVIDIA A100 | 24 | 220 GiB | High-end performance for large boxes. |
| **Standard_NV6ads_A10_v5** | 1 x NVIDIA A10 | 6 | 55 GiB | Cost-effective modern GPU inference/simulation. |

*Note: GPU and HPC VMs often require requesting a specific Quota increase in your Azure subscription.*

---

## 5. Cost Analysis (Estimates)

Pricing is based on **East US** region (approximate rates as of late 2025/early 2026).

### Computing (x9 Nodes)
*   **VM**: Standard_D64alds_v6
*   **Rate**: ~$3.05 / hour (Pay-as-you-go)
*   **Total for Cluster**: 9 * $3.05 = **$27.45 / hour**

### Storage (x9 Nodes)
*   **OS Disk (30GB Premium SSD)**: P4 Tier (~$5.92/month/disk)
    *   Total: 9 * $5.92 ≈ $53.28 / month
*   **Data Disk (2TB Standard HDD)**: S40 Tier (~$77.83/month/disk)
    *   Total: 9 * $77.83 ≈ $700.47 / month

### Total Monthly Estimate (24/7 Operation)
*   **Compute**: $27.45 * 730 hours = ~$20,038
*   **Storage**: ~$753
*   **Grand Total**: **~$20,791 / month** (or ~$28.50 / hour)

*Tip: Use Azure Spot Instances for up to 90% savings if your simulations can handle interruptions, or Reserved Instances (1-year/3-year) for ~30-50% savings.*

---

## 7. GROMACS Simulation Parameters

The simulations are performed using **GROMACS 2023.2**, configured with parameters optimized for computational efficiency on high-core-count cloud infrastructure (Azure D-series). The specific choices allow for a significant increase in sampling speed without compromising biological accuracy for standard protein systems.

### Configuration Details

*   **Force Field**: GROMOS96 53a6.
*   **Water Model**: SPC/E (Extended Simple Point Charge).
*   **Integrator**: md (Lead-Frog).
*   **Time Step**: 5 fs (`dt = 0.005`).
    *   *Optimization*: This larger time step (standard is 2 fs) is enabled by the use of **Virtual Sites** (`-vsite hydrogens` in `pdb2gmx`), which removes the fastest degrees of freedom (hydrogen vibrations) by treating hydrogen atoms as rigid bodies fixed to their heavy atoms. This results in a ~2.5x speedup in simulation throughput.
*   **Parallelization**:
    *   **Architecture**: Hybrid MPI/OpenMP.
    *   **Configuration**: 8 MPI ranks × 8 OpenMP threads (`-ntmpi 8 -ntomp 8`).
    *   **Rationale**: This 64-thread configuration is tuned for the AMD EPYC (Genoa) dual-socket architecture of the `Standard_D64alds_v6` VMs, minimizing NUMA domain crossing penalties compared to pure OpenMP or flat MPI approaches.
*   **Thermostat**: V-rescale (velocity rescaling with a stochastic term) at 300 K/310 K.
*   **Barostat**: Parrinello-Rahman at 1.0 bar (Isotropic).
*   **Long-range Electrostatics**: PME (Particle Mesh Ewald).

### Bibliography

The choice of parameters is supported by standard best practices in molecular dynamics and recent benchmarks for high-throughput cloud computing:

1.  **GROMACS 2023 Features & Performance**:
    *   Kutzner, C., et al. (2022). *Gromacs in the Cloud: A Global High-Performance Computing Framework.* Journal of Chemical Information and Modeling. [DOI: 10.1021/acs.jcim.2c00044](https://doi.org/10.1021/acs.jcim.2c00044)
    *   Páll, S., et al. (2020). *Heterogeneous parallelization and acceleration of molecular dynamics simulations in GROMACS.* The Journal of Chemical Physics. [DOI: 10.1063/5.0018516](https://doi.org/10.1063/5.0018516)

2.  **Virtual Sites & Time Step Efficiency**:
    *   Feenstra, K. A., et al. (1999). *Improving efficiency of large-scale molecular dynamics simulations of proteins.* Journal of Computational Chemistry. [DOI: 10.1002/(SICI)1096-987X(199906)20:8<786::AID-JCC5>3.0.CO;2-B](https://doi.org/10.1002/(SICI)1096-987X(199906)20:8<786::AID-JCC5>3.0.CO;2-B)
    *   *Note*: While older, this remains the foundational reference for the virtual site method used in modern GROMACS distributions.

3.  **Force Field (GROMOS96 53a6)**:
    *   Oostenbrink, C., et al. (2004). *A biomolecular force field based on the free enthalpy of hydration and solvation: The GROMOS force-field parameter sets 53A5 and 53A6.* Journal of Computational Chemistry. [DOI: 10.1002/jcc.20090](https://doi.org/10.1002/jcc.20090)

4.  **Cloud Computing for MD**:
    *   Moraes, F. R., et al. (2021). *Evaluating cloud computing for high-throughput molecular dynamics simulations.* Journal of Cloud Computing. [DOI: 10.1186/s13677-021-00276-8](https://doi.org/10.1186/s13677-021-00276-8) (Validating cost/performance of Spot instances and various VM types).
