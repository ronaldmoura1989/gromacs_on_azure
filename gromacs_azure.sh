#!/bin/bash
set -e

# Arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <BASENAME> <PDB_FILE>"
    exit 1
fi

BASENAME="$1"
PDB_FILE="$2"

echo "=== Processing $BASENAME ==="
echo "Input PDB: $PDB_FILE"
echo "Running on $(hostname) with 64 vCPUs optimization"

# GMX Flags for Standard_D64alds_v6 (64 vCPUs)
# Using 8 MPI ranks with 8 OpenMP threads each = 64 threads total.
# 8x8 is often optimal for NUMA nodes on 64-core EPYC.
# -dlb yes: Enable dynamic load balancing for better CPU utilization
# -tunepme yes: Auto-tune PME grid vs direct space split for optimal performance
GMX_FLAGS="-ntmpi 8 -ntomp 8 -nb cpu -dlb yes"
GMX_FLAGS_TUNE="-ntmpi 8 -ntomp 8 -nb cpu -dlb yes -tunepme yes"

# 1. Topology Generation
echo "-> Generating topology..."
# Force field 13: GROMOS96 53a6
# Added -vsite hydrogens to allow 5fs time step (requires update in mdp file dt=0.005)
printf '13\n' | gmx pdb2gmx -f $PDB_FILE -o ${BASENAME}_processed.gro -water spce -ignh

# 2. Box Definition
echo "-> Defining box..."
gmx editconf -f ${BASENAME}_processed.gro -o ${BASENAME}_newbox.gro -c -d 1.0 -bt cubic

# 3. Solvation
echo "-> Solvating..."
gmx solvate -cp ${BASENAME}_newbox.gro -cs spc216.gro -o ${BASENAME}_solv.gro -p topol.top

# 4. Ions
echo "-> Adding ions..."
gmx grompp -f ions.mdp -c ${BASENAME}_solv.gro -p topol.top -o ions.tpr -maxwarn 10
printf '13\n' | gmx genion -s ions.tpr -o ${BASENAME}_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# 5. Energy Minimization
echo "-> Energy Minimization..."
gmx grompp -f em.mdp -c ${BASENAME}_solv_ions.gro -p topol.top -o em.tpr -maxwarn 10
gmx mdrun -deffnm em $GMX_FLAGS -v

# 6. NVT Equilibration
echo "-> NVT Equilibration..."
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 10
gmx mdrun -deffnm nvt $GMX_FLAGS -v

# 7. NPT Equilibration
echo "-> NPT Equilibration..."
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 10
gmx mdrun -deffnm npt $GMX_FLAGS -v

# 8. Production MD
echo "-> Production MD..."
gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 10
gmx mdrun -deffnm md $GMX_FLAGS_TUNE -v

# 9. Analysis
echo "-> Analysis..."
# PBC correction
printf '0\n' | gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc nojump -ur compact

# RMSD (Backbone)
printf '4\n4\n' | gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns

# Gyrate
printf '1\n' | gmx gyrate -s md.tpr -f md_noPBC.xtc -o ${BASENAME}_gyrate.xvg

# SASA
printf '1\n' | gmx sasa -s md.tpr -f md_noPBC.xtc -o ${BASENAME}_sasa.xvg -tu ns

# RMSF
printf '1\n' | gmx rmsf -s md.tpr -f md_noPBC.xtc -res -o ${BASENAME}_rmsf.xvg

echo "=== Simulation Complete for $BASENAME ==="
