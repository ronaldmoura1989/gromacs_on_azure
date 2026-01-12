#!/bin/bash
# Resume NVT Equilibration step and continue workflow

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <SIMULATION_DIR>"
    exit 1
fi

SIM_DIR="$1"
cd "$SIM_DIR"

# Source GROMACS environment
source /usr/local/gromacs/bin/GMXRC

# GMX Flags for Standard_D64alds_v6 (64 vCPUs)
GMX_FLAGS="-ntmpi 8 -ntomp 8 -nb cpu -dlb yes"
GMX_FLAGS_TUNE="-ntmpi 8 -ntomp 8 -nb cpu -dlb yes -tunepme yes"

# Run continuation workflow in background
nohup bash -c "
    set -e
    cd '$SIM_DIR'
    source /usr/local/gromacs/bin/GMXRC

    echo '[$(date)] Resuming NVT equilibration...' >> workflow.log
    gmx mdrun -deffnm nvt -cpi nvt.cpt $GMX_FLAGS -v >> nvt_resume.log 2>&1

    echo '[$(date)] NVT completed, starting NPT equilibration...' >> workflow.log
    gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 10 >> workflow.log 2>&1
    gmx mdrun -deffnm npt $GMX_FLAGS -v >> npt.log 2>&1

    echo '[$(date)] NPT completed, starting Production MD...' >> workflow.log
    gmx grompp -f md.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 10 >> workflow.log 2>&1
    gmx mdrun -deffnm md $GMX_FLAGS_TUNE -v >> md.log 2>&1

    echo '[$(date)] MD completed, starting analysis...' >> workflow.log
    # Get basename from directory name
    BASENAME=\$(basename '$SIM_DIR')

    # PBC correction
    printf '0\n' | gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc nojump -ur compact >> workflow.log 2>&1

    # RMSD (Backbone)
    printf '4\n4\n' | gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns >> workflow.log 2>&1

    # Gyrate
    printf '1\n' | gmx gyrate -s md.tpr -f md_noPBC.xtc -o \${BASENAME}_gyrate.xvg >> workflow.log 2>&1

    # SASA
    printf '1\n' | gmx sasa -s md.tpr -f md_noPBC.xtc -o \${BASENAME}_sasa.xvg -tu ns >> workflow.log 2>&1

    # RMSF
    printf '1\n' | gmx rmsf -s md.tpr -f md_noPBC.xtc -res -o \${BASENAME}_rmsf.xvg >> workflow.log 2>&1

    echo '[$(date)] ✓ Complete simulation pipeline finished!' >> workflow.log
" >> workflow.log 2>&1 &

echo "NVT resumed with PID $! - Full workflow will continue automatically"
