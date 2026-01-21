#!/bin/bash

#rm -rf *.tpr *.xtc \#*

source /DAMM/software/gromacs/trixie/AVX_512/gromacs-2024.5/bin/GMXRC

gmx grompp -f min.mdp -p topol.top -c start.gro -o min.tpr -maxwarn 1 -r start.gro
gmx  mdrun -nt 8  -v -deffnm min 
gmx  grompp -f eq.mdp -p topol.top -c min.gro    -maxwarn 1 -o eq.tpr -r start.gro
gmx  mdrun -nt 8  -v -deffnm eq 
gmx grompp -f md.mdp -c eq.gro -p topol.top  -maxwarn 1
gmx mdrun -v -nt 8
