"pmemd.sh" 

#!/bin/bash

pmemd16='mpirun -np 4  pmemd.MPI'

$pmemd16 -O -i min.in -o min.out -p complex_vdw_bonded.prmtop -c complex_vdw_bonded.inpcrd -r complex_min.rst -ref complex_vdw_bonded.inpcrd

$pmemd16 -O -i heat.in -o heat.out -p complex_vdw_bonded.prmtop -c complex_min.rst -r complex_heat.rst -x complex_heat.mdcrd -ref complex_min.rst

$pmemd16 -O -i density.in -o density.out -p complex_vdw_bonded.prmtop -c complex_heat.rst -r complex_density.rst -x complex_density.mdcrd -ref complex_heat.rst

pmemd.cuda  -O -i equil.in -o equil.out -p complex_vdw_bonded.prmtop  -c complex_density.rst -r complex_equil.rst -x complex_equil.mdcrd

pmemd.cuda  -O -i equil2.in -o equil2.out -p complex_vdw_bonded.prmtop  -c complex_equil.rst -r complex_equil2.rst -x complex_equil2.mdcrd

pmemd.cuda  -O -i equil3.in -o equil3.out -p complex_vdw_bonded.prmtop  -c complex_equil2.rst -r complex_equil3.rst -x complex_equil3.mdcrd
