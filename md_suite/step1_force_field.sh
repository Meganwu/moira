#!/bin/bash
### generate force field given initial ligand and protein structures ###

pdblist=`cat list_00`  # the folder path for each molecular dyanmic task

for id in $pdblist
do
     j=`echo ${id:0:4}`
     md_name2=${j}_2a
     cd   $md_name2
     lig_name=${id}
#     pro_name='clean.pdb'  # protein file name
#     reduce -Trim $pro_name  > clean.pdb  # remove hydrogens of protein
     antechamber -i $lig_name  -fi pdb -o lig.mol2 -fo mol2 -c bcc -s 2
     parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod
     tleap -f leap1.in  #generate the force field of ligand
     tleap -f leap2.in  #generate force fields of protein and complex
     cd  /data/wunian2/neutral_collection/neutral_2a
done
