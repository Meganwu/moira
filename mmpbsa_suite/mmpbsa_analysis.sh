#! /bin/bash
filename=`cat  md_222.txt `
cd    /data/mdml/neutral_collection/neutral_2a 
for  i in $filename
do
  cd ${i}MD_neu_2a

  cp   /data/mdml/neutral_collection/result_sh_all/mmpbsa/*       .

  
  ligfirst=`sed -n '3p' complex_vdw_bonded.pdb | awk '{print $4}'`
  ligsec=`grep $ligfirst complex_vdw_bonded.pdb | tail -n 1|  awk '{print $2}'`
  let prostart=$ligsec+1
  iscl=`grep Cl-  complex_vdw_bonded.pdb`
  isna=`grep Na+  complex_vdw_bonded.pdb`
  if [ -n "$iscl" ]; then
     noprostart=`grep "Cl-"  complex_vdw_bonded.pdb | head -n 1 | awk '{print $2}'`
     let proend=$noprostart-1

  elif  [ -n "$isna" ]; then
     noprostart=`grep "Na+"  complex_vdw_bonded.pdb | head -n 1 | awk '{print $2}'`
     let proend=$noprostart-1

  else
     noprostart=`grep "WAT"  complex_vdw_bonded.pdb | head -n 1 | awk '{print $2}'`
     let proend=$noprostart-1
  fi

  comall=`grep WAT complex_vdw_bonded.pdb  | tail -n 1 | awk '{print $2}'`
  comall=$(echo $comall)
  ligsec=$(echo $ligsec)
  prostart=$(echo $prostart)
  proend=$(echo $proend)
  sed -i 's/num_total/'"$comall"'/g' extract_coords.mmpbsa
  sed -i 's/lig_stop/'"$ligsec"'/g' extract_coords.mmpbsa
  sed -i 's/pro_start/'"$prostart"'/g' extract_coords.mmpbsa
  sed -i 's/pro_end/'"$proend"'/g' extract_coords.mmpbsa
  
  cpptraj   script_traj
  mm_pbsa.pl extract_coords.mmpbsa > extract_coords.log
  mm_pbsa.pl binding_energy.mmpbsa > binding_energy.log


  cd   /data/mdml/neutral_collection/neutral_2a 
done 

