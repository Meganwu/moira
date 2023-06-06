#! /bin/bash
filename=`cat  /data/mdml/neutral_collection/md_222.txt`
cd   /data/mdml/neutral_collection/neutral_2a
for  i in $filename
do
  cd ${i}MD_neu_2a
  cp   /data/mdml/neutral_collection/result_sh_all/total_interaction_input/total.cpptraj   .
  iscl=`grep Cl-  complex_vdw_bonded.pdb`
  isna=`grep Na+  complex_vdw_bonded.pdb`

  if [ -n "$iscl" ]; then
     noprostart=`grep "Cl-"  complex_vdw_bonded.pdb | head -n 1 | awk '{print $5}'`
     let proend=$noprostart-1
  elif  [ -n "$isna" ]; then
     noprostart=`grep Na+  complex_vdw_bonded.pdb  | head -1  |  awk '{print $5}'`
     let proend=$noprostart-1
  else
     noprostart=`grep "WAT"  complex_vdw_bonded.pdb | head -n 1 | awk '{print $5}'`
     let proend=$noprostart-1
  fi

  proend=$(echo $proend)
  sed -i 's/:1/:1-'"$proend"'/g'   total.cpptraj
  sed -i 's/:3com/:1-'"$proend"'/g'   total.cpptraj
  sed -i 's/:2pro/:2-'"$proend"'/g'   total.cpptraj
  sed -i 's/:4lig/:1/g'   total.cpptraj

  cpptraj total.cpptraj
  cd   /data/mdml/neutral_collection/neutral_2a
done
