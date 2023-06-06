#! /bin/bash
filename=`cat  fine.txt`
cd    /data/mdml/Data/2021-9-10-amber 
for  i in $filename
do
  cd ${i}MD_5a

  cp  ../total.cpptraj     .  

  iscl=`grep Cl-  complex_vdw_bonded.pdb`
  isna=`grep Na+  complex_vdw_bonded.pdb`
  if [ -n "$iscl" ]; then
     noprostart=`grep "Cl-"  complex_vdw_bonded.pdb | head -n 1 | awk '{print $5}'`
     let proend=$noprostart-1    
  elif  [ -n "$isna" ]; then
     noprostart=`grep "Na+"  complex_vdw_bonded.pdb | head -n 1 | awk '{print $5}'`
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
  cd   /data/mdml/Data/2021-9-10-amber  
done 

