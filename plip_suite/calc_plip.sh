#! /bin/bash
filename=`cat  fine_2a.txt `
cd  /data/mdml/neutral_collection/neutral_2a
for  i in $filename
do
  cd ${i}MD_neu_2a
  cp  /data/mdml/neutral_collection/result_sh_all/save_traj_plip/save_traj_5      .
  sed -i 's/names/'"$i"'_native.pdb/g'   save_traj_5 
  mkdir  trajs
  cpptraj   save_traj_5

  cd   /data/mdml/neutral_collection/neutral_2a
done
