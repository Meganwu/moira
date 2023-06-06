id_num=`cat confs_good_ls.txt | awk '{print $1}'`

pdb_list=`cat confs_good_ls.txt | awk '{print $2}'`
cd step2_ligand_mol2
for id in $pdb_list
do
dockname=${id}dock
ref_ligand=${id}_ligand.mol2
cd $dockname
  for cid in $(seq 0 1 9)
  do
  num_dock=${id}${cid}
  cd $num_dock
  cd result-vina-final
      for cid_dock in $(seq 1 1 10)
      do
        rmsd_ligand=${num_dock}-cmd-${cid_dock}.pdb
        a=`isoRMSD.py  -r ../../$ref_ligand  -p  $rmsd_ligand`
        align=`echo $a | awk '{print $13}'`
        noalign=`echo $a | awk '{print $14}'`
        echo $id  $num_dock  $rmsd_ligand   $align  $noalign >>  ../../../../dock_all_rmsd.txt
      done
  cd ../../
  done
cd ..
done

