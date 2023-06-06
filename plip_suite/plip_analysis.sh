#! /bin/bash
filename=`cat /data/mdml/neutral_collection/result_sh_all/md_list_split/list_00 `
cd  /data/mdml/neutral_collection/neutral_2a
for  i in $filename
do
cd ${i}MD_neu_2a/trajs

for j in $(seq  1    2550)
do
        sed   -i   '/Na+/d'   ${i}_native.pdb.${j}
        plip  -f   ${i}_native.pdb.${j}    -x      -o  ${i}_native_f_${j}    --name  ${i}_native_f_${j}
        rm   ${i}_native_f_${j}/*protonated.pdb  ${i}_native_f_${j}/plipfixe*pdb
done

cd   /data/mdml/neutral_collection/neutral_2a
done
