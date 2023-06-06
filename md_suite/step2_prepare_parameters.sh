pdblist=`cat /data/wunian2/new-bags/list.txt`
cd  /data/wunian2/new-bags

for  id  in $pdblist
do

md_name=${id}MD_new

cd   $md_name
cp   /data/wunian2/new-bags/md_input/*     .
tleap -f leap2.in

iscl=`grep Cl-  complex_vdw_bonded.pdb`
isna=`grep Na+  complex_vdw_bonded.pdb`

if [ -n "$iscl" ]; then
      ion_cl=`grep Cl-  complex_vdw_bonded.pdb  | head -1  |  awk '{print $5}'`
      let line=$ion_cl-1
elif  [ -n "$isna" ]; then
      ion_na=`grep Na+  complex_vdw_bonded.pdb  | head -1  |  awk '{print $5}'`
      let line=$ion_na-1
else
      ion_wat=`grep WAT  complex_vdw_bonded.pdb  | head -1  |  awk '{print $5}'`
      let line=$ion_wat-1
fi


sed  -i  "s/1-num_resid/1-${line}/g"    min.in
sed  -i  "s/1-num_resid/1-${line}/g"    heat.in
sed  -i  "s/1-num_resid/1-${line}/g"    density.in

cd  /data/wunian2/new-bags
done
