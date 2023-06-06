#!/bin/bash

pdblist=`cat  md_222.txt` 

for id in $pdblist
do
     j=`echo ${id:0:4}`
     md_name2=${j}MD2_neu
     cd   native_neutral/$md_name2
     cp  /data/wunian2/neutral_collection/inputs/*     .
     ./pmemd.sh
     cd   /data/wunian2/neutral_collection
done
