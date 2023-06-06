import os
import numpy as np
from pymol import *
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem
from rdkit.Chem import Draw 
from rdkit import DataStructs
import pandas as pd
from Bio import PDB

from collections import Counter

def keep_chain(pdbid,dist=4):
#    pro=pdbid+'_protein.pdb'
#    pocket=pdbid+'_pocket.pdb'
    lig=pdbid+'_ligand.mol2'
#    cmd.load(pro)
#    cmd.load(pocket)
    cmd.load(lig)
    #aa=cmd.select('cha1','4tmn_ligand around 3')
 #   aa=cmd.select('cha1',pdbid+'_ligand '+'around '+str(dist))
#    bb=cmd.get_model('cha1')
#    att=[]
#    for i in bb.atom:
#        att.append(i.chain)
#    att=[i for i in att if i!='']
#    cha_sel=Counter(att).most_common(1)[0][0]
#    cmd.save(pdbid+'_chain.pdb',pdbid+'_protein '+'and chain '+cha_sel)
    cmd.save(pdbid+'_ligand.pdb',pdbid+'_ligand')
    cmd.remove('all')
#    return cha_sel

def read_pdbid(filename='refine_name.txt',col_name=0):
    list=pd.read_table(filename,header=None)[col_name].to_list()
    return list


if __name__ == '__main__':
    outfile=open('result3.txt','w')
    pdb_list=read_pdbid()
#    list=pd.read_table('result1.txt',sep='\s+', header=None)
#    pdb_list=list[list[2]==0][1].to_list()
    os.chdir('refined-set')
    for j,i in enumerate(pdb_list,start=1):
        os.chdir(i)
        try:
           keep_chain(i)
           os.chdir('../')
           outfile.write('%s  %s  %s\n' % (str(j),str(i),str(1)))
        except:
           os.chdir('../') 
           outfile.write('%s  %s  %s\n' % (str(j),str(i),str(0)))

    outfile.close()
