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

def get_center(pymol_model):
    model=cmd.get_model(pymol_model)
    a=[]
    for i in range(len(model.atom)):
        coord=model.atom[i].coord
        a.append(coord)
    a=np.array(a)
    center=np.mean(a, axis=0)
    return center

def write_center(pdbid):
    lig=pdbid+'_ligand.pdb'
    cmd.load(lig)
    pymol_model=pdbid+'_ligand'
    center=get_center(pymol_model)
    cmd.remove('all')
    return center


def read_pdbid(filename='refine_name.txt',col_name=0):
    list=pd.read_table(filename,header=None)[col_name].to_list()
    return list


if __name__ == '__main__':
    outfile_good=open('center_good.txt','w')
    outfile_bad=open('center_bad.txt','w')
    pdb_list=read_pdbid()
    os.chdir('step1_prepared_pdb')
    for j,i in enumerate(pdb_list,start=1):
        try:
           center=write_center(i)
           outfile_good.write('%s  %s  %s \n' % (str(j),str(i),str(center)))
        except:
           outfile_bad.write('%s  %s  %s \n' % (str(j),str(i),str(0)))

    outfile_good.close()
    outfile_bad.close()
