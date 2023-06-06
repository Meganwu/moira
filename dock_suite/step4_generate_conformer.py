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
from rdkit.Chem import rdMolAlign,rdForceFieldHelpers,rdDistGeom,rdMolDescriptors

import shutil



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


def py3D_show(mol1,mol2,confId_1=-1,confId_2=-1):
    p_O3A = py3Dmol.view()
    p_O3A.addModel(Chem.MolToMolBlock(mol1, confId=confId_1), 'sdf')
    p_O3A.addModel(Chem.MolToMolBlock(mol2, confId=confId_2), 'sdf')
    p_O3A.setStyle({'stick':{}})
    p_O3A.zoomTo()
    p_O3A.show()

def get_confs(mol2,confs=10,save_prefix='id'):
    #mol=Chem.MolFromPDBFile(pdb,removeHs=False)
    mol=Chem.MolFromMol2File(mol2)
    smi=Chem.MolToSmiles(mol,allBondsExplicit=True, allHsExplicit=True)
    #inchi=Chem.MolToInchi(mol)
    #mol1=Chem.MolFromInchi(inchi,sanitize=False)
    mol_t=Chem.MolFromSmiles(smi,sanitize=False)
    m3d=Chem.AddHs(mol_t)
    AllChem.EmbedMultipleConfs(m3d, numConfs=confs,pruneRmsThresh=1)
    #AllChem.MMFFOptimizeMolecule(m3d,mmffVariant='MMFF94')
    rmsd=[]
    for k in range(confs):
        pyO3A=rdMolAlign.GetO3A(m3d,mol,prbCid=k,refCid=-1)
        rmsd_0=pyO3A.Align()
        rmsd.append(rmsd_0)
        filename=save_prefix+str(k)+'.pdb'
        Chem.MolToPDBFile(m3d,filename, confId=k)


def get_confs2(mol2,confs=10,save_prefix='id'):
    mol=Chem.MolFromMol2File(mol2,removeHs=False)
    smi=Chem.MolToSmiles(mol,allBondsExplicit=True, allHsExplicit=True)
    #inchi=Chem.MolToInchi(mol)
    #mol1=Chem.MolFromInchi(inchi,sanitize=False)
    mol_t=Chem.MolFromSmiles(smi,sanitize=False)
    m3d=Chem.AddHs(mol_t)
    AllChem.EmbedMultipleConfs(m3d, numConfs=confs,maxAttempts=500,pruneRmsThresh=0.8,randomSeed=1,numThreads=8)
    crip1=rdMolDescriptors._CalcCrippenContribs(mol)
    crip2=rdMolDescriptors._CalcCrippenContribs(m3d)
    #AllChem.MMFFOptimizeMolecule(m3d,mmffVariant='MMFF94')
    tmp=open('rmsd.txt','w')
    rmsd=[]
    for k in range(confs):
        # pyO3A=rdMolAlign.GetO3A(m3d,mol,mmff_param,mmff_param,prbCid=k,refCid=-1)
        crippenO3A = rdMolAlign.GetCrippenO3A(m3d, mol, crip2, crip1, k, -1)
        rmsd_0=crippenO3A.Align()
        rmsd.append(rmsd_0)
        tmp.write('cid:  %s rmsd: %s \n' % (str(k),rmsd_0))
        filename=save_prefix+str(k)+'.pdb'
        Chem.MolToPDBFile(m3d,filename, confId=k)
    tmp.close() 



if __name__ == '__main__':
    outfile_good=open('confs_good.txt','w')
    outfile_bad=open('confs_bad.txt','w')
    pdb_list=read_pdbid()
    os.chdir('step2_ligand_mol2')
    for j,i in enumerate(pdb_list,start=1):
        dock_name=i+'dock'
        lig_name=i+'_ligand.mol2'
        os.mkdir(dock_name)
        shutil.move(lig_name, dock_name)
        os.chdir(dock_name)
        try:
           get_confs2(lig_name,save_prefix=i)
           outfile_good.write('%s  %s  %s \n' % (str(j),str(i),str(1)))
        except:
           outfile_bad.write('%s  %s  %s \n' % (str(j),str(i),str(0)))
        os.chdir('../')
    outfile_good.close()
    outfile_bad.close()
