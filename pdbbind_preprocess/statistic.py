#!/usr/bin/env python
# coding: utf-8

### analyze all pdbid in pdbbind database 2019

import pandas as pd
import os
import numpy as np
from pandas import Series, DataFrame
import matplotlib.pyplot as plt
import matplotlib

import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score


from sklearn.feature_extraction.text import CountVectorizer
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score,confusion_matrix,classification_report
pd.set_option('display.max_columns', None)

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem
from rdkit.Chem import Draw 
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdFMCS

from rdkit.Chem import MCS
from rdkit.Chem import rdFMCS

from rdkit.Chem import rdRGroupDecomposition as rdRGD


import math 
import re
import time

from pymol import cmd
import pymol2

from tqdm import tqdm



amino={'G':'Gly','A':'Ala','V':'Val','L':'Leu','I':'Ile','P':'Pro','F':'Phe','Y':'Tyr','W':'Trp','S':'Ser','T':'Thr',
       'C':'Cys','M':'Met','N':'Asn','Q':'Gln','D':'Asp','E':'Glu','K':'Lys','R':'Arg','H':'His'}
amino_v=list(amino.values())
amino_std=[amino_v[i].upper()  for   i in range(len(amino_v))]+list(amino.keys())
amino_std=amino_std+['HOH']


p=pymol2.PyMOL()
p.start()
cmd=p.cmd


class pdbid_properties:
    def __init__(self, pdbid):
        self.path='/data/wunian/master-github/md-ml-nian/v2019-other-PL-refine'
        self.pdbid=pdbid
        self.lig=self.pdbid+'_ligand'
        self.pro=self.pdbid+'_protein'
        self.lig_file=self.pdbid+'_ligand.mol2'
        self.pro_file=self.pdbid+'_protein.pdb'
        self.lig_file_path=os.path.join(self.path, self.pdbid, self.lig_file)
        self.pro_file_path=os.path.join(self.path, self.pdbid, self.pro_file)

    def get_smiles(self):
        mol=Chem.MolFromMol2File(self.lig_file_path)
        smi=Chem.MolToSmiles(mol,isomericSmiles=True,canonical=True)
        return smi

    def get_lig_name(self):
        cmd.load(self.lig_file_path)
        lig_name=cmd.get_model(self.lig).atom[0].resn
        cmd.remove("all")
        return lig_name

    def get_ions(self):
        cmd.load(self.pro_file_path)      
        cmd.select('sel_metal',self.pro+" and metal")     
        metals=cmd.get_model('sel_metal').atom
        metals=[i.name for i in metals]
        metals_num=len(metals)
        cmd.remove("all")
        return metals, metals_num
    
    def get_non_aa_chains_seq(self):
        cmd.load(self.pro_file_path)
        gg=cmd.get_model(self.pro)
        chs=cmd.get_chains(self.pro)
        pattern = re.compile('[A-Za-z]')
        if pattern.findall(str(chs)):
             try:
                chs.remove('')
             except:
                pass
             seqs=cmd.get_fastastr(self.pro+' and chain %s' % chs[0])
             seqs=''.join([s0 for s0 in seqs.strip().split()[1:]])
        else:
             seqs=cmd.get_fastastr(self.pro)
             seqs=''.join([s0 for s0 in seqs.strip().split()[1:]])
        aa_abc=[]
        for atom in gg.atom:
            if atom.resn not in amino_std: 
                #print(atom.resi, atom.resn)
                aa_cdd=atom.resi+'_'+atom.resn
                aa_abc.append(aa_cdd)
        cmd.delete('all') 
        return set(aa_abc), chs, seqs, len(seqs)

def write_information():   
    content=open('2019-results-final-chs-seq-2-aa.txt','a')
    content.write('pdbid \t lig_name \t smiles \t chains \t chains_num \t ion num \t ion  \t aa num  \t aa  \n')
    content.close()
    os.chdir('../v2019-other-PL-refine')
    pdblist=os.listdir()
    for i in tqdm(pdblist):
        with open('/data/wunian/master-github/md-ml-nian/PL-refine-analysis/2019-results-final-chs-seq-2-aa.txt','a') as  content:
            try:
                pp=pdbid_properties(i)
                # pp_lig_name=pp.get_lig_name()
                # pp_smis=pp.get_smiles()     
                # pp_chains,pp_chains_num=pp.get_chains()
                # ions,ions_num=pp.get_ions()
                aa,bb,cc,dd=pp.get_non_aa_chains_seq()
#                aa.remove('')
                cmd.delete('all')
                # content.write('%s \t %s \t %s \t %s \t  %s \t %s \t %s \t %s \t  %s \t %s \n' % (i,str(pp_lig_name), str(pp_smis), str(pp_chains), str(pp_chains_num), str(ions_num),str(ions),str(aa_num),str(aa), str(bb)))
                content.write('%s \t  %s  \t %s  \t %s \t %s  \n' % (i, str(aa), str(bb), str(cc), dd))
            except:
                content.write('%s \n' % i )
        # os.chdir('../')
    #    content.write('pdbid: %s \t ion num: %s \t ion: %s \t aa num: %s \t aa: %s \n' % (i,str(ions_num),str(ions),str(aa_num),str(aa)))
    os.chdir('/data/wunian/master-github/md-ml-nian/PL-refine-analysis')
    

write_information()




