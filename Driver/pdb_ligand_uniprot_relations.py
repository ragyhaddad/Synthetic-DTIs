#!/bin/bash/python 
# Author: Ragy Haddad 
import sys,os,json 
import pandas as pd 
import numpy as np 
sys.path.append('../')
from Utils.pdb_utils import pdb_ligands 
from Utils.sequence_utils import uniprotToSequenceLocal 

# Usage: python3 pdb_ligand_uniprot_relations.py ../Data/formatted_data/pdb_uniprot_homo_reviewed_formatted.tsv 

columns = ['pdb_id','uniprot_id','uniprot_name','status','protein_names','gene_names','organism','length'] 
input_df = pd.read_csv(sys.argv[1],sep='\t') 
output_df = []

# Load Ligand Info 
with open('../HashMaps/pdb_ligands_map.json',"r") as f:
    ligand_map = json.load(f)  
with open('../HashMaps/pdb_ligands_relations.json',"r") as f:
    pdb_ligand_map = json.load(f)

skipped_ligands = []
new_columns = ['smiles','target_seq'] + columns + ['ligand_molecular_weight','ligand_inchi_key','ligand_name','ligand_id']
for idx,v in enumerate(input_df.values):
    print('-- Progress: [%i/%i]' % (idx,len(input_df.values)))
    uniprot_id = v[1]
    pdb_id = v[0]
    target_seq = str(uniprotToSequenceLocal(v[1]))
    ligands = pdb_ligand_map[pdb_id] 
    for ligand in ligands:
        new_v = []
        try:
            ligand_info = ligand_map[ligand]
        except:
            skipped_ligands.append(ligand)
            continue 
        try:
            smiles = ligand_info['smiles']
        except:
            skipped_ligands.append(ligand) 
            continue 
        if smiles == '*':
            continue 
        try:
            inchi_key = ligand_info['inchi'] 
        except:
            skipped_ligands.append(ligand) 
            continue 
        new_v.append(smiles)
        new_v.append(target_seq)
        new_v += list(v) 
        new_v.append(ligand_info['molecular_weight'])
        new_v.append(inchi_key)
        new_v.append(ligand_info['chemical_name']) 
        new_v.append(ligand_info['ligand_id'])
        output_df.append(new_v) 
output_df = pd.DataFrame(output_df) 
output_df.columns = new_columns 
output_df.to_csv('pdb_ligands_uniprot_relations_home_reviewed.tsv',index=False,sep='\t')

stats_file = open('skipped_ligands.txt',"w") 
for x in skipped_ligands:
    stats_file.write('%s\n' % x)