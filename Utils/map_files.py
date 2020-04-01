#!/bin/bash/python 
# Author: Ragy Haddad 
import sys,os,json 
import pandas as pd  
import numpy  
from rdkit.Chem.Descriptors import ExactMolWt  
from rdkit import Chem
# Usage: python3 map_files.py ../Data/supporting_tables/Components-smiles-oe.smi 
# Parse All PDB Ligand Information to a JSON File for data retrieval 


final_json = {}

inchi_map = sys.argv[2] 
smiles_map = sys.argv[1] 


malformed_count = 0
# Parse Smiles 
with open(smiles_map,"r") as f:
    for line in f:
        line = line.strip() 
        cols = line.split('\t') 
        if len(cols) == 1:
            continue
        smiles = cols[0] 
        ligand_id = cols[1] 
        try:
            ligand_name = cols[2] 
        except:
            ligand_name = "NA"
        if ligand_id not in final_json:
            final_json[ligand_id] = {}
            final_json[ligand_id]["smiles"] = smiles
            final_json[ligand_id]["chemical_name"] = ligand_name
            final_json[ligand_id]["ligand_id"] = ligand_id
            try: 
                mol  = Chem.MolFromSmiles(smiles) 
                mw = ExactMolWt(mol)  
                final_json[ligand_id]['molecular_weight'] = mw 
                final_json[ligand_id]['explicit_valence'] = True
                
            except:
                malformed_count += 1
                final_json[ligand_id]['molecular_weight'] = 0.0
                final_json[ligand_id]['explicit_valence'] = False
 
with open(inchi_map,"r") as f:
    for line in f:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) == 1:
            continue
        inchi = cols[0] 
        ligand_id = cols[1]
        if ligand_id not in final_json:
            final_json[ligand_id] = {}
            final_json[ligand_id]['inchi'] = inchi 
        else:
            final_json[ligand_id]['inchi'] = inchi

# Save file 
print('-- Saving JSON pdb_ligands_map.json')
with open('pdb_ligands_map.json', 'w') as outfile:
    json.dump(final_json, outfile)