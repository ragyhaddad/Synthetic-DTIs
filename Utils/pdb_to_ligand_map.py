#!/bin/bash/python 
import sys,os,json  
import pandas as pd 
import numpy as np 

# Usage: python3 pdb_to_ligand_map.py ../Data/supporting_tables/cc-to-pdb.tdd 

infile = sys.argv[1] 
pdb_ligand_map = {}
with open(infile,"r") as f:
    for line in f:
        line = line.strip() 
        cols = line.split('\t') 
        ligand_id = cols[0] 
        pdb_ids = cols[1].split(' ') 
        for pdb in pdb_ids:
            if pdb not in pdb_ligand_map:
                pdb_ligand_map[pdb] = [] 
                pdb_ligand_map[pdb].append(ligand_id) 
            else:
                pdb_ligand_map[pdb].append(ligand_id)

# Save file 
print('-- Saving JSON pdb_ligands_relations.json')
with open('pdb_ligands_relations.json', 'w') as outfile:
    json.dump(pdb_ligand_map, outfile)


