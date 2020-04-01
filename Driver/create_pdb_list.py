#!/bin/bash/python
#Author: Ragy Haddad
import sys,os 
import pandas as pd 
import numpy as np 
sys.path.append('../') 
from Utils.pdb_utils import pdb_mapping,pdb_ligands

# Usage: python3 create_pdb_list.py ../Data/pdb_ids.txt

# Create a list of Human PDB Chains their associated Uniprot IDs and their information 
pdb_ids = open(sys.argv[1],"r").readlines() 
pdb_ids = [i.strip() for i in pdb_ids]  

# Make The List Unique 
pdb_ids = list(set(pdb_ids))
n_ids = len(pdb_ids) 
missing_uniprot = 0 # Number of PDB IDs without a uniprot accession 
human_pdb_chains = 0 

# Select Only Human PDB Complexes
human_only = True

# Outfile 
outfile = open('PDB_UNIPROT_HOMO.tsv',"w")
outfile.write('pdb_id\ttaxonomy\tuniprot_id\tdescription\tmacromolecule_name\n')
def write_to_file(output_json):
    pdb_id = output_json['pdb_id'] 
    taxonomy = output_json['taxonomy'] 
    uniprot_id = output_json['uniprot_id']
    description = output_json['description'] 
    macromolecule_name = output_json['macromolecule_name'] 
    outfile.write('%s\t%s\t%s\t%s\t%s\n' % (pdb_id,taxonomy,uniprot_id,description,macromolecule_name))

# Get Uniprot Mapping + Info for each PDB ID 
for idx,pdb_id in enumerate(pdb_ids):
    print('-- Progress [%i/%i]' % (idx,n_ids)) 
    # If an ID is malformed skip
    if len(pdb_id) != 4:
        continue 
    # Get PDB to Uniprot Mapping
    pdb_info = pdb_mapping(pdb_id)  
    # If there is no mapping count and skip
    if len(pdb_info['uniprot_id']) == 0:
        missing_uniprot += 1  
        continue 
    # Select only records with human taxonomy
    if human_only and pdb_info['taxonomy'].lower() == 'homo sapiens':
        write_to_file(pdb_info) 

outfile.close()


