#!/bin/bash/python 
import sys,os,json 
import pandas as pd 
import numpy as np 
sys.path.append('../')
from Similarity.utils import find_negative_dti
# Usage: python3 create_synthetic_dataset.py ../Data/formatted_data/pdb_ligands_uniprot_relations_home_reviewed.tsv

# This Dataframe Contains All Human Positive Interactions from PDB 
input_df = pd.read_csv(sys.argv[1],sep='\t') 

# Filter Interactions by MW So we Can remove ION to Target Interaction
input_df = input_df[input_df['ligand_molecular_weight'] > 179.0]  

# 55,000 Human Interaction 
input_df = input_df.sort_values(by=['uniprot_id']) 

all_targets = set(input_df['uniprot_id'].values.tolist())
all_targets = list(all_targets)

columns = ['smiles','target_seq','ligand_id','target_id','target_description','interaction']  

# Keep Track of interactions that are already in dataset
in_dataset = {}

output_df = []

count = 0 
for index,uniprot_id in enumerate(all_targets): 
    print('-- Progress [%i/%i]' % (index,len(all_targets)))

    # Get All Entries For A Given Target
    values = input_df[input_df['uniprot_id'] == uniprot_id]
    # print(values)
    # exit()
    target_seq = values['target_seq'].values[0] 
    # Discard Protein Longer Than 2500
    if len(target_seq) > 2500:
        continue  
    target_desc = values['protein_names'].values[0]
    true_positive_ligands = list(set(values['ligand_id'].values.tolist())) 

    # Set True Negatives to be same number as true pos
    number_of_synthetic_negatives = len(true_positive_ligands) # Imbalance test - have more negatives than positives 

    synthetic_negatives,synthetic_positives = find_negative_dti(true_positive_ligands,negative_count=number_of_synthetic_negatives) 
    
    # True Positives
    for row in values.values:
        smiles = row[0]
        ligand_id = row[-1]
        value = [smiles,target_seq,ligand_id,uniprot_id,target_desc,1.0]
        interaction_key = uniprot_id + '-' + ligand_id 
        # Make sure there is no duplicate interactions
        if interaction_key not in in_dataset:
            output_df.append(value) 
            in_dataset[interaction_key] = 1
    # Synthetic Negatives 
    for l in synthetic_negatives:
        smiles = l['smiles']
        ligand_id = l['ligand_id']
        value = [smiles,target_seq,ligand_id,uniprot_id,target_desc,0.0] 
        interaction_key = uniprot_id + '-' + ligand_id 
         # Make sure there is no duplicate interactions
        if interaction_key not in in_dataset:
            output_df.append(value) 
            in_dataset[interaction_key] = 1
        
    # # Synthetic Positives  
    # for l in synthetic_positives:
    #     smiles = l['smiles']
    #     ligand_id = l['ligand_id']
    #     value = [smiles,target_seq,uniprot_id,ligand_id,target_desc,1.0]
    #     output_df.append(value)  


# Save output
df = pd.DataFrame(output_df) 
df.columns = columns 
df.to_csv('synthetic_dataset_v1.tsv',index=False,sep='\t')
exit()



    