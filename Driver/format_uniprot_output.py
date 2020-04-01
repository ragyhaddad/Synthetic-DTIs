#!/bin/bash/python 
# Author: Ragy Haddad 
import sys,os 
import pandas as pd 

df = pd.read_csv(sys.argv[1],sep='\t') 

# columns = ['pdb_id','Entry','Entry_name','Status','Protein_names','Gene_names','Organism','Length']
columns = ['pdb_id','uniprot_id','uniprot_name','status','protein_names','gene_names','organism','length'] 
df.columns = columns 


final_df = [] 

for idx,v in enumerate(df.values):
    pdb_ids = v[0] 
    if ',' in pdb_ids:
        split_ids = pdb_ids.split(',')
        for id in split_ids:
            new_v = [] 
            new_v.append(id) 
            new_v += list(v[1:])
            final_df.append(new_v) 
    else:
        final_df.append(v) 

df_out = pd.DataFrame(final_df) 
df_out.columns = columns 

# Save Df 
df_out.to_csv('pdb_uniprot_homo_reviewed_formatted.tsv',index=False,sep='\t')