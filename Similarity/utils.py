#!/bin/bash/python 
# Author: Ragy Haddad
import sys,os,json 
import pandas as pd 
import numpy as np 
from rdkit.Chem import MACCSkeys 
from rdkit import DataStructs 
from rdkit import Chem
sys.path.append('../') 
from Descriptor.ligand_descriptor import getMorganFP  
from random import randrange
# Load Ligand Info 
with open('../HashMaps/pdb_ligands_map.json',"r") as f:
    ligand_map = json.load(f)  

def check_mols(smiles):
    try:
        getMorganFP(smiles)
        return True 
    except:
        return False
n_mols = len(ligand_map)
all_ligands = list(ligand_map.keys())
# Takes all True Positives for a Target
def find_negative_dti(true_positives,negative_count=4):
    proportion_similarity = 0.0 
    true_positive_fps = [] 
    smiles = [ligand_map[mol_smile]['smiles'] for mol_smile in true_positives] 
    true_fps = [getMorganFP(smi,string=False,bit=2048,radius=2) for smi in smiles] 

    # Scan Through Availble Compounds in PDB and
    # Select a Molecule that is dissimilar to all the true positives
    # Associated with a given target
    minimum_similarity = 9999 
    hq_negatives = [] 
    hq_positives = []
    n_tries = 0
    while True:
        average_similarity = 0.0 # The average similarity detected by comparing one potential true negative with all positives
        random_idx = randrange(n_mols) 
        n_tries += 1
        ligand = all_ligands[random_idx] 
        if ligand in true_positives:
            continue 
        try:
            neg_smiles = ligand_map[ligand]['smiles']
        except:
            continue
        # A potential Negative Fingerprint  
        if check_mols(neg_smiles) == False:
            continue
        neg_fp = getMorganFP(neg_smiles,string=False,bit=2048,radius=2)
        
        # Compare the potential negative to all the true positives
        for fp in true_fps:
            instance_similarity = DataStructs.FingerprintSimilarity(neg_fp,fp)
            average_similarity += instance_similarity  
        
        # Calculate the Average Similarity Across True Positives
        average_similarity = average_similarity/len(true_positives)
        # Check Rules:
        if average_similarity <= 0.1:
            hq_negatives.append(ligand_map[ligand]) 
        if average_similarity >= 0.90 and ligand not in true_positives: 
            hq_positives.append(ligand_map[ligand])
        # if average_similarity >= 0.9:
        if len(hq_negatives) == negative_count or n_tries == 30000:
            break 
    return hq_negatives,hq_positives



        


    