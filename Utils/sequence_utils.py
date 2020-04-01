#!/bin/bash/python
# Author: Ragy Haddad 
import sys,os
import pandas as pd 
import numpy as np 
import requests 
from pyfaidx import Fasta

print(os.getcwd())
ROOT_DIR = os.path.dirname(os.path.abspath(__file__)) # This is your Project Root 

load_dud_csv = False
dud_gene_map = {}
genes = None 
local_seqs = {}
# Kegg to sequence
def keggToSequence(kegg_target_id):
    formatted_id = ''
    if ':' in kegg_target_id:
        formatted_id = kegg_target_id 
    else:
        formatted_id = kegg_target_id[0:3] + ':' + kegg_target_id[3:]
    url = 'http://rest.kegg.jp/get/' + formatted_id +'/aaseq' 
    r = requests.get(url)
    print(r) 
    r = r.text.split('>')[1].split('\n')[1:]
    r = "".join(r)
    seq = r.strip()
    return seq 

def uniprotToSequence(uniprot_id):
    url = "https://www.uniprot.org/uniprot/" + uniprot_id + '.fasta' 
    r = requests.get(url) 
    r = r.text.split('>')[1].split('\n')[1:]
    r = "".join(r)
    seq = r.strip()
    return seq
def uniprotToSequenceLocal(uniprot_id):
    global genes 
    if genes == None:
        print('-- Indexing Uniprot database')
        genes = Fasta(ROOT_DIR + '/uniprot_acc.fasta') 
    try:
        seq = genes[uniprot_id]
    except:
        seq = uniprotToSequence(uniprot_id)
    return seq
# Kegg IF to Smiles
def keggToSmiles(kegg_drug_id):
    url = 'http://rest.kegg.jp/get/' + kegg_drug_id 
    r = requests.get(url)
    r = r.text.split('\n')
    smiles = ''
    d_name = ''
    for i in r:
        if 'NAME' in i:
            smiles = i.replace('NAME','')
            smiles = smiles.strip().replace(';','').split(' ')[0]
            d_name = smiles 
            try:
                smiles = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+ smiles + '/property/CanonicalSMILES/json')
                smiles = smiles.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            except:
                return '','' 
            break 
    return smiles,d_name

# Gene Name to Sequence
def geneToSequence(gene_name):
    pass 


# PDB To Sequence
def pdbToSequence(pdb_id):
    url = 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=' + str(pdb_id).upper() + '&compressionType=uncompressed' 
    r = requests.get(url) 
    r = r.text.split('>')[1].split('\n')[1:]
    r = "".join(r)
    seq = r.strip()
    return seq

    

