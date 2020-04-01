#!/bin/bash 
# Author: Ragy Haddad
import sys,os 
from rdkit import Chem  
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys 
import numpy as np 

# Encode Molecule to RDKit ECFP - Default Radius 2 
def getMorganFP(smiles,bit=1024,radius=2,string=True):  
    m = Chem.MolFromSmiles(smiles)
    if string == False:
        morgan = AllChem.GetMorganFingerprintAsBitVect(m,radius,nBits=bit)
        return morgan
    morgan = AllChem.GetMorganFingerprintAsBitVect(m,radius,nBits=bit).ToBitString()
    bit_array = np.array([int(c) for c in morgan],dtype=np.float32)
    return bit_array 