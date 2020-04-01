#!/bin/bash/python 
# Author: Ragy Haddad
import sys,os 
import pandas as pd 
import xml.etree.ElementTree as ET
import requests 

# Format PDB ID To Uniprot
def pdb_mapping(pdb_id): 
    url = "https://www.rcsb.org/pdb/rest/describeMol?structureId=" + str(pdb_id).lower() + "&outputFormat=JSON" 
    r = requests.get(url) 
    response_xml_as_string = r.text 
    responseXml = ET.fromstring(response_xml_as_string) 
    pdb_id,taxonomy,uniprot_id,description,macromolecule_name = "","","","",""
    for child in responseXml:
        pdb_id = child.get('id') 
        for c in child:
            for j in c:
                if j.tag == 'Taxonomy':
                    taxonomy = j.get('name')
                if j.tag == 'macroMolecule':
                    macromolecule_name = j.get('name')  
                if j.tag == 'polymerDescription':
                    description = j.get("description")
                for k in j:
                    if k.tag == 'accession':
                        uniprot_id = k.get('id') 
    out = {"pdb_id":pdb_id,"taxonomy":taxonomy,"uniprot_id":uniprot_id,"description":description,"macromolecule_name":macromolecule_name} 
    return out  

# Get Ligands Associated with a PDB ID and all their information    
def pdb_ligands(pdb_id):
    url = "https://www.rcsb.org/pdb/rest/ligandInfo?structureId=" + str(pdb_id).lower() 
    r = requests.get(url) 
    response_xml_as_string = r.text    
    responseXml = ET.fromstring(response_xml_as_string)  
    ligands = []
    for i in responseXml:
        ligand = {}
        ligand_id,chemical_name,smiles,inchl_key,inchi,formula,molecular_weight = "","","","","","",0.0
        for j in i:
            if j.tag == 'ligand':
                molecular_weight = float(j.get('molecularWeight'))
                ligand_id = j.get("chemicalID").strip()
            for k in j:
                if k.tag == 'chemicalName':
                    chemical_name = k.text.strip() 
                if k.tag == 'formula':
                    formula = k.text.strip() 
                if k.tag == 'InChI':
                    inchi = k.text.strip() 
                if k.tag == 'InChIKey':
                    inchl_key = k.text.strip() 
                if k.tag == 'smiles':
                    smiles = k.text.strip() 
            ligand = {"ligand_id":str(ligand_id),"chemical_name":str(chemical_name),"smiles":str(smiles),"inchl_key":str(inchl_key),"inchi":str(inchi),"formula":str(formula),"molecular_weight":float(molecular_weight)}
            ligands.append(ligand) 
    return ligands


def pdbLigandsLocal(pdb_id):
    pass