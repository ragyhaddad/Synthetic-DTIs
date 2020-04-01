#!/bin/bash/python 
# Author: Ragy Haddad
import sys,os 
from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors, viewVoxelFeatures
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
from moleculekit.smallmol.smallmol import SmallMol
from moleculekit.home import home
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np 
import Bio  
import subprocess
from moleculekit.smallmol.smallmol import SmallMol 
import time 
## MOLECULE KIT VERSION NEEDS TO BE VERSION 0.1.16 <----- HAS FIXED BUGS

# Create Temp 
if os.path.isdir('tmp') == False:
    os.mkdir('tmp') 

non_prot_mols = []

# Remove HetAtoms for Protonation Step - Required for voxelization
def removeHETAtoms(pdb_file,pdb_id):
    class NonHetSelect(Select):
        def accept_residue(self, residue):
            return 1 if residue.id[0] == " " else 0
    class RemoveErroneous(Select):
        def accept_residue(self, residue):
            if residue.resname == 'MG5' or residue.resname == 'HOH' or residue.resname == 'ZN3' or residue.resname == 'ZN4':
                return 0 
            else:
                return 1
    pdb = PDBParser().get_structure(pdb_file.split('.')[0], pdb_file)
    io = PDBIO()
    io.set_structure(pdb)
    out_file = os.path.join('./tmp',pdb_id +'_nonhet.pdb')
    io.save(out_file, NonHetSelect()) 
    io.save(out_file, RemoveErroneous()) 
    return out_file 

# Generate Descriptor using the Moleculekit - Cited from Moleculekit tutorial 
def generateProteinDescriptor(pdb_file,voxel_size=64,mode='standard'): 
    # Remove HetAtoms
    pdb_id = os.path.basename(pdb_file).split('.')[0]
    # If your structure is fully protonated and contains all bond information in prot.bonds skip this step!
    if mode != 'scPDB':
        pdb_file_nonhet = removeHETAtoms(pdb_file,pdb_id)
        prot = Molecule(pdb_file_nonhet)
        prot = prepareProteinForAtomtyping(prot)
        prot.center()
        prot.write(pdb_file_nonhet)
        prot = Molecule(pdb_file_nonhet)
        prot.view(guessBonds=False)
    else:
        prot = SmallMol(pdb_file)     
    try:
        prot_vox, prot_centers, prot_N = getVoxelDescriptors(prot,voxelsize=1, buffer=1,center=(0,0,0),boxsize=(voxel_size,voxel_size,voxel_size))
    except:
        return False
    nchannels = prot_vox.shape[1]
    # Reshape Voxels
    prot_vox_t = prot_vox.transpose().reshape([1, nchannels, prot_N[0], prot_N[1], prot_N[2]]) 
    return prot_vox_t


# scPDB binding site file
def fuzCavDescriptor(bs_mol):
    os.system('rm Count.txt')
    fuzcav_dir = '../Packages/FuzCav'
    cmd_1 = 'perl %s/utils/CaTagger.pl %s > site1_Tagged.mol2' % (fuzcav_dir,bs_mol)
    l_file = open('listCavTagged',"w")
    l_file.write('site1_Tagged.mol2\n')
    os.system(cmd_1)
    cmd_2 = "java -jar %s/dist/3pointPharCav.jar -d %s/utils/resDef/tableDefCA.txt -t %s/utils/triplCav/interval.txt -l listCavTagged -o Count.txt -c" % (fuzcav_dir,fuzcav_dir,fuzcav_dir)
    p = subprocess.Popen([cmd_2],shell=True,stdout=subprocess.PIPE)



fuzCavDescriptor(sys.argv[1])