#!/usr/bin/env python
import sys
import numpy as np
import rdkit
import rdkit.Chem as Chem
from rdkit.Chem import Draw
import rdkit.Chem.rdchem
import rdkit.Chem.rdmolfiles
from rdkit.Chem.rdchem import BondType,HybridizationType
from typification import *
import xml.etree.ElementTree as ET

def read_ligpargen_lib(path):
    lib_xml_naive = {}
    lib_pdb_naive = {}
    with open(path+'/file_index.csv', 'r') as ifile:
        for line in ifile:
            words = line.split(',')
            smiles = words[0]
            name = words[1].strip()
            if words[1] != 'None':
                lib_xml_naive[smiles] = '/'.join([path, name, name+'.xml'])
                lib_pdb_naive[smiles] = '/'.join([path, name, name+'.pdb'])
    return lib_xml_naive, lib_pdb_naive

def sanity_check(mol1, mol2, o1, o2):
    elem1 = np.array([a.GetSymbol() for a in mol1.GetAtoms()])
    connect1 = np.array(Chem.GetAdjacencyMatrix(mol1))
    elem2 = np.array([a.GetSymbol() for a in mol1.GetAtoms()])
    connect2 = np.array(Chem.GetAdjacencyMatrix(mol2))
    if np.all(elem1[o1] == elem2[o2]) and np.all(connect1[o1,:][:,o1] == connect2[o2,:][:,o2]):
        return True
    else:
        return False

if __name__ == '__main__':
    path = 'ligpargen_files'
    lib_xml_naive, lib_pdb_naive = read_ligpargen_lib('ligpargen_files')

    # template 
    smiles = 'CC1COC(=O)O1'
    # mol_template = Chem.MolFromSmiles(smiles)
    mol1_template, atoms, atypes = typify_mol(smiles, depth=1, withH=True)

    # structure from pdb
    mol1 = Chem.MolFromPDBFile(lib_pdb_naive[smiles], removeHs=False)
    if not sanity_check(mol1_template, mol1, range(mol1_template.GetNumAtoms()), range(mol1.GetNumAtoms())):
        print('ERROR: connectivity mismatch')

    # sanity check passed for heavy atoms
    xmltree = ET.parse(lib_xml_naive[smiles])
