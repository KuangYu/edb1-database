#!/usr/bin/env python
import sys
import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem import Draw, AddHs
import re

smiles = "[B-]12(OC(=O)C(=O)O1)OC(=O)C(=O)O2" #sys.argv[1]
oname = "test"

mol = Chem.MolFromSmiles(smiles)
# add hydrogen atoms
mol = Chem.AddHs(mol)

# add conformer
params = AllChem.ETKDGv3()
AllChem.EmbedMolecule(mol, params)

Chem.MolToPDBFile(mol, '%s.pdb'%oname)
