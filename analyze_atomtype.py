#!/usr/bin/env python
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw

ifile = open('res1.xvg', 'r')
ifile.readline()
ifile.readline()

class Molecule:

    def __init__(self, smiles, freq=0, pricetag='N', is_salt=False):
        self.smiles = smiles
        if '+' in smiles:
            is_salt = True
        self.is_salt = is_salt
        self.freq = freq
        self.pricetag = pricetag
        self.rdkitmol = Chem.MolFromSmiles('CCOCC')
        self.count = 0

    def __repr__(self):
        return 'Molecule with smiles: %s'%self.smiles

    def draw_2d_struct(self, ofn):
        Draw.MolToFile(self.rdkitmol, ofn)
        return

class Formula:

    def __init__(self, line, index):
        words = line.split(',')
        solv_list = []
        for w in words[0:4]:
            if len(w) > 0:
                solv_list.append(w)
        self.solv_list = np.array(solv_list)
        self.salt = words[4]
        self.solv_ratio = []
        for w in words[5:9]:
            if len(w) > 0:
                self.solv_ratio.append(float(w))
        self.solv_ratio = np.array(self.solv_ratio)
        self.salt_c = float(words[9])
        self.T = float(words[10])
        self.conductivity = float(words[11])
        self.index = index

def read_mol_lib(ifn):
    ifile = open(ifn, 'r')
    ifile.readline()
    ifile.readline()
    mol_lib = {}
    i = 0
    for line in ifile:
        words = line.split()
        if len(words) == 3:
            mol_lib[i] = Molecule(words[0], int(words[1]), words[2])
            i += 1
    ifile.close()
    return mol_lib

def read_formula_lib(ifn):
    ifile = open(ifn, 'r')
    ifile.readline()
    formulas = {}
    i = 0
    for line in ifile:
        formulas[i] = Formula(line, i)
        i += 1
    ifile.close()
    return formulas

def find_formual_with_component(mol, lib_formula):
    formulas = []
    for i_formula in lib_formula:
        formula = lib_formula[i_formula]
        if mol.smiles in formula.solv_list or mol.smiles == formula.salt:
            formulas.append(formula)
    return formulas

def find_by_smiles(smiles, lib_mols):
    for index in lib_mols:
        if lib_mols[index].smiles == smiles:
            return index, lib_mols[index]
    return None

if __name__ == '__main__':
    lib_mols = read_mol_lib('res1.xvg')
    lib_formula = read_formula_lib('edb1.csv')
    i_mol, mol = find_by_smiles('COCOC', lib_mols)
    # for formula in find_formual_with_component(mol, lib_formula):
    i_formulas = []
    ifile = open('index_formula.txt', 'r')
    for line in ifile:
        i_formulas.append(int(line))
    ifile.close()

    for index in i_formulas:
        formula = lib_formula[index]
        i, salt = find_by_smiles(formula.salt, lib_mols)
        salt.count += 1

    count = 0
    for i_mol in lib_mols:
        mol = lib_mols[i_mol]
        if mol.is_salt:
            print(mol.smiles, mol.freq, mol.count)
            count += mol.count
    print(count)

