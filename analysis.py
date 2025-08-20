#!/usr/bin/env python
import sys
import numpy as np
from typification import mk_dir

ifile = open('edb1.csv', 'r')
ifile.readline()
molecules = {}
salts = {}
formulas = {}

i_formula = 0
for line in ifile:
    words = line.split(',')
    formula = []
    for word in words[0:4]:
        if len(word) == 0:
            continue
        if word in molecules:
            molecules[word] += 1
        else:
            molecules[word] = 1
        formula.append(word)

        # if '[Li+]' in word or 'Na' in word:
    salt = words[4]
    if salt in salts:
        salts[salt] +=1
    else:
        salts[salt] = 1
    formula.append(salt)
    formulas[i_formula] = np.array(formula)
    i_formula += 1
ifile.close()

# sort the molecules and salts according to frequencies
mol_list = np.array(list(molecules.keys()))
n_formula = np.array([molecules[m] for m in mol_list])
indices = np.argsort(n_formula)
mol_list = mol_list[indices[::-1]]
n_formula_mol = n_formula[indices[::-1]]

salt_list = np.array(list(salts.keys()))
n_formula = np.array([salts[m] for m in salt_list])
indices = np.argsort(n_formula)
salt_list = salt_list[indices[::-1]]
n_formula_salt = n_formula[indices[::-1]]

n_mol = len(mol_list)
n_salt = len(salt_list)
# for i in range(n_mol):
#     print(mol_list[i], n_formula_mol[i])

# print('-------')
# for i in range(n_salt):
#     print(salt_list[i], n_formula_salt[i])

# filter

n_thresh_mol = 10
n_thresh_salt = 20

# filter_mol = (n_formula_mol>n_thresh_mol) #* np.array(['Si' not in w for w in mol_list])
mol_excl = ['CS(N)(=O)=O', 'COC=O']
filter_mol = (n_formula_mol>n_thresh_mol) * np.array([w not in mol_excl for w in mol_list])
filter_salt = np.array((n_formula_salt>n_thresh_salt)) * np.array(['Na' not in w for w in salt_list])

mol_list = mol_list[filter_mol]
n_formula_mol = n_formula_mol[filter_mol]
salt_list = salt_list[filter_salt]
salt_list = np.array([
'[Li+].F[P-](F)(F)(F)(F)F',
'[Li+].[B-]12(OC(=O)C(=O)O1)OC(=O)C(=O)O2',
'[Li+].C(F)(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F',
'[Li+].F[As-](F)(F)(F)(F)F',
'CC[N+](CC)(CC)CC.F[P-](F)(F)(F)(F)F',
'[B-](F)(F)(F)F.CC[N+](CC)(CC)CC',
'[Li+].[B-](F)(F)(F)F',
'[Li+].[O-]Cl(=O)(=O)=O',
'[Li+].C(F)(F)(F)S(=O)(=O)[O-]',
'[Li+].[N-](S(=O)(=O)F)S(=O)(=O)F',
# '[Li+].FC(F)(F)C1=NC2=C([N-]1)C=C(C#N)C(=C2)C#N',
# '[Li+].C1=CC=C2C(=C1)O[P-]34(O2)(OC5=CC=CC=C5O3)OC6=CC=CC=C6O4',
'[Li+].[B-]1(OC(=O)C(=O)O1)(F)F',
'[Li+].C(C(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(C(F)(F)F)(F)F)(F)(F)F',
# '[Li+].c1ccc2O[B-]3(Oc2c1)Oc1ccccc1O3',
'[B-](F)(F)(F)F.CC[N+](C)(CC)CC',
'[O-]S(=O)(=O)C(F)(F)F.CCCC[N+](CCCC)(CCCC)CCCC',
'CCCC[N+](CCCC)(CCCC)CCCC.[O-]Cl(=O)(=O)=O',
# '[Li+].FC(F)(F)S(=O)(=O)[C-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F',
# '[Li+].[B-]12(OC3=CC4=CC=CC=C4C=C3O1)OC5=CC6=CC=CC=C6C=C5O2',
# '[B-](C1=CC=CC=C1)(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4.CCCC[N+](CCCC)(CCCC)CCCC',
# 'CCN1C=C[N+](=C1)C.[B-]12(OC(=O)C(=O)O1)OC(=O)C(=O)O2',
# 'CCN1C=C[N+](=C1)C.O=C1O[B-](F)(F)OC1=O',
# 'CC[N+](CC)(CC)CC.[B-]12(OC(=O)C(=O)O1)OC(=O)C(=O)O2',
# 'CC[N+](CC)(CC)CC.O=C1O[B-](F)(F)OC1=O',
])

n_formula_salt = n_formula_salt[filter_salt]

ofile = open('index_formula.txt', 'w')
i = 0
for i_formula in formulas:
    f = formulas[i_formula]
    flag = 0
    for comp in f:
        if comp not in mol_list and comp not in salt_list:
            flag = 1
    if flag == 0:
        print(i_formula, file=ofile)
        i += 1
ofile.close()

# price and availability information
prices = {}
ifile = open('prices.csv', 'r')
for line in ifile:
    if '------' in line:
        continue
    words = line.split(',')
    mol = words[1]
    if len(words[5]) > 0 and words[5].startswith('Y'):
        price = float(words[5][1:])
    else:
        price = -1.000
    prices[mol] = price
ifile.close()


# Write the targeted molecules
ofile = open('mols_target.dat', 'w')
print('Covered # formula:', i, file=ofile)
print('Target # mols and # salts:', len(mol_list), len(salt_list), file=ofile)

print('-----------', file=ofile)
for i_mol in range(len(mol_list)):
    if prices[mol_list[i_mol]] > 1000:
        flag = 'YE'
    elif prices[mol_list[i_mol]] > 0:
        flag = 'Y'
    else:
        flag = 'N'
    print(mol_list[i_mol], n_formula_mol[i_mol], flag, file=ofile)

print('-----------', file=ofile)
for i_salt in range(len(salt_list)):
    if prices[salt_list[i_salt]] > 1000:
        flag = 'YE'
    elif prices[salt_list[i_salt]] > 0:
        flag = 'Y'
    else:
        flag = 'N'
    print(salt_list[i_salt], n_formula_salt[i_salt], flag, file=ofile)

ofile.close()

# draw structure
from rdkit import Chem
from rdkit.Chem import Draw

mk_dir('salts_target')
mk_dir('mols_target')
for i_salt in range(len(salt_list)):
    mol = Chem.MolFromSmiles(salt_list[i_salt])
    smiles = salt_list[i_salt]
    if prices[smiles] > 0:
        flag = 'Y'
        if prices[smiles] > 1000:
            flag += 'E'
    else:
        flag = 'N'
    Draw.MolToFile(mol, 'salts_target/%d_%d_%s.png'%(i_salt, n_formula_salt[i_salt], flag), size=(300, 300))

for i_mol in range(len(mol_list)):
    mol = Chem.MolFromSmiles(mol_list[i_mol])
    smiles = mol_list[i_mol]
    if prices[smiles] > 0:
        flag = 'Y'
        if prices[smiles] > 1000:
            flag += 'E'
    else:
        flag = 'N'
    Draw.MolToFile(mol, 'mols_target/%d_%d_%s.png'%(i_mol, n_formula_mol[i_mol], flag), size=(300, 300))
