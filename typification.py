#!/usr/bin/env python
import jax
# import dmff
import rdkit
import rdkit.Chem as Chem
from rdkit.Chem import Draw
import rdkit.Chem.rdchem
from rdkit.Chem.rdchem import BondType,HybridizationType
import numpy as np
import sys
import subprocess

BOND_LABELS = {BondType.SINGLE: '-',
               BondType.DOUBLE: '=',
               BondType.TRIPLE: '#',
               BondType.AROMATIC:'@'
        }

HYB_LABELS = {HybridizationType.UNSPECIFIED: 'U',
        HybridizationType.S: 'S',
        HybridizationType.SP: 'SP',
        HybridizationType.SP2: 'SP2',
        HybridizationType.SP3: 'SP3',
        HybridizationType.SP2D: 'SP2D',
        HybridizationType.SP3D: 'SP3D',
        HybridizationType.SP3D2: 'SP3D2',
        HybridizationType.OTHER: 'O'}

def mk_dir(dirname):
    subprocess.call(['rm', '-r', dirname])
    subprocess.call(['mkdir', dirname])
    return

def typify_atom(atom, atypes, depth=0, excl=None):
    i = atom.GetIdx()
    if depth == 0:
        return atom.GetSymbol() + '/%s'%atom.GetHybridization().name
    else:
        atype = atypes[i]
        atype_nbs = []
        btype_nbs = []
        htype_nbs = []
        # for neighbor in atom.GetNeighbors():
        for b in atom.GetBonds():
            if b.GetBeginAtomIdx() == i:
                neighbor = b.GetEndAtom()
            else:
                neighbor = b.GetBeginAtom()
            j = neighbor.GetIdx()
            if j != excl:
                atype_nbs.append(typify_atom(neighbor, atypes, depth=depth-1, excl=i))
                btype_nbs.append(BOND_LABELS[b.GetBondType()])
                htype_nbs.append(HYB_LABELS[neighbor.GetHybridization()])
        atype_nbs = np.array(atype_nbs)
        btype_nbs = np.array(btype_nbs)
        htype_nbs = np.array(htype_nbs)
        # atype_nbs = np.array([btype_nbs[i] + atype_nbs[i] + '/' + htype_nbs[i] for i in range(len(atype_nbs))])
        atype_nbs = np.array([btype_nbs[i] + atype_nbs[i] for i in range(len(atype_nbs))])
        order = np.argsort(atype_nbs)
        atype_nbs = atype_nbs[order]
        if len(atype_nbs) == 0:
            return atype
        else:
            atype = atype + '(' + ','.join(atype_nbs) + ')'
            return atype
    #     atypes = np.array([a.GetSymbol()+'(%s)'%a.GetHybridization().name for a in mol.GetAtoms()])
    #     atype_nbs = []
    #     for 
    #     for j in np.where(self.connectivity[i] == 1)[0]:
    return



def typify_mol(smiles, depth=1, withH=False):
    mol = Chem.MolFromSmiles(smiles)
    atom_types = []
    atypes = np.array([a.GetSymbol() for a in mol.GetAtoms()])
    for a in mol.GetAtoms():
        atom_types.append(typify_atom(a, atypes, depth=depth))
    if withH:
        mol1 = Chem.AddHs(mol)
        atoms = list(mol1.GetAtoms())
        # loop over H
        for ia in range(mol.GetNumAtoms(), mol1.GetNumAtoms()):
            a = atoms[ia]
            ihost = a.GetNeighbors()[0].GetIdx()
            atom_types.append('H-'+atom_types[ihost])
        return mol1, list(mol1.GetAtoms()), np.array(atom_types)
    return mol, list(mol.GetAtoms()), np.array(atom_types)


if __name__ == '__main__':

    depth = 1

    ifile = open('mols_target.dat', 'r')
    ifile.readline()
    ifile.readline()
    ifile.readline()
    lib_mols = {}
    list_mol_smiles = []
    list_salt_smiles = []
    lib_atomtypes = []
    flag = 'mol'
    for line in ifile:
        if '-----' in line:
            flag = 'salt'
            continue
        words = line.split()
        lib_mols[words[0]] = {'freq': int(words[1]), 'availability': words[2]}
        if flag == 'mol':
            list_mol_smiles.append(words[0])
        else:
            list_salt_smiles.append(words[0])
    ifile.close()
    # build the full molecule library
    for smiles in lib_mols:
        _, atoms, atypes = typify_mol(smiles, depth=depth)
        for atype in atypes:
            if atype not in lib_atomtypes:
                lib_atomtypes.append(atype)
    # print out all atomtypes
    ofile = open('atomtypes.dat', 'w')
    for itype in range(len(lib_atomtypes)):
        print(itype, lib_atomtypes[itype], file=ofile)
    ofile.close()

    print(len(lib_atomtypes))
    print('-----')

    print(len(list_mol_smiles), len(list_salt_smiles))
    # double check coverage
    atypes_covered = []
    mols_selected = []
    for smiles in list_salt_smiles:
        _, atoms, atypes = typify_mol(smiles, depth=depth)
        for atype in atypes:
            if atype not in atypes_covered:
                atypes_covered.append(atype)
    print(len(atypes_covered))

    n_mol_selected = 0
    for smiles in list_mol_smiles:
        if lib_mols[smiles]['availability'] == 'N':
            continue
        _, atoms, atypes = typify_mol(smiles, depth=depth)
        flag = 0
        for atype in atypes:
            if atype not in atypes_covered:
                atypes_covered.append(atype)
                flag = 1
        # a molecule that contains new atomtypes
        if flag == 1:
            mols_selected.append(smiles)
            n_mol_selected += 1
            print(n_mol_selected, len(atypes_covered))

    ofile = open('mols_selected.dat', 'w')
    mk_dir('mols_selected')
    for i, smiles in enumerate(mols_selected):
        Draw.MolToFile(rdkit.Chem.MolFromSmiles(smiles), 'mols_selected/%d_%d.png'%(i, lib_mols[smiles]['freq']), size=((300, 300)))
        print(smiles, file=ofile)
    ofile.close()
    ofile = open('salts_selected.dat', 'w')
    mk_dir('salts_selected')
    for i, smiles in enumerate(list_salt_smiles):
        Draw.MolToFile(rdkit.Chem.MolFromSmiles(smiles), 'salts_selected/%d_%d.png'%(i, lib_mols[smiles]['freq']), size=((300, 300)))
        print(smiles, file=ofile)
    ofile.close()


    # analyze which molecules are not covered?
    ofile = open('mols_not_covered.dat', 'w')
    mk_dir('mols_not_covered')
    i = 0
    for smiles in list_mol_smiles:
        _, atoms, atypes = typify_mol(smiles, depth=depth)
        flag = 0
        atypes_not_covered = []
        for atype in atypes:
            if atype not in atypes_covered:
                flag = 1
                if atype not in atypes_not_covered:
                    atypes_not_covered.append(atype)
        if flag == 1:
            print('%-50s'%smiles, ' /// '.join(atypes_not_covered), file=ofile)
            Draw.MolToFile(rdkit.Chem.MolFromSmiles(smiles), 'mols_not_covered/%d_%d.png'%(i, lib_mols[smiles]['freq']), size=((300, 300)))
            i += 1
    ofile.close()
