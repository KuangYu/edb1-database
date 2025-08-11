#!/usr/bin/env python
import jax
# import dmff
import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.rdchem
from rdkit.Chem.rdchem import BondType,HybridizationType
import numpy as np

BOND_LABELS = {BondType.SINGLE: '-',
               BondType.DOUBLE: '=',
               BondType.TRIPLE: '#',
               BondType.AROMATIC:'@'
        }


def typify_atom(atom, atypes, depth=0, excl=None):
    if depth == 0:
        return atom.GetSymbol() + '(%s)'%atom.GetHybridization().name
    else:
    #     atypes = np.array([a.GetSymbol()+'(%s)'%a.GetHybridization().name for a in mol.GetAtoms()])
    #     atype_nbs = []
    #     for 
    #     for j in np.where(self.connectivity[i] == 1)[0]:
    return



def typify_mol(smiles, depth=0):
    mol = Chem.MolFromSmiles(smiles)
    atom_types = []
    atypes = np.array([a.GetSymbol()+'(%s)'%a.GetHybridization().name for a in mol.GetAtoms()])
    for a in mol.GetAtoms():
        atom_types.append(typify_atom(a, atypes, depth=depth))
    return np.array(atom_types)


if __name__ == '__main__':

    # smiles = 'CCOC(=O)OCC'
    smiles = 'Cc1ccccc1'
    print(typify_mol(smiles))
