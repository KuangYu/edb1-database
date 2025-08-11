#!/usr/bin/env python
import sys

list_smiles = []

with open('mols_selected.dat', 'r') as ifile:
    for line in ifile:
        words = line.split()
        list_smiles.append(words[0])
with open('salts_selected.dat', 'r') as ifile:
    for line in ifile:
        words = line.split()
        list_smiles.append(words[0])

order_info = {}

with open('order_info.csv', 'r') as ifile:
    for line in ifile:
        words = line.split(',')
        if len(words) < 8:
            continue
        smiles = words[1]
        val = words[7]
        if len(smiles) > 0:
            order_info[smiles] = line

for s in list_smiles:
    print(order_info[s])
