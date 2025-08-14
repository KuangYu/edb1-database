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
# import xml.etree.ElementTree as etree
from lxml import etree
import re
import math
import copy

MASSES = {'H': 1.00784,
        'C': 12.011,
        'O': 15.999,
        'N': 14.0067,
        'S': 32.065,
        'P': 30.973762,
        'F': 18.9984,
        'As': 74.9216,
        'B': 10.811,
        'Si': 28.0855
        }

THRESH = 1e-5


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


class ForceField:

    # create an empty fftree
    def __init__(self, xml=None, lib_atypes=None):
        if xml is None:
            root = etree.Element('ForceField')
            self.atomtypes = etree.SubElement(root, 'AtomTypes')
            self.residues = etree.SubElement(root, 'Residues')
            self.bonds = etree.SubElement(root,'HarmonicBondForce')
            self.angles = etree.SubElement(root,'HarmonicAngleForce')
            self.torsions = etree.SubElement(root,'PeriodicTorsionForce')
            self.nonbonded = etree.SubElement(root,'NonbondedForce', coulomb14scale="0.5", lj14scale="0.5")
            etree.SubElement(self.nonbonded, 'UseAttributeFromResidue', name="charge")
            self.fftree = etree.ElementTree(root)
            self.root = root
            self.lib_atypes = read_atypes_defs(lib_atypes)
        else:
            self.fftree = etree.parse(xml)
            fftree = self.fftree
            self.root = fftree.getroot()
            self.atomtypes = self.root.find('AtomTypes')
            self.residues = self.root.find('Residues')
            self.bonds = self.root.find('HarmonicBondForce')
            self.angles = self.root.find('HarmonicAngleForce')
            self.torsions = self.root.find('PeriodicTorsionForce')
            self.nonbonded = self.root.find('NonbondedForce')
            etree.SubElement(self.nonbonded, 'UseAttributeFromResidue', name="charge")
            self.lib_atypes = {}
            for atype in self.atomtypes:
                name = atype.attrib['name']
                aclass = atype.attrib['class']
                element = atype.attrib['element']
                toplabel = ''
                idx = int(re.sub('opls_', '', name)) - 800
                self.lib_atypes[idx] = {
                    'idx': idx,
                    'toplabel': toplabel,
                    'element': element,
                    'name': name,
                    'class': aclass
                    }


    def add_atype(self, idx, sigma=-1.0, epsilon=-1.0):
        lib_atypes = self.lib_atypes
        # use -1.0 as placeholder
        # check consistency with existing atom types
        element = lib_atypes[idx]['element']
        name = lib_atypes[idx]['name']
        aclass = lib_atypes[idx]['class']
        atom = self.get_nbtype_node(idx)
        if atom is not None:
            # a hit, check parameter consistency
            if atom.attrib['type'] == name:
                flag = 1
                if sigma is not None:
                    if abs(float(atom.attrib['sigma']) - float(sigma)) > THRESH:
                        sys.exit('Conficting sigma for %s: %s %f'%(name, atom.attrib['sigma'], sigma))
                if epsilon is not None:
                    if abs(float(atom.attrib['epsilon']) - float(epsilon)) > THRESH:
                        sys.exit('Conficting sigma for %s: %s %f'%(name, atom.attrib['epsilon'], epsilon))
        # if this is a unseen type, add it
        else:
            # add atypes
            etree.SubElement(self.nonbonded, 'Atom', type=name, sigma="%.6f"%sigma, epsilon="%.6f"%epsilon)
            node = etree.SubElement(self.atomtypes, 'Type', name=name, mass="%.6f"%MASSES[element])
            node.set('class', aclass)
        return

    def get_lj_params(self, idx):
        name = 'opls_' + '%d'%(idx+800)
        for atom in self.nonbonded.findall('Atom'):
            if atom.attrib['type'] == name:
                return {'sigma':float(atom.attrib['sigma']), 'epsilon':float(atom.attrib['epsilon'])}

    # find if atype is in folder?
    def get_nbtype_node(self, idx):
        name = self.lib_atypes[idx]['name']
        for atom in self.nonbonded.findall('Atom'):
            if atom.attrib['type'] == name:
                return atom
        return None

    def get_bondtype_node(self, idx1, idx2):
        class1 = self.lib_atypes[idx1]['class']
        class2 = self.lib_atypes[idx2]['class']
        for b in self.bonds:
            indices = set([b.attrib['class1'], b.attrib['class2']])
            if indices == set([class1, class2]):
                return b
        return None

    def get_angletype_node(self, idx1, idx2, idx3):
        class1 = self.lib_atypes[idx1]['class']
        class2 = self.lib_atypes[idx2]['class']
        class3 = self.lib_atypes[idx3]['class']
        for a in self.angles:
            indices = np.array([a.attrib['class1'], a.attrib['class2'], a.attrib['class3']])
            if np.all(indices == np.array([class1, class2, class3])) or np.all(indices == np.array([class3, class2, class1])):
                return a
        return None

    def get_torsiontype_node(self, idx1, idx2, idx3, idx4, torsiontype):
        class1 = self.lib_atypes[idx1]['class']
        class2 = self.lib_atypes[idx2]['class']
        class3 = self.lib_atypes[idx3]['class']
        class4 = self.lib_atypes[idx4]['class']
        for d in self.torsions:
            indices = np.array([d.attrib['class1'], d.attrib['class2'], d.attrib['class3'], d.attrib['class4']])
            if torsiontype == 'Improper' and indices[0] == class1 \
                    and np.all(np.sort(indices[1:]) == np.sort(np.array([class2, class3, class4]))) \
                    and torsiontype == d.tag:
                return d
            if torsiontype == 'Proper' and (np.all(indices == np.array([class1, class2, class3, class4])) \
                    or np.all(indices == np.array([class4, class3, class2, class1]))):
                return d
        return None

    def add_bondtype(self, idx1, idx2, length=-1.0, k=-1.0):
        lib_atypes = self.lib_atypes
        # check existing bond types
        bond = self.get_bondtype_node(idx1, idx2)
        class1 = lib_atypes[idx1]['class']
        class2 = lib_atypes[idx2]['class']
        if bond is not None:
            # check consistency
            length0 = float(bond.attrib['length'])
            k0 = float(bond.attrib['k'])
            if not math.isclose(length0, length, abs_tol=THRESH):
                sys.exit('Conflicting bond length: %s %s - %f %f'%(class1, class2, length0, length))
            if not math.isclose(k0, k, abs_tol=THRESH):
                sys.exit('Conflicting bond k: %s %s - %f %f'%(class1, class2, k0, k))
        else:
            # new bond type
            etree.SubElement(self.bonds, 'Bond', class1=class1, class2=class2, length='%.6f'%length, k='%.6f'%k)
        return 

    def add_angletype(self, idx1, idx2, idx3, angle=-1.0, k=-1.0):
        lib_atypes = self.lib_atypes
        angletype = self.get_angletype_node(idx1, idx2, idx3)
        class1 = lib_atypes[idx1]['class']
        class2 = lib_atypes[idx2]['class']
        class3 = lib_atypes[idx3]['class']
        if angletype is not None:
            # check consistency
            angle0 = float(angletype.attrib['angle'])
            k0 = float(angletype.attrib['k'])
            if not math.isclose(angle0, angle, abs_tol=THRESH):
                sys.exit('Conflicting angle value: %s %s %s - %f %f'%(class1, class2, class3, angle0, angle))
            if not math.isclose(k0, k, abs_tol=THRESH):
                sys.exit('Conflicting angle k: %s %s %s - %f %f'%(class1, class2, class3, k0, k))
        else:
            etree.SubElement(self.angles, 'Angle', class1=class1, class2=class2, class3=class3, angle='%.6f'%angle, k='%.6f'%k)
        return

    def add_torsiontype(self, idx1, idx2, idx3, idx4, attrib, torsiontype):
        lib_atypes = self.lib_atypes
        torsion = self.get_torsiontype_node(idx1, idx2, idx3, idx4, torsiontype)
        class1 = lib_atypes[idx1]['class']
        class2 = lib_atypes[idx2]['class']
        class3 = lib_atypes[idx3]['class']
        class4 = lib_atypes[idx4]['class']
        attrib['class1'] = class1
        attrib['class2'] = class2
        attrib['class3'] = class3
        attrib['class4'] = class4
        if torsion is not None:
            for k in torsion.attrib.keys():
                if 'class' in k:
                    continue
                if 'periodicity' in k:
                    val0 = int(torsion.attrib[k])
                    val = int(attrib[k])
                    if val0 != val:
                        sys.exit('Conflicting torsion %s: %s %s %s %s - %d %d'%(k, class1, class2, class3, class4, val0, val))
                else:
                    val0 = float(torsion.attrib[k])
                    val = float(attrib[k])
                    if not math.isclose(val0, val, abs_tol=THRESH):
                        sys.exit('Conflicting torsion %s: %s %s %s %s - %f %f'%(k, class1, class2, class3, class4, val0, val))
        else:
            etree.SubElement(self.torsions, torsiontype.lower().capitalize(), attrib=attrib)
        return

    # write out file
    def write(self, ofn='forcefield.xml'):
        self.fftree.write(ofn, pretty_print=True)

    # def add_atype(self, index, element, charge=None, sigma=None, )


def read_atypes_defs(ifn):
    lib_atypes = {}
    with open(ifn) as ifile:
        for line in ifile:
            # deal with heavy atoms
            words = line.split()
            idx = int(words[0])
            toplabel = words[1]
            element = toplabel[0]
            name = 'opls_' + '%d'%(idx+800)
            aclass = element + '%d'%(idx+800)
            lib_atypes[idx] = {
                    'idx': idx,
                    'toplabel': toplabel,
                    'element': element,
                    'name': name,
                    'class': aclass
                    }
            # deal with Hs
            idx = idx + 100
            toplabel = 'H-' + toplabel
            element = 'H'
            name = 'opls_' + '%d'%(idx+800)
            aclass = element + '%d'%(idx+800)
            lib_atypes[idx] = {
                    'idx': idx,
                    'toplabel': toplabel,
                    'element': element,
                    'name': name,
                    'class': aclass
                    }
    return lib_atypes

def name2idx(name):
    return int(re.sub('opls_', '', name)) - 800

def class2idx(name):
    return int(re.sub('[A-Za-z]+', '', name)) - 800

if __name__ == '__main__':
    path = 'ligpargen_files'
    lib_xml_naive, lib_pdb_naive = read_ligpargen_lib('ligpargen_files')
    # Create an empty FF
    forcefield = ForceField(lib_atypes='../atomtypes.dat')
    map_toplabel2atypes = {}
    for idx in forcefield.lib_atypes:
        label = forcefield.lib_atypes[idx]['toplabel']
        map_toplabel2atypes[label] = forcefield.lib_atypes[idx]

    # template 
    list_smiles = ['CC1COC(=O)O1', 'O=C1OCCO1', 'CCOC(=O)OCC']
    # list_smiles = ['CC1COC(=O)O1']
    for i_mol, smiles in enumerate(list_smiles):
        resname = '%d'%i_mol
        if len(resname) == 1:
            resname = 'R0' + resname
        else:
            resname = 'R' + resname
        # template
        mol1_template, atoms, atypes = typify_mol(smiles, depth=1, withH=True)
        # structure from pdb
        mol1 = Chem.MolFromPDBFile(lib_pdb_naive[smiles], removeHs=False)
        # sanity check with template
        if not sanity_check(mol1_template, mol1, range(mol1_template.GetNumAtoms()), range(mol1.GetNumAtoms())):
            print('ERROR: connectivity mismatch')

        ff_naive = ForceField(xml=lib_xml_naive[smiles])

        # deal with atomtypes
        # map from idx in naive xml to new atomtypes
        charges = np.zeros(mol1.GetNumAtoms())
        map_atypes = {}
        for nbtype in ff_naive.nonbonded.findall('Atom'):
            name0 = nbtype.attrib['type']
            sigma = float(nbtype.attrib['sigma'])
            epsilon = float(nbtype.attrib['epsilon'])
            idx0 = name2idx(name0)
            toplabel = atypes[idx0]
            atype = map_toplabel2atypes[toplabel]
            map_atypes[idx0] = atype
            # add atomtypes
            forcefield.add_atype(atype['idx'], sigma=sigma, epsilon=epsilon)
            charges[idx0] = float(nbtype.attrib['charge'])

        # add bonding terms
        for btype in ff_naive.bonds.findall('Bond'):
            idx01 = class2idx(btype.attrib['class1'])
            idx02 = class2idx(btype.attrib['class2'])
            atype1 = map_atypes[idx01]
            atype2 = map_atypes[idx02]
            # class1 = atype1['class']
            # class2 = atype2['class']
            idx1 = atype1['idx']
            idx2 = atype2['idx']
            length = float(btype.attrib['length'])
            k = float(btype.attrib['k'])
            forcefield.add_bondtype(idx1, idx2, length=length, k=k)

        # add angle terms
        for angletype in ff_naive.angles.findall('Angle'):
            idx01 = class2idx(angletype.attrib['class1'])
            idx02 = class2idx(angletype.attrib['class2'])
            idx03 = class2idx(angletype.attrib['class3'])
            atype1 = map_atypes[idx01]
            atype2 = map_atypes[idx02]
            atype3 = map_atypes[idx03]
            idx1 = atype1['idx']
            idx2 = atype2['idx']
            idx3 = atype3['idx']
            angle = float(angletype.attrib['angle'])
            k = float(angletype.attrib['k'])
            forcefield.add_angletype(idx1, idx2, idx3, angle=angle, k=k)

        # add dihedral torsions
        for torsiontype in ff_naive.torsions:
            idx01 = class2idx(torsiontype.attrib['class1'])
            idx02 = class2idx(torsiontype.attrib['class2'])
            idx03 = class2idx(torsiontype.attrib['class3'])
            idx04 = class2idx(torsiontype.attrib['class4'])
            atype1 = map_atypes[idx01]
            atype2 = map_atypes[idx02]
            atype3 = map_atypes[idx03]
            atype4 = map_atypes[idx04]
            idx1 = atype1['idx']
            idx2 = atype2['idx']
            idx3 = atype3['idx']
            idx4 = atype4['idx']
            forcefield.add_torsiontype(idx1, idx2, idx3, idx4, torsiontype.attrib, torsiontype.tag)

        # deal with charge and residue
        top_group_labels = np.array(list(Chem.CanonicalRankAtoms(mol1, breakTies=False)))
        set_labels = set(top_group_labels)
        for label in set_labels:
            idxs = np.where(top_group_labels == label)[0]
            # do average over charge
            stdev = np.std(charges[idxs])
            if stdev > 1e-3:
                print('WARNING: charge variations may be too large: ', charges[idxs])
                print(smiles)
                print(idxs)
            charges[idxs] = np.average(charges[idxs])

        residue = copy.deepcopy(ff_naive.residues.find('Residue'))
        residue.attrib['name'] = resname
        for atom in residue.findall('Atom'):
            name0 = atom.attrib['type']
            idx0 = name2idx(name0)
            atype = map_atypes[idx0]
            name = atype['name']
            atom.attrib['type'] = name
            atom.attrib['charge'] = '%.6f'%charges[idx0]
        forcefield.residues.append(residue)

    # Write out force field
    forcefield.write('forcefield.xml')
