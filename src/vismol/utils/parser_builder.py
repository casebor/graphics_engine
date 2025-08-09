#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from logging import getLogger
from vismol.model.atom import Atom
from vismol.model.residue import Residue
from vismol.model.chain import Chain
from vismol.model.molecule import Molecule


logger = getLogger(__name__)


def add_atom_to_vm_object(vm_object: "VismolObject", coords: np.array, atom_unique_id: int) -> None:
    molecule = Molecule(vm_object)
    vm_object.molecule = molecule
    chain = Chain(vm_object, name="BX", molecule=molecule)
    molecule.chains[chain.name] = chain
    residue = Residue(vm_object, name="RX", index=1, chain=chain)
    chain.residues[residue.index] = residue
    atom = Atom(vm_object, name="X", index=0, residue=residue, chain=chain,
                pos=coords, atom_id=0, molecule=molecule)
    atom.spheres = True
    atom.unique_id = atom_unique_id
    residue.atoms[0] = atom
    molecule.atoms[0] = atom
    molecule.frames = np.empty([1,1,3], dtype=np.float32)
    molecule.frames[0,0] = atom.pos
    molecule.build_bonded_and_nonbonded_atoms()

def add_new_atom(vm_object: "VismolObject", coords: np.array, atom_unique_id: int) -> None:
    pass
