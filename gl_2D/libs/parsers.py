#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import libs.objects as vmos

def parse_pdb(pdb_file):
    """ Function doc """
    atoms = []
    gl_id = 0
    with open(pdb_file, "r") as pdbin:
        for line in pdbin:
            if line.startswith("ATOM"):
                atom_id = int(line[6:11])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                name = line[12:16].strip()
                atom = vmos.VMAtom(gl_id, atom_id, x, y, z, name=name)
                atoms.append(atom)
                gl_id += 1
    return vmos.GLContainer(atoms)

def parse_obj(obj_file):
    """ Function doc """
    i, j = 0, 0
    nodes = []
    indexes = []
    with open(obj_file, "r") as objin:
        for line in objin:
            if line.startswith("v "):
                vc = [float(v) for v in line.split()[1:]]
                node = vmos.VMNode(i, *vc)
                nodes.append(node)
                i += 1
            elif line.startswith("vn "):
                vn = [float(v) for v in line.split()[1:]]
                nodes[j].set_normal(vn)
                j += 1
            elif line.startswith("f "):
                chunks = line.split()[1:]
                face = [int(chunk.split("/")[0])-1 for chunk in chunks]
                indexes.extend(face)
    return vmos.GLContainer(nodes, indexes)
