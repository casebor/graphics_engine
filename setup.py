#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension("VisMol.glCore.matrix_operations",
              sources=["VisMol/glCore/matrix_operations.pyx"],
              include_dirs=[np.get_include()]),
    Extension("VisMol.glCore.sphere_representation",
              sources=["VisMol/glCore/sphere_representation.pyx"],
              include_dirs=[np.get_include()]),
    Extension("VisMol.vBabel.PDBFiles",
              sources=["VisMol/vBabel/PDBFiles.pyx"],
              include_dirs=[np.get_include()]),
    Extension("VisMol.vBabel.GROFiles",
              sources=["VisMol/vBabel/GROFiles.pyx"],
              include_dirs=[np.get_include()]),
    Extension("VisMol.vModel.cDistances",
              sources=["VisMol/vModel/cDistances.pyx"],
              include_dirs=[np.get_include()]),
    # Extension("VisMol.vModel.VismolObject",
    #           sources=["VisMol/vModel/VismolObject.pyx"],
    #           include_dirs=[np.get_include()]),
    # Extension("VisMol.vModel.Atom",
    #           sources=["VisMol/vModel/Atom.pyx"],
    #           include_dirs=[np.get_include()])
]

setup(
    name="VisMol",
    ext_modules=cythonize(extensions),
)
