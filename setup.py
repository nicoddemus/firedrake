from distutils.core import setup
from distutils.extension import Extension
from glob import glob
from os import environ as env
from pyop2.utils import get_petsc_dir
import numpy as np
import petsc4py

import versioneer
versioneer.versionfile_source = 'firedrake/_version.py'
versioneer.versionfile_build = 'firedrake/_version.py'
versioneer.tag_prefix = 'v'
versioneer.parentdir_prefix = 'firedrake-'
versioneer.VCS = "git"

cmdclass = versioneer.get_cmdclass()

try:
    from Cython.Distutils import build_ext
    cmdclass['build_ext'] = build_ext
    dmplex_sources = ["firedrake/dmplex.pyx"]
    evtk_sources = ['evtk/cevtk.pyx']
except ImportError:
    # No cython, dmplex.c must be generated in distributions.
    dmplex_sources = ["firedrake/dmplex.c"]
    evtk_sources = ['evtk/cevtk.c']

if 'CC' not in env:
    env['CC'] = "mpicc"

petsc_dirs = get_petsc_dir()
include_dirs = [np.get_include(), petsc4py.get_include()]
include_dirs += ["%s/include" % d for d in petsc_dirs]

setup(name='firedrake',
      version=versioneer.get_version(),
      cmdclass=cmdclass,
      description="""Firedrake is an automated system for the portable solution
          of partial differential equations using the finite element method
          (FEM)""",
      author="Imperial College London and others",
      author_email="firedrake@imperial.ac.uk",
      url="http://firedrakeproject.org",
      packages=["firedrake", "evtk"],
      package_data={"firedrake": ["firedrake_geometry.h"]},
      scripts=glob('scripts/*'),
      ext_modules=[Extension('firedrake.dmplex',
                             sources=dmplex_sources,
                             include_dirs=include_dirs,
                             libraries=["petsc"],
                             extra_link_args=["-L%s/lib" % d for d in petsc_dirs] +
                             ["-Wl,-rpath,%s/lib" % d for d in petsc_dirs]),
                   Extension('evtk.cevtk', evtk_sources,
                             include_dirs=[np.get_include()])])
