################################################################################
# 24.09.2024
# - This Python file holds all the global variables to manipulated by the other
#   Python code files qmmm.py and cmd_line.py
################################################################################

# Import CCDC modules (CSD Python API)
from ccdc          import io                    # read crystal structures
csd_reader = io.EntryReader('CSD')

import ccdc.search
from ccdc.search   import TextNumericSearch     # search by text or numbers
from ccdc.search   import SubstructureSearch    # search for substructures

from ccdc.molecule import Molecule              # Build a molecule
from ccdc.molecule import Atom                  # Atomic data and properties
from ccdc.molecule import Bond                  # Bond properties

# Import OpenBabel to assist with atom classification
from openbabel import openbabel as ob           # basic OpenBabel

def init():
  global VerboseFlag
  VerboseFlag = 0
  global mol, obmol
  mol   = Molecule()
  obmol = ob.OBMol()
  global qm, mm, dnts
  qm = []
  mm = []
  dnts = {}
  