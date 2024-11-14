################################################################################
# 12.11.2024 (Q17_CCDCTY, adjustments for Ting Yee's Au cluster)
# This Python file holds all the global variables to manipulated by the other
# Python code files qmmm.py and cmd_line.py
# - add a switch vor distant neighbours via the coordination number
# - add a version string
################################################################################

# set version string
Version='gb18-20241113'

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
  global conju_twist_angle
  conju_twist_angle = 30.0  #  twist angle in degrees
  global met_min_lig_num
  met_min_lig_num = 6       # minimum number of ligands attached to metal atoms
  global del_H_from_list
  del_H_from_list=bool(True)