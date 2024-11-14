################################################################################
# 24.09.2024 (Q15_CCDC)
# - This Python file holds all the global variables to manipulated by the other
#   Python code files qmmm.py and cmd_line.py
# 31.10.2024 (Q16_CCDC)
# - fix the problem with atoms in a spiro substructure
# - make the abgle arument for the twist in comjugated links
#   a) tunable -> change the switch-off angle
#   b) switchable -> on/off sttings for the whole check
#   angle : -1  => test is turned off
#   angle :  0  => use the default of 30 degrees
#   angle : >0  => use this angle to switch
#   add the global variable conju_twist_angle with a default value of 30 deg
# - code added fix missing H atoms
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
  global conju_twist_angle
  conju_twist_angle = 30.0  #  twist angle in degrees
  