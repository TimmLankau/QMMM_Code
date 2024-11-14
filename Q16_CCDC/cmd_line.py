################################################################################
# 31.10.2024 (Q15_CCDC)
# - start building a cmd line version of the QM/MM code using the Python import
#   methods 
# 31.10.2024 (Q16_CCDC)
# - fix the problem with atoms in a spiro substructure
# - make the abgle arument for the twist in comjugated links
#   a) tunable -> change the switch-off angle
#   b) switchable -> on/off sttings for the whole check
#   angle : -1  => test is turned off
#   angle :  0  => use the default of 30 degrees
#   angle : >0  => use this angle to switch
# - code added fix missing H atoms
################################################################################

# Import basic Python modules
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import math                     # mathematical  functions

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

# Import my modules
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q16_CCDC')
import globals as     gb      # short for globals
import qmmm    as     qms     # short for qm/mm separation
from   qmmm    import printf
from   qmmm    import fprintf

################################################################################
# general purpose function                                                     #
# function to delete global variables to limit the scope of these variables    #
# input  vlist  a list with the names (strings) of the variables to be         #
#               deleted.                                                       #
# output anz    integer with the number of deleted variables                   #
################################################################################

def clean_glob_var(vlist):
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  anz = 0
  g = globals()
  for var in vlist:
    if var in globals():
      if gb.VerboseFlag==2: printf("%s in globals\n", var)
      anz += 1
      del g[var]
    else:
      if gb.VerboseFlag==2: printf("%s NOT in globals\n", var)
  if gb.VerboseFlag==2: printf("%i global variables deleted\n", anz)
  return(anz)

################################################################################

# VERY IMPORTANT, initialize globals
gb.init()

# Managing the cmd line is not part of the Python import
my_arglist = sys.argv.copy() # get command line argument
my_arglist.pop(0) # remove the command name

# look for optional verbose parameter
argpos = -1
for argcnt, my_arg in enumerate(my_arglist):
  if my_arg == "-v":
    argpos = argcnt
    break
my_arglist.pop(argpos)
my_arg = my_arglist.pop(argpos)
if my_arg in ["0", "1", "2"]:
  gb.VerboseFlag = int(my_arg)
else:
  printf("bad command line parameter\n")
  printf("use %s [-v 0,1,2] file name\n", sys.argv[0])
  exit()

# the last surviving argument is the file name
if len(my_arglist) != 1:
    printf("bad command line parameter\n")
    printf("use %s [-v 0,1,2] file name\n", sys.argv[0])
    exit()
inp_name = my_arglist[0]

if gb.VerboseFlag > 0:
  printf("Summarize command line parameter\n")
  printf("  VerboseFlag = %i", gb.VerboseFlag)
  if gb.VerboseFlag == 1:
    printf(" (default or cmd line)\n")
  else:
    printf("\n")
  printf("  Input File  = %s\n", inp_name)

# delete variables I don't need anay more
clean_glob_var(['my_arg', 'my_arglist', 'argpos', 'argcnt'])

# turn of OpenBabel warnings for no verbose runs
if gb.VerboseFlag==0:
  ob.obErrorLog.SetOutputLevel(0)

# clean up old files
if gb.VerboseFlag > 0: printf("Clean up old files\n")
for file in glob.glob("stp*xyz"):
  os.remove(file)
if os.path.exists("./combi.xyz"): os.remove("./combi.xyz")

################################################################################

# read input file
gb.mol, gb.obmol = qms.read_input_file(inp_name)

# set the twist angle for conjugated chains
qms.set_conju_angle(50)
# use qms.get_conju_angle() to get the current value
print(qms.get_conju_angle())

# standard QM/MM separation (deviation from the old code)
qms.std_qmmm_sep()

# write output file with the QM/MM separation data
stats = qms.write_output_file("combi.xyz")
print(stats)