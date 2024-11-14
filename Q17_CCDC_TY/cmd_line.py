################################################################################
# 12.11.2024 (Q17_CCDCTY, adjustments for Ting Yee's Au cluster)
# cmd line version of the QM/MM code using the Python import methods 
# - add a switch vor distant neighbours via the coordination number
# - add a version string
################################################################################

# set version string
Version='cm17-20241112'

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

# add version numbers to the output
if gb.VerboseFlag > 0:
  print("Versions used for thisQM/MM separation")
  print("cmd line    :",     Version) 
  print("QM/MM module:", qms.Version) 
  print("qlobal vars :" ,  gb.Version) 
  print()

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
# This is the critical part of the QM/MM separation and can be copied into
# other Python scripts
################################################################################

# read input file
gb.mol, gb.obmol = qms.read_input_file(inp_name)

# optional part ofthe setup
if gb.VerboseFlag>0: print()
# modification of the default settings for the QM/MM separation
# set the twist angle for conjugated chains
# -1: turn the check for twisted chains off
#  0: set the default value of 30 deg
qms.set_conju_angle(0) # request default value
# set the minimum coordination number for metal atoms
# -1: turn the check for distant neighbors off
#  0: set the default value of 6 (octahedral coordination)
qms.set_metal_coor(3)  # match it to reality
# verify the changed settings
if gb.VerboseFlag > 0:
  print('Parameter which can be changed by the user')
  print('Twist angle for conjugated chains', qms.get_conju_angle())
  print('Minimum metal coordiantion number', qms.get_metal_coor())

# standard QM/MM separation (deviation from the old code)
qms.std_qmmm_sep()

# write output file with the QM/MM separation data
stats = qms.write_output_file("combi.xyz")
if gb.VerboseFlag>0:
  print("List with a summary of the results", stats)

# end the code
exit()