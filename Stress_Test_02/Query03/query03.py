# 19.09.2024 - Building a data base for the stress tests
# - Starting with ccdc_t01, but use my experience from ccdc_q13.py

# Load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import math                     # mathematical  functions

# Load modules and prepare to read the CCS data base
import ccdc.search

from ccdc          import io                    # read crystal structures
from ccdc.search   import Search
from ccdc.search   import TextNumericSearch     # search by text or numbers
from ccdc.search   import SubstructureSearch    # search for substructures

from ccdc.molecule import Molecule              # Build a molecule
from ccdc.molecule import Atom                  # Atomic data and properties
from ccdc.molecule import Bond                  # Bond properties

# Fake printf and fprintf as used in C
def printf(format, *args):
  sys.stdout.write(format % args)
def fprintf(stream, format, *args):
  stream.write(format % args)

################################################################################
# Set up environment for the data base search                                  #
################################################################################

basedir = os.path.dirname(os.path.realpath(__file__))
resdir = basedir + "/Results"
printf("%-16s: ", "Base directory")
printf("%s\n",  basedir)
printf("%-16s: ", "Res directory")
printf("%s\n",  resdir)

max_read = 2000000
max_hits =   10000
printf("%-16s: %7i (about 1.5 mio CCDC entries)\n", "max reads", max_read)
printf("%-16s: %7i\n", "max hits", max_hits)
printf("\n")

if os.path.isdir(resdir):
  print("Directory for results exists")
else:
  print("Directory for results does NOT exist -> Create")
  os.makedirs(resdir)

################################################################################
# Search the CCDC database                                                     #
################################################################################

# Define target transition metal
target="Fe"
# list of unwanted transtion metals
tm=["Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y" , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg"]
# remove Fe from the list of transition metals
tm.remove(target)
csd_reader = io.EntryReader('CSD')
cnt=0
csd_cnt=0
no_3d_struc = 0
disorder = 0
too_small = 0
not_target = 0
tm_metal = 0
organo_metallic = 0
polymeric = 0
coordination = 0
target_cnt = 0
for entry in csd_reader:
  # savefty break
  if csd_cnt==max_read: break
  csd_cnt += 1
  # using entry information to preselect
  if not entry.is_organometallic: continue
  organo_metallic += 1
  if entry.is_polymeric:
    polymeric += 1
    continue
  if not entry.has_3d_structure:
    no_3d_struc += 1
    continue
  if entry.has_disorder:
    disorder += 1
    continue
  
  # select based on molecular information
  mol = entry.molecule 
  # Skipping all unwanted entries
  if len(mol.atoms) < 75:
    too_small += 1
    continue
  if (target+"1") not in mol.formula:
    not_target += 1
    continue
  if any(x in mol.formula for x in tm):
    tm_metal += 1
    continue

  # counting the target atoms
  tcnt = 0
  for at in mol.atoms:
    if at.atomic_symbol == target:
      tcnt += 1
  if tcnt > 1:
    target_cnt += 1
    continue

  # testing the coordination of the metal atoms
  skip=bool(False)
  for at in mol.atoms:
    if at.atomic_symbol == target:
      koor=len(at.neighbours)
      # printf("%s %i\n", at.atomic_symbol, koor)
      if koor==6: skip=bool(False)
      break
  if skip:
    coordination += 1
    continue

  # Analyze survivors
  cnt+=1
  printf("%7i  ", cnt)
  printf("%-10s  ", entry.identifier)
  fname=entry.identifier
  fname=fname.lower()
  fname=resdir+"/"+fname+".mol2"
  printf("%4i  ", len(mol.atoms))
  printf("%s  ", mol.formula)
  printf("\n")
  # print(fname)
  mol_writer = io.MoleculeWriter(fname, format="mol2")
  mol_writer.write(mol)

  # Check break conditions
  if cnt==max_hits: break

printf("\n")
printf("Summary of the search\n")
printf("%10i CCDC entries checked\n", csd_cnt)
printf("%10i organometallic entries\n", organo_metallic)
printf("%10i skipped entries - is_polymeric\n", polymeric)
printf("%10i skipped entries - no 3D struc\n", no_3d_struc)
printf("%10i skipped entries - disorder\n", disorder)
printf("%10i skipped entries - too small\n", too_small)
printf("%10i skipped entries - no target TM ion\n", not_target)
printf("%10i skipped entries - 2nd transition metal\n", tm_metal)
printf("%10i skipped entries - target counter\n", target_cnt)
printf("%10i skipped entries - wrong coordination\n", coordination)

exit()
