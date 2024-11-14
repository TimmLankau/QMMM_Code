# 23.09.2024
# - picking random files without repetions

import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import random                   # create random number
import time                     # gather timing info

# Import OpenBabel to assist with atom classification
from openbabel import openbabel as ob           # basic OpenBabel

# Import my modules
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q15_CCDC')
import globals as     gb      # short for globals
import qmmm    as     qms     # short for qm/mm separation
from   qmmm    import printf
from   qmmm    import fprintf

# VERY IMPORTANT, initialize globals
gb.init()

# define vital variables
file_dir = "/dicos_ui_home/tlankau/QMMM/Query02/QRes-134592/Results"
stress_anz = 100

def mol2_num_atoms(file_name):
  name=""
  anz=0
  f = open(file_name,"r")
  linelist = f.readlines()
  cnt=0
  for line in linelist:
    if "@<TRIPOS>MOLECULE" in line : break
    cnt += 1
  name = linelist[cnt+1].strip()
  anz = linelist[cnt+2].strip()
  anz = anz.split()
  anz = int(anz[0])
  return(name, anz)  

total_time = 0
job_list = glob.glob(file_dir+'/*.mol2')
for cnt in range(1, stress_anz+1, 1):
  printf("%3i  ", cnt)
  # printf("%4i  ", len(job_list))
  file = random.choice(job_list)
  job_list.remove(file)
  name, numat = mol2_num_atoms(file)
  # printf("%-8s  %3i  ", name, numat)
  printf("%-8s  ", name)

  ##############################################################################
  # Do the QMM/MM separation                                                   #
  ##############################################################################

  # set output verbosity
  gb.VerboseFlag = 0

  if gb.VerboseFlag==0:
    ob.obErrorLog.SetOutputLevel(0)
    
  # read input file
  gb.mol, gb.obmol = qms.read_input_file(file)

  # standard QM/MM separation (deviation from the old code)
  start_time = time.time()
  qms.std_qmmm_sep()
  end_time   = time.time()

  # write output file with the QM/MM separation data
  outfile = file[:-5]
  outfile += ".xyz"
  outfile = outfile.replace(file_dir, ".")
  stats = qms.write_output_file(outfile)

  # printf("(%3i)", numat)
  for number in stats:
    printf("  %4i", number)
  
  if stats[0]!=numat:
    printf("Atom numbes don't match")
    exit()

  ##############################################################################
  # Leave the QMM/MM separation                                                #
  ##############################################################################
   
  time_diff = end_time - start_time
  total_time += time_diff
  printf("  %5.1f", time_diff)
  printf("\n")
printf("Avergae time per QM/MM separation %.1f\n", total_time/stress_anz)  