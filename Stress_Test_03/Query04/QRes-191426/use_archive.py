# 23.10.2024 - Building a data base for the stress tests based on the archive file
# - Starting with query03.py

# load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import math                     # mathematical  functions
import random                   # create random number


# Load modules and prepare to read the CCS data base
import ccdc.search

from ccdc          import io                    # read crystal structures
from ccdc.search   import Search
from ccdc.search   import TextNumericSearch     # search by text or numbers
from ccdc.search   import SubstructureSearch    # search for substructures

from ccdc.molecule import Molecule              # Build a molecule
from ccdc.molecule import Atom                  # Atomic data and properties
from ccdc.molecule import Bond                  # Bond properties

# load rdkit module for a quick chek of the molecule
from rdkit import Chem

# Fake printf and fprintf as used in C
def printf(format, *args):
  sys.stdout.write(format % args)
def fprintf(stream, format, *args):
  stream.write(format % args)

# define the archive file
archive = '/dicos_ui_home/tlankau/QMMM/Query04/OM_non_polymer.csv'

# define the number of hits for the stress test
max_hits = 750

# a list to hold the CCDC identifier
CCDC_ident = []

# list of main group metals
mm =["Li", "Na", " K", "Rb", "Cs",
     "Be", "Mg", "Ca", "Sr", "Ba",
           "Al", "Ga", "In", "Tl",
                 "Ge", "Sn", "Pb",
                       "Sb", "Bi",
                             "Po"]
# list of transtion metals
tm=["Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y" , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg"]
# list of rare earth elements
re = ["Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
# list of radiocative elements
ra = ["Tc", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
# only mono-nuclear Fe, Mn cluster
bad_tm = tm.copy()
bad_tm.remove('Fe')
bad_tm.remove('Mn')

# initialize dictionary for counting stats
cnt = {}
def inc_dict_cnt(key):
    if key not in cnt:
        cnt[key] = 0
    cnt[key] += 1
    return()
def list_dict(kw, iw):
    for key, value in cnt.items():
        fstr = "%-"+str(kw)+"s: %"+str(iw)+"i\n"
        sys.stdout.write(fstr % (key, value))
    return()

# test atom count by rdkit
# rdk_mol = Chem.MolFromSmiles('C')
# rdk_mol = Chem.MolFromSmiles('CCO')
# rdk_mol = Chem.MolFromSmiles('c1ccccc1')
# rdk_mol = Chem.AddHs(rdk_mol)
# print(rdk_mol.GetNumAtoms())
# exit()

# test for the archive csv file
if os.path.isfile(archive):
    print('Archive csv file found, continue')
else:
    print('Archive csv file NOT found, exit code')
    exit()

# open archive abd readit line by line in a gigantic while loop
with open(archive, 'r', encoding='UTF-8') as file:
    for line in file:
        token = line.rstrip().split(',')
        # make sure it's organometallic
        if token[3] != 'True':
            continue
        # make sure it's not a polymer
        if token[4] != 'False':
            continue
        inc_dict_cnt('file')
        SMILES=token[2]
        # use the SMILES string to check for main group metals
        if any(x in SMILES for x in mm):
            inc_dict_cnt('mg_metal')
            continue
        # use the SMILES string to check for rare earth metals
        if any(x in SMILES for x in re):
            inc_dict_cnt('rare_earth')
            continue
        # use the SMILES string to check for radioactive elements
        if any(x in SMILES for x in ra):
            inc_dict_cnt('radio_active')
            continue
        # use the SMILES string to check for the correct transition metal ions
        if any(x in SMILES for x in bad_tm):
            inc_dict_cnt('bad_transition_metal')
            continue
        if ('Fe' in SMILES) and ('Mn' in SMILES):
            inc_dict_cnt('fe_and_mn')
            continue
        if SMILES.count('Fe') > 1:
            inc_dict_cnt('not_fe1')
            continue
        if SMILES.count('Mn') > 1:
            inc_dict_cnt('not_mn1')
            continue
        # identify substructructures in the SMILES string
        if bool(True):
            if '.' in SMILES:
                inc_dict_cnt('smiles_dot')
                continue
        # use rdkit to count the number of atoms
        # The translation of SMILES strings to molecules proved to be unreliable
        # 43% failure rate, because dative bonds are not indicated in the SMILES
        # string with a '>' or a '<'. An undeclared dative bond results in a too
        # large valence count for a ligand atom (typically N and O) and consequently
        # to a failure for the translation and an undefined molecule object.
        if bool(False):
            rdk_mol = Chem.MolFromSmiles(SMILES)
            if rdk_mol == None:
                inc_dict_cnt('bad_SMILES')
                print(SMILES)
                continue
            rdk_mol = Chem.AddHs(rdk_mol)
            anz = rdk_mol.GetNumAtoms()
            if anz<75:
                inc_dict_cnt('too_small')
                continue
        # add the good CCDC identifier to a list
        CCDC_ident.append(token[1])

# stats form the csv file
list_dict(21, 6)
print(f"%-21s: %6i" % ("survivors", len(CCDC_ident)))

# end the data base build using the csv file
# recycle the code from pervious querry (query03.py)

print()
print("Continue with the CCSD data base")

# reset dictionaries for a set of new counter
cnt = {}
cnt['list_len'] = len(CCDC_ident)

# get the directories and 'Results'
basedir = os.path.dirname(os.path.realpath(__file__))
resdir = basedir + "/Results"
printf("%-16s: ", "Base directory")
printf("%s\n",  basedir)
printf("%-16s: ", "Res directory")
printf("%s\n",  resdir)
if os.path.isdir(resdir):
  print("Directory for results exists")
else:
  print("Directory for results does NOT exist -> Create")
  os.makedirs(resdir)

# sort list entries
print("Sort list entries")
CCDC_ident.sort()
# loop the list CCDC_ident
print(f"Pick random entries until I have %i hits" % (max_hits))
csd_reader = io.EntryReader('CSD')

#  file = 
#  job_list.remove(file)

cnt['hits'] = 0
while cnt['hits'] < max_hits:
    inc_dict_cnt('trials')
    CCDC_trial = random.choice(CCDC_ident)
    CCDC_ident.remove(CCDC_trial)
    entry = csd_reader.entry(CCDC_trial)
    # check the entry
    if not entry.has_3d_structure:
        inc_dict_cnt('3d_struc')
        continue
    if entry.has_disorder:
        inc_dict_cnt('disorder')
        continue
    # select based on molecular information
    mol = entry.molecule 
    if len(mol.atoms) < 75:
        inc_dict_cnt('too_small')
        continue
    # testing the coordination of the metal atoms
    skip=bool(True)
    for at in mol.atoms:
        if at.atomic_symbol in ['Fe', 'Mn']:
            koor=len(at.neighbours)
            # printf("%-12s %s %i\n", CCDC_trial, at.atomic_symbol, koor)
            if koor==6: skip=bool(False)
            break
    if skip:
        inc_dict_cnt('coord')
        continue
    print(f"%5i %2s %-12s %4i %i" % (cnt['hits'], at.atomic_symbol, CCDC_trial, len(mol.atoms), koor))
    if at.atomic_symbol == 'Fe':
        inc_dict_cnt('good_Fe')
    if at.atomic_symbol == 'Mn':
        inc_dict_cnt('good_Mn')
    inc_dict_cnt('hits')
    if bool(True):
        fname=entry.identifier
        fname=fname.lower()
        fname=resdir+"/"+fname+".mol2"
        mol_writer = io.MoleculeWriter(fname, format="mol2")
        mol_writer.write(mol)
list_dict(21, 6)
# exit code for good
exit()