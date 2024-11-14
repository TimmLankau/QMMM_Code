################################################################################
# 24.09.2024
# - This Puthon code holds all the functions to be imported
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
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q15_CCDC')
import globals as gb  # short for globals

################################################################################
# general purpose function                                                     #
# fake printf and fprintf as used in C                                         #
################################################################################

def printf(format, *args):
  sys.stdout.write(format % args)
def fprintf(stream, format, *args):
  stream.write(format % args)

################################################################################
# general purpose function                                                     #
# write a list of atoms to stdout                                              #
# input  ml  list with CSD API atoms to be listed on stdout                    #
#        indent  opt. arg. to define the beginning of the output lines         #
# output None                                                                  #
################################################################################

def my_atom_list(ml, indent=""):
  for at in ml:
    printf("%s%3i  ",  indent, at.index)
    printf("%2s  ",    at.atomic_symbol)
    printf("%10.6f  ", at.coordinates.x)
    printf("%10.6f  ", at.coordinates.y)
    printf("%10.6f\n", at.coordinates.z)
  return()

################################################################################
# general purpose function                                                     #
# writing the xyz file from a list of atoms                                    #
# input  alist  list with CSD API atoms containing the geometry information    #
#        name   file name (with suffix)                                        #
#        remark 2nd line (comment) in the xyz file                             #
# output None                                                                  #
################################################################################

def list_xyz(alist, name, remark):
  out = open(name, "w")
  fprintf(out, "%s\n", len(alist))
  fprintf(out, "%s\n", remark)
  for n in range(len(alist)):
    fprintf(out, "%-2s  ", alist[n].atomic_symbol)
    fprintf(out, "%10.6f  ", alist[n].coordinates.x)
    fprintf(out, "%10.6f  ", alist[n].coordinates.y)
    fprintf(out, "%10.6f\n", alist[n].coordinates.z)
  out.close()
  return()

################################################################################
# general purpose function                                                     #
# create a custom label for an atom                                            #
# input  at     atom as defined by the CSD API                                 #
# output mlabel custom label for the atom                                      #
################################################################################

def my_label(at):
  labelstrg = "%s-%i" % (at.atomic_symbol, at.index)
  return(labelstrg)

################################################################################
# general purpose function                                                     #
# create a string to label a triad of atoms                                    #
# input  trip a list of 3 atoms (no sanity control)                            #
# output mstr string with the label for the triad                              #
################################################################################

def triad_string(trip):
  mstr  = "%s%i" % (trip[0].atomic_symbol, trip[0].index)
  mstr += "-"
  mstr += "%s%i" % (trip[1].atomic_symbol, trip[1].index)
  mstr += "-"
  mstr += "%s%i" % (trip[2].atomic_symbol, trip[2].index)
  return(mstr)

################################################################################
# functions to make it easier with the API                                     #
# these two come together to sort and clean up (unique elements) lists         #
# my_sort_by_index - sort an atom list by the index of the atom                #
#   the function to process the search key in embedded in the function         #
#   input  ml  list to be sorted                                               #
#   output ml  returns the sorted list                                         #
# my_clean_list - making the list elements unique and sort by index            #
#   input  ml  list to be cleaned up                                           #
#   output ml  returns the cleaned up list                                     #
################################################################################

# sorting
def my_sort_by_index(ml):
  def my_key(ma):
    return(ma.index)
  ml.sort(key=my_key)
  return(ml)

# clean up
def my_clean_list(ml):
  ml=list(set(ml))     # remove dublicates
  ml=my_sort_by_index(ml) # sort the atom list by atom index
  return(ml)

################################################################################
# functions to make it easier with the API                                     #
# function in_mol_bonds to test if two atoms are joined by registered bond     #
# input  At1       index of the first atom                                     #
#        At2       index of the second atom                                    #
# output -1        the bond is not MyMol.bonds                                 #
#        integer   index number in the bond in the list                        #
################################################################################

def in_mol_bonds(At1, At2):
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  SearchResult = -1  # default value
  for mb in gb.mol.bonds:
    if bool(False):
      printf("  %3i", gb.mol.bonds.index(mb))
      MyLabel = "%s-%i" % (mb.atoms[0].atomic_symbol, mb.atoms[0].index)
      printf("  %6s", MyLabel)
      MyLabel = "%s-%i" % (mb.atoms[1].atomic_symbol, mb.atoms[1].index)
      printf("  %6s\n", MyLabel)
    if mb.atoms[0].index in [At1, At2] and mb.atoms[1].index in [At1, At2]:
      SearchResult = gb.mol.bonds.index(mb)
  return(SearchResult)

################################################################################
# functions to make it easier with the API                                     #
# calculate distance between atoms, no bond needed                             #
# input  at1  1st CSD atom                                                     #
#        at2  2nd CSD atom                                                     #
# output dist float with the distance between the atoms                        #
#             (The unit depends on the unit of the cartesian coordinates)      #
################################################################################

def atom_dist(at1, at2):
  dist = 0.0
  dist += (at1.coordinates.x-at2.coordinates.x)**2
  dist += (at1.coordinates.y-at2.coordinates.y)**2
  dist += (at1.coordinates.z-at2.coordinates.z)**2
  dist = math.sqrt(dist)
  return(dist)

################################################################################
# functions to to help with the QM/MM separation                               #
# can theese tree atoms be a unit of conjugated chain?                         #
# this code is looking only for standard 2nd row bond types                    #
# A doublebond folled by a single bond is not enough to detect a unit, because #
# all atoms need to  have p-orbitals ready to engage in pi-bonding.            #
# input  mymol molecule containing the atoms                                   #
#        at1   1st CSD atom                                                    #
#        at2   2nd CSD atom                                                    #
#        at3   3rd CSD atom                                                    #
# output utype interger to indicate the type of 3-atom unit                    #
#        -3  triple-single                                                     #
#        -2  double-single                                                     #
#         0  not a viable candidate (default)                                  #
#         2  single-double                                                     #
#         3  single-triple                                                     #
# bond type definitions form                                                   #
# downloads.ccdc.camd.ac.uk/documentation/API/descriptive_docs/molecule.html   #
################################################################################

def twisted(dihedral, threshold):
  ang = abs(dihedral)
  if (ang >= threshold) and (ang <= 180.0-threshold):
    return(bool(True))
  else:
    return(bool(False))

def con_unit(at1, at2, at3):
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # No H atoms allowed
  if 1 in [at1.atomic_number, at2.atomic_number, at3.atomic_number]:
    # print("No H atoms allowd")
    return(0)
  if gb.VerboseFlag==2:
    printf("Testing triad %s\n", triad_string([at1, at2, at3]))
  if bool(True):
    # Use the number of neighbours to estimate the hybridization state
    # Each atom should have less then 4 direct neighbours
    if len(at1.neighbours)>3 or len(at1.neighbours)<1: return(0)  # not a candidate
    if len(at2.neighbours)>3 or len(at2.neighbours)<1: return(0)  # not a candidate
    if len(at3.neighbours)>3 or len(at3.neighbours)<1: return(0)  # not a candidate
  else:
    # I try to read the OpenBabel hybridization information 
    # 1 for sp, 2 for sp2, 3 for sp3, 4 for sq. planar, 5 for trig. bipy, 6 for octahedral
    # https://openbabel.org/api/3.0/classOpenBabel_1_1OBAtom.shtml
    printf("  Hyb")
    obatom=gb.obmol.GetAtom(at1.index+1)
    printf("  %i", obatom.GetHyb())
    if obatom.GetHyb() not in [1, 2]:
      return(0)
    obatom=gb.obmol.GetAtom(at2.index+1)
    printf("  %i", obatom.GetHyb())
    if obatom.GetHyb() not in [1, 2]:
      return(0)
    obatom=gb.obmol.GetAtom(at3.index+1)
    printf("  %i\n", obatom.GetHyb())
    if obatom.GetHyb() not in [1, 2]:
      return(0)
  # conjugated chains can start or end at an aromatic rind, but should not be
  # part of an aromatic ring
  obatom=gb.obmol.GetAtom(at2.index+1)
  if obatom.IsAromatic(): return(0)  # not a candidate
  # get the indices of the bonds linking the atoms and check whether the atoms
  # are bonded
  utype = 0 # set default, not a candidate
  ind1 = in_mol_bonds(at1.index, at2.index)
  ind2 = in_mol_bonds(at2.index, at3.index)
  if gb.VerboseFlag==2:
    printf("  ind1: %i, ind2 %i\n", ind1, ind2)
  if ind1 == -1: return(0) # no bond between at1 and at2
  if ind2 == -1: return(0) # no bond between at2 and at3
  # Reply by bond type
  # How likely is it that I catch bond types 7 and 9 ?
  if gb.VerboseFlag==2:
    printf("  %s", my_label(at1))
    printf(" %s", gb.mol.bonds[ind1].bond_type)
    printf(" %s  ", my_label(at2))
    printf(" %s", gb.mol.bonds[ind2].bond_type)
    printf(" %s", my_label(at3))
  if gb.mol.bonds[ind1].bond_type == 3 and gb.mol.bonds[ind2].bond_type == 1:
    if gb.VerboseFlag==2: printf("\n")
    return(-3)
  if gb.mol.bonds[ind1].bond_type == 2 and gb.mol.bonds[ind2].bond_type == 1:
    if len(at3.neighbours)>0:
      if gb.VerboseFlag==2: printf("  %2i", len(at3.neighbours))
      for pat in at3.neighbours:
        if pat !=  at2: break
      dihedral = gb.obmol.GetTorsion(at1.index+1, at2.index+1, at3.index+1, pat.index+1)
      if gb.VerboseFlag==2:
        printf("  %s", my_label(pat))
        printf("  %f", dihedral)
        if twisted(dihedral, 30.0): printf("  twisted")
        else: printf("  conjugated")
    if gb.VerboseFlag==2: printf("\n")
    if twisted(dihedral, 30.0): return(0)
    else: return(-2)
  if gb.mol.bonds[ind1].bond_type == 1 and gb.mol.bonds[ind2].bond_type == 3:
    if gb.VerboseFlag==2: printf("\n")
    return(3)
  if gb.mol.bonds[ind1].bond_type == 1 and gb.mol.bonds[ind2].bond_type == 2:
    if len(at1.neighbours)>0:
      if gb.VerboseFlag==2: printf("  %2i", len(at1.neighbours))
      for pat in at1.neighbours:
        if pat !=  at2: break
      dihedral = gb.obmol.GetTorsion(pat.index+1, at1.index+1, at2.index+1, at3.index+1)
      if gb.VerboseFlag==2:
        printf("  %s", my_label(pat))
        printf("  %f", dihedral)
        if twisted(dihedral, 30.0): printf("  twisted")
        else: printf("  conjugated")
    if gb.VerboseFlag==2: printf("\n")
    if twisted(dihedral, 30.0): return(0)
    else: return(2)
  return(utype)

################################################################################
# functions to to help with the QM/MM separation                               #
# test wether a conjugated 3-atom unit starts at the given atom                #
# The search does not detct all possible units. The serach stops after the one #
# valid unit has been found.                                                   #
# input  at1        1st atom (CDC atoms type) of the triad                     #
# output con_level  integer describing the bond pattern of the unit            #
#        triad      list with the 3 atoms of the conjugated unit               #
################################################################################

def check_conjugation(at1):
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual code
  triad = []
  for at2 in at1.neighbours:
    # if at2 in gb.qm:
    #   continue
    for at3 in at2.neighbours:
      if at3 in gb.qm:
        continue
      if at1 == at3:
        continue
      con_level = con_unit(at1, at2, at3)
      if con_level != 0:
        triad += [at1, at2, at3]
        return(con_level, triad)
  return(0, triad)

################################################################################
# functions for the compartimensation of the  QM/MM separation                 #
# read input file                                                              #
# global VerboseFlag  controll the output to stdout                            #
# input  inp_name     name of the input file with the path to it               #
# output mol          CSD molecule object                                      #
#        obmol        OpenBabel molecule object                                #
################################################################################

def read_input_file(inp_name):

  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()

  # split input names into useful token
  root, ext = os.path.splitext(inp_name)
  path, file = os.path.split(inp_name)
  path += '/'
  base = root.replace(path, '')
  if gb.VerboseFlag == 2:
    print("ext  ", ext)
    print("path ", path)
    print("base ", base)

  # check if the file exists
  if not os.path.isfile(inp_name):
    printf("file %s not found\n", inp_name)
    exit()

  # handle the input file by extention
  if gb.VerboseFlag > 0: printf("Process input file %s\n", inp_name)
  input_type = ext.lower()
  if input_type not in ['.xyz', '.mol2']:
    printf("Wrong input file type\n")
    exit()

  # handel xyz file by converting them to mol2
  # If a mol2 file with same name is available, I will read the available 
  # file. If no substitute is available, the xyz file will be converted
  # into a mol2 file. Finally, the name of the input file will be updated
  # so that the input can be processed in the standard way.
  if input_type==".xyz":
    if gb.VerboseFlag > 0: printf("  process xyz input file (depreciated)\n")
    mol2_inp_name = inp_name.replace(ext, ".mol2")
    if os.path.isfile(mol2_inp_name):
      if gb.VerboseFlag > 0:
        printf("  %s (easier to use file format) found\n", mol2_inp_name)
      # update input information
      inp_name = mol2_inp_name
      ext = ".mol2"
      input_type = ".mol2"
    else:
      if gb.VerboseFlag > 0:
        printf("  use Python module OpenBabel to convert %s%s\n", base, ext)
      obConversion = ob.OBConversion()
      obConversion.SetInAndOutFormats("xyz", "mol2")
      gb.obmol = ob.OBMol()
      obConversion.ReadFile(gb.obmol, inp_name)
      if gb.VerboseFlag > 0:
        printf("  write converted file %s\n", mol2_inp_name)
      obConversion.WriteFile(gb.obmol, mol2_inp_name)
      if gb.VerboseFlag > 0:
        printf("  switch to %s\n", mol2_inp_name)
      if gb.VerboseFlag == 2:
        printf("  delete ObenBabel object (rebuild later)\n")
      del obConversion, gb.obmol
      # update input information
      inp_name = mol2_inp_name
      ext = ".mol2"
      input_type = ".mol2"

  # process the preferred mol2 input file
  # This part is in a conditional block so that the number of input
  # file types can be extended later.
  if input_type==".mol2":
    if gb.VerboseFlag > 0:
      printf("  process mol2 input file (preferred)\n")
      printf("  read input file directly into a CSD object\n")
    mol_reader = io.MoleculeReader(inp_name)
    gb.mol = mol_reader[0]
    if gb.VerboseFlag > 0:
      printf("  standarize bonds\n")
    gb.mol.assign_bond_types()
    gb.mol.standardise_aromatic_bonds()
    gb.mol.standardise_delocalised_bonds()

  # At this point I should have a CSD molecule object.
  # build auxiliary OpenBabel molecule object by reading the mol2 file
  if gb.VerboseFlag > 0:
    printf("  build supplementary OpenBabel object from mol2 file\n")
  obConversion = ob.OBConversion()
  obConversion.SetInAndOutFormats("mol2", "mol2")
  gb.obmol = ob.OBMol()
  obConversion.ReadFile(gb.obmol, inp_name)
  if gb.VerboseFlag > 0:
    printf("  compare the molecular geometries of both objects\n")
  csdanz = len(gb.mol.atoms)
  obanz  = gb.obmol.NumAtoms()
  if gb.VerboseFlag==2:
    printf("  Number of atoms: CSD %i   OB %i\n", csdanz, obanz)
  if csdanz != obanz:
    print("The number of atoms in both molecule objects (CSD:  %i, OB: %i) don't match",
          csdanz, obanz)
    exit()
  else:
    if gb.VerboseFlag>0:
      printf("    the number of atoms %i in both objects match\n", csdanz)
  if gb.VerboseFlag>0:
    printf("    compare geometries line by line ...\n")
  for n in range(csdanz):
    obatom = gb.obmol.GetAtom(n+1)
    csdnum = gb.mol.atoms[n].atomic_number
    obnum  = obatom.GetAtomicNum()
    csdX = gb.mol.atoms[n].coordinates.x
    csdY = gb.mol.atoms[n].coordinates.y
    csdZ = gb.mol.atoms[n].coordinates.z
    obX  = obatom.GetX()
    obY  = obatom.GetY()
    obZ  = obatom.GetZ()
    deltaX = abs(csdX-obX)
    deltaY = abs(csdY-obY)
    deltaZ = abs(csdZ-obZ)
    if gb.VerboseFlag==2:
      printf("    %3i", n)
      printf(  "  %3i", csdnum)
      printf(  "  %3i", obnum)
      printf(  "  %10.6f", deltaX)
      printf(  "  %10.6f", deltaY)
      printf(  "  %10.6f", deltaZ)
      printf("\n")
    if csdnum != obnum:
      printf("The atomic numbers (CSD %i, OB %i) for entry #%i don't match\n",
            csdnum, obnum, n)
      exit()
    if deltaX >= 0.0001:
      printf("The X coordinates (CSD %.6f, OB %.6f) for entry #%i don't match\n",
            csdX, obX, n)
      exit()
    if deltaY >= 0.0001:
      printf("The Y coordinates (CSD %.6f, OB %.6f) for entry #%i don't match\n",
            csdY, obY, n)
      exit()
    if deltaZ >= 0.0001:
      printf("The Z coordinates (CSD %.6f, OB %.6f) for entry #%i don't match\n",
            csdZ, obZ, n)
      exit()
  if gb.VerboseFlag>0: printf("    the geometries of both objects match\n")
  # return both molecule objects
  return(gb.mol, gb.obmol)

################################################################################
# functions for the compartimensation of the  QM/MM separation                 #
# read input file                                                              #
# global VerboseFlag  controll the output to stdout                            #
# input  inp_name     name of the input file with the path to it               #
# output mol          CSD molecule object                                      #
#        obmol        OpenBabel molecule object                                #
################################################################################

def write_output_file(out_name):
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # create combined output file
  combi=[]       # list holding the complete molecule
  link_list = [] # list hoding the atom pairs of broken bonds
  # building the link list
  for qat in gb.qm:
    for mat in qat.neighbours:
      if mat in gb.mm:
        pair=[]
        pair.append(qat)
        pair.append(mat)
        link_list.append(pair)
  combi = gb.qm + gb.mm
  tat = len(combi)
  qat = len(gb.qm)
  mat = len(gb.mm)
  qel = 0
  for at in gb.qm: qel += at.atomic_number
  mel = 0
  for at in gb.mm: mel += at.atomic_number
  tel = qel + mel
  li = len(link_list)
  stats = [tat, qat, mat, tel, qel, mel, li]
  if gb.VerboseFlag > 0:
    printf("Combined output:\n")
    printf("  All atoms     %4i\n", tat)
    printf("  QM atoms      %4i (%4.1f%%)\n", qat, 100.0*qat/tat)
    printf("  MM atoms      %4i (%4.1f%%)\n", mat, 100.0*mat/tat)
    printf("  All electrons %4i (based on Znuc)\n", tel)
    printf("  QM electrons  %4i (%4.1f%%)\n", qel, 100.0*qel/tel)
    printf("  MM electrons  %4i (%4.1f%%)\n", mel, 100.0*mel/tel)
    printf("  Links         %4i\n", len(link_list))
  comment = "QM: %i at, MM: %i at, Links: %i" % (len(gb.qm), len(gb.mm), len(link_list))
  list_xyz(combi, out_name, comment)
  # append xyz file with the link list
  out = open(out_name, "a")
  fprintf(out, "\n")
  fprintf(out, "QM/MM links (counting starts at zero)\n")
  for pair in link_list:
    fprintf(out, "%3i  ", combi.index(pair[0]))
    fprintf(out, "%3i\n", combi.index(pair[1]))
  out.close()
  return(stats)

################################################################################
# functions for the compartimensatiom of the  QM/MM separation                 #
# this function summarizes an individual step of the search                    #
#   global var VerboseFlag will be used to controll the output                 #
# input  cnt step counter                                                      #
#        nlist new atom list                                                   #
#        flist atom list send to file                                          #
#        mytxt string for the description                                      #
################################################################################

def summarize_step(cnt, nlist, flist, mytxt):
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # output to stdout
  if gb.VerboseFlag == 0:
    return()
  printf("Step %02i - %s\n", cnt, mytxt)
  printf("  %i new atoms found\n", len(nlist))
  if gb.VerboseFlag == 2 and len(nlist)>0:
    my_atom_list(nlist, "  ")
  # write atom list to disk
  fname = "stp%02i.xyz" % cnt
  if len(flist)>0:
    printf("  Write %s\n", fname)
    list_xyz(flist, fname, mytxt)
  else:
    printf("  Skip  %s\n", fname)
  return()

################################################################################
# functions to to help with the QM/MM separation                               #
# find metal center                                                            #
# input   at_num       atomic number of the lightest metal atom                #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
################################################################################

def metal_center(at_num):
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  new_atoms = []  # empty list for serach results
  for at in gb.mol.atoms:
    if at.is_metal and at.atomic_number >= at_num:
      new_atoms.append(at)
  new_atoms = my_clean_list(new_atoms)
  gb.qm += new_atoms
  gb.qm = my_clean_list(gb.qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find close neighbors                                                         #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def close_neighbors():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  new_atoms = []  # empty list for serach results
  for mc in gb.dnts['cent']:
    for at in mc.neighbours:
      new_atoms.append(at)
      gb.qm.append(at)
  new_atoms = my_clean_list(new_atoms)
  gb.qm = my_clean_list(gb.qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# distant neighbors                                                            #
# input   max          maximum distance in Angs                                #
#         fac          scaling factor for the sum of VdW radii                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def distant_neighbors(max, vdw_fac):
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  new_atoms = []  # empty list for serach results
  bnd_cnt = 0      # counter for new bonds
  vdw_sum = 0.0    # sum of the VdW radii of a long distance contact

  if gb.VerboseFlag == 2:
    printf("Distant neighbour search\n")
    printf("  <= %3.1f A and <= %.0f%% of the VdW radius\n", max, vdw_fac*100)
  for a1 in gb.dnts['cent']: # loop over all center atoms
    for a2 in gb.mol.atoms: # loop over all atoms in the molecule
      if a1 == a2 : continue
      dist = atom_dist(a1, a2)
      vdw_sum =  (1.0 * a1.vdw_radius)  # The type of the vdw_property appears to be
      vdw_sum += (1.0 * a2.vdw_radius)  # context dependent. I force it to a float.
      vdw_sum *= vdw_fac                # Apply cut-off factor
      if dist < max and a2 not in gb.qm:
        if gb.VerboseFlag == 2:
          printf("  %s%i",    a1.atomic_symbol, a1.index)
          printf( "-%s%i",    a2.atomic_symbol, a2.index)
          printf("  %6.4f",   dist)
          printf("  %6.4f",   vdw_sum)
          if dist<=vdw_sum: printf("  Passed\n")
          else: printf("  too long\n")
        if dist<=vdw_sum:
          new_atoms.append(a2)
          # Adding new bonds messes with the CDS assignment of rings. The new 
          # bonds can create small 3-member rings looping back to the metal 
          # atoms. To reactivate the option for new bonds change False to True
          if in_mol_bonds(a1, a2)==-1 and bool(True):
            gb.mol.add_bond(1, a1, a2)
            bnd_cnt += 1
  new_atoms = my_clean_list(new_atoms)
  gb.qm += new_atoms
  gb.qm = my_clean_list(gb.qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find the rings containing neighbours to the metal center                     #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def neighbour_rings():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  new_atoms = []  # empty list for serach results
  # create an intermediate test list for atoms to be tested
  itl = gb.dnts['close'] + gb.dnts['dist']
  itl = my_clean_list(itl)
  if gb.VerboseFlag==2:
    printf("Rings containing neighbours atoms\n")
    printf("    atom  r  s\n")
  for ne in itl:
    if gb.VerboseFlag==2:
      printf("  %6s", my_label(ne))
      if len(ne.rings)>0:
        printf("  %i", len(ne.rings[0]))
      else:
        printf("  -")
      if ne.is_spiro:
        printf("  %i\n", len(ne.rings[1]))
      else:
        printf("  -\n")
    if len(ne.rings)>0 and len(ne.rings[0])<=10 : # Is the atom member of a ring?
      for at in ne.rings[0].atoms: # focus on the smallest ring
        if at not in gb.qm:
          new_atoms.append(at)
    if ne.is_spiro and len(ne.rings)>1: # Is the atom member a spiro atom?
      for at in ne.rings[1].atoms: # focus on the smallest ring
        if at not in gb.qm:
          new_atoms.append(at)
  new_atoms = my_clean_list(new_atoms)
  gb.qm += new_atoms
  gb.qm  = my_clean_list(gb.qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# detect the inner aromatic rings                                              #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def inner_aromatic_rings():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  new_atoms = []  # empty list for search results
  # create an intermediate test list for atoms to be tested
  itl = gb.dnts['close'] + gb.dnts['dist'] + gb.dnts['nrings']
  itl = my_clean_list(itl)
  # find new MM aromatic ring atoms
  for sa in itl:
    for ri in sa.rings: # test all rings
      if ri.is_aromatic:
        for at in ri.atoms: # add all aromatic ring atoms not in qm
          if at not in gb.qm:
            new_atoms.append(at)
  new_atoms = my_clean_list(new_atoms)
  gb.qm += new_atoms
  gb.qm = my_clean_list(gb.qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find chains of conjugated 3-atom units                                       #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
################################################################################

def conjugated_chains():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  new_atoms = []
  while bool(True): # infinite loop to detect chains of units
    cl1 = 0
    trip = []
    gotcha = 0
    for at1 in my_clean_list(gb.qm): # Loop over all QM atoms
      if gb.VerboseFlag==2:
        printf("Testing %s-%i\n", at1.atomic_symbol, at1.index)
      cl1, trip = check_conjugation(at1)
      if cl1 != 0:
        if gb.VerboseFlag==2:
          printf("  %-10s  %2i\n", triad_string(trip), cl1)
        for nat in trip:
          if nat not in gb.qm:
            gotcha = 1
            new_atoms.append(nat)
    new_atoms = my_clean_list(new_atoms)
    gb.qm += new_atoms
    gb.qm = my_clean_list(gb.qm)
    if gotcha == 0:
      break
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# outward growth of aromatic rings                                             #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def grow_aromatic_rings():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  new_atoms = []  # empty list for search results
  it_cnt = 0      # counter for iterative rounds          
  # create an intermediate test list for atoms to be tested
  # itl = dnts['close']+dnts['dist']+dnts['nrings']+dnts['cArings']+dnts['conju']
  itl = list(set(gb.qm) - set(gb.dnts['cent']))
  itl = my_clean_list(itl)
  if gb.VerboseFlag == 2:
    printf("Details for the iterative search for aromatic rings\n")
  while bool(True):   # infinite loop for the iterative search
    new_at = []       # atoms found by each loop of the iterative search
    it_cnt += 1       # increase counter
    # check whether the new atoms are part of an aromatic ring and that these atoms
    # are not in the QM core
    for at in itl:
      for ri in at.rings:
        if ri.is_aromatic:
          for tat in ri.atoms:
            if tat not in gb.qm:
              new_at.append(tat)
    if len(new_at) == 0:
      break
    else:
      if gb.VerboseFlag == 2:
        printf("  %3i atoms found in step %2i\n", len(new_at), it_cnt)
      new_atoms += new_at
      new_atoms = my_clean_list(new_atoms)
      itl = my_clean_list(itl + new_atoms)
      new_at = []
      gb.qm += new_atoms
      gb.qm = my_clean_list(gb.qm)
  if gb.VerboseFlag == 2:
    printf("  %i iterative search rounds needed\n", it_cnt-1)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# looking for double bonds exo to the quantum core                             #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def double_check():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  results = []
  while bool(True):
    new_atoms = []  # empty list for search results
    if gb.VerboseFlag==2:
      print("  qm      mm      bond    action")
    for qa in gb.qm: # loop over quantum core
      for na in qa.neighbours:
        if na in gb.qm:
          continue
        if gb.VerboseFlag==2:
          printf("  %-6s", my_label(qa))
          printf("  %-6s", my_label(na))
        bidx = in_mol_bonds(qa.index, na.index)
        btype = gb.mol.bonds[bidx].bond_type
        if gb.VerboseFlag==2:
          printf("  %-6s", gb.mol.bonds[bidx].bond_type)
        if gb.mol.bonds[bidx].bond_type != 1:
          new_atoms.append(na)
        if gb.VerboseFlag==2:
          if gb.mol.bonds[bidx].bond_type != 1:
            printf("  add %s", my_label(na))
          printf("\n")
    if len(new_atoms) == 0:
      break
    results += my_clean_list(new_atoms)
    gb.qm += new_atoms
    gb.qm = my_clean_list(gb.qm)
  return(results)

################################################################################
# functions to to help with the QM/MM separation                               #
# looking for a single bond connected to an atom with a lone pair              #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def my_trig_plan(at):
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  if len(at.neighbours) != 3:
    return(bool(False))
  n1 = at.neighbours[0]
  n2 = at.neighbours[1]
  n3 = at.neighbours[2]
  dihedral = gb.obmol.GetTorsion(n1.index+1, n2.index+1, at.index+1, n3.index+1)
  if (dihedral>= -10.0) and (dihedral<=  10.0):
    return(bool(True))
  if (dihedral>=-190.0) and (dihedral<=-170.0):
    return(bool(True))
  if (dihedral>= 170.0) and (dihedral<= 190.0):
    return(bool(True))
  return(bool(False))
  
def lp_num(at):
  # dictionary with the valence electons of main group elements
  ve_dic= { 1: 1,                                            2: 2,
            3: 1,  4: 2,  5: 3,  6: 4,  7: 5,  8: 6,  9: 7, 10: 8,
           11: 1, 12: 2, 13: 3, 14: 4, 15: 5, 16: 6, 17: 7, 18: 8,
           19: 1, 20: 2, 31: 3, 32: 4, 33: 5, 34: 6, 35: 7, 36: 8,
           37: 1, 38: 2, 49: 3, 50: 4, 51: 5, 52: 6, 53: 7, 54: 8,
           55: 1, 56: 2, 81: 3, 82: 4, 83: 5}
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual code
  if at.atomic_number not in ve_dic:
    print("Function 'lp_mum' called for an atom not defined in 've_dic'.")
    exit()
  lp = 0
  elec = ve_dic[at.atomic_number]
  formal = at.formal_charge
  bond=0
  for na in at.neighbours:
    bidx = in_mol_bonds(at.index, na.index)
    btype = gb.mol.bonds[bidx].bond_type
    # print(na, btype)
    if btype=="Single":    bond = bond+1
    if btype=="Double":    bond = bond+2
    if btype=="Triple":    bond = bond+3
    if btype=="Quadruple": bond = bond+4
    if btype=="Aromatic":  bond = bond+1.5
  lp = (elec-formal-bond)/2.0  
  return(lp)

def lp_check():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  results = []
  new_atoms = []
  for qa in gb.qm:
    for na in qa.neighbours:
      dihedral = 999
      if na in gb.qm: continue
      if na.atomic_number==1: continue
      bidx = in_mol_bonds(qa.index, na.index)
      btype = gb.mol.bonds[bidx].bond_type
      if btype == "Single":
        lp = lp_num(na)
        if my_trig_plan(qa) and my_trig_plan(na):
          for qn in qa.neighbours:
            if qn != na: break
          for nn in na.neighbours:
            if nn != qa: break
          dihedral = gb.obmol.GetTorsion(qn.index+1, qa.index+1, na.index+1, nn.index+1)
        if gb.VerboseFlag==2:
          printf("  %-4s", my_label(qa))
          printf("  %s", btype)
          printf("  %-5s", my_label(na))
          printf("  %i", lp)
          if my_trig_plan(na): printf("  plan")
          else: printf("  not")
          if my_trig_plan(qa) and my_trig_plan(na):
            printf("  %.1f", dihedral)
            if twisted(dihedral, 30.0): printf("  twisted")
            else: printf("  planar ")
          printf("\n")
        if lp>0:
          if dihedral != 999:
            if twisted(dihedral, 30.0):
              continue
            else:
              new_atoms.append(na)
          else:
            new_atoms.append(na)
  new_atoms = my_clean_list(new_atoms)
  gb.qm += new_atoms
  gb.qm += my_clean_list(gb.qm)
  results += new_atoms
  return(results)

################################################################################
# functions to to help with the QM/MM separation                               #
# find one atom bridges                                                        #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
################################################################################

def one_atom_links():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  new_atoms = []
  if gb.VerboseFlag == 2:
    printf("List of all found 1-atom links\n")
  for sa in gb.qm:
    for ea in gb.qm:
      if sa.index == ea.index:
        continue
      for la in sa.neighbours:
        if la in ea.neighbours and la not in gb.qm:
          new_atoms.append(la)
          if gb.VerboseFlag == 2:
            printf("  %s%i-", sa.atomic_symbol, sa.index)
            printf("%s%i-", la.atomic_symbol, la.index)
            printf("%s%i\n", ea.atomic_symbol, ea.index)
  new_atoms = my_clean_list(new_atoms)
  gb.qm += new_atoms
  gb.qm = my_clean_list(gb.qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find two atom bridges                                                        #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
################################################################################

def two_atom_links():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  new_atoms = []
  if gb.VerboseFlag == 2:
    printf("List of all found 2-atom links\n")
  for b0 in gb.qm:
    for b1 in b0.neighbours:
      for b2 in b1.neighbours:
        for b3 in b2.neighbours:
          if gb.VerboseFlag == 2:
            printf("  %s%i-", b0.atomic_symbol, b0.index)
            printf("%s%i-",   b1.atomic_symbol, b1.index)
            printf("%s%i-",   b2.atomic_symbol, b2.index)
            printf("%s%i  ",  b3.atomic_symbol, b3.index)
            if (b1 not in gb.qm and b2 not in gb.qm and b3 in gb.qm):
              printf("Add\n")
            else:
              printf("Del\n")
          if (b1 not in gb.qm and b2 not in gb.qm and b3 in gb.qm):
            new_atoms.append(b1)
            new_atoms.append(b2)
  new_atoms = my_clean_list(new_atoms)
  gb.qm += new_atoms
  gb.qm = my_clean_list(gb.qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find H atoms attached to the QM core                                         #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
################################################################################

def lone_H_atoms():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  new_atoms = []
  for at in gb.qm:
    for na in at.neighbours:
      if na.atomic_number == 1 and na not in gb.qm:
        new_atoms.append(na)
  new_atoms = my_clean_list(new_atoms)
  gb.qm += new_atoms
  gb.qm = my_clean_list(gb.qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# build the MM layer by subtracting the QM atoms from the molecule             #
# input   none                                                                 #
# output  new_atoms    list of the MM atoms                                    #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         mm           list with all atoms in the mechanical layer             #
################################################################################

def create_mm_layer():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()
  # the actual search
  gb.mm = []
  for at in gb.mol.atoms:
    if at in gb.qm:
      pass
    else:
      gb.mm.append(at)
  gb.mm = my_clean_list(gb.mm)
  return(gb.mm)

################################################################################
################################################################################

def std_qmmm_sep():
  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()

  # test if both molecular objects can be read
  if bool(False):
    print("CSD:", len(gb.mol.atoms), "at")
    print("OB :", gb.obmol.NumAtoms(), "at")

  # initialize output
  if gb.VerboseFlag > 0:
    printf("\n")
    printf("Start the actual QM/MM partioning process\n")

  # start with a clean slate
  gb.qm = []
  gb.mm = []
  gb.dnts = {}

  # find metal centers
  stp_cnt = 1
  gb.dnts['cent'] = metal_center(19) # start with potassium
  gb.dnts.update()
  summarize_step(stp_cnt, gb.dnts['cent'], gb.qm, "metal center")

  # find the direct neighbours to the metal center
  stp_cnt += 1
  gb.dnts['close'] = close_neighbors()
  gb.dnts.update()
  summarize_step(stp_cnt, gb.dnts['close'], gb.qm, "direct neighbours")

  # additional search based on distance, use the VdW radii for validity
  stp_cnt += 1
  # are there enough binding ligand atoms (6-fold coordination preferred)
  if len(gb.dnts['close']) < 4:
    if gb.VerboseFlag>0:
      printf("{len(dnts['close'])<4}, do distant neightbour search\n")
    gb.dnts['dist'] = distant_neighbors(3.0, 0.8)
  else:
    if gb.VerboseFlag>0:
      printf("{len(dnts['close'])>=4}, skip distant neightbour search\n")
    gb.dnts['dist'] = []
  gb.dnts.update()
  summarize_step(stp_cnt, gb.dnts['dist'], gb.qm, "distant neighbours")

  # add the atoms of rings containing direct neighbours
  stp_cnt += 1
  gb.dnts['nrings'] = neighbour_rings()
  gb.dnts.update()
  summarize_step(stp_cnt, gb.dnts['nrings'], gb.qm, "rings containing neighbours atoms")

  # first round of aromatic ring joined to inner rings
  stp_cnt   += 1    # increase QM file counter
  gb.dnts['cArings'] = inner_aromatic_rings()
  gb.dnts.update()
  summarize_step(stp_cnt, gb.dnts['cArings'], gb.qm, "inner aromatic rings")

  # enter the loop for the coninous growth of the QM core
  if gb.VerboseFlag>0: printf("\nEntering search loop\n")
  gb.dnts['conju']   = []
  gb.dnts['gArings'] = []
  gb.dnts['dBonds'] = []
  gb.dnts['lpBonds'] = []
  gb.dnts.update()
  hc = []  # help/intermediate variable for conjugated chains
  ha = []  # help/intermediate variable for aromatic rings
  hd = []  # help/intermediate variable for exo double bonds from the QM core
  hl = []  # help/intermediate variable for sigma bonds to atoms with lone pairs
  cnt = 0
  while bool(True):
    cnt+=1
    num_new_at = 0
    # Test step to check for conjugated chanins
    stp_cnt += 1
    hc = conjugated_chains()
    num_new_at += len(hc)
    gb.dnts['conju'] = my_clean_list(gb.dnts['conju'] + hc)
    gb.dnts.update()
    summarize_step(stp_cnt, hc, gb.qm, "conjugated chains")
    # grow the number of joined aromatic ring
    stp_cnt   += 1
    ha = grow_aromatic_rings()
    num_new_at += len(ha)
    gb.dnts['gArings'] += ha
    gb.dnts.update()
    summarize_step(stp_cnt, ha, gb.qm, "aromatic rings")
    # check for exo-double bonds
    stp_cnt   += 1
    hd = double_check()
    num_new_at += len(hd)
    gb.dnts['dBonds'] += hd
    gb.dnts.update()
    summarize_step(stp_cnt, hd, gb.qm, "exo double bonds")
    # check for sigma bonds to ligands with lone pairs
    stp_cnt   += 1
    hl = lp_check()
    num_new_at += len(hl)
    gb.dnts['lpBonds'] += hl
    gb.dnts.update()
    summarize_step(stp_cnt, hl, gb.qm, "sigma bonds to Lewis bases")
    # summarize results from the current round of the loop
    if num_new_at==0:
      if gb.VerboseFlag>0: printf("Round %i:  no new atoms\n", cnt)
      break
    else:
      if gb.VerboseFlag>0: printf("Round %i:  %2i new atoms\n\n", cnt, num_new_at)
  if gb.VerboseFlag>0: printf("Leaving search loop\n\n")
  del hc, ha, hd, cnt
  # leave the loop for the coninous growth of the QM core

  # looking for 1 atom links between QM atoms
  stp_cnt   += 1
  gb.dnts['aLink'] = one_atom_links()
  gb.dnts.update()
  summarize_step(stp_cnt, gb.dnts['aLink'], gb.qm, "1 atom links")

  # looking for 2 atom links between the rings
  stp_cnt   += 1
  gb.dnts['aaLink'] = two_atom_links()
  gb.dnts.update()
  summarize_step(stp_cnt, gb.dnts['aaLink'], gb.qm, "2 atom links")

  # adding single hydrogen atoms
  stp_cnt   += 1
  gb.dnts['hatom'] = lone_H_atoms()
  gb.dnts.update()
  summarize_step(stp_cnt, gb.dnts['hatom'], gb.qm, "H atoms")

  # create mm layer
  stp_cnt   += 1
  gb.dnts['mmatom'] = create_mm_layer()
  gb.dnts.update()
  summarize_step(stp_cnt, gb.dnts['mmatom'], gb.mm, "only MM atoms")

  return()