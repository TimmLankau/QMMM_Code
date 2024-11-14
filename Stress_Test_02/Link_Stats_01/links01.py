# Import basic Python modules
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import math                     # mathematical  functions
import operator                 # used to sort dictionaries
import pickle                   # save Python objects to file

# Import OpenBabel to assist with atom classification
from openbabel import openbabel as ob           # basic OpenBabel
ob.obErrorLog.SetOutputLevel(0)           # suppress warnings

# Import my modules
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q15_CCDC')
from qmmm import printf
from qmmm import fprintf

# Define paths
res_dir="/dicos_ui_home/tlankau/QMMM/Stress02/Test-137829/"

# Loop over all result files
res_list = glob.glob(res_dir+'/*.xyz')
printf("%4i Result files found\n", len(res_list))

# Get info from the result file
link_dict={} # dictionary for all results
link_cnt=0   # total number of links analyzed
no_links=0   # cluster without links
for res in res_list:
  # print(res)
  obConversion = ob.OBConversion()
  obConversion.SetInAndOutFormats("xyz", "xyz")
  obmol = ob.OBMol()
  obConversion.ReadFile(obmol, res)
  atom_number = obmol.NumAtoms()
  # printf("%i Atoms read\n", atom_number)
  with open(res) as file:
    xyz_lines = [line.rstrip() for line in file]
  links_num = xyz_lines[1]
  links_num = links_num.split()
  links_num = int(links_num[-1])
  # printf("%i links in the file\n", links_num)
  if links_num == 0:
    no_links += 1
    continue
  link_list = []
  for n in range(0, links_num, 1):
    link_cnt += 1
    pair = xyz_lines[atom_number+n+4]
    pair = pair.split()
    pair[0] = int(pair[0])
    pair[1] = int(pair[1])
    # printf("%i %i\n", pair[0], pair[1])
    link_list.append(pair)
  # print(link_list)

  # Use the link list to analyze bonds
  for link in link_list:
    # print(pair)
    obatom = obmol.GetAtom(pair[0]+1)
    label1=ob.GetSymbol(obatom.GetAtomicNum())
    label1=label1.lower()
    label1+=str(obatom.GetTotalValence())
    # print(label1)
    obatom = obmol.GetAtom(pair[1]+1)
    label2=ob.GetSymbol(obatom.GetAtomicNum())
    label2=label2.lower()
    label2+=str(obatom.GetTotalValence())
    # print(label2)
    label=label1+label2
    # print(label)
    # add the link to the dictionary
    if label not in link_dict.keys():
      link_dict[label] = 1
      # print(link_dict)
    else:
      link_dict[label] += 1
      # print(link_dict)

printf("%4i Structures without liks\n", no_links)
printf("%4i Links analyzed\n", link_cnt)
printf("%4i Bond types found\n", len(link_dict))
# print(link_dict)
# Sorting the dictionary by value in descending order
sorted_dict = dict(sorted(link_dict.items(), key=lambda item: item[1], reverse=True))
print(sorted_dict)

exit()

# Save the dictionary
with open('link_dict.pkl', 'wb') as fp:
  pickle.dump(sorted_dict, fp)
  print('dictionary saved successfully to file')

exit()