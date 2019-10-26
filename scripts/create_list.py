import numpy as np
import sys 
import os

#qm path directory
qmpath = "/mnt/raidc3/kbelfon/18.cnx/02.QM/CNX_cavern/CNX/" 
# mm path directory 
mmpath = "/mnt/raidc3/kbelfon/18.cnx/04.MM_after_QM_mod3"

struct_list = np.arange(1, 501, 1) 

t = open("not_used_for_fitting", "w") 
k = open("used_for_fitting", "w") 
backbone = ['alpha', 'opt']
crd_mm = []; crd_qm = []
for bb in backbone:
  suspect = []
  suspect = np.genfromtxt("%s/filter/%s.suspect_list_sorted" %(qmpath,bb), usecols=0, dtype=str)
  rej_list=np.genfromtxt('%s/%s_reject.dat' %(qmpath,bb), usecols=0, dtype=str)
  print (suspect)
  for om in struct_list:
    # if the structure not in the suspect list for that dataset e.g alpha or opt
    if "struct%s.rst7" %(om) not in suspect and "struct%s.rst7" %(om) not in rej_list: 
      string_qm = "%s/%s_struct%s/Energy.log" %(qmpath, bb, om)
      string_mm = "%s/%s_struct%s/struct%s.zero" %(mmpath, bb, om, om)
      crd_mm.append(string_mm)
      crd_qm.append(string_qm)
      k.write("%s_struct%s\n" %(bb,om))
    elif "struct%s.rst7" %(om) in suspect: 
      t.write("%s_struct%s\n" %(bb,om))

  crd_mm.append("break")
  crd_qm.append("break")

k1 = open("list_mmstructure.dat", "w")
for string in crd_mm:
  #os.system("cp %s.r %s.rst7" %(string, string))
  k1.write("%s\n" %(string))

k2 = open("list_qmstructure.dat", "w")
for string in crd_qm:
  k2.write("%s\n" %(string))

#for debugging 
debug = open("all_suspects", "w")
for bb in backbone:
  suspect = np.genfromtxt("%s/filter/%s.suspect_list_sorted" %(qmpath, bb), usecols=0, dtype=str)
  for struct in suspect:
    debug.write("%s\n" %(struct)) 
