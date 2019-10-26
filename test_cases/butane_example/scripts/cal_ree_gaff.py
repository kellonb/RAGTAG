import numpy as np
import sys


# use col = 2 can change based on dihedrals
QM = np.genfromtxt("genA_input", skip_header=2, skip_footer=1, usecols=1)
#MM0 = np.genfromtxt("genA_input", skip_header=2, skip_footer=1, usecols=2)
MM = np.genfromtxt("MM_gaff_energy.dat", usecols=0)

#print (MM0)

#QM = np.arange(4,8,1)

f = open("GAFF_REE_for_pairs", "w")
f.write("pair   pair conf     QMi          QMj           MMi        MMj      REE\n")  

counter = -1  # begin the counter
num_conf = len(QM) # number of conformation 
QM = np.insert(QM, num_conf, 0) # insert 0 at the end of the QM array, to loop all the way to the end of the array
MM = np.insert(MM, num_conf, 0) # insert 0 at the end of the MM array, 
print (QM)
print (num_conf)
weight = 2/(num_conf * (num_conf-1)) # weight
REE_SUM = 0 # array that keep summing the REEs for each pairs
pair_counter = 0 # to count the total number of pairs
indx = np.arange(1, len(QM), 1) # hold the index for the conformations, e.g conformation 1

# do for each pair of i and j 
for i, im, i_counter in zip(QM,MM,indx):
  counter +=1
  QMp1 = np.delete(QM, np.s_[0:counter+1]) # delete the element from 0 to i+1 i n the array
  MMp1 = np.delete(MM, np.s_[0:counter+1]) # delete the element from 0 to i+1 i n the array
  indxp1 = np.delete(indx, np.s_[0:counter+1]) # delete the element from 0 to i+1 i n the array
  for j, jm, j_counter in zip(QMp1,MMp1,indxp1):
    pair_counter += 1
    QMpair = i - j
    MMpair = im - jm
    REE_pair = abs(QMpair - MMpair)
    f.write("%4d (%4d,%4d)  %8.4f  %8.4f  %4.4f  %4.4f  %4.4f\n" %(pair_counter, i_counter, j_counter, i, j, im, jm, REE_pair))
    REE_SUM += REE_pair

AAE = weight * REE_SUM
#print ("REE sum")
#print (REE_SUM)
f.write("Final AAE is %s" %(AAE))
print ("FInal AAE:")
print (AAE)
