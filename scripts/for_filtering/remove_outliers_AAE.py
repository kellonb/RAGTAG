import numpy as np
import matplotlib.pylab as plt
import seaborn as sns; sns.set(style="ticks")
import sys

def plot_bar(AAE, title):
    fig = plt.figure(figsize=(16,8))
    sns.set_style({"xtick.direction": "out","ytick.direction": "out"})
    ax = fig.add_subplot(111)
    struct = np.arange(len(AAE))
    # plot data points 
    #plt.hist((eng14,eng14new,eng18), bins=np.arange(0, 3, 0.1), label = ['14SB', '14SBnew', '18SBnew'], color = ['blue', 'red', 'purple']) #alpha = 0.8, color='b')
    plt.bar(struct, AAE, align='center')
    #plt.legend()
    plt.xlabel('Structures')
    plt.ylabel('AAE')
    #plt.xlim(0,3)
    #plt.xticks((np.arange(-180,181,30)))
    #plt.yticks((np.arange(-180,181,30)))
    #plt.xlim(-210,210)
    return plt.savefig('%s.png' % title)  


backbone = ['alpha', 'opt']
for bb in backbone:
  # QM energies 
  QM = np.genfromtxt("%s_energies.dat" %(bb), usecols=1)
  # MM0 + MM_RAGTAG 
  MM = np.genfromtxt("%s_energies.dat" %(bb), usecols=2)
  
  f = open("Final_%s_AAE_for_structure" %(bb), "w")
  f.write("structure   AAE\n")  

  counter = -1  # begin the counter
  num_conf = len(QM) # number of conformation 
  weight = 2/(num_conf * (num_conf-1)) # weight
  pair_counter = 0 # to count the total number of pairs
  indx = np.arange(1, len(QM), 1) # hold the index for the conformations, e.g conformation 1

  ALL_AAE = []; struct = []
  # for a give QM and MM energy (ensuring all conformation is used as a reference) 
  for qmref, mmref in zip(QM,MM):
    REE_SUM = 0 # array that keep summing the REEs for each pairs
    counter +=1
    # loop through all QM and MM energies
    for qm, mm in zip(QM,MM):
      # calculate E_QM - E_QMref
      QMpair = qm - qmref
      # calculate E_MM - E_MM ref 
      MMpair = mm - mmref
      # Now calculate the relative energy of the pair (pairwise energy to a given reference)
      REE_pair = abs(QMpair - MMpair)
      # SUmm over the REE to the get the Average absolute error after
      REE_SUM += REE_pair
    # Absolute Average error is the 
    AAE = weight * REE_SUM
    f.write("%8d %4.10f\n" %(counter, AAE))
    ALL_AAE.append(AAE)
    struct.append(counter)
  
  # set up a dictionary
  AAE_dict = dict(zip(ALL_AAE, struct))
  plot_bar(ALL_AAE, "%s_structAAEnew" %(bb))
  print ("minimum AAE for %s is %s (struct%s)" %(bb, min(ALL_AAE), AAE_dict[min(ALL_AAE)]))
  f.write("minimum AAE for %s is %s (struct%s)" %(bb, min(ALL_AAE),AAE_dict[min(ALL_AAE)] ))
  print ("Final AAE for %s is %s" %(bb, (sum(ALL_AAE))/2))
  f.write("Final AAE for %s is %s" %(bb, (sum(ALL_AAE))/2))

  # Here setting up the filtering by IQR
  mean = np.mean(ALL_AAE)
  std = np.std(ALL_AAE)
  maxd = np.max(ALL_AAE)
  mind = np.min(ALL_AAE)
  median = np.median(ALL_AAE)
  
  # sort the data
  sort = sorted(ALL_AAE)
  # find 1st and 3rd quartile
  q1, q3 = np.percentile(sort, [25,75])
  # find the IQR
  iqr = q3 - q1 
  # find the lower and upper bound
  lbound = q1 - (1.5 * iqr)
  ubound = q3 + (1.5 * iqr)

  AAE_IQR_pass=[]
  sus = open("%s_suspect.dat" %(bb), "w")
  nosus = open("%s_nonsuspect.dat" %(bb), "w")
  for ae in (sort):
    if (ae <= lbound) or (ae >= ubound):
       z = (ae - mean) / std
       sus.write("struct%s  %s\n" %(AAE_dict[ae], ae)) 
    else: 
       AAE_IQR_pass.append(ae)
       nosus.write("struct%s  %s\n" %(AAE_dict[ae], ae))

#

  
