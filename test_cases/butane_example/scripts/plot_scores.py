import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import pandas as pd
import pytraj as pt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import os


# python program_name ncp peng ngen filename
#python plot_scores.py 5 5 10000 100 scores1.dat

ncp = sys.argv[1] #number of chromosome print 
peng = sys.argv[2] # print every n generation 
ngen = sys.argv[3] #number of generation 
max_x = sys.argv[4] # maxiumum number of generation u want to plot

filename = np.genfromtxt(sys.argv[5])
scores = np.genfromtxt(sys.argv[5], usecols=2) # get the scores
# but scores are printed based on ncp 

# need to split  scores to give an array with ncp for each peng
scores = np.split(scores, (len(scores))/int(peng))

def plot_avg_score(scores, ncp, peng, ngen, max_x):
  Avg = [];Std = []
  for i in scores:
    sum = np.sum(i)
    std = np.std(i)
    Std.append(std)
    avg = float(sum) / int(ncp)
    Avg.append(avg)
  Avg = np.array(Avg)
  Std = np.array(Std)

  import seaborn as sns; sns.set(style="ticks")
  fig = plt.figure(figsize=(8,4))
  Yarray = Avg
    
  #plot my data, +ncp, because of initial scores (-1)
  x=np.arange(0, int(ngen)+1, int(peng))
  plt.plot(x, Yarray, 'k-', color="blue", label='Butane')
  plt.fill_between(x, Yarray-Std, Yarray+Std)
  plt.scatter(x, Yarray, s=15, color="black")
  #plt.axhline(y=2, color="black", linewidth=1)
 
  plt.ylabel("Score (AAE) ", fontsize=32)
  plt.xlabel('nGen', fontsize=32)
  plt.ylim(0, max(Yarray)) 
  plt.legend(loc=1, prop={'size': 28})
  plt.yticks(np.arange(0,0.31,0.1))
  #plt.xlim(0,max(x))
  plt.xlim(0, 30)
  b = np.arange(-5, 31, 5)
  plt.xticks(b, ['Init', '0', '5', '10', '15', '20', '25', '30'])
  plt.grid(b=True, which='both', color='gray', linewidth=0.13, linestyle='--')
  plt.tick_params(axis='both', which='major', labelsize=30, length=12, width=2, colors='black')
  plt.tight_layout()
  plt.savefig('Average_score_vs_ngen.png')




plot_avg_score(scores, ncp, peng, ngen, max_x) 
#print (scores)
