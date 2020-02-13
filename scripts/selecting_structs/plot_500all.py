import warnings
import pytraj as pt
import seaborn as sns
from pytraj import matrix
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os 

#usage /mnt/raidc2/haichit/anaconda3/bin/python
#home/haichit/anaconda3/lib/python3.5/

#top = "../CNX.ff14SB.MM.parm7" 
#crd = "output.alltraj.nc"
#traj = pt.iterload(crd, top)

#function to plot ramachadaran:
def ram(x,y,z,tit,l1,l2):
    #alpha set transparency
    import seaborn as sns; sns.set(style="ticks")
    import matplotlib as mpl
    fig = plt.figure(figsize=(8,8))#width, height in  inches
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    g = sns.JointGrid(x, y)
    col = np.where(z<=43.00, 'b', np.where(z<60,'r', 'k'))
    g = g.plot_joint(plt.scatter, color=col, s=5) #, alpha=1, edgecolor="white")
    #g = g.plot_joint(sns.kdeplot, cmap='Spectral', shade=False, shade_lowest=True)
    _ = g.ax_marg_y.hist(y, color="cornflowerblue", orientation="horizontal", bins=np.arange(-180,180, 2))
    _ = g.ax_marg_x.hist(x, color="palevioletred", bins=np.arange(-180, 180, 2)) #bins is the range from (min max interval)
    #g.set_axis_labels(r"%s ($\deg$)", r"%s ($\deg$)")
    plt.xlim(-180,180)
    plt.ylim(-180,180)
    plt.xticks(np.arange(-180, 180, 30))
    plt.yticks(np.arange(-180, 180, 30))
    plt.xlabel(r"%s $\degree$" %l1, fontsize=20)
    plt.ylabel(r"%s $\degree$" %l2, fontsize=20)
    plt.grid(b=True, which='both', color='gray', linewidth=0.13, linestyle='--')
    plt.tight_layout()
    return plt.savefig('%s.png' % tit)

# calculate and plot dihedral
#chi1 = pt.dihedral(traj, ':2@N :2@CA :2@CB :2@CG') 
#chi2 = pt.dihedral(traj, ':2@CA :2@CB :2@CG :2@CD1') 

chi1 = np.genfromtxt('./prod.dih.all', usecols=0, dtype=float)
chi2 = np.genfromtxt('./prod.dih.all', usecols=1, dtype=float)
chi3 = np.genfromtxt('./prod.dih.all', usecols=2, dtype=float)
chi4 = np.genfromtxt('./prod.dih.all', usecols=3, dtype=float)
chi5 = np.genfromtxt('./prod.dih.all', usecols=4, dtype=float)
z = np.genfromtxt('./prod.E.all', usecols=0, dtype=float)


ram(chi1, chi2, z, 'chi1_vs_chi2.all', 'chi1', 'chi2')
ram(chi5, chi4, z, 'chi4_vs_chi5.all', 'chi5', 'chi4')

ram(chi1, chi3, z, 'chi1_vs_chi3.all', 'chi1', 'chi3')
ram(chi2, chi3, z, 'chi2_vs_chi3.all', 'chi2', 'chi3')
ram(chi4, chi3, z, 'chi4_vs_chi3.all', 'chi4', 'chi3')
ram(chi5, chi3, z, 'chi5_vs_chi3.all', 'chi5', 'chi3')

