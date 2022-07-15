import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import os
from . import data 
from pathlib import Path

hep.style.use(hep.style.ROOT)

DATA_PATH = Path(data.__path__[0])

def plot_group_structure(root_dir, save=False): 
    os.chdir(root_dir)
    energy = np.loadtxt(f"{DATA_PATH}/162_group.txt")
    energy = energy / 1e6
    upper = energy[:,1]
    lower = energy[:,2]
    mean_162 = (upper+lower)/2
    error_162 = mean_162 - lower 
    a_162 = np.empty((162,))
    a_162[::2] = 1
    a_162[1::2] = 1
    energy = np.loadtxt(f"{DATA_PATH}/709_group.txt") 
    energy = energy / 1e6
    upper = energy[:,1]
    lower = energy[:,2]
    mean_709 = (upper+lower)/2
    error_709 = mean_709 - lower 
    a_709 = np.empty((709,))
    a_709[::2] = -1
    a_709[1::2] = -1
    fig = plt.figure(figsize=(10,3))
    fig.patch.set_facecolor('white')
    plt.errorbar(mean_709, a_709,  xerr=error_709, capsize=15, lw=0, elinewidth=2, label='CCFE-709')
    plt.errorbar(mean_162, a_162,  xerr=error_162, capsize=15, lw=0, elinewidth=2, label='CCFE-162')
    plt.ylim(-2, 2)
    plt.xlim(1, 20)
    plt.xlabel("Energy (MeV)")
    #plt.yticks([-1,1])
    plt.yticks([-1, 1], ['CCFE-709', 'CCFE-162'],
           rotation=0)
    #plt.legend()
    plt.ylabel("Group structure")
    plt.tick_params('both', length=0, width=2, which='major')
    plt.tick_params('both', length=0, width=1, which='minor')
    if save: 
        plt.savefig(save)
    return

def plotting_projected_range(energy, projected_range):
    fig = plt.figure(figsize=(10,10)) 
    plt.xlabel("Projected Range(mm)")
    plt.ylabel("Energy (MeV)")
    plt.plot(projected_range, energy, label="SRIM")
    plt.legend()
    return