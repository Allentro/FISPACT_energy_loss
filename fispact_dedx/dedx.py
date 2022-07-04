import numpy as np
from scipy import interpolate
import srim

def find_nearest(array, energy):
    idx = (np.abs(array - energy)).argmin()
    return idx

def projectile_range(reaction, density, top_energy=40):
    energy, stopping_power, projected_range = pysrim.run_srim_wrap(reaction, density)
    energy = np.insert(energy, 0, 0)
    projected_range = np.insert(projected_range, 0, 0)
    f = interpolate.interp1d(energy, projected_range, kind='quadratic')
    steps = 0.0001
    energy = np.arange(0, top_energy+steps, steps)
    projected_range = f(energy)
    return energy, projected_range

def arb_flux_thickness(energy, projected_range, thickness, number_of_slices, top_energy):
    energy_drop = top_energy / number_of_slices
    energy_array = np.arange(0,top_energy+energy_drop, energy_drop)
    energy_list = list(energy_array)
    range_list = np.array([projected_range[find_nearest(energy,  x)] for x in energy_list])
    thickness_array = np.diff(range_list)
    return thickness_array, energy_array[:-1], energy_array[1:]

def thickness_between_energies(energy, projected_range, top_energy, bottom_energy): 
    top_range = projected_range[find_nearest(energy, top_energy)]
    bottom_range = projected_range[find_nearest(energy, bottom_energy)]
    return top_range - bottom_range

def group_structure_thickness(energy, projected_range, group_structure, top_energy, flux, min_value=1):
    energy = energy * 1e6
    top_energy = top_energy * 1e6 
    min_value = min_value * 1e6
    array = np.loadtxt(f"{group_structure}energy.txt")
    idx_top_energy = int(array[find_nearest(array[:,1],top_energy),0])
    if top_energy != array[find_nearest(array[:,1],top_energy),1]: 
        print(f"The top energy of {top_energy} is not a boundary in the group structure CCFE-{group_structure}. A new top_energy of {array[find_nearest(array[:,1],top_energy),1]} has been used")
    idx_min_value = int(array[find_nearest(array[:,1],min_value),0])
    if min_value != array[find_nearest(array[:,1],min_value),1]: 
        print(f"The 'min_value' of {min_value} is not a boundary in the group structure CCFE-{group_structure}. A new top_energy of {array[find_nearest(array[:,1],top_energy),1]} has been used")
    empty = np.arange(0, int(group_structure), 1)
    coutner = 0
    lower_energy = []
    upper_energy = []
    thickness = []
    flux_dic = {}
    for value in range(idx_top_energy-1, idx_min_value): 
        if value+1 == idx_min_value:
            new_array = (empty >= value).astype(int)*float(flux)
            b_energy = 0
        else: 
            new_array = (empty == value).astype(int)*float(flux) # +1 is used as the array goes from 
            b_energy = array[value+1, 1]
        lower_energy.append(b_energy)
        upper_energy.append(array[value, 1])
        flux_dic[array[value, 1]] = new_array
        thickness.append(thickness_between_energies(energy, projected_range, array[value,1], array[value+1,1]))
    
    return np.array(thickness), np.array(lower_energy), np.array(upper_energy), flux_dic