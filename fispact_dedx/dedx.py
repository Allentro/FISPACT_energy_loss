import numpy as np
from . import data
from scipy import interpolate
from srim import TRIM, Ion, Layer, Target, SR
from srim.output import Results
from pathlib import Path


DATA_PATH = Path(data.__path__[0])

# pySRIM wrapper ===================================================================== 
# ====================================================================================

def split_projectile(reaction):
    '''
    Function to split the target material and projectile
    '''
    try: 
        return reaction.split('-')
    except: 
        raise ValueError('Entered target must be in the form of <TARGET>-<PROJECTILE>')

def single_dictionary(element, dictionary): 
    '''
    For any individual component of a target, 
    the dictionary required by srim is appended to
    '''
    if '(' in element:
        print(element)
        print(element.split(')')[0].split('(')[1])
        try: 
            stoci = int(element.split(')')[0].split('(')[1])
            element = element.split(')')[1]
        except: 
            raise ValueError('With greater than one stochiometry, the required format is an integar as: (<Stochiometry>)<Element>')
    else: 
        stoci = 1
    dictionary[element]= {
            "stoich": stoci,
            "E_d": 35.0,  # Displacement Energy
            "lattice": 0.0,
            "surface": 3.0,
        }
    return dictionary

def dictionary_construct(target):  
    '''
    This function checks if the target is composed of multiple elements 
    If so, the target is broken up and a dictionary is constructed for each.
    '''
    target_dic = {}
    if ':' in target: 
        sub_target = target.split(':') 
        for item in sub_target:
            target_dic = single_dictionary(item, target_dic)
    else: 
        target_dic = single_dictionary(target, target_dic)
    return target_dic

def check_density(reaction, density):
    ''' 
    Checks for user input density
    For single elements if the density isn't input standard NIST data will be used
    '''
    density_dict = {'H': 0.0708, 'He': 0.147, 'Li': 0.534, 'Be': 1.848, 'B': 2.34, 'C': 2.25, 'N': 0.808, 'O': 1.149, 
                'F': 1.108, 'Ne': 1.204, 'Na': 0.971, 'Mg': 1.738, 'Al': 2.6989, 'Si': 2.33, 'P': 1.82, 'S': 2.07, 
                'Cl': 1.56, 'Ar': 1.4, 'K': 0.856, 'Ca': 1.55, 'Sc': 2.99, 'Ti': 4.54, 'V': 6.11, 'Cr': 7.18, 
                'Mn': 7.21, 'Fe': 7.874, 'Co': 8.9, 'Ni': 8.902, 'Cu': 8.96, 'Zn': 7.133, 'Ga': 5.91, 'Ge': 5.323, 
                'As': 5.73, 'Se': 4.79, 'Br': 3.12, 'Kr': 2.155, 'Rb': 1.532, 'Sr': 2.54, 'Y': 4.47, 'Zr': 6.506, 
                'Nb': 8.57, 'Mo': 10.22, 'Tc': 11.5, 'Ru': 12.41, 'Rh': 12.41, 'Pd': 12.02, 'Ag': 10.5, 'Cd': 8.65, 
                'In': 7.31, 'Sn': 7.31, 'Sb': 6.691, 'Te': 6.24, 'I': 4.93, 'Xe': 3.52, 'Cs': 1.873, 'Ba': 3.5, 
                'La': 6.15, 'Ce': 6.757, 'Pr': 6.773, 'Nd': 7.007, 'Pm': 7.2, 'Sm': 7.52, 'Eu': 5.243, 'Gd': 7.9, 
                'Tb': 8.229, 'Dy': 8.55, 'Ho': 8.795, 'Er': 9.06, 'Tm': 9.321, 'Yb': 6.9654, 'Lu': 9.8404, 'Hf': 13.31, 
                'Ta': 16.654, 'W': 19.3, 'Re': 21.02, 'Os': 22.57, 'Ir': 22.42, 'Pt': 21.45, 'Au': 19.3, 'Hg': 13.546, 
                'Tl': 11.85, 'Pb': 11.35, 'Bi': 9.747, 'Po': 9.32, 'Rn': 4.4, 'Ra': 5.5, 'Th': 11.78, 'Pa': 15.37, 
                'U': 19.05}
    if density != None: 
        return density 
    elif ':' in reaction: 
        raise ValueError('The density must be stated in g/cm3 for all compounds and molecules')
    else: 
        if '(' in reaction: 
            element = reactin.split(')')[1].split('-')[0]
        else: 
            element = reaction.split('-')[0]
        density = density_dict[element]
    return density 

def convert_projectile(projectile):
    '''
    Converts the nuclear projectile notation to element and mass
    '''
    projectile_dic = {'p':['H',1], 'd':['H',2], 't':['H', 3], 'He3':['He',3], 'a':['He', 4]}
    projectile = projectile_dic[projectile]
    return projectile[0], projectile[1]


def run_srim(reaction, srim_executable_directory, density=None):
    '''
    The main function run srim. 
    Contains densities for all elements. 
    If a compound or molecule is given this must be input by the user
    =================================================================
    Inputs
    ------
    reaction: <element>-<projectile>. For molecules <elem1>:<elem2>-<projectile> 
              brackets are used for stochiometry (stochi)<element>-<projectile>
              projectiles: p,d,t,he3,he4
              example: (2)Gd:(3)O-p
    srim_executable_directory: 
    
    density: target density 
             units of g/cm3
    ==================================================================
    Outputs: 
    --------
    energy_array, stopping_power, projected_range
    
    energy array: array of the energy corresponding to the two other outputs
                  units of MeV 
    stopping_power: array of total stopping power (total = electric + nuclear)
                    units of MeV/mm
    projected_range: array of projected range at given energies 
                     units of mm 
    
    '''
    target, projectile = split_projectile(reaction)
    target_dictionary = dictionary_construct(target)
    if density == None: 
        density = check_density(reaction, density)
    
    projectile, mass = convert_projectile(projectile)
    proj = Ion(projectile, energy=1e8, mass=mass)
    
    layer = Layer(
        target_dictionary,
        density=density, 
        width=10,)
    
    srim = SR(layer, proj, number_ions=1000, calculation=1)
    
    results = srim.run(srim_executable_directory)

    (energy_kev,
        electronic_stopping_power,
        nuclear_stopping_power,
        projected_range, #microns
        long_straggle,
        lang_straggle,
    ) = results.data
    projected_range = projected_range/1000
    total_stopping_power =100*(nuclear_stopping_power + electronic_stopping_power)/density
    energy = energy_kev / 1000
    return energy, total_stopping_power, projected_range


# ====================================================================================
# ====================================================================================

def find_nearest(array, energy):
    idx = (np.abs(array - energy)).argmin()
    return idx

def projectile_range(reaction, density, srim_executable_directory, top_energy=100):
    energy, stopping_power, projected_range = run_srim(reaction, srim_executable_directory, density)
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
    array = np.loadtxt(f"{DATA_PATH}/{group_structure}energy.txt")
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