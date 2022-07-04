import numpy as np
from srim import TRIM, Ion, Layer, Target, SR
from srim.output import Results
srim_executable_directory = "/home/alletro/SRIM"

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


def run_srim(reaction, density=None):
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

def bottom_energy(top_energy, energy, projected_range, thickness, prnt=True): 
    '''
    Calculating the energy at the rear of a target given a beam energy and target thickness 
    top_energy = energy of the beam entering the target (MeV)
    energy = array of energy output from SRIM (MeV)
    projected_range = array of projected_range output from SRIM (MeV/mm)
    thickness = target thickness (mm) 
    '''
    def inner_energy(top_energy, energy, projected_range, thickness, prnt):
        bottom_range = np.interp(float(top_energy), energy, projected_range) - float(thickness)
        if bottom_range < 0: 
            if prnt == True:
                print(f'{top_energy}MeV --> {0}MeV (stopped in target)')
            return 0
        #print(np.interp(bottom_range, range_y, energy_x))
        else:
            if prnt == True: 
                print(f'{top_energy}MeV --> {round(np.interp(bottom_range, projected_range, energy), 2)}MeV')
            return round(np.interp(bottom_range, projected_range, energy), 2)
    
    if isinstance(top_energy, list): 
        for item in top_energy:
            bottom_energy = inner_energy(item, energy, projected_range, thickness, prnt)
    else: 
        bottom_energy = inner_energy(top_energy, energy, projected_range, thickness, prnt)
    return bottom_energy

def extract_xs(xs_path,product):
    data = pd.read_csv(xs_path)
    if 'index' in data.columns: 
        del data['index']
    xs_array = data.loc[data['product'] == product]
    energy = data.columns[2:].astype(float)
    xs = xs_array.values[0][2:].astype(float)
    return energy, xs

def plotting_range(energy, stopping_power, projected_range=None, xscale='lin', yscale='lin', top_energy=None, thickness=None, xmax=None, ti=0, al=0):
    fig, ax1 = plt.subplots()
    
    color = 'tab:red'
    #ax1.set_title('Stopping Power')
    ax1.set_xlabel('Energy (MeV)')
    xmin=0
    if xscale == 'log':
        ax1.set_xscale('log')
        xmin = 0.001
    ax1.plot(energy, stopping_power, color=color)
    ax1.set_xlim(left=xmin, right=xmax)
    ax1.set_ylabel('Stopping Power (MeV/mm)', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    if yscale =='log':
        ax1.set_yscale('log')
    if isinstance(projected_range, np.ndarray):
        ax2 = ax1.twinx()  

        color = 'tab:blue'
        if top_energy!= None and thickness != None:
            e,dedx, pr = run_srim('Ti-d')
            lower_ti = bottom_energy(top_energy, e, pr, ti, prnt=False)
            ax2.axvspan(lower_ti, top_energy, alpha=0.5, color='grey')
            e,dedx, pr = run_srim('Al-d')
            lower_al = bottom_energy(lower_ti, e, pr, al, prnt=False)
            ax2.axvspan(lower_al, lower_ti, alpha=0.5, color='black')
            ax2.axvspan(bottom_energy(lower_al, energy, projected_range, thickness, prnt=False), lower_al, alpha=0.5, color='green')
            ax2.axvline(x=bottom_energy(lower_al, energy, projected_range, thickness/2, prnt=False), ls='--', color='green')
        ax2.set_ylabel('Projected range (mm)', color=color)  
        ax2.plot(energy, projected_range, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        if yscale =='log':
            ax2.set_yscale('log')
    return 

def plotting_xs(energy, stopping_power, xs=None, xscale='lin', yscale='lin', top_energy=None, thickness=None, xmax=None, product=None, ti=0, al=0):
    fig, ax1 = plt.subplots()
    
    color = 'tab:red'
    #ax1.set_title('Stopping Power')
    ax1.set_xlabel('Energy (MeV)')
    xmin=0
    if xscale == 'log':
        ax1.set_xscale('log')
        xmin = 0.001
    ax1.plot(energy, stopping_power, color=color)
    ax1.set_xlim(left=xmin, right=xmax)
    ax1.set_ylabel('Stopping Power (MeV/mm)', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    if yscale =='log':
        ax1.set_yscale('log')
    if isinstance(xs, str):
        ax2 = ax1.twinx()  

        color = 'tab:blue'
        if top_energy!= None and thickness != None:
            e,dedx, pr = run_srim('Ti-d')
            lower_ti = bottom_energy(top_energy, e, pr, ti, prnt=False)
            ax2.axvspan(lower_ti, top_energy, alpha=0.5, color='grey')
            e,dedx, pr = run_srim('Al-d')
            lower_al = bottom_energy(lower_ti, e, pr, al, prnt=False)
            ax2.axvspan(lower_al, lower_ti, alpha=0.5, color='black')
            ax2.axvspan(bottom_energy(lower_al, energy, projected_range, thickness, prnt=False), lower_al, alpha=0.5, color='green')
            ax2.axvline(x=bottom_energy(lower_al, energy, projected_range, thickness/2, prnt=False), ls='--', color='green')
        ax2.set_ylabel('Cross-section (mb)', color=color)  
        xs_energy, xs = extract_xs(xs, product)
        ax2.plot(xs_energy, xs, color=color, label=f'XS({product})')
        ax2.tick_params(axis='y', labelcolor=color)
        if yscale =='log':
            ax2.set_yscale('log')
    return 