import os 
from stat import S_IEXEC 
from .utilities import split_reaction
import subprocess 
import time 
import numpy
import tqdm

def projectile_conversion(projectile): 
    if projectile == 'n': 
        projectile = '1' 
    elif projectile == 'd': 
        projectile = '2' 
    elif projectile == 'p': 
        projectile = '3'
    elif projectile == 'a': 
        projectile = '4'
    elif projectile == 'g': 
        projectile = '5' 
    elif projectile == 't': 
        projectile = '6'
    elif projectile == 'h': 
        projectile = '7'
    else: 
        raise ValueError(f"Projectile must be either 'n', 'p', 'd', 'g', 'h', 't' or 'a' not '{projectile}'")
    return projectile

def create_directory(parent_name, upper_energy):
    os.makedirs(f"{parent_name}", exist_ok=True)
    for energy in upper_energy:
        os.makedirs(f"{parent_name}/{round(energy,4)}", exist_ok=True)
    return


def files_convert(parent_name, energy): 
    with open(f"{parent_name}/{round(energy,4)}/files.convert", "w") as file: 
        file.write("# index of nuclides to be included \n")
        file.write("ind_nuc /opt/fispact/nuclear_data/decay/decay_2012_index_2012 \n")
        file.write("\n")
        file.write("# fluxes \n")
        file.write("fluxes  fluxes \n")
        file.write("arb_flux  arb_flux \n")
    return


def arb_flux(parent_name, upper_energy, lower_energy, flux): 
    with open(f"{parent_name}/{round(upper_energy,4)}/arb_flux", "w") as file: 
        energy_diff = (upper_energy - lower_energy)/6
        file.write(f"{(upper_energy+energy_diff*2)*1e6} {(upper_energy+energy_diff)*1e6} {(upper_energy)*1e6} {(lower_energy+energy_diff*5)*1e6} {(lower_energy+energy_diff*4)*1e6} {(lower_energy+energy_diff*3)*1e6} {(lower_energy+energy_diff*2)*1e6} {(lower_energy+energy_diff)*1e6} {(lower_energy)*1e6} {(lower_energy-energy_diff)*1e6} {(lower_energy-energy_diff*2)*1e6}\n")
        file.write(f"0 0 {flux} {flux} {flux} {flux} {flux} {flux} 0 0\n")
        file.write("1 \n")
        file.write(f"# Flux spectrum for {lower_energy}-{upper_energy} energy range")
    return

def create_flux_file(parent_name, upper_energy, flux): 
    if len(flux) == 162:
        flux = flux.reshape((-1,6))
        with open(f"{parent_name}/{upper_energy}/fluxes", "w") as file:
            for i in range(flux.shape[0]): 
                file.write(f"  {' '.join(map(str, flux[i,:]))} \n")
            file.write(" 0.01 \n")
            file.write(f"# Flux spectrum for energy up to {upper_energy} \n")
    elif len(flux) == 709: 
        v = 0
    else:
        raise ValueError("The length of the flux file does not match either of group structures 162 or 709.")
    return flux


def files(parent_name, energy, projectile, groupset, arb_flux):
    with open(f"{parent_name}/{round(energy,4)}/files", "w") as file: 
        file.write(f"# gamma attenuation data \n")
        file.write(f"absorp  /opt/fispact/nuclear_data/decay/abs_2012 \n")
        file.write(f" \n")
        file.write(f"# index of nuclides to be included \n")
        file.write(f"ind_nuc  /opt/fispact/nuclear_data/TENDL2017data/tendl17_decay12_index \n")
        file.write(f" \n")
        file.write(f"# Library cross section data \n")
        file.write(f"xs_endf  /opt/fispact/nuclear_data/TENDL2017data/tal2017-{projectile}/gxs-{groupset} \n")
        file.write(f" \n")
        file.write(f"# Library of spallation cross-section data \n")
        file.write(f"sp_endf  /opt/fispact/nuclear_data/HEIR01data/heir01-{projectile}/gxs-{groupset} \n")
        #file.write(f" \n")
        #file.write(f"# Additional cross section data \n")
        ##file.write(f"xs_extra ./xs_extra \n")
        file.write(f"  \n")
        file.write(f"# Library probability tables for self-shielding \n")
        file.write(f"prob_tab  /opt/fispact/nuclear_data/TENDL2017data/tal2017-n/tp-709-294 \n")
        file.write(f" \n")
        if arb_flux:
            file.write(f"# abs fluxes \n")
            file.write(f"arb_flux  arb_fluxes \n")
            file.write(f" \n")
        file.write(f"# fluxes \n")
        file.write(f"fluxes  fluxes \n")
        file.write(f" \n")
        file.write(f"# Library decay data \n")
        file.write(f"dk_endf /opt/fispact/nuclear_data/decay/decay_2012 \n")
        file.write(f" \n")
        file.write(f"# Library fission  data \n")
        file.write(f"fy_endf /opt/fispact/nuclear_data/UKFY41data/ukfy4_1p \n")
        file.write(f" \n")
        file.write(f"# Spontaneous fission data \n")
        file.write(f"sf_endf /opt/fispact/nuclear_data/GEFY61data/gefy61_sfy \n")
        file.write(f" \n")
        file.write(f"# collapsed cross section data (in and out) \n")
        file.write(f"collapxi  COLLAPX \n")
        file.write(f"collapxo  COLLAPX \n")
        file.write(f" \n")
        file.write(f"# condensed decay and fission data (in and out) \n")
        file.write(f"arrayx  ARRAYX \n")
    return


def convert_i(parent_name, energy, groupset):
    with open(f"{parent_name}/{round(energy,4)}/convert.i", "w") as file: 
        file.write(f"<< convert flux to {groupset} group structure>> \n")
        file.write("CLOBBER \n")
        file.write(f"GRPCONVERT 10 {groupset} \n")
        file.write("FISPACT \n")
        file.write("* SPECTRAL MODIFICATION \n")
        file.write("END \n")
        file.write("* END \n")
    return
        
    
def fisprun_sh(parent_name, energy, inventory, arb_flux): 
    with open(f"{parent_name}/{round(energy,4)}/fisprun.sh", "w") as file: 
        file.write("#!/bin/bash \n")
        file.write("# getting_started/proton_HEIR \n")
        file.write(" \n")
        file.write("if [ -z ${FISPACT+x} ]; then \n")
        file.write("    echo -e '\033[1;91mEnvironment variable $FISPACT is not set, you need to set this to point to $FISPACT executable. For example $$FISPACT=/path/to/$FISPACT\033[0m' \n")
        file.write("exit \n")
        file.write("    fi \n")
        file.write(" \n")
        file.write(f"echo 'Simulation TENDL-2017 to calculate {energy}MeV yield' \n")
        file.write(" \n")
        file.write("rm *.log ARRAYX* COLLAPX* *.out *.gra *.tab* *.plt &>/dev/null \n")
        if arb_flux:
            file.write("$FISPACT convert files.convert \n")
        file.write("$FISPACT collapse \n")
        if inventory:
            file.write("$FISPACT inventory \n")
    os.chmod(f"{parent_name}/{round(energy,4)}/fisprun.sh", S_IEXEC | os.stat(f"{parent_name}/{round(energy,4)}/fisprun.sh").st_mode)
    return


def collapse(parent_name, energy, groupset, projectile): 
    projectile = projectile_conversion(projectile)
    with open(f"{parent_name}/{round(energy,4)}/collapse.i", "w") as file: 
        file.write("<< -----collapse cross section data----- >>\n")
        file.write("CLOBBER\n")
        file.write(f"PROJ {projectile} \n")
        file.write(f"GETXS 1 {groupset}\n")
        file.write("GETDECAY 1\n")
        file.write("FISPACT\n")
        file.write("* COLLAPSE tal2017-a/gxs-709 decay12 index\n")
        file.write("PRINTLIB 4\n")
        file.write("END\n")
        file.write("* END OF RUN\n")
    return

def calculate_mass(thickness): 
    thickness_m = thickness / 1000 
    area_m = 1 * 0.0001
    volume = thickness_m * area_m
    return volume * (10.28*1000)

def inventory_i(parent_name, energy, thickness, element, projectile, density, flux):
    mass = calculate_mass(thickness)
    projectile = projectile_conversion(projectile)
    with open(f"{parent_name}/{round(energy,4)}/inventory.i", "w") as file: 
        file.write("<< -----get nuclear data----- >> \n")
        #file.write("MONITOR 1 \n")
        file.write(f"PROJ {projectile} \n")
        file.write("NOERROR \n")
        file.write("GETXS 0 \n")
        file.write("GETDECAY 0 \n")
        file.write("FISPACT \n")
        file.write("* IRRADIATION OF Pb \n")
        file.write("<< -----set initial conditions----- >> \n")
        file.write(f"MASS {mass:e} 1 \n")
        file.write(f"{element} 100.0 \n")
        file.write(f"DENSITY {density} \n")
        file.write("MIND 1.E3 \n")
        file.write("GRAPH 1 2 1 1 \n")
        file.write(f"FLUX {round(float(flux), 5)} \n")
        file.write("ATOMS \n")
        file.write("HALF \n")
        file.write("<< -----irradiation phase----- >> \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("TIME 30 MINS \n")
        file.write("ATOMS \n")
        file.write("<< -----cooling phase----- >> \n")
        file.write("FLUX 0. \n")
        file.write("ZERO \n")
        file.write("TIME 1 HOURS ATOMS \n")
        file.write("TIME 1 HOURS ATOMS \n")
        file.write("TIME 1 HOURS ATOMS \n")
        file.write("TIME 1 HOURS ATOMS \n")
        file.write("TIME 1 HOURS ATOMS \n")
        file.write("TIME 1 HOURS ATOMS \n")
        file.write("TIME 6 HOURS ATOMS \n")
        file.write("TIME 6 HOURS ATOMS \n")
        file.write("TIME 6 HOURS ATOMS \n")
        file.write("TIME 1 DAYS ATOMS \n")
        file.write("TIME 1 DAYS ATOMS \n")
        file.write("TIME 1 DAYS ATOMS \n")
        file.write("TIME 1 DAYS ATOMS \n")
        file.write("TIME 1 DAYS ATOMS \n")
        file.write("TIME 1 DAYS ATOMS \n")
        file.write("TIME 7 DAYS ATOMS \n")
        file.write("END \n")
        file.write("* END \n")
    return 

def create_files(parent_name, reaction, density, groupset, upper_energy, lower_energy, thickness, flux, inventory, arb_flux=False, flux_dic=False):
    element, projectile = split_reaction(reaction)
    create_directory(parent_name, upper_energy)
    for idx in range(len(upper_energy)):
        files(parent_name, upper_energy[idx], projectile, groupset, arb_flux)
        if arb_flux:
            files_convert(parent_name, upper_energy[idx])
            arb_flux(parent_name, upper_energy[idx], lower_energy[idx], flux)
            convert_i(parent_name, upper_energy[idx], groupset)
        else: 
            create_flux_file(parent_name, upper_energy[idx], flux_dic[upper_energy[idx]])
        collapse(parent_name, upper_energy[idx], groupset, projectile)
        fisprun_sh(parent_name, upper_energy[idx], inventory, arb_flux)
        inventory_i(parent_name, upper_energy[idx], thickness[idx], element, projectile, density, flux)
    return

def execute_file(root_dir, energy, output=False):
    os.chdir(f"{root_dir}/{energy}")
    if output == False: 
        subprocess.call(f"{root_dir}/{energy}/fisprun.sh", shell=True, stdout=subprocess.DEVNULL)
    else: 
        subprocess.call(f"{root_dir}/{energy}/fisprun.sh", shell=True)
    return

def execute_run(root_dir, upper_energy=None, output=False, progress=True): 
    energy_list = next(os.walk(root_dir))[1]
    if progress: 
        for idx in tdqm(range(len(energy_list)), desc='Run progress'):
            execute_file(root_dir, energy_list[idx], output)  
    else: 
        for idx in range(len(energy_list)):
            execute_file(root_dir, energy_list[idx], output)
    return