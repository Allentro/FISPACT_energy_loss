from . import create_files
from . import convert
from . import dedx

def run_dedx_fispact(root_dir, flux_type='group'):
    energy, projected_range = dedx.projectile_range(reaction, density)
    element, projectile = convert.convert_reaction(reaction)
    if flux_type == 'group': 
        if projectile == 'n': 
            group_structure == '709'
        else: 
            group_structure == '162'
        arb_flux=False
        thickness, lower_energy, upper_energy, flux_dic = dedx.group_structure_thickness(energy, projected_range, group_structure, top_energy, flux)
    else: 
        arb_flux=True
        thickness, lower_energy, upper_energy = dedx.arb_flux_thickness(energy, projected_range, thickness, number_of_slices, top_energy)
    create_files.create_files(element, projectile, density, groupset, parent_bane, upper_energy, lower_energy, thickness, flux, inventory, arb_flux=False, flux_dic=False)
    create_files.execute_files(upper_energy)
    return