def current_to_flux(I, charge=1): 
    c_alpha = charge * (1.6022*10**-19)
    return I / c_alpha 

def split_reaction(reaction): 
    target, projectile = reaction.split('-')
    return target, projectile