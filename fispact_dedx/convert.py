def current_to_flux(I, charge=2): 
    c_alpha = charge * (1.6022*10**-19)
    return I / c_alpha 

def convert_reaction(reaction): 
    e, p = reaction.split('-')
    return e, p