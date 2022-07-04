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