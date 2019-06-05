water = {
    "concentration": {
        "Na":   340,   # mg/L
        "Ca":   223,    # mg/L
        "Mg":   107,   # mg/L
        "K":    7,    # mg/L
        "Cl":   702,  # mg/L
        "SO4":  176,   # mg/L
        "CO3":  0,      # mg/L
        "HCO3": 242,    # mg/L
        "NO3":  416       # mg/L
    },
    "temperature": 29,  # celcius
    "pH": 7.1,
    "water type": "BW"  # "SW": seawater well or "BW": brackish well
}

membrane = {
    "surface area": 7.9,         # m2
    "module diameter": 0.2,     # m
    "spacer thickness": 0.0008, # m
    "module length": 1,         # m
    "a1": 0.802,                # temperatur
    "a2": -0.549,               # feec concentration
    "a3": 0,                    # scaling 
    "b1": 0.970,                # temperature
    "b2": 0.432,                # feec conc. 
    "b3": 1.212,                # scaling 
    "lw0": 24.478,              # solvent permeability
    "ls0": 1.879,               # solute permeability
    "x": 0.000233,              # friction drop 
    "y": 1.735,                 # friction drop
    "membrane type": "BWRO"     # "SWRO": seawater RO or "BWRO": brackish RO or "NF": nanofiltration
}

system = {
    "module recovery": 0.19, 
    "system recovery": 0.20,
    "flow rate": 25,            # m3/h
    "erd efficiency": 0.8,
    "pump efficiency": 0.8            
}