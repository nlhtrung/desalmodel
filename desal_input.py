water = {
    "concentration": {
        "Na":   4658,   # mg/L
        "Ca":   546,    # mg/L
        "Mg":   1130,   # mg/L
        "K":    100,    # mg/L
        "Cl":   10645,  # mg/L
        "SO4":  1057,   # mg/L
        "CO3":  0,      # mg/L
        "HCO3": 154,    # mg/L
        "NO3":  0       # mg/L
    },
    "temperature": 21,  # celcius
    "pH": 7.2,
    "water type": "SW"  # "SW": seawater well or "BW": brackish well
}

membrane = {
    "surface area": 37,         # m2
    "module diameter": 0.1,     # m
    "spacer thickness": 0.0007, # m
    "module length": 1,         # m
    "a1": 0.916,
    "a2": -0.566,
    "a3": -1.790,
    "b1": 1.316,
    "b2": 0,
    "b3": 5.709,
    "lw0": 128.442,
    "ls0": 0.00486,
    "x": 12,
    "y": 0,
    "membrane type": "SWRO"     # "SWRO": seawater RO or "BWRO": brackish RO or "NF": nanofiltration
}

system = {
    "module recovery": 0.15, 
    "system recovery": 0.15,
    "flow rate": 25,            # m3/h
    "erd efficiency": 0.8,
    "pump efficiency": 0.8            
}