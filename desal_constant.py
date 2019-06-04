# equivalent weight
eq_w = {
    "Na": 22.99,
    "Ca": 20.04,
    "Mg": 12.16,
    "K": 39.10,
    "Cl": 35.45,
    "SO4": 48.03,
    "CO3": 30.00,
    "HCO3": 61.02,
    "NO3": 62.00,
}

# ion valency
val = {
    "Na": 1,
    "Ca": 2,
    "Mg": 2,
    "K": 1,
    "Cl": 1,
    "SO4": 2,
    "CO3": 2,
    "HCO3": 1,
    "NO3": 1,
}

# diffusion coefficient
diff_coef = {
    "Na": 1.334,
    "Ca": 0.792,
    "Mg": 0.706,
    "K": 1.957,
    "Cl": 2.032,
    "SO4": 1.065,
    "CO3": 0.923,
    "HCO3": 1.185,
    "NO3": 1.902,
}

# rejection coefficient
rjec_coef = { 
    "BWRO": {
        "Na": 1.00,
        "Ca": 0.25,
        "Mg": 0.25,
        "K": 2.50,
        "Cl": 1.00,
        "SO4": 0.65,
        "CO3": 0.005,
        "HCO3": 0.90,
        "NO3": 1.20,
    },
    "SWRO": {
        "Na": 1.00,
        "Ca": 0.65,
        "Mg": 0.65,
        "K": 4.50,
        "Cl": 1.00,
        "SO4": 1.50,
        "CO3": 0.005,
        "HCO3": 0.50,
        "NO3": 0.80,
    },
    "NF": {
        "Na": 1.00,
        "Ca": 0.06,
        "Mg": 0.06,
        "K": 2.50,
        "Cl": 1.00,
        "SO4": 0.35,
        "CO3": 0.005,
        "HCO3": 0.25,
        "NO3": 4.00,
    },
}

# design guideline: maximum flux
max_flux = {
    "BW": 0.034,
    "SW": 0.020,
}