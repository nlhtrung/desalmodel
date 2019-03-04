import math
import operator
import phreeqpython # https://github.com/Vitens/phreeqpython

# equivalent weight
eq_w = {
    "Na": 22.99,
    "Ca": 20.04,
    "Mg": 12.16,
    "K": 39.1,
    "Cl": 35.45,
    "SO4": 48.03,
    "CO3": 30,
    "HCO3": 61.02,
    "NO3": 62,
}

# valency
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
    "Na": 1.33,
    "Ca": 0.793,
    "Mg": 0.705,
    "K": 1.96,
    "Cl": 2.03,
    "SO4": 1.07,
    "CO3": 0.955,
    "HCO3": 1.18,
    "NO3": 1.9,
}

# rejection coefficient
rjec_coef = {
    "Na": 1,
    "Ca": 0.35,
    "Mg": 0.35,
    "K": 1.5,
    "Cl": 1,
    "SO4": 0.2,
    "CO3": 0.05,
    "HCO3": 1,
    "NO3": 8,
}

# charge balance (concentration in eq/m3) #
def do_charge_balance(conc):
    c = 0
    for n in range (4): # cations
        c += conc[n]
    a = 0
    for n in range (4, 9): # anions
        a += conc[n]
    if a > c:
        conc[0] += a - c # add Na
    else:
        conc[4] += c - a # add Cl
    return conc
# end #

# phreeqc water analysis (concentration in eq/m3, temperature, pH, dissolved CO2) #
def do_phreeqc_input(conc, temp, pH):
    pp = phreeqpython.PhreeqPython()
    solution = pp.add_solution({
        "units": "mg/L",
        "pH": pH,
        "temp": temp,
        "Na": conc[0],
        "Ca": conc[1],
        "Mg": conc[2],
        "K": conc[3],
        "Cl": str(conc[4]) + " as Cl",
        "S(6)": str(conc[5]) + " as SO4",
        "Alkalinity": str(conc[7]) + " as HCO3",
        "N(5)": str(conc[8]) + " as NO3"
    })
    conc[6] = solution.total("CO3", units="mg")
    co2 = solution.total("CO2", units="mg")
    return (conc, solution.phases, co2)

def do_phreeqc_output(conc, temp, co2):
    pp = phreeqpython.PhreeqPython()
    solution = pp.add_solution({
        "units": "mg/L",
        "temp": temp,
        "Na": conc[0],
        "Ca": conc[1],
        "Mg": conc[2],
        "K": conc[3],
        "Cl": str(conc[4]) + " as Cl",
        "S(6)": str(conc[5]) + " as SO4",
        "Alkalinity": str(conc[6]*2 + conc[7]) + " as HCO3",
        "C(4)": str(conc[6]*12/30 + conc[7]*12/61 + co2*12/44) + " as C",
        "N(5)": str(conc[8]) + " as NO3"
    })
    pH = solution.pH
    return (pH, solution.phases)
# end #

# scaling index (Calcium concentration, TDS, Alkalinity, pH, temperature, ionic strength) # 
# https://www.usbr.gov/tsc/techreferences/mands/mands-pdfs/WQeval_documentation.pdf
def calculate_scaling_index(ca, tds, alk, pH, temp, istr):
    pca = -math.log10(ca/1000)
    palk = -math.log10(alk/1000)
    temp = temp*9/5 + 32
    c = 3.26*math.exp(-0.005*temp) - 0.0116*math.log10(tds**3) + 0.0905*math.log10(tds**2) - 0.133*math.log10(tds) - 0.02
    lsi = pH - (pca + palk + c)
    if istr < 1.2:
        k = 2.022*math.exp((math.log(istr) + 7.544)**2/102.60) - 0.0002*temp**2 + 0.00097*temp + 0.262
    else:
        k = -0.1*istr - 0.0002*temp**2 -0.00097*temp + 3.887
    sdsi = pH - (pca + palk + k)
    return (lsi, sdsi)
# end #

# density kg/m3 (temperature, TDS) #
def calculate_density(temp, tds):
    sal = tds/1000
    a = (2*temp - 200)/160
    b = (2*sal - 150)/150
    f1 = 0.5
    f2 = a
    f3 = 2*a**2 - 1
    f4 = 4*a**3 - 3*a
    g1 = 0.5
    g2 = b
    g3 = 2*b**2 -1
    a1 = 4.032219*g1 + 0.115313*g2 + 3.26*10**-4*g3
    a2 = -0.108199*g1 + 1.571*10**-3*g2 - 4.23*10**-4*g3
    a3 = -0.012247*g1 + 1.74*10**-3*g2 - 9*10**-6*g3
    a4 = 6.92*10**-4*g1 - 8.7*10**-5*g2 -5.3*10**-5*g3
    dens = (a1*f1 + a2*f2 + a3*f3 + a4*f4)*10**3
    return dens
# end #

# kinematic viscosity m2/s (temperature, TDS) #
def calculate_viscosity(temp, tds):
    sal = tds/1000
    a = 1.474*10**-3 + 1.5*10**-5*temp - 3.927*10**-8*temp**2
    b = 1.0734*10**-5 - 8.5*10**-8*temp + 2.23*10**-10*temp**2
    dvis_w = math.exp(-3.79418 + 604.129/(139.18 + temp))
    dvis_s = 1 + a*sal + b*sal**2
    dvis = dvis_w*dvis_s*10**-3
    kvis = dvis/calculate_density(temp, tds)
    return kvis
# end #

# osmotic coefficient (temperature, ionic strength) #
def calculate_osmotic_coefficient(temp, istr):
    a = 20.661 - 432.5797/(temp + 273) - 3.712*math.log(temp + 273) + 8.638*10**-3*(temp + 273)
    b = 2.303/istr**1.5*((1 + istr**0.5) - 1/(1 + istr**0.5) - 2*math.log(1 + istr**0.5))
    c = -831.659 + 17022.399/(temp + 273) + 157.653*math.log(temp + 273) - 0.493*(temp + 273) + 2.595*10**-4*(temp + 273)**2
    d = 553.905 - 11200.445/(temp + 273) - 105.239*math.log(temp + 273) + 0.333*(temp + 273) - 1.774*10**-4*(temp + 273)**2
    e = -0.15112
    return 1 - (a*b*istr**0.5 + c*istr + d*istr**1.5 + e*istr**2)
# end #

# ionic strength (concentration in eq/m3)
def calculate_ionic_strength(conc):
    return sum(list(map(operator.mul, conc, list(val.values()))))/2/1000

# diffusivity (cation name, anion name)
def calculate_diffusivity(cation, anion):
    return (val[cation] + val[anion])*diff_coef[cation]*diff_coef[anion]/ \
        (val[cation]*diff_coef[cation] + val[anion]*diff_coef[anion])

# concentration fraction (concentration)
def calculate_fraction(conc):
    return [x/sum(conc) for x in conc]

# flow (feed flow, recovery rate)
def calculate_flow(feed, rec):
    return (feed*rec, feed*(1-rec))

# concentrate species (feed/permeate/concentrate flow, feed/permeate concentration)
def calculate_concentrate_species(flow_f, flow_p, flow_c, conc_f, conc_p):
    f = (flow_f*x for x in conc_f)
    p = (flow_p*x for x in conc_p)
    c = list(map(operator.sub, f, p))
    return [x/flow_c for x in c]

# reynolds (feed/concentrate flow, kinematic viscosity, spacer thickness, effective area)
def calculate_reynolds_number(feed, conc, kvis, thik, effa):
    return (feed + conc)/2/effa*thik*2/kvis

# schmidt (kinematic viscosity, diffusivity)
def calculate_schmidt_number(kvis, diff):
    return kvis/diff*10**9

# transport coefficient (diffusivity, reynolds, schmidt number, spacer thickness)
def calculate_transport_coefficient(diff, rey, scm, thik):
    return 0.065*rey**0.875*scm**0.25*diff*10**-9/(thik*2)

# flux (permeate flow, total surface area)
def calculate_permeate_flux(flow, area):
    return flow/area*1000

# concentration polarisation factor (permeate flux, transport coefficient)
def calculate_cp_factor(flux, k):
    return math.exp(flux/k/1000)

# solvent permeability (membrane intrinsic permeabilty, temperature, feed concentration, cp factor, model parameters)
def calculate_solvent_permeabilty(lw0, temp, conc, cpf, a1, a2, a3):
    return lw0*math.exp(a1*math.log(temp/273) + a2*conc/1000 + a3*cpf)

# solute permeability (membrane intrinsic permeability, temperature, feed concentration, cp factor, model parameters)
def calculate_solute_permeability(ls0, temp, conc, cpf, b1, b2, b3):
    return ls0*math.exp(b1*math.log(temp/273) + b2*math.log(conc/1000) + b3*cpf)

# friction drop (reynolds number, model parameters)
def calculate_friction_drop(rey, x, y):
    return x*rey**y/100

# concentrate pressure (feed/concentrate/permeate concentration, temperature, osmotic coefficient, solvent permeability, permeate flux)
def calculate_concentrate_pressure(conc_f, conc_c, conc_p, temp, osm_coef, lw, flux):
    dpi = ((conc_f + conc_c)/2 - conc_p)*8.3144598*(temp + 273)*osm_coef*10**-5
    return flux*3600/lw + dpi

# permeate species (feed concentration, rejection coefficient, solute permeability, permeate flux) #
def calculate_permeate_species(conc, rjec_coef, solute, flux):
    fr_eq = calculate_fraction(conc)
    rjec_coef = list(rjec_coef.values())
    a = []
    i = 0
    for x in fr_eq:
        if x == 0:
            a.append(0)
        elif i < 4: # cations
            z = 0
            for n in range(4): 
                z += fr_eq[n]*rjec_coef[n]
            y = rjec_coef[i]/z/2
            a.append(y)
        else: # anions
            z = 0
            for n in range(4, 9): 
                z += fr_eq[n]*rjec_coef[n]
            y = rjec_coef[i]/z/2
            a.append(y)
        i += 1
    b = sum(list(map(operator.mul, fr_eq, rjec_coef)))  
    c = list(map(operator.mul, [solute*sum(conc)*x*b/(flux*3600) for x in fr_eq], a))
    # return permeate concentration in eq/m3, mg/L and mol/m3
    return (c, list(map(operator.mul, c, list(eq_w.values()))), list(map(operator.truediv, c, list(val.values()))))    
# end #


