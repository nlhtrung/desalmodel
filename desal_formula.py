import math
import operator
import phreeqpython # https://github.com/Vitens/phreeqpython
import desal_constant


# system design #
def do_system_design(flow, flux, area, rec, blend_fr):
    el_flow = flux*area
    flow = flow*blend_fr
    no_of_pass = math.ceil(flow*rec/el_flow)
    flow = flow/no_of_pass
    return (flow, no_of_pass)
# end #

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
def calculate_scaling_index(ca, tds, hco3, co3, pH, temp, istr):
    alk = hco3 + 2*co3
    if alk != 0:
        pca = -math.log10(ca/1000)
        palk = -math.log10(alk/1000)
        temp = temp*9/5 + 32
        c = 3.26*math.exp(-0.005*temp) - 0.0116*math.log10(tds**3) + 0.0905*math.log10(tds**2) - 0.133*math.log10(tds) - 0.02
        lsi = round(pH - (pca + palk + c), 3)
        if istr < 1.2:
            k = 2.022*math.exp((math.log(istr) + 7.544)**2/102.60) - 0.0002*temp**2 + 0.00097*temp + 0.262
        else:
            k = -0.1*istr - 0.0002*temp**2 -0.00097*temp + 3.887
        sdsi = round(pH - (pca + palk + k), 3)
    else:
        lsi = "n/a"
        sdsi = "n/a"   
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
    a4 = 6.92*10**-4*g1 - 8.7*10**-5*g2 - 5.3*10**-5*g3
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
    return sum(list(map(operator.mul, conc, list(desal_constant.val.values()))))/2/1000

# diffusivity (cation name, anion name)
def calculate_diffusivity(cation, anion):
    return (desal_constant.val[cation] + desal_constant.val[anion])*desal_constant.diff_coef[cation]*desal_constant.diff_coef[anion]/ \
        (desal_constant.val[cation]*desal_constant.diff_coef[cation] + desal_constant.val[anion]*desal_constant.diff_coef[anion])

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
    return 0.023*rey**0.875*scm**0.25*diff*10**-9/(thik*2)

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

# transmembrane pressure (feed/concentrate/permeate concentration, temperature, osmotic coefficient, solvent permeability, permeate flux)
def calculate_transmembrane_pressure(conc_f, conc_c, conc_p, temp, osm_coef, lw, flux):
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
    c = list(map(operator.mul, [solute*x*b/(flux*3600 + solute) for x in conc], a))
    # return permeate concentration in eq/m3, mg/L and mol/m3
    return (c, list(map(operator.mul, c, list(desal_constant.eq_w.values()))), list(map(operator.truediv, c, list(desal_constant.eq_w.values()))))    
# end #

# salt rejection (feed TDS, permeate TDS)
def calculate_salt_rejection(feed, permeate):
    sr = 1 - permeate/feed
    return sr

# blending fraction (module recovery, system recovery)
def calculate_blending_fraction(mod_rec, sys_rec):
    blend_fr = (1 - sys_rec)/(1 - mod_rec)
    return blend_fr

# specific energy consumption (blending fraction, ERD efficiency, module recovery, system recovery, TMP)
def calculate_energy_consumption(blend_fr, erd_eff, pump_eff, mod_rec, tmp):
    sec = blend_fr*(1 - erd_eff*(1 - mod_rec))/(1 - blend_fr*(1 - mod_rec))/pump_eff*tmp/36
    return sec

# solution blending #
def do_phreeqc_blending(conc_f, temp_f, pH_f, conc_p, temp_p, pH_p, blend_fr, mod_rec):
    pp = phreeqpython.PhreeqPython()
    solution_f = pp.add_solution({
        "units": "mg/L",
        "pH": pH_f,
        "temp": temp_f,
        "Na": conc_f[0],
        "Ca": conc_f[1],
        "Mg": conc_f[2],
        "K": conc_f[3],
        "Cl": str(conc_f[4]) + " as Cl",
        "S(6)": str(conc_f[5]) + " as SO4",
        "Alkalinity": str(conc_f[6]*2 + conc_f[7]) + " as HCO3",
        "N(5)": str(conc_f[8]) + " as NO3"
    })
    solution_p = pp.add_solution({
        "units": "mg/L",
        "pH": pH_p,
        "temp": temp_p,
        "Na": conc_p[0],
        "Ca": conc_p[1],
        "Mg": conc_p[2],
        "K": conc_p[3],
        "Cl": str(conc_p[4]) + " as Cl",
        "S(6)": str(conc_p[5]) + " as SO4",
        "Alkalinity": str(conc_p[6]*2 + conc_p[7]) + " as HCO3",
        "N(5)": str(conc_p[8]) + " as NO3"
    })
    a = (1 - blend_fr)/(1 - blend_fr + blend_fr*mod_rec)
    solution_b = solution_f*a + solution_p*(1 - a)
    conc_b = [
        solution_b.total("Na", units="mg"),
        solution_b.total("Ca", units="mg"),
        solution_b.total("Mg", units="mg"),
        solution_b.total("K", units="mg"),
        solution_b.total("Cl", units="mg"),
        solution_b.total("SO4", units="mg"),
        solution_b.total("CO3", units="mg"),
        solution_b.total("HCO3", units="mg"),
        solution_b.total("NO3", units="mg"),
    ]
    pH_b = solution_b.pH
    temp_b = solution_b.temperature
    return (conc_b, pH_b, temp_b)
# end #

