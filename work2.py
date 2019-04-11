import data
import formula

input_water = {
    "concentration": {
        "Na": 340,
        "Ca": 223,
        "Mg": 107,
        "K": 7,
        "Cl": 702,
        "SO4": 176,
        "CO3": 0,
        "HCO3": 242,
        "NO3": 416},
    "temperature": 29,
    "pH": 7.1,
    "module recovery": 0.19,
    "system recovery": 0.20,
    "flow rate": 25,
    "water type": "BW"}

input_membrane = {
    "surface area": 7.9,
    "module diameter": 0.2,
    "spacer thickness": 0.0008,
    "module length": 1,
    "a1": 0.802,
    "a2": -0.549,
    "a3": 0,
    "b1": 0.970,
    "b2": 0.432,
    "b3": 1.212,
    "lw0": 24.478,
    "ls0": 1.879,
    "x": 0.000233,
    "y": 1.735,
    "membrane type": "BW"
}



input_water["concentration"] = formula.do_phreeqc_input(list(input_water["concentration"].values()), input_water["temperature"], input_water["pH"])[0]
co2 = formula.do_phreeqc_input(input_water["concentration"], input_water["temperature"], input_water["pH"])[2]

blend_fr = formula.calculate_blending_fraction(input_water["module recovery"], input_water["system recovery"])
feed_flow = formula.do_system_design(input_water["flow rate"], formula.max_flux[input_water["water type"]], input_membrane["surface area"], input_water["module recovery"], blend_fr)[0]
no_of_pass = formula.do_system_design(input_water["flow rate"], formula.max_flux[input_water["water type"]], input_membrane["surface area"], input_water["module recovery"], blend_fr)[1]

feed1 = data.Feed(input_water["concentration"], input_water["temperature"], input_water["pH"], input_water["module recovery"], feed_flow)
permeate1 = data.Permeate()
concentrate1 = data.Concentrate()
membrane1 = data.Membrane(input_membrane["surface area"], input_membrane["module diameter"], input_membrane["spacer thickness"], input_membrane["module length"], \
    input_membrane["a1"], input_membrane["a2"], input_membrane["a3"], input_membrane["b1"], input_membrane["b2"], input_membrane["b3"], \
        input_membrane["lw0"], input_membrane["ls0"], input_membrane["x"], input_membrane["y"])

permeate1.fl_m3s = formula.calculate_flow(feed1.fl_m3s, feed1.mod_rec)[0]
concentrate1.fl_m3s = formula.calculate_flow(feed1.fl_m3s, feed1.mod_rec)[1]

rey = formula.calculate_reynolds_number(feed1.fl_m3s, concentrate1.fl_m3s, feed1.kvis, membrane1.thik, membrane1.effa)
scm = formula.calculate_schmidt_number(feed1.kvis, feed1.diff)
k = formula.calculate_transport_coefficient(feed1.diff, rey, scm, membrane1.thik)
flux = formula.calculate_permeate_flux(permeate1.fl_m3s, membrane1.area)
cpf = formula.calculate_cp_factor(flux, k)
fric = formula.calculate_friction_drop(rey, membrane1.x, membrane1.y)
lw = formula.calculate_solvent_permeabilty(membrane1.lw0, feed1.temp, feed1.tds_eq, cpf, membrane1.a1, membrane1.a2, membrane1.a3)
ls = formula.calculate_solute_permeability(membrane1.ls0, feed1.temp, feed1.tds_eq, cpf, membrane1.b1, membrane1.b2, membrane1.b3)


# need shorten!!!!!!!!!!!!!!!!
permeate1.conc_eq = formula.calculate_permeate_species(feed1.conc_eq, formula.rjec_coef[input_membrane["membrane type"]], ls, flux)[0]
permeate1.conc_mg = formula.calculate_permeate_species(feed1.conc_eq, formula.rjec_coef[input_membrane["membrane type"]], ls, flux)[1]
permeate1.conc_mol = formula.calculate_permeate_species(feed1.conc_eq, formula.rjec_coef[input_membrane["membrane type"]], ls, flux)[2]
permeate1.tds_eq = sum(permeate1.conc_eq)
permeate1.tds_mg = sum(permeate1.conc_mg)
permeate1.tds_mol = sum(permeate1.conc_mol)
concentrate1.conc_eq = formula.calculate_concentrate_species(feed1.fl_m3s, permeate1.fl_m3s, concentrate1.fl_m3s, feed1.conc_eq, permeate1.conc_eq)
concentrate1.conc_mg = formula.calculate_concentrate_species(feed1.fl_m3s, permeate1.fl_m3s, concentrate1.fl_m3s, feed1.conc_mg, permeate1.conc_mg)
concentrate1.conc_mol = formula.calculate_concentrate_species(feed1.fl_m3s, permeate1.fl_m3s, concentrate1.fl_m3s, feed1.conc_mol, permeate1.conc_mol)
concentrate1.tds_eq = sum(concentrate1.conc_eq)
concentrate1.tds_mg = sum(concentrate1.conc_mg)
concentrate1.tds_mol = sum(concentrate1.conc_mol)

tmp = formula.calculate_transmembrane_pressure(feed1.tds_mol, concentrate1.tds_mol, permeate1.tds_mol, feed1.temp, feed1.osm_coef, lw, flux) 
feed1.pres = tmp + fric/2
feed1.si = formula.do_phreeqc_input(input_water["concentration"], input_water["temperature"], input_water["pH"])[1]
permeate1.temp = concentrate1.temp = feed1.temp
permeate1.pH = formula.do_phreeqc_output(permeate1.conc_mg, permeate1.temp, co2)[0]
concentrate1.pH = formula.do_phreeqc_output(concentrate1.conc_mg, concentrate1.temp, co2)[0]
permeate1.si = formula.do_phreeqc_output(permeate1.conc_mg, feed1.temp, co2)[1]
concentrate1.si = formula.do_phreeqc_output(concentrate1.conc_mg, feed1.temp, co2)[1]
concentrate1.istr = formula.calculate_ionic_strength(concentrate1.conc_eq)
permeate1.istr = formula.calculate_ionic_strength(permeate1.conc_eq)
feed1.scale = formula.calculate_scaling_index(feed1.conc_mol[1], feed1.tds_mg, feed1.conc_mol[7], feed1.conc_mol[6], feed1.pH, feed1.temp, feed1.istr)
concentrate1.scale = formula.calculate_scaling_index(concentrate1.conc_mol[1], concentrate1.tds_mg, concentrate1.conc_mol[7], concentrate1.conc_mol[6], \
    concentrate1.pH, concentrate1.temp, concentrate1.istr)
permeate1.scale = formula.calculate_scaling_index(permeate1.conc_mol[1], permeate1.tds_mg, permeate1.conc_mol[7], permeate1.conc_mol[6], \
    permeate1.pH, permeate1.temp, permeate1.istr)
sr = formula.calculate_salt_rejection(feed1.tds_mg, permeate1.tds_mg)
sec = formula.calculate_energy_consumption(blend_fr, formula.erd_eff, formula.pump_eff, feed1.mod_rec, tmp)

out1 = formula.do_phreeqc_blending(feed1.conc_mg, feed1.temp, feed1.pH, permeate1.conc_mg, permeate1.temp, permeate1.pH, blend_fr, feed1.mod_rec)[0]


out1 = {
    "Na": out1[0],
    "Ca": out1[1],
    "Mg": out1[2],
    "K": out1[3],
    "Cl": out1[4],
    "SO4": out1[5],
    "CO3": out1[6],
    "HCO3": out1[7],
    "NO3": out1[8],
    "TDS": sum(out1)
}

out1 = {x: round(out1[x], 1) for x in out1}
    
out2 = round(tmp, 2)

out3 = {x: round(concentrate1.si[x], 3) for x in concentrate1.si}

out4 = round(sec*no_of_pass, 2)

out5 = round(flux*3600, 1)

print("\nNumber of pass:")
print(no_of_pass)

print("\nPermeate flux:")
print(out5, "LMH")

print("\nPermeate concentration:")
for x in out1:
    print(x, out1[x], "mg/L")

print("\nTransmembrane pressure:", out2, "bar")

print("\nSpecific energy consumption:", out4, "kWs/m3")

print("\nConcentrate SI:")
for x in out3:
    if out3[x] > 0:
        print(x, out3[x])
print("LSI", concentrate1.scale[0])
print("SDSI", concentrate1.scale[1])
print("\n")