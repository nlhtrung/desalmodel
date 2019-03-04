import data
import formula

in1 = [{
    "Na": 82.76,
    "Ca": 176.41,
    "Mg": 33.31,
    "K": 1.75,
    "Cl": 269.66,
    "SO4": 95.97,
    "CO3": 0,
    "HCO3": 116.49,
    "NO3": 227.56},
    25,
    7.4,
    0.2,
    1]

membrane1 = data.Membrane(7.9, 0.1, 0.0008, 0.964, \
    0.782360338895866, -0.717007003850838, -0.275798741964664, \
        0.969656473719167, 0.424364559615452, 4.46574778569048, \
            31.0518275461852, 0.06918357125442, \
                0.000203655002873562, 1.74547844949712)




in1[0] = formula.do_phreeqc_input(list(in1[0].values()), in1[1], in1[2])[0]
co2 = formula.do_phreeqc_input(in1[0], in1[1], in1[2])[2]

feed1 = data.Feed(in1[0], in1[1], in1[2], in1[3], in1[4])
permeate1 = data.Permeate()
concentrate1 = data.Concentrate()

permeate1.fl_m3s = formula.calculate_flow(feed1.fl_m3s, feed1.rec)[0]
concentrate1.fl_m3s = formula.calculate_flow(feed1.fl_m3s, feed1.rec)[1]

rey = formula.calculate_reynolds_number(feed1.fl_m3s, concentrate1.fl_m3s, feed1.kvis, membrane1.thik, membrane1.effa)
scm = formula.calculate_schmidt_number(feed1.kvis, feed1.diff)
k = formula.calculate_transport_coefficient(feed1.diff, rey, scm, membrane1.thik)
flux = formula.calculate_permeate_flux(permeate1.fl_m3s, membrane1.area)
cpf = formula.calculate_cp_factor(flux, k)
fric = formula.calculate_friction_drop(rey, membrane1.x, membrane1.y)
lw = formula.calculate_solvent_permeabilty(membrane1.lw0, feed1.temp, feed1.tds_eq, cpf, membrane1.a1, membrane1.a2, membrane1.a3)
ls = formula.calculate_solute_permeability(membrane1.ls0, feed1.temp, feed1.tds_eq, cpf, membrane1.b1, membrane1.b2, membrane1.b3)


# need shorten!!!!!!!!!!!!!!!!
permeate1.conc_eq = formula.calculate_permeate_species(feed1.conc_eq, formula.rjec_coef, ls, flux)[0]
permeate1.conc_mg = formula.calculate_permeate_species(feed1.conc_eq, formula.rjec_coef, ls, flux)[1]
permeate1.conc_mol = formula.calculate_permeate_species(feed1.conc_eq, formula.rjec_coef, ls, flux)[2]
permeate1.tds_eq = sum(permeate1.conc_eq)
permeate1.tds_mg = sum(permeate1.conc_mg)
permeate1.tds_mol = sum(permeate1.conc_mol)
concentrate1.conc_eq = formula.calculate_concentrate_species(feed1.fl_m3s, permeate1.fl_m3s, concentrate1.fl_m3s, feed1.conc_eq, permeate1.conc_eq)
concentrate1.conc_mg = formula.calculate_concentrate_species(feed1.fl_m3s, permeate1.fl_m3s, concentrate1.fl_m3s, feed1.conc_mg, permeate1.conc_mg)
concentrate1.conc_mol = formula.calculate_concentrate_species(feed1.fl_m3s, permeate1.fl_m3s, concentrate1.fl_m3s, feed1.conc_mol, permeate1.conc_mol)
concentrate1.tds_eq = sum(concentrate1.conc_eq)
concentrate1.tds_mg = sum(concentrate1.conc_mg)
concentrate1.tds_mol = sum(concentrate1.conc_mol)

concentrate1.pres = formula.calculate_concentrate_pressure(feed1.tds_mol, concentrate1.tds_mol, permeate1.tds_mol, feed1.temp, feed1.osm_coef, lw, flux) 
feed1.pres = concentrate1.pres + fric
feed1.si = formula.do_phreeqc_input(in1[0], in1[1], in1[2])[1]
permeate1.temp = concentrate1.temp = feed1.temp
permeate1.pH = formula.do_phreeqc_output(permeate1.conc_mg, permeate1.temp, co2)[0]
concentrate1.pH = formula.do_phreeqc_output(concentrate1.conc_mg, concentrate1.temp, co2)[0]
permeate1.si = formula.do_phreeqc_output(permeate1.conc_mg, feed1.temp, co2)[1]
concentrate1.si = formula.do_phreeqc_output(concentrate1.conc_mg, feed1.temp, co2)[1]
concentrate1.istr = formula.calculate_ionic_strength(concentrate1.conc_eq)
permeate1.istr = formula.calculate_ionic_strength(permeate1.conc_eq)
feed1.foul = formula.calculate_scaling_index(feed1.conc_mol[1], feed1.tds_mg, feed1.conc_mol[7], feed1.pH, feed1.temp, feed1.istr)
concentrate1.foul = formula.calculate_scaling_index(concentrate1.conc_mol[1], concentrate1.tds_mg, concentrate1.conc_mol[7], concentrate1.pH, concentrate1.temp, concentrate1.istr)
permeate1.foul = formula.calculate_scaling_index(permeate1.conc_mol[1], permeate1.tds_mg, permeate1.conc_mol[7], permeate1.pH, permeate1.temp, permeate1.istr)

out1 = {
    "Na": permeate1.conc_mg[0],
    "Ca": permeate1.conc_mg[1],
    "Mg": permeate1.conc_mg[2],
    "K": permeate1.conc_mg[3],
    "Cl": permeate1.conc_mg[4],
    "SO4": permeate1.conc_mg[5],
    "CO3": permeate1.conc_mg[6],
    "HCO3": permeate1.conc_mg[7],
    "NO3": permeate1.conc_mg[8],
}

out2 = feed1.pres

print("\nPermeate concentration:")
for x in out1:
    print(x, round(out1[x], 4), "mg/L")

print("\nFeed pressure:", round(out2*100, 2), "kPa")

print("\nConcentrate SI:")
for x in concentrate1.si:
    if concentrate1.si[x] > 0:
        print(x, round(concentrate1.si[x], 4))
print("LSI", round(concentrate1.foul[0], 4))
print("SDSI", round(concentrate1.foul[1], 4))
print("\n")