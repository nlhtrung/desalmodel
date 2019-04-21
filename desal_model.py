import desal_data
import desal_formula
import desal_input
import desal_constant

# fetch data from input file
input_water = desal_input.water
input_membrane = desal_input.membrane
input_system = desal_input.system

# water analysis with phreeqc
input_water["concentration"] = desal_formula.do_phreeqc_input(list(input_water["concentration"].values()), input_water["temperature"], input_water["pH"])[0]
co2 = desal_formula.do_phreeqc_input(input_water["concentration"], input_water["temperature"], input_water["pH"])[2]

# system design
blend_fr = desal_formula.calculate_blending_fraction(input_system["module recovery"], input_system["system recovery"])
feed_flow = desal_formula.do_system_design(input_system["flow rate"], desal_constant.max_flux[input_water["water type"]], input_membrane["surface area"], \
    input_system["module recovery"], blend_fr)[0]
no_of_pass = desal_formula.do_system_design(input_system["flow rate"], desal_constant.max_flux[input_water["water type"]], input_membrane["surface area"], \
    input_system["module recovery"], blend_fr)[1]

# create main objects feed, permeate, concentrate, membrane
feed1 = desal_data.Feed(input_water["concentration"], input_water["temperature"], input_water["pH"], input_system["module recovery"], feed_flow)
permeate1 = desal_data.Permeate()
concentrate1 = desal_data.Concentrate()
membrane1 = desal_data.Membrane(input_membrane["surface area"], input_membrane["module diameter"], input_membrane["spacer thickness"], input_membrane["module length"], \
    input_membrane["a1"], input_membrane["a2"], input_membrane["a3"], input_membrane["b1"], input_membrane["b2"], input_membrane["b3"], \
        input_membrane["lw0"], input_membrane["ls0"], input_membrane["x"], input_membrane["y"])

# flow balance
permeate1.fl_m3s = desal_formula.calculate_flow(feed1.fl_m3s, feed1.mod_rec)[0]
concentrate1.fl_m3s = desal_formula.calculate_flow(feed1.fl_m3s, feed1.mod_rec)[1]

# fluid characteristcs
rey = desal_formula.calculate_reynolds_number(feed1.fl_m3s, concentrate1.fl_m3s, feed1.kvis, membrane1.thik, membrane1.effa)
scm = desal_formula.calculate_schmidt_number(feed1.kvis, feed1.diff)
k = desal_formula.calculate_transport_coefficient(feed1.diff, rey, scm, membrane1.thik)
flux = desal_formula.calculate_permeate_flux(permeate1.fl_m3s, membrane1.area)
cpf = desal_formula.calculate_cp_factor(flux, k)
fric = desal_formula.calculate_friction_drop(rey, membrane1.x, membrane1.y)

# permeability
lw = desal_formula.calculate_solvent_permeabilty(membrane1.lw0, feed1.temp, feed1.tds_eq, cpf, membrane1.a1, membrane1.a2, membrane1.a3)
ls = desal_formula.calculate_solute_permeability(membrane1.ls0, feed1.temp, feed1.tds_eq, cpf, membrane1.b1, membrane1.b2, membrane1.b3)

# concentrations estimation
permeate1.conc_eq = desal_formula.calculate_permeate_species(feed1.conc_eq, desal_constant.rjec_coef[input_membrane["membrane type"]], ls, flux)[0]
permeate1.conc_mg = desal_formula.calculate_permeate_species(feed1.conc_eq, desal_constant.rjec_coef[input_membrane["membrane type"]], ls, flux)[1]
permeate1.conc_mol = desal_formula.calculate_permeate_species(feed1.conc_eq, desal_constant.rjec_coef[input_membrane["membrane type"]], ls, flux)[2]
permeate1.tds_eq = sum(permeate1.conc_eq)
permeate1.tds_mg = sum(permeate1.conc_mg)
permeate1.tds_mol = sum(permeate1.conc_mol)
concentrate1.conc_eq = desal_formula.calculate_concentrate_species(feed1.fl_m3s, permeate1.fl_m3s, concentrate1.fl_m3s, feed1.conc_eq, permeate1.conc_eq)
concentrate1.conc_mg = desal_formula.calculate_concentrate_species(feed1.fl_m3s, permeate1.fl_m3s, concentrate1.fl_m3s, feed1.conc_mg, permeate1.conc_mg)
concentrate1.conc_mol = desal_formula.calculate_concentrate_species(feed1.fl_m3s, permeate1.fl_m3s, concentrate1.fl_m3s, feed1.conc_mol, permeate1.conc_mol)
concentrate1.tds_eq = sum(concentrate1.conc_eq)
concentrate1.tds_mg = sum(concentrate1.conc_mg)
concentrate1.tds_mol = sum(concentrate1.conc_mol)

# pressures
tmp = desal_formula.calculate_transmembrane_pressure(feed1.tds_mol, concentrate1.tds_mol, permeate1.tds_mol, feed1.temp, feed1.osm_coef, lw, flux) 
feed1.pres = tmp + fric/2

# temperature, pH
permeate1.temp = concentrate1.temp = feed1.temp
permeate1.pH = desal_formula.do_phreeqc_output(permeate1.conc_mg, permeate1.temp, co2)[0]
concentrate1.pH = desal_formula.do_phreeqc_output(concentrate1.conc_mg, concentrate1.temp, co2)[0]

# scaling indices
feed1.si = desal_formula.do_phreeqc_input(input_water["concentration"], input_water["temperature"], input_water["pH"])[1]

permeate1.si = desal_formula.do_phreeqc_output(permeate1.conc_mg, permeate1.temp, co2)[1]
concentrate1.si = desal_formula.do_phreeqc_output(concentrate1.conc_mg, concentrate1.temp, co2)[1]
concentrate1.istr = desal_formula.calculate_ionic_strength(concentrate1.conc_eq)
permeate1.istr = desal_formula.calculate_ionic_strength(permeate1.conc_eq)
feed1.scale = desal_formula.calculate_scaling_index(feed1.conc_mol[1], feed1.tds_mg, feed1.conc_mol[7], feed1.conc_mol[6], feed1.pH, feed1.temp, feed1.istr)
concentrate1.scale = desal_formula.calculate_scaling_index(concentrate1.conc_mol[1], concentrate1.tds_mg, concentrate1.conc_mol[7], concentrate1.conc_mol[6], \
    concentrate1.pH, concentrate1.temp, concentrate1.istr)
permeate1.scale = desal_formula.calculate_scaling_index(permeate1.conc_mol[1], permeate1.tds_mg, permeate1.conc_mol[7], permeate1.conc_mol[6], \
    permeate1.pH, permeate1.temp, permeate1.istr)

# salt rejection
sr = desal_formula.calculate_salt_rejection(feed1.tds_mg, permeate1.tds_mg)

# specific energy consumption
sec = desal_formula.calculate_energy_consumption(blend_fr, input_system["erd efficiency"], input_system["pump efficiency"], input_system["module recovery"], tmp)

# output blending
output_water = desal_formula.do_phreeqc_blending(feed1.conc_mg, feed1.temp, feed1.pH, permeate1.conc_mg, permeate1.temp, permeate1.pH, \
    blend_fr, input_system["module recovery"])[0]

