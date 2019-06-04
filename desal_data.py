import operator
import desal_formula
import desal_constant


# define feed (concentrations, temperature, pH, module recovery, system recovery, flow rate in m3/h)
class Feed:
    def __init__ (self, conc_mg, temp, pH, mod_rec, fl_m3h):   
        conc_eq = list(map(operator.truediv, conc_mg, list(desal_constant.eq_w.values()))) # calculate eq concentration
        conc_eq = desal_formula.do_charge_balance(conc_eq) # ion balance
        self.conc_eq = conc_eq
        conc_mg = list(map(operator.mul, conc_eq, list(desal_constant.eq_w.values()))) # calculate mass concentration
        self.conc_mg = conc_mg
        conc_mol = list(map(operator.truediv, conc_eq, list(desal_constant.val.values()))) # calculate mol concentration
        self.conc_mol = conc_mol
        tds_mg = sum(conc_mg) # calculate tds mg/L
        self.tds_mg = tds_mg
        tds_eq = sum(conc_eq) # calculate tds eq/m3
        self.tds_eq = tds_eq
        tds_mol = sum(conc_mol) # calculate tds mol/m3
        self.tds_mol = tds_mol 
        istr = desal_formula.calculate_ionic_strength(conc_eq) # calculate ionic strength
        self.istr = istr
        osm_coef = desal_formula.calculate_osmotic_coefficient(temp, istr) # calulate osmotic coefficient
        self.osm_coef = osm_coef
        fr_mol = desal_formula.calculate_fraction(conc_mol) # calculate mol fraction
        self.fr_mol = fr_mol
        fr_eq = desal_formula.calculate_fraction(conc_eq) # calculate equivalent fraction
        self.fr_eq = fr_eq
        diff = desal_formula.calculate_diffusivity("Na", "Cl") # calculate diffusivity
        self.diff = diff
        self.temp = temp
        self.pH = pH
        self.mod_rec = mod_rec
        self.fl_m3h = fl_m3h
        fl_m3s = fl_m3h/3600 # calculate flow m3/s
        self.fl_m3s = fl_m3s
        self.dens = desal_formula.calculate_density(temp, tds_mg) # calculate water density
        self.kvis = desal_formula.calculate_viscosity(temp, tds_mg) # calculate water viscosity
        self.pres = 0
        self.si = 0
        self.scale = []

# define membrane (area, module diameter, spacer thickness, module length, effective area, a1,2,3, b1,2,3, Lw0, Ls0, x, y)
class Membrane:
    def __init__ (self, area, diam, thik, leng, a1, a2, a3, b1, b2, b3, lw0, ls0, x, y):
        self.area = area
        self.diam = diam
        self.thik = thik
        self.leng = leng
        effa = thik*leng # calculate effective area
        self.effa = effa
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3
        self.lw0 = lw0
        self.ls0 = ls0
        self.x = x
        self.y = y

# define permeate (placeholder)
class Permeate:
    def __init__ (self):
        pass

# define permeate (placeholder)
class Concentrate:
    def __init__ (self):
        pass