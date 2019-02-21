import operator
import formula



# define feed (concentrations, TDS, ionic strength, osmotic coefficient, diffusivity, temperature, pH, recovery, flow)
class feed:
    def __init__ (self, conc_mg, temp, pH, rec, fl_m3h):
        self.conc_mg = conc_mg
        conc_eq = list(map(operator.truediv, conc_mg, list(formula.eq_w.values()))) # calculate eq concentration
        conc_eq = formula.balance(conc_eq) # ion balance
        self.conc_eq = conc_eq
        conc_mol = list(map(operator.truediv, conc_eq, list(formula.val.values()))) # calculate mol concentration
        self.conc_mol = conc_mol
        tds_mg = sum(conc_mg) # calculate tds mg/L
        self.tds_mg = tds_mg
        tds_eq = sum(conc_eq) # calculate tds eq/m3
        self.tds_eq = tds_eq
        tds_mol = sum(conc_mol) # calculate tds mol/m3
        self.tds_mol = tds_mol 
        istr = formula.ionic(conc_eq) # calculate ionic strength
        self.istr = istr
        osm_coef = formula.osmotic(temp, istr) # calulate osmotic coefficient
        self.osm_coef = osm_coef
        fr_mol = formula.fraction(conc_mol) # calculate mol fraction
        self.fr_mol = fr_mol
        fr_eq = formula.fraction(conc_eq) # calculate equivalent fraction
        self.fr_eq = fr_eq
        diff = sum(list(map(operator.mul, fr_mol, list(formula.diff_coef.values())))) # calculate diffusivity
        self.diff = diff
        self.temp = temp
        self.pH = pH
        self.rec = rec
        self.fl_m3h = fl_m3h
        fl_m3s = fl_m3h/3600 # calculate flow m3/s
        self.fl_m3s = fl_m3s
        self.dens = formula.density(temp, tds_mg)
        self.kvis = formula.viscosity(temp, tds_mg)
        self.pres = 0
        self.si = 0
        self.foul = []




# define membrane (area, module diameter, spacer thickness, modul length, effective area, a1,2,3, b1,2,3, Lw0, Ls0, x, y)
class membrane:
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

class permeate:
    def __init__ (self):
        pass


class concentrate:
    def __init__ (self):
        pass