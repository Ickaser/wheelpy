#MUC: modules, units, constants (and properties, now)
import pint
import numpy as np

uReg = pint.UnitRegistry()
uReg.define("ppm=mg/kg")
uReg.define("cm1 = 1/cm")
uReg.ndim = uReg.dimensionless

c = 2.998E8 * uReg.m/uReg.s
h = 6.626E-34 *uReg.m*uReg.m/uReg.s*uReg.kg
h_ = h/2/np.pi
kB = 1.38064852E-23 *uReg.J/uReg.K
NA = 6.02214086E23 / uReg.mol
R = kB * NA
m_e = 9.10938356E-31 * uReg.kg
m_p = 1*uReg.amu
ϵ0 = 8.8541878128E-12 *uReg.C**2/uReg.J/uReg.m
e = 1.602176634E-19*uReg.C
sigma = 5.670374419E-8*uReg.W/uReg.m**2/uReg.K**4

def CtoK(temp):
    return uReg.Quantity(temp, uReg.degC).to(uReg.K)

def wgtInterp(x, xa, xb, ya, yb):
    w = (x - xa)/(xb-xa)
    return w*yb + (1-w)*ya

def valprint(name, var, unit=False, fmt='.3f'):
    if unit:
        var = var.to(unit)
    print(f'{name} = {var:{fmt}}')

def list_unit(lis):
    unit = lis[0].units
    nondim = [l.to(unit).magnitude for l in lis]
    return np.array(nondim)*unit

def ddx(fun, x, dx):
    """
    Does a two-sided finite difference on fun at x, using dx as a spacing.
    """
    dy = fun(x+dx) - fun(x-dx)
    return dy/2/dx

def EE(fun, t_arr, y0, vec=True):
    """
    Explicit Euler using fun(t, y).
    t_arr array of times to use; dt equal throughout.
    vec: if True, fun takes a list of y and returns a list of dydt.
    """
    if vec:
        y_all = [y0]
        dt = t_arr[1]-t_arr[0]
        for t in t_arr:
            y_old = y_all[-1]
            dydt = fun(t, y_old)
            y_new = [y_o + dt*dy for y_o, dy in zip(y_old, dydt)]
            y_all.append(y_new)
            #transpose the lists
            y_vars = [list(i) for i in zip(*y_all[:-1])]
            for i in range(len(y_vars)):
                unit = y_vars[i][0].units
                y_vars[i] = np.array([y.to(unit).magnitude for y in y_vars[i]])*unit
    else:
        y_all = [y0]
        dt = t_arr[1]-t_arr[0]
        for t in t_arr:
            y_old = y_all[-1]
            dydt = fun(t, y_old)
            y_new = y_old + dt*dydt
            y_all.append(y_new)
            y_vars = list_unit(y_all)[:-1]
    return y_vars

def ivp_wrapper(ivp, y_dim, t_unit):
    y_unit = [y.units for y in y_dim]
    ret1 = ivp(0*t_unit, y_dim)
    def new_ivp(t, y_ndim):
        y_dim = [n * u for n, u in zip(y_ndim, y_unit)]
        step = ivp(t, y_dim)
        ret = [s.to(u/t_unit).magnitude for s, u in zip(step, y_unit)]
        return ret
    y_ndim = [y.magnitude for y in y_dim]
    return new_ivp, y_ndim


# -----------------------------
# Dr. Knott's property files

import wheelpy.props.waterproperties as wat
import wheelpy.props.airproperties as air
import wheelpy.props.benzeneproperties as ben

class wwat:
    """
    Wrapper on Dr. Knotts' file.
    tc          critical temperature in K               
    pc          critical pressure in Pa                 
    vc          critical volume in m**3/mol             
    zc          critical compressibility factor         
    mw          molecular weight in kg/mol              
    acen        acentric factor                         
    ldn(t)      liquid density in kg/m**3               
    lcp(t)      liquid heat capacity in J/mol/K         
    ltc(t)      liquid thermal conductivity in W/m/K    
    vp(t)       liquid vapor pressure in Pa             
    hvp(t)      heat of vaporization in J/mol           
    Values evaluated at 1 atm:
    pr(t)       Prandtl number                          
    lvs(t)      liquid viscosity in Pa*s                
    nu(t)       liquid kinematic viscosity in m**2/s    
    tsat(p)     temperature at saturation in K          
    vvs(t)      vapor (steam) viscosity in Pa*s         
    vtc(t)      vapor (steam) therm. conductiv. in W/m/K
    vdnsat(t)   vapor (steam) density at saturation     
    dabair(t)   D_AB for H2O in air
    """
    # Values
    tc = wat.tc*uReg.K
    pc = wat.pc*uReg.Pa
    vc = wat.vc*uReg.m**3/uReg.mol
    zc = wat.zc
    mw = wat.mw*uReg.kg/uReg.mol
    acen = wat.acen

    # Functions
    ldn = uReg.wraps("kg/m**3", "K")(wat.ldn)
    lcp = uReg.wraps("J/mol/K", "K")(wat.lcp)
    ltc = uReg.wraps("W/m/K", "K")(wat.ltc)
    vp = uReg.wraps("Pa", "K")(wat.vp)
    hvp = uReg.wraps("J/mol", "K")(wat.hvp)
    pr = uReg.wraps("", "K")(wat.pr)
    lvs = uReg.wraps("Pa*s", "K")(wat.lvs)
    nu = uReg.wraps("m**2/s", "K")(wat.nu)
    tsat = uReg.wraps("K", "Pa")(wat.tsat)
    vvs = uReg.wraps("Pa*s", "K")(wat.vvs)
    vtc = uReg.wraps("W/m/K", "K")(wat.vtc)
    vdnsat = uReg.wraps("kg/m**3", "K")(wat.vdnsat)
    def dab_air(T, P=1*uReg.atm):
        Tref = 298*uReg.K
        Pref = 1*uReg.atm
        dab_ref = .26E-4*uReg.m**2/uReg.s
        dab = (T/Tref)**1.5 * (Pref/P) * dab_ref
        return dab

class wair:
    """
    Wrapper on Dr. Knotts' file.
    Values:
    tc          critical temperature in K                                                                
    pc          critical pressure in Pa                                                                  
    vc          critical volume in m**3/mol                                                              
    zc          critical compressibility factor                                                          
    mw          molecular weight in kg/mol                                                               
    acen        acentric factor                                                                          
    Functions:
    icp(t)      ideal gas heat capacity in J/mol/K                                                      
    vtc(t)      vapor thermal conductivity in W/m/K                                                      
    vvs(t)      vapor viscosity in Pa*s                                                                  
    hvp(t)      heat of vaporization in J/mol                                                            
    Functions which assume 1 atm:
    rho(t)  density at 1 atm in kg/m**3                                                              
    nu(t)   kinematic viscosity at 1 atm in m**2/s                                                   
    alp(t)thermal diffusivity at 1atm in m**2/s                                                    
    pr(t)   Prandtl number at 1 atm                                                                  
    """
    # Values
    tc   = air.tc*uReg.K
    pc   = air.pc*uReg.Pa
    vc   = air.vc*uReg.m**3/uReg.mol
    zc   = air.zc
    mw   = air.mw*uReg.kg/uReg.mol
    acen = air.acen

    # Functions
    icp = lambda t: air.icp(t.to("K").magnitude) * uReg.J/uReg.mol/uReg.K
    vtc = lambda t: air.vtc(t.to("K").magnitude) * uReg.W/uReg.m/uReg.K
    vvs = lambda t: air.vvs(t.to("K").magnitude) * uReg.Pa*uReg.s
    hvp = lambda t: air.hvp(t.to("K").magnitude) * uReg.J/uReg.mol

    # Properties at 1 atm
    rho = lambda t: air.rho1atm(t.to("K").magnitude) * uReg.kg/uReg.m**3
    nu  = lambda t: air.nu1atm(t.to("K").magnitude)   * uReg.m**2/uReg.s
    alp = lambda t: air.alpha1atm(t.to("K").magnitude) * uReg.m**2/uReg.s
    pr  = lambda t: air.pr1atm(t.to("K").magnitude)   


class wben:
    """
    Wrapper on Dr. Knotts' file.
    tc          critical temperature in K               
    pc          critical pressure in Pa                 
    mw          molecular weight in kg/mol              
    ldn(t)      liquid density in kg/m**3               
    lcp(t)      liquid heat capacity in J/mol/K         
    ltc(t)      liquid thermal conductivity in W/m/K    
    vp(t)       liquid vapor pressure in Pa             
    hvp(t)      heat of vaporization in J/mol           
    lpr(t)      liquid Prandtl number                          
    lvs(t)      liquid viscosity in Pa*s                
    nu(t)       liquid kinematic viscosity in m**2/s    
    tsat(p)     temperature at saturation in K          
    vvs(t)      vapor viscosity in Pa*s         
    vtc(t)      vapor therm. conductiv. in W/m/K
    vPr(t)      vapor Prandtl number     
    """
    # Values
    tc   = ben.tc*uReg.K
    pc   = ben.pc*uReg.Pa
    mw   = ben.mw*uReg.kg/uReg.mol

    # Functions
    ldn  = lambda t: ben.ldn(t.to("K").magnitude) * uReg.kg/uReg.m**3
    lcp  = lambda t: ben.lcp(t.to("K").magnitude) * uReg.J/uReg.mol/uReg.K
    ltc  = lambda t: ben.ltc(t.to("K").magnitude) * uReg.W/uReg.m/uReg.K
    vp   = lambda t: ben.vp(t.to("K").magnitude)   * uReg.Pa
    hvp  = lambda t: ben.hvp(t.to("K").magnitude) * uReg.J/uReg.mol
    lpr  = lambda t: ben.lpr(t.to("K").magnitude)   
    lvs  = lambda t: ben.lvs(t.to("K").magnitude) * uReg.Pa*uReg.s
    nu   = lambda t: ben.nu(t.to("K").magnitude)   * uReg.m**2/uReg.s
    tsat = lambda p: ben.tsat(p.to("Pa").magnitude)*uReg.K
    vvs  = lambda t: ben.vvs(t.to("K").magnitude) * uReg.Pa*uReg.s
    vtc  = lambda t: ben.vtc(t.to("K").magnitude) * uReg.W/uReg.m/uReg.K
    vPr  = lambda t: ben.vPr(t.to("K").magnitude)
    
    def dab_air(T, P=1*uReg.atm):
        dab_ref = .88E-5*uReg.m**2/uReg.s
        Tref = 298*uReg.K
        Pref = 1*uReg.atm
        dab = (T/Tref)**1.5 * (Pref/P) * dab_ref
        return dab

# from iapws import IAPWS97
# 
# class Steam:
#     """
#     Wrapper by Nathan Barrett for iapws.
#     See https://pypi.org/project/iapws/ for documentation
#     """
#     def __init__(self,T="NotGiven",P="NotGiven",x="NotGiven",Print=False):
#         """
#         T    : Temperature
#         P    : Pressure
#         x    : Steam Quality
#         H    : Enthalpy
#         S    : Entropy
#         l    : wavelength of light emmited
#         V    : Specific volume
#         """
#         PossibleVars = [T,P,x]
#         varsGiven = []
#         for i in range(len(PossibleVars)):
#             if PossibleVars[i] != "NotGiven":
#                 varsGiven.append(i)
#         if len(varsGiven) != 2:
#             raise Exception("Error! Steam Class definition can only take two arguments. " + str(len(varsGiven)) + " were provided.")
#         
#         if varsGiven == [0,1]:
#             T = T.to(uReg.K).magnitude
#             P = P.to(uReg.MPa).magnitude
#             self.steam = IAPWS97(T=T,P=P)
# 
#         elif varsGiven == [0,2]:
#             T = T.to(uReg.K).magnitude
#             self.steam = IAPWS97(T=T,x=x)
#         elif varsGiven == [1,2]:
#             P = P.to(uReg.MPa).magnitude
#             self.steam = IAPWS97(P=P,x=x)
#             
#         propertyFound = False
#         
#         try:
#             self.H = (self.steam.h) * uReg.kJ / uReg.kg
#             propertyFound = True
#         except:
#             self.H = "NotFound"
#         try:
#             self.S = (self.steam.s) * uReg.kJ / uReg.kg / uReg.K
#             propertyFound = True
#         except:
#             self.S = "NotFound"
#         try:
#             self.P = (self.steam.P) * uReg.MPa
#             propertyFound = True
#         except:
#             self.P = "NotFound"
#         try:
#             self.T = (self.steam.T) * uReg.K
#             propertyFound = True
#         except:
#             self.T = "NotFound"
#         try:
#             self.x = self.steam.x
#             propertyFound = True
#             if Print:
#                 if self.x > 1:
#                     print("This is Saturated Vapor.")
#                 elif self.x > 0:
#                     print("This is in vapor liquid equilibirum.")
#                 elif self.x <= 0:
#                     print("This is Saturated Liquid.")
#         except:
#             self.x = "NotFound"
#             print("WARNING! Steam condition could not be determined.")
#         try:
#             self.l = (self.steam.l) * uReg.nm
#             propertyFound = True
#         except:
#             self.l = "NotFound"
#         try:
#             self.V = (self.steam.v) * uReg.m**3 / uReg.kg
#             propertyFound = True
#         except:
#             self.V = "NotFound"
#             
#         if not propertyFound:
#             raise Exception("There is no data for the conditions given.")
