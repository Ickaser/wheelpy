# Possibly to move into thermo.py ?
import numpy as np
import wheelpy.muc as muc
un = muc.uReg

def calc_Psat_arb(T, coeff):
    """
    DIPPr style Psat calculation.
    """
    T = T.to("K").magnitude
    A, B, C, D, E = coeff
    Y = np.exp(A + B/T + C*np.log(T) + D*T**E)
    return Y*un.Pa
pen_VP_coeff = (78.741, -5420.3, -8.8253, 9.6171e-6, 2)
hex_VP_coeff = (104.65, -6995.5, -12.702, 1.2381e-5, 2)
pen_Psat = lambda t: calc_Psat_arb(t, pen_VP_coeff)
hex_Psat = lambda t: calc_Psat_arb(t, hex_VP_coeff)
def calc_Psat_ant(T, coeff, unit="degC"):
    """
    Antoine vapor pressure calculation.
    """
    A, B, C = coeff
    try:
        T = T.to("degC").magnitude
    except:
        pass
    Y = A - B/(T+C)
    return 10**Y *un.mmHg
def calc_Ki_arb(T, P, coeff):
    T = T.to("degR").magnitude
    P = P.to("psi").magnitude
    at1, at2, at6, ap1, ap2, ap3 = coeff
    Y = at1/T/T + at2/T + at6 + ap1*np.log(P) + ap2/P/P + ap3/P
    K = np.exp(Y)
    return K
met_Ki_coeff = [-292860, 0, 8.2445, -.8951, 59.8465, 0]
eth_Ki_coeff = [-600076.875, 0, 7.90595, -.84677, 42.94594, 0]
pro_Ki_coeff = [-923484.6875, 0, 7.71725, -0.87871, 47.67624, 0]
but_Ki_coeff = [-1280557, 0, 7.94986, -0.96455, 0, 0]
pen_Ki_coeff = [-1524891, 0, 7.33129, -.89143, 0, 0]
hex_Ki_coeff = [-1778901, 0, 6.96783, -.84634, 0, 0]
hep_Ki_coeff = [-2013803, 0, 6.52914, -.79543, 0, 0]
oct_Ki_coeff = [0, -7646.81641, 12.48457, -.73152, 0, 0]
calc_K_met = lambda t, p: calc_Ki_arb(t, p, met_Ki_coeff)
calc_K_eth = lambda t, p: calc_Ki_arb(t, p, eth_Ki_coeff)
calc_K_pro = lambda t, p: calc_Ki_arb(t, p, pro_Ki_coeff)
calc_K_but = lambda t, p: calc_Ki_arb(t, p, but_Ki_coeff)
calc_K_pen = lambda t, p: calc_Ki_arb(t, p, pen_Ki_coeff)
calc_K_hex = lambda t, p: calc_Ki_arb(t, p, hex_Ki_coeff)
calc_K_hep = lambda t, p: calc_Ki_arb(t, p, hep_Ki_coeff)
calc_K_oct = lambda t, p: calc_Ki_arb(t, p, oct_Ki_coeff)

def calc_RachRice(VF, zi, Ki):
    zi = np.array(zi)
    Ki = np.array(Ki)
    terms = (Ki-1)*zi/(1 + (Ki-1)*VF)
    return np.sum(terms)
def calc_RR_xi(VF, zi, Ki):
    xi = zi/(1+(Ki-1)*VF)
    return xi
def calc_RR_yi(VF, zi, Ki):
    yi = Ki*zi/(1+(Ki-1)*VF)
    return yi
