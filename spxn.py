# Possibly to move into thermo.py ?
import numpy as np
import wheelpy.muc as muc
from scipy.optimize import fsolve
un = muc.uReg

def calc_Psat_arb(T, coeff):
    """
    DIPPr style Psat calculation.
    """
    T = T.to("K").magnitude
    A, B, C, D, E = coeff
    Y = np.exp(A + B/T + C*np.log(T) + D*T**E)
    return Y*un.Pa
methanol_VP_coeff = (82.718, -6904.5, -8.8622, 7.4664e-6, 2)
water_VP_coeff = (73.649, -7258.2, -7.3037, 4.1653e-6, 2)
pen_VP_coeff = (78.741, -5420.3, -8.8253, 9.6171e-6, 2)
hex_VP_coeff = (104.65, -6995.5, -12.702, 1.2381e-5, 2)
methanol_Psat = lambda t: calc_Psat_arb(t, methanol_VP_coeff)
water_Psat = lambda t: calc_Psat_arb(t, water_VP_coeff)
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

methwat_WilsonLL = [.41830958, 1.0717958]

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
non_Ki_coeff = [-2551040, 0, 5.69313, -.67818, 0, 0]
calc_K_met = lambda t, p: calc_Ki_arb(t, p, met_Ki_coeff)
calc_K_eth = lambda t, p: calc_Ki_arb(t, p, eth_Ki_coeff)
calc_K_pro = lambda t, p: calc_Ki_arb(t, p, pro_Ki_coeff)
calc_K_but = lambda t, p: calc_Ki_arb(t, p, but_Ki_coeff)
calc_K_pen = lambda t, p: calc_Ki_arb(t, p, pen_Ki_coeff)
calc_K_hex = lambda t, p: calc_Ki_arb(t, p, hex_Ki_coeff)
calc_K_hep = lambda t, p: calc_Ki_arb(t, p, hep_Ki_coeff)
calc_K_oct = lambda t, p: calc_Ki_arb(t, p, oct_Ki_coeff)
calc_K_non = lambda t, p: calc_Ki_arb(t, p, non_Ki_coeff)

def calc_RachRice(VF, zi, Ki):
    # Simply provides RHS for Rachford Rice equation; use in a solver loop.
    # Ki should be passed as values, not functions.
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

def bin_xy(const_alpha, numPoints = 101):
    x_arr = np.linspace(0, 1, numPoints)
    y_arr = (x_arr*const_alpha)/(1+x_arr*(const_alpha-1))
    return x_arr, y_arr

def fenske(xAD, xAR, αAB, xBD=None, xBR=None):
    if xBR == None:
        xBR = 1-xAR
    if xBD == None:
        xBD = 1-xAD
    Nmin = np.log(xAD/xBD / (xAR/xBR))/np.log(αAB)
    return Nmin

def underwood(α_arr, F, z_arr, DxD_arr, dVf):
    """
    Takes α_arr, F, z as arguments.
    F, z_arr feed flow rate and composition.
    dVf = V - Vb = change in V at feed
    α_arr an array relative volatilities; reference species has α=1. Array is sorted internally.
    """
    def calc_f(ϕ):
        return np.sum(α_arr*F*z_arr/(α_arr-ϕ)).magnitude
    α_sort = np.sort(α_arr)
    guess_vals = (α_sort[1:]+α_sort[:-1])/2
#     guess_vals = α_sort+.1
    ϕ_arr = np.array([fsolve(calc_f, guess)[0] for guess in guess_vals])
#     print(ϕ_arr)
    Vmin = np.sum([α_arr*DxD_arr/(α_arr-ϕ) for ϕ in ϕ_arr], axis=1)*DxD_arr.units
    Vmin = np.max(Vmin)
    D = np.sum(DxD_arr)
    Lmin = Vmin-D
    LDmin = Lmin/D
    return LDmin

def gilliland(LDmin, M, Nmin, Nfmin,):
    """
    Arguments: LDmin, M, Nmin, Nfmin. M= (L/D) / (L/D)_min
    Nmin is minimum stages from Fenske, Nfmin is minimum stages from feed to distillate from Fenske
    Uses Liddle's 1986 fit to Gilliland's correlation.
    """
    X = (M-1)/(1/LDmin + M)
    if 0 <= X and X<= .01:
        Y = 1 - 18.5715 * X
    elif .01 < X and X< .90:
        Y = .545827 - .591422*X + .002743/X
    elif .90 < X and X<= 1.0:
        Y = 1 - 18.5715 * X
    else:
        print("Invalid abscissa for Gilliland correlation.")
        raise
    # Y*N + Y = N - Nmin
    # N * (1-Y) = Y + Nmin
    N = (Y+Nmin)/(1-Y)
    Nf = Nfmin/Nmin * N
    return N, Nf

def FUG(xDi, xBi, zi, alpha_i, F,dVf, M, LKi = 0, HKi = 1):

    """
    Takes arrays of: xD, xB, z, alpha; F, dVf.
    F and dVf feed rate and change in vapor flow at feed.
    M the factor by which reflux is increased above minimum.
    LKi and HKi optional arguments, indicating index of LK and HK components respectively. Default to 0 and 1.
    
    Returns N, Nf, LDmin, Nmin
    """
    lf = (xDi[LKi]-zi[LKi])/(xBi[LKi]-xDi[LKi])
    vf = 1-lf
    D = F*vf
    B = F*lf
    DxD = np.array(xDi)*D
    
    Nmin = fenske(xDi[LKi], xBi[LKi], alpha_i[LKi]/alpha_i[HKi], xBD=xDi[HKi], xBR=xBi[HKi])
    Nfmin = fenske(xDi[LKi], zi[LKi], alpha_i[LKi]/alpha_i[HKi], xBD=xDi[HKi], xBR=zi[HKi])
                   
    LDmin = underwood(alpha_i, F, zi, DxD, dVf)
    N, Nf = gilliland(LDmin, M, Nmin, Nfmin)
    
    return N, Nf, LDmin, Nmin
                

