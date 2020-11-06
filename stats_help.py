import numpy as np
import wheelpy.muc as muc
from scipy.optimize import curve_fit
import scipy.stats as st

def calc_confint_xb(X, α):
    """
    Simple confidence interval for the mean, using an array of values X and (1-alpha) % confidence.
    """
    n = len(X)
    xb = np.mean(X)
    s = np.std(X, ddof=1)
    serr = s/np.sqrt(n)
    t_s = st.t.ppf(1-α/2, n-1)
    return xb, t_s*serr


def curve_deriv(curve, x, ph, abs_der_width=True):
    # Takes partial derivatives of curve with respect to parameters
    # Args: curve function, x-value of interest, estimated parameters
    der_list = []
    for i in range(len(ph)):
        def listSub(l,i, u):
            l = l.copy()
            l[i]=u
            return l
        partFun = lambda bi: curve(x, *listSub(ph, i, bi))
        der = muc.ddx(partFun, x, .01 if abs_der_width else ph[i]/10)
        der_list.append(der)
    return np.array(der_list)
def calc_curve_var(x, curve, covar, fit):
    Gp = curve_deriv(curve, x, fit)
    c_var = np.dot(Gp, np.dot(covar, Gp))
    return c_var
def calc_conf_band(x_arr, y_arr, curve, covar, fit, alpha=.05):
    n=len(x_arr)
    c_var = np.array([calc_curve_var(x, curve, covar, fit) for x in x_arr])
    yh = curve(x_arr, *fit)
    ts = st.t.ppf(1-alpha/2, n-1)
    wid = np.sqrt(c_var)*ts
    return yh+wid, yh-wid
def calc_pred_band(x_arr, y_arr, curve, covar, fit, alpha=.05):
    n=len(x_arr)
    yh = curve(x_arr, *fit)
    e_var = np.var(y_arr - yh, ddof=1)
    c_var = np.array([calc_curve_var(x, curve, covar, fit) for x in x_arr]) + e_var
    ts = st.t.ppf(1-alpha/2, n-1)
    wid = np.sqrt(c_var)*ts
    return yh+wid, yh-wid
def calc_conf_wid(x, x_arr, y_arr, curve, covar, fit, alpha=.05):
    n = len(x_arr)
    try:
        c_var = np.array([calc_curve_var(xi, curve, covar, fit) for xi in x])
    except:
        c_var = calc_curve_var(x, curve, covar, fit)
    ts = st.t.ppf(1-alpha/2, n-1)
    wid = np.sqrt(c_var)*ts
    return wid
def calc_pred_wid(x, x_arr, y_arr, curve, covar, fit, alpha=.05):
    n=len(x_arr)
    yh = curve(x_arr, *fit)
    try:
        c_var = np.array([calc_curve_var(xi, curve, covar, fit) for xi in x])
    except:
        c_var = calc_curve_var(x, curve, covar, fit)
    e_var = np.var(y_arr - yh, ddof=1)
    c_var += e_var
    ts = st.t.ppf(1-alpha/2, n-1)
    wid = np.sqrt(c_var)*ts
    return wid
def nonlin_analysis(x_arr, y_arr, fun, p0, alpha=.05):
    """
    Pass in x_arr, y_arr, curve_func, parameter guess p0
    Returns [yhat, conf_upper, conf_lower, pred_upper, pred_lower], [fit, covar]
    
    If confidence or prediction bands seem off, the most likely culprit is the curve_deriv function used inside here.
    To make curve_deriv a little more well-behaved, you can try making abs_der_width False, which will then use .1*fit parameter to calculate a finite difference instead of a fixed .01. ALternatively, you can take the gradient by hand of your curve function with respect to fitting parameters and use that instead of the finite differences.
    """
    n = len(x_arr)
    fit, covar = curve_fit(fun, x_arr, y_arr, p0)
    yh_arr = fun(x_arr, *fit)
    wid_args = (x_arr, y_arr, fun, covar, fit)
    conf_upper, conf_lower = calc_conf_band(*wid_args, alpha)
    pred_upper, pred_lower = calc_pred_band(*wid_args, alpha)
    return [yh_arr, conf_upper, conf_lower, pred_upper, pred_lower], [fit, covar]

# ----------------------------------------
# Provided here is sample code for how to use the above functions.
# curve_func is the curve you fit, which should have the form y = curve_func(x, b0, b1, ..., bn)
# p0 is an estimate of parameters. This only needs to be within a couple orders of magnitude, the solver is very stable for that.

# lines, fitstuff = nonlin_analysis(x_dat, y_dat, curve_func, p0) 
# yh, conf_u, conf_l, pred_u, pred_l = lines
# y_fit, y_covar = fitstuff

# # Sample code to make smooth prediction bands, instead of just evaluating at x_dat
# wid_args = (x_dat, y_dat, curve_func, y_covar, y_fit)
# x_smooth = np.linspace(x_dat[0], x_dat[-1])
# confwd_sm = calc_conf_wid(x_smooth, *wid_args)
# predwd_sm = calc_pred_wid(x_smooth, *wid_args)
# y_smooth = y_curve(x_smooth, *y_fit)
# conf_u_sm = y_smooth + confwd_sm
# conf_l_sm = y_smooth - confwd_sm
# pred_u_sm = y_smooth + predwd_sm
# pred_l_sm = y_smooth - predwd_sm

# plt.plot(x_dat, y_dat, label="data")
# plt.plot(x_smooth, y_smooth, label="curve fit")
# plt.plot(x_smooth, conf_u_sm, label="upper confidence band")
# plt.plot(x_smooth, conf_l_sm, label="lower confidence band")
# plt.plot(x_smooth, pred_u_sm, label="upper prediction band")
# plt.plot(x_smooth, pred_l_sm, label="lower prediction band")
