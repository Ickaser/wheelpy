import numpy as np
import wheelpy.muc as muc
from scipy.optimize import curve_fit, minimize
from inspect import signature
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

class line_estimates:
    def __init__(self, x_arr, y_arr):
        sf = self
        sf.x_arr = x_arr
        sf.y_arr = y_arr
        sf.n = len(x_arr)
        sf.xb = np.mean(x_arr)
        sf.yb = np.mean(y_arr)
        
        sf.Sxx = np.sum((sf.x_arr-sf.xb)**2)
        sf.Syy = np.sum((sf.y_arr-sf.yb)**2)
        sf.Sxy = np.sum((sf.x_arr-sf.xb)*(sf.y_arr-sf.yb))
        
        sf.b1 = sf.Sxy/sf.Sxx
        sf.b0 = self.yb - self.b1*self.xb
        sf.yh_arr = self.b0 + self.b1*self.x_arr

        # Correlation coefficient between parameters
        sf.rho = -sf.xb/np.sqrt(sf.Sxx/sf.n + sf.xb**2)

        sf.SS = np.sum((sf.y_arr-sf.yh_arr)**2)
        sf.df = sf.n-2
        sf.se = np.sqrt(sf.SS/sf.df)

        b0_var = self.se**2 * (1/self.n + self.xb**2/self.Sxx)
        b1_var = self.se**2 / self.Sxx
        self.b0_var = b0_var
        self.b1_var = b1_var
    def pred_conf(self, x, alpha=.05):
        yh = self.b0 + self.b1*x
        ts = st.t.ppf(1-alpha/2, self.n)
        pred = self.se*np.sqrt(1/self.n+(x-self.xb)**2/self.Sxx)*ts
        return yh+pred, yh-pred
    def onep_conf(self, x, alpha=.05):
        yh = self.b0 + self.b1*x
        ts = st.t.ppf(1-alpha/2, self.n)
        onep = self.se*np.sqrt(1/1 + 1/self.n+(x-self.xb)**2/self.Sxx)*ts
        return yh+onep, yh-onep
    def conf_ints_wid(self, alpha=.05):
        ts = st.t.ppf(1-alpha/2, self.n-2)
        b0_wid = ts*np.sqrt(self.b0_var)
        b1_wid = ts*np.sqrt(self.b1_var)
        b0, b1 = self.b0, self.b1
        self.conf_rectx = [b0-b0_wid, b0-b0_wid, b0+b0_wid, b0+b0_wid, b0-b0_wid]
        self.conf_recty = [b1-b1_wid, b1+b1_wid, b1+b1_wid, b1-b1_wid, b1-b1_wid]
        return b0_wid, b1_wid
    def conf_regi_xyy(self, alpha=.05):
        S_xsq = np.sum(self.x_arr)**2
        S_sqx = np.sum(self.x_arr**2)
        S_x = np.sum(self.x_arr)
        alpha1 = 1 - alpha # .95
        alpha2 = alpha1**2 # .9025
        F = st.f.ppf(alpha2, 2, self.n-2)
        β0_wid = np.sqrt(F * 2 * self.se**2 * S_sqx / (self.n*S_sqx-S_xsq))
        β0_lo = self.b0 - β0_wid
        β0_hi = self.b0 + β0_wid
        β0_arr = np.linspace(β0_lo, β0_hi)
        root = np.sqrt((self.b0 - β0_arr)**2*(S_xsq - self.n*S_sqx) + 2*self.se**2*F*S_sqx)
        root = np.nan_to_num(root)
        β1_p = self.b1 + ((self.b0-β0_arr) * S_x + root) / S_sqx
        β1_m = self.b1 + ((self.b0-β0_arr) * S_x - root) / S_sqx
        return β0_arr, β1_p, β1_m

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
        der = muc.ddx(partFun, x, 1e-6 if abs_der_width else ph[i]*1-3)
        der_list.append(der)
    return np.array(der_list)

class trad_nonlin_estimates:
    def __init__(self, xdat, ydat, model, p0):
        """
        Args:
            x_arr, y_arr: data to fit
            model: function, with signature yh = model(x, a0, a1, ..., an)
            p0: guess values for model parameters
        After setting up, parameter best estimates are available as obj.fit;
        call obj.calc_Fmat, and obj.F, obj.FTF, and obj.FTF1 become available;
        call obj.pred_conf(x) and obj.onep_conf(x) to get confidence bands.
        """
        sf = self
        sf.xdat = xdat
        sf.ydat = ydat
        sf.n = len(xdat)
        sf.p = len(p0)
        sf.fit, sf.covar = curve_fit(fun, xdat, ydat, p0)
        sf.yh_arr = fun(xdat, *fit)
        sf.s2e = sf.calc_S(sf.fit)
        sf.se = np.sqrt(s2e)
        # sf.wid_args = (x_arr, y_arr, fun, covar, fit)
        # conf_upper, conf_lower = calc_conf_band(*wid_args, alpha)
        # pred_upper, pred_lower = calc_pred_band(*wid_args, alpha)
    def calc_S(sf, params):
        md = self.model(sf.xdat, *params)
        sumsq = np.sum((md - sf.ydat)**2)
        return sumsq
    def calc_Fmat(sf, getJac=None):
        """
        Calculate the F, F'F, and F'F-1 matrices. 
        Pass in a function getJac(x, b0, b1,...,bn) to use analytical derivatives;
        takes numerical differences if getJac is not supplied.
        """
        sf.getJac = getJac
        if getJac == None:
            sf.F = np.array([curve_deriv(sf.model, xi, sf.params) for xi in sf.xdat])
        else:
            sf.F = np.array([getJac(xi, *sf.params) for xi in sf.xdat])
        sf.FTF = np.dot(sf.F.T, sf.F)
        sf.FTF1 = np.linalg.inv(sf.FTF)
        return sf.F
    def param_conf_int(self, alpha=.05):
        """
        Returns the width of the confidence intervals for each parameter.
        """
        try:
            sf.Vmat = sf.FTF1*sf.s2e
        except:
            sf.calc_Fmat()
            sf.Vmat = sf.FTF1*sf.s2e
        ts = st.t.ppf(1-alpha/2, sf.n-sf.p)
        widths = ts*np.diag(sf.Vmat)
        return widths
    # def jcr_2param(self, alpha=.05):
    #     """
    #     If model only has two parameters, calculate the joint confidence region from F'F.
    #     """
    #     # b0_var = self.se**2 * (1/self.n + self.xb**2/self.Sxx)
    #     # b1_var = self.se**2 / self.Sxx
    #     # ts = st.t.ppf(1-alpha/2, self.n-2)
    #     # b0_wid = ts*np.sqrt(b0_var)
    #     # b1_wid = ts*np.sqrt(b1_var)
    #     # b0, b1 = self.b0, self.b1
    #     # self.conf_rectx = [b0-b0_wid, b0-b0_wid, b0+b0_wid, b0+b0_wid, b0-b0_wid]
    #     # self.conf_recty = [b1-b1_wid, b1+b1_wid, b1+b1_wid, b1-b1_wid, b1-b1_wid]
    #     # return b0_wid, b1_wid

    def pred_conf(sf, x_arr, alpha=.05):
        F0_arr = np.array([calc_curve_deriv(sf.model, x, sf.fit) for x in x_arr])
        c_var = [np.linalg.multi_dot((F0, sf.FTF1, F0)) for F0 in F0_arr]
        yh = sf.model(x_arr, *fit)
        ts = st.t.ppf(1-alpha/2, sf.n-sf.p)
        wid = sf.se*ts * np.sqrt(c_var)
        return yh+wid, yh-wid

    def onep_conf(self, x, alpha=.05):
        F0_arr = np.array([calc_curve_deriv(sf.model, x, sf.fit) for x in x_arr])
        c_var = [np.linalg.multi_dot((F0, sf.FTF1, F0)) for F0 in F0_arr]
        yh = sf.model(x_arr, *fit)
        ts = st.t.ppf(1-alpha/2, sf.n-sf.p)
        wid = sf.se*ts * np.sqrt(1 + c_var)
        return yh+wid, yh-wid
    # def conf_regi_xyy(self, alpha=.05):
    #     """
    #     Not currently working.
    #     """
    #     S_xsq = np.sum(self.x_arr)**2
    #     S_sqx = np.sum(self.x_arr**2)
    #     S_x = np.sum(self.x_arr)
    #     alpha1 = 1 - alpha # .95
    #     alpha2 = alpha1**2 # .9025
    #     F = st.f.ppf(alpha2, 2, self.n-2)
    #     β0_wid = np.sqrt(F * 2 * self.se**2 * S_sqx / (self.n*S_sqx-S_xsq))
    #     β0_lo = self.b0 - β0_wid
    #     β0_hi = self.b0 + β0_wid
    #     β0_arr = np.linspace(β0_lo, β0_hi)
    #     root = np.sqrt((self.b0 - β0_arr)**2*(S_xsq - self.n*S_sqx) + 2*self.se**2*F*S_sqx)
    #     root = np.nan_to_num(root)
    #     β1_p = self.b1 + ((self.b0-β0_arr) * S_x + root) / S_sqx
    #     β1_m = self.b1 + ((self.b0-β0_arr) * S_x - root) / S_sqx
    #     return β0_arr, β1_p, β1_m


# ----------------------------------------
# Provided here is sample code for how to use the above functions.
# curve_func is the curve you fit, which should have the form y = curve_func(x, b0, b1, ..., bn)
# p0 is an estimate of parameters. This only needs to be within a couple orders of magnitude, the solver is very stable for that.

# nonlin_obj = trad_nonlin_estimates(x_dat, y_dat, curve_func, p0) 
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

# ---------------------- Obsolete functions I didn't dare delete
# def calc_curve_var(x, curve, covar, fit):
#     Gp = curve_deriv(curve, x, fit)
#     c_var = np.dot(Gp, np.dot(covar, Gp))
#     return c_var
# def calc_pred_band(x_arr, y_arr, curve, covar, fit, alpha=.05):
#     n=len(x_arr)
#     yh = curve(x_arr, *fit)
#     e_var = np.var(y_arr - yh, ddof=1)
#     c_var = np.array([calc_curve_var(x, curve, covar, fit) for x in x_arr]) + e_var
#     ts = st.t.ppf(1-alpha/2, n-1)
#     wid = np.sqrt(c_var)*ts
#     return yh+wid, yh-wid
# def calc_conf_wid(x, x_arr, y_arr, curve, covar, fit, alpha=.05):
#     n = len(x_arr)
#     try:
#         c_var = np.array([calc_curve_var(xi, curve, covar, fit) for xi in x])
#     except:
#         c_var = calc_curve_var(x, curve, covar, fit)
#     ts = st.t.ppf(1-alpha/2, n-1)
#     wid = np.sqrt(c_var)*ts
#     return wid
# def calc_pred_wid(x, x_arr, y_arr, curve, covar, fit, alpha=.05):
#     n=len(x_arr)
#     yh = curve(x_arr, *fit)
#     try:
#         c_var = np.array([calc_curve_var(xi, curve, covar, fit) for xi in x])
#     except:
#         c_var = calc_curve_var(x, curve, covar, fit)
#     e_var = np.var(y_arr - yh, ddof=1)
#     c_var += e_var
#     ts = st.t.ppf(1-alpha/2, n-1)
#     wid = np.sqrt(c_var)*ts
#     return wid
# def nonlin_analysis(x_arr, y_arr, fun, p0, alpha=.05):
#     """
#     Pass in x_arr, y_arr, curve_func, parameter guess p0
#     Returns [yhat, conf_upper, conf_lower, pred_upper, pred_lower], [fit, covar]
    
#     If confidence or prediction bands seem off, the most likely culprit is the curve_deriv function used inside here.
#     To make curve_deriv a little more well-behaved, you can try making abs_der_width False, which will then use .1*fit parameter to calculate a finite difference instead of a fixed .01. ALternatively, you can take the gradient by hand of your curve function with respect to fitting parameters and use that instead of the finite differences.
#     """
#     n = len(x_arr)
#     fit, covar = curve_fit(fun, x_arr, y_arr, p0)
#     yh_arr = fun(x_arr, *fit)
#     wid_args = (x_arr, y_arr, fun, covar, fit)
#     conf_upper, conf_lower = calc_conf_band(*wid_args, alpha)
#     pred_upper, pred_lower = calc_pred_band(*wid_args, alpha)
#     return [yh_arr, conf_upper, conf_lower, pred_upper, pred_lower], [fit, covar]

