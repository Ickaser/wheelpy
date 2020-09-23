import wheelpy.muc as muc
un = muc.uReg
R = muc.R
from scipy.optimize import fsolve, curve_fit, root
from scipy.integrate import quad
import numpy as np

class EOS:
    """
    Initialize the object:
    obj = EOS(kind, T, P, Tc, Pc, omega, pint=True)
    Then, call:
    obj.calc_Z()
    Afterwards, other functions are available: (if root is an argument, pass "liq" or "vap")
    obj.calc_V(root="all")
    obj.calc_fugacity(root="vap")
    obj.calc_residG(root="liq")
    obj.calc_residSH()

    """

    def __init__(self, kind, T, P, Tc, Pc, omega, pint=True):
        """
        Cubic kinds: ["vdWg", "RK", "SRK", "PR"] (vdWg is generalized)
        Other kinds: ["Pitzer", "vdWn", "LeeKesler"] (vdWn is normal form)
        pint = True: sets whether to account for using pint module. If True, removes units before using fsolve.
            If pint=False, uses SI units for R.
        """
        self.kind = kind
        self.pint = pint
        self.P = P
        self.T = T
        self.Pc = Pc
        self.Tc = Tc
        self.omega = omega

        # Used to avoid recalculating Z unnecessarily
        self.calc_finished = False

        self.Pr = P/Pc
        self.Tr = T/Tc
        if pint:
            self.Pr.ito(un.dimensionless)
            self.Tr.ito(un.dimensionless)

        self.cubic = ["vdWg", "RK", "SRK", "PR"]
        self.other = ["Pitzer", "LeeKesler"]
        if kind in self.cubic:
            dat = np.array([
            #    sig           eps            Omega   Psi       Zc
                [0           , 0           ,  1/8   , 27/64  ,  3/8   ],
                [1           , 0           ,  .08664, .42748,  1/3   ],
                [1           , 0           ,  .08664, .42748,  1/3   ],
                [1+np.sqrt(2), 1-np.sqrt(2),  .07780, .45724,  .30740],
            ])
            alphas = [
                    lambda omega, Tr: 1,
                    lambda omega, Tr: np.sqrt(Tr),
                    lambda omega, Tr: (1 + (.480 + 1.574*omega - .176*omega**2)*(1- np.sqrt(Tr)))**2,
                    lambda omega, Tr: (1 + (.37464 + 1.54226*omega - .26992*omega**2)*(1- np.sqrt(Tr)))**2,
                    ]
            if kind == "vdWg":
                self.cnst = dat[0]
                self.alph = alphas[0]
            elif kind == "RK":
                self.cnst = dat[1]
                self.alph = alphas[1]
            elif kind == "SRK":
                self.cnst = dat[2]
                self.alph = alphas[2]
            elif kind == "PR":
                self.cnst = dat[3]
                self.alph = alphas[3]
            self.sig, self.eps, self.Omega, self.Psi, self.Zc = self.cnst
        # TODO  
        elif kind in self.other:
            if kind == "LeeKesler":
                "nothing to see here"
                #Plan to iterate!
            elif kind == "Pitzer":
                self.B0 = 0.083 - .422/(self.Tr**1.6)
                self.B1 = 0.139 - .172/(self.Tr**4.2)
        else:
            raise ValueError("Passed an invalid kind of EOS.")

    # Intermediate functions for cubics
    @staticmethod
    def calc_beta(Omega, Pr, Tr):
        return Omega*Pr/Tr
    # @staticmethod
    # def calc_q(Psi, alpha, omega, Tr, Omega):
    #     # alpha is a callable, alpha(omega, Tr)
    #     return Psi * alpha(omega, Tr)/Omega/Tr
    def calc_q(self):
        # alpha is a callable, alpha(omega, Tr)
        return self.Psi * self.alph(self.omega, self.Tr)/self.Omega/self.Tr

    # This has been tested and proven for vdWg, and gives similar answers for other methods.
    def calc_Z(self, guess = np.logspace(-5, 2, 7)):
        """
        Takes no inputs.
        Returns a list of solutions.
        Uses values defined when EOS object is initialized.
        For cubic EOS:
        Runs a function solver on the EOS to compute Z
        Uses a default of 7 guess values from 10^-5 to 10^2
        For Pitzer:
        Returns the result of the calculation.
        For Lee-Kesler:
        Uses Tr, Pr to prompt user for values from the table. If has been calculated already, does not repeat. 
        """
        if self.kind in self.cubic:
            cnst = self.cnst
            # Define two separate functions, depending on if pint is used. Nearly identical.
            if self.pint:
                def z_solve(Z):
                    beta = self.calc_beta(self.cnst[2], self.Pr.magnitude, self.Tr.magnitude)
                    q = self.calc_q()
                    r1 = 1 + beta - q*beta*(Z-beta)/(Z+cnst[1]*beta)/(Z+cnst[0]*beta) - Z
                    return r1.magnitude
            else:
                def z_solve(Z):
                    beta = self.calc_beta(self.cnst[2], self.Pr, self.Tr)
                    q = self.calc_q()
                    r1 = 1 + beta - q*beta*(Z-beta)/(Z+cnst[1]*beta)/(Z+cnst[0]*beta) - Z
                    return r1

            Z_guess = guess
            Z_check = [fsolve(z_solve, z_g)[0] for z_g in Z_guess]
            Z_sol = []
            for z in Z_check:
                add = True
                for zj in Z_sol:
                    if np.abs(z-zj) < .001 :
                        add = False
                if z < 0:
                    add = False
                if add:
                    Z_sol.append(z)
            Z_sol = np.sort(Z_sol)
            if len(Z_sol) == 3:
                self.Z_liq = np.min(Z_sol)
                self.Z_vap = np.max(Z_sol)
            self.Z_sol = Z_sol

            b = self.Omega*muc.R*self.Tc/self.Pc
            self.b = b
            self.q = self.calc_q()

        elif self.kind == "Pitzer":
            self.Z_sol = [1 + self.Pr/self.Tr * (self.B0 + self.omega*self.B1)]

        elif self.kind == "LeeKesler":
            if self.calc_finished:
                return self.Z_sol
            refTr = np.array([.30, .35, .40, .45, .50, .55, .60, .65, .70, .75, .80, .85, .90,
                              .93, .95, .97, .98, .99, 1.00, 1.01, 1.02, 1.05, 1.10, 1.15, 1.20,
                              1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60,
                              2.80, 3.00, 3.50, 4.00])
            refPr = np.array([.0100, .0500, .1000, .2000, .4000, .6000, .8000, 1.0000,
                              1.2000, 1.5000, 2.0000, 3.0000, 5.0000, 7.0000, 10.0000])
            iTr = 0
            jPr = 0
            if (self.Tr < refTr[0] or self.Tr > refTr[-1]):
                raise ValueError("Tr out of range for LeeKesler table.")
            elif (self.Pr < refPr[0] or self.Pr > refPr[-1]):
                raise ValueError("Pr out of range for Lee-Kesler table.")
            for i, t in enumerate(refTr):
                if t > self.Tr:
                    iTr = i-1
                    break
            for j, p in enumerate(refPr):
                if p > self.Pr:
                    jPr = j-1
                    break
            Z0_vals = np.zeros((2,2))
            for i in range(0,2):
                for j in range(0,2):
                    Z0_vals[i, j] = float(input(f"Z0: Tr = {refTr[iTr + i]}, Pr = {refPr[jPr + j]}"))
            Z0 = self.interp2d(self.Tr, self.Pr, refTr[iTr:iTr+2], refPr[jPr:jPr+2], Z0_vals )
            Z1_vals = np.zeros((2,2))
            for i in range(0,2):
                for j in range(0,2):
                    Z1_vals[i, j] = float(input(f"Z1: Tr = {refTr[iTr + i]}, Pr = {refPr[jPr + j]}"))
            Z1 = self.interp2d(self.Tr, self.Pr, refTr[iTr:iTr+2], refPr[jPr:jPr+2], Z1_vals )
            self.Z_sol = [Z0 + self.omega*Z1]

        else:
            raise ValueError("Requested EOS not programmed yet.")
        self.calc_finished = True
        return self.Z_sol

    def calc_Z_alt(self, V):
        """
        Uses stored value of T and a given value of V to work out Z, which is explicit.
        """
        if not self.calc_finished:
            self.calc_Z()
        if self.kind in self.cubic:
            rho = 1/V
            b = self.b
            q = self.q
            eps = self.eps
            sig = self.sig
            def rho_Z(rho):
                r1 = 1/(1-rho*b)
                r2 = q*rho*b / (1 + eps*rho*b) / (1 + sig*rho*b)
                return (r1 - r2)
            self.z_alt = rho_Z(rho)
            return self.z_alt
        else:
            raise ValueError("EOS kind not implemented yet")


    @staticmethod
    def interp2d(x, y, xmp, ymp, vals):
        """
        x: actual value
        y: actual value
        xmp: (lower x ref, upper x ref)
        ymp: (lower y ref, upper y ref)
        vals: ((z(xm, ym), z(xp, ym)), (z(xm, yp), z(xp, yp)))
        """
        xm, xp = xmp
        ym, yp = ymp
        xs = (x-xm)/(xp-xm)
        ys = (y-ym)/(yp-ym)
        r1a = xs*vals[0][0] + (1-xs)*vals[1][0]
        r1b = xs*vals[0][1] + (1-xs)*vals[1][1]
        r2  = ys*r1a + (1-ys)*r1b
        return r2
    
    def calc_V(self, root="all"):
        """"
        Takes "liq" or "vap" as an argument, to determine which root; default returns all roots.
        Uses values of Z computed previously if available; otherwise runs a Z calculation.
        """
        if not self.calc_finished:
            self.calc_Z()
        if self.pint:
            self.V = np.array(self.Z_sol)*muc.R*self.T/self.P
        else:
            self.V = np.array(self.Z_sol)*8.3145*self.T/self.P
        if root == "vap":
            return self.V[-1]
        elif root == "liq":
            return self.V[0]
        else:
            return self.V

    def calc_residSH(self):
        """
        Uses data as already passed to work out S_resid/R and H_resid/RT.
        Uses the Z value which should correspond to vapor phase.
        """
        if not self.calc_finished:
            self.calc_Z()
        if self.kind in self.cubic:
            sig = self.cnst[0]
            eps = self.cnst[1]
            beta = self.calc_beta(self.cnst[2], self.Pr.magnitude, self.Tr.magnitude)
            q = self.calc_q()
            Z = self.Z_sol[-1]

            if self.kind == "RK":
                dlnaTr = -.5
            else:
                raise ValueError("EOS kind residuals not implemented yet")
            # TODO
            if self.kind == "vdW":
                I = beta/(Z+eps*beta)
            else:
                I = 1/(sig+eps)*np.log((Z+sig*beta)/(Z+eps*beta))
            self.HRRT = Z - 1 + (dlnaTr -1)*q*I
            self.SRR  = np.log(Z-beta) + dlnaTr*q*I
            
        else:
            raise ValueError("EOS kind residuals not implemented yet")
        if self.pint:
            self.HR = self.HRRT*muc.R*self.T
            self.SR = self.SRR *muc.R

    def calc_residG(self, root="liq"):
        """
        Uses data as already passed to work out G_resid/RT.
        Reference: Lecture 16, Eq. 6.58, Eq. 3.51, Eq. 3.44, Z equation not in book
        """
        if not self.calc_finished:
            self.calc_Z()
        if root == "liq":
            Z = np.min(self.Z_sol)
        elif root == "vap":
            Z = np.max(self.Z_sol)
        else:
            raise ValueError("Invalid root type given.")
        if self.kind in self.cubic:
            sig = self.cnst[0]
            eps = self.cnst[1]
            q = self.calc_q()
            b = self.Omega*muc.R*self.Tc/self.Pc
            self.b = b
            rho = self.P/Z/muc.R/self.T
            rho.ito("mol/m**3")
            
            def rho_Z_int(rho):
                r1 = 1/(1-rho*b)
                r2 = q*rho*b / (1 + eps*rho*b) / (1 + sig*rho*b)
                return (r1 - r2)

            # def rho_Z_integrated(rho):
            #     r1 = -np.log(1-rho*b)
            #     r2 = q*(np.log(1+eps*rho*b)-np.log(1+sig*rho*b))/(eps-sig)
            #     return r1-r2

            
            rho_unit = rho.units
            rho_int, err = quad(lambda rh: (rho_Z_int(rh*rho_unit)-1)/rh, 0, rho.magnitude)
            # rho_int = rho_Z_integrated(rho) - rho_Z_integrated(0*rho_unit)
            # print(rho_int)
            self.GRRT = rho_int + Z - 1 - np.log(Z)
        else:
            raise ValueError("EOS kind residuals not implemented yet")
        if self.pint:
            self.GR = self.GRRT*muc.R*self.T
        else: 
            self.GR = self.GRRT*muc.R*self.T
        return self.GR
        
    def calc_fugacity(self, root="vap"):
        self.calc_residG(root=root)
        phi = np.exp(self.GRRT)
        self.f = phi * self.P
        return self.f

        #-----------------------------------
    # Leftovers from ChEn 273

    

class calc:

    @staticmethod
    def book_CpR(T, ABCD, pint=True):
        """
        Provide T in K.
        Takes coefficients ABCD in an iterable and unpacks them; should all be dimensionless.
        Equation form: C_p/R = A + BT + CT^2 + DT^-2 with B*10^-3, C*10^-6, D*10^5 evaluated internally
        Optional bool pint=True: if True, divides by muc.uReg.K to enforce dimensional consistency
        Returns C_p/R, to simplify units problems
        """
        A, B, C, D = ABCD
        if pint:
            return A + B*T/un.K*1E-3 + C*T*T/un.K/un.K*1E-6 + D/T/T*un.K*un.K*1E5
        else:
            return A + B*T*1E-3 + C*T*T*1E-6 + D/T/T*1E5

    @staticmethod
    def DS_ig(T1, T2, P1, P2, Cp, pint=True):
        """
        Cp a function, takes T and returns Cp
        pint: tells whether to strip units or not
        """
        if pint:
            t_part = quad(lambda t: (Cp(t*un.K)/t).magnitude, T1.to("K").magnitude, T2.to("K").magnitude)[0] * Cp(T1).units
            p_part = muc.R*np.log(P2/P1)
        else:
            t_part = quad(lambda t: Cp(t)/t, T1, T2)[0]
            p_part = 8.3145**np.log(P2/P1)
        return t_part - p_part

    @staticmethod
    def DH_ig(T1, T2, Cp, pint=True):
        if pint:
            h_unit = (Cp(T1)*T1).units
            t_part = quad(lambda t: Cp(t*un.K).to(Cp(T1).units).magnitude, T1.to("K").magnitude, T2.to("K").magnitude)[0]*(h_unit)
        else:
            t_part = quad(lambda t: Cp(t), T1, T2)
        return t_part 

    @staticmethod
    def B(T, Tc, omega):
        """
        Computes second virial coefficient B.
        Takes T, Tc, Pc, omega, as arguments; all should be in absolute units
        Optional argument: R, to match units being used. Defaults to m^3*Pa/K/mol
        Tr: reference temperature
        Tc, Pc: critical temperature, pressure
        omega: "Pitzer acentric factor", deals with geometry and polarity
        """
        Tr = T / Tc
        B0 = 0.083 - .422/(Tr**1.6)
        B1 = 0.139 - .172/(Tr**4.2)
        return B0 + omega*B1

    @staticmethod
    def dBdTr(T, Tc, omega):
        """
        Computes second virial coefficient B.
        Takes T, Tc, Pc, omega, as arguments; all should be in absolute units
        Optional argument: R, to match units being used. Defaults to m^3*Pa/K/mol
        Tr: reference temperature
        Tc, Pc: critical temperature, pressure
        omega: "Pitzer acentric factor", deals with geometry and polarity
        """
        Tr = T / Tc
        dB0dTr =  .675/(Tr**2.6)
        dB1dTr =  .722/(Tr**5.2)
        return dB0dTr + omega*dB1dTr

    @staticmethod
    def B_resid(T, Tc, P, Pc, omega):
        Tr = T/Tc
        Pr = P/Pc
        HRRT = Pr * (calc.B(T, Tc, omega) - Tr*calc.dBdTr(T, Tc, omega))
        SRR = -Pr * (calc.dBdTr(T, Tc, omega))
        return HRRT, SRR

    @staticmethod
    def Tsat(T_guess, P, Tc, Pc, omega, pint=True):
        mat_prop = (Tc, Pc, omega)
        if pint:
            def sol_T(T):
                T = un.Quantity(T, un.K)
                state = EOS("PR", T, P, *mat_prop, pint)
                state.calc_Z()
                vap = state.calc_fugacity(root="vap")
                liq = state.calc_fugacity(root="liq")
                return (vap - liq).magnitude
            T_sol = fsolve(sol_T, T_guess.magnitude)[0]*T_guess.units
        else:
            def sol_T(T):
                state = EOS("PR", T, P, *mat_prop, pint)
                state.calc_Z(Z_guess)
                vap = state.calc_fugacity(root="vap")
                liq = state.calc_fugacity(root="liq")
                return (vap - liq)
            T_sol = fsolve(sol_T, T_guess)[0]
        return T_sol

    @staticmethod
    def Psat(T, P_guess, Tc, Pc, omega, pint=True):
        mat_prop = (Tc, Pc, omega)
        if pint:
            def sol_P(P):
                P = un.Quantity(P, Pc.units)
                state = EOS("PR", T, P, *mat_prop, pint)
                state.calc_Z()
                vap = state.calc_fugacity(root="vap")
                liq = state.calc_fugacity(root="liq")
                return (vap - liq).magnitude
            P_sol = fsolve(sol_P, P_guess.magnitude)[0]*P_guess.units
        else:
            def sol_P(P):
                state = EOS("PR", T, P, *mat_prop, pint)
                state.calc_Z(Z_guess)
                vap = state.calc_fugacity(root="vap")
                liq = state.calc_fugacity(root="liq")
                return (vap - liq)
            P_sol = fsolve(sol_P, P_guess)[0]
        return P_sol

    @staticmethod
    def DIPPR_Psat(T, coeff):
        A, B, C, D, E = coeff
        T = T.to(un.K).magnitude
        Y = np.exp(A + B/T + C*np.log(T) + D*T**E)
        return Y*un.Pa
    
    @staticmethod
    def Bmix1(nFracs, Bvals):
        """
        Probably outdated; from a ChEn 273 problem.
        Mixing rule: sum_i(sum_j(yi*yj*Bij)), with Bii pure species, Bij = .5(Bii+Bjj)
        Takes two same-length iterables as arguments: nFracs, Bvals
        nFracs: list of mole fractions by species
        Bvals: list of B by pure species (second virial coefficient), matched with nFracs
        """
        if len(nFracs) != len(Bvals):
            raise Exception("Number of mole fractions does not equal number of B coefficients")
        num = len(nFracs) # number of species, effectively
        # List comprehensions!
        Bij = [[.5 * (Bvals[i] + Bvals[j]) for j in range(num)] for i in range(num)]
        Bmix = 0
        for i in range(len(nFracs)):
            for j in range(len(nFracs)):
                Bmix += nFracs[i]*nFracs[j]*Bij[i][j]
        return Bmix

    @staticmethod
    def multi_reac_equil(rxn_list, Ka_list, ext_guess, n0):
        """
        Function for calculating equilibrium of multiple reactions.
        Arguments: rxn_list, Ka_list, ext_guess, n0
        rxn_list is a list of reac_equil objects, which already have everything set (phase, activity) and have the same species in common (just different activities).
        Ka_list is a list of Ka values for each reaction, at the appropriate temperature.
        ext_guess is a list of guess values for the extent of reaction for each reaction. Pass without units; moles assumed.
        n0 is the overall inlet feed, as a dictionary by species with mole units.
        """
        def sol_ee(ext_guess):
            nn = {}
            names = rxn_list[0].names
            for nm in names:
                nn[nm] = n0[nm] + sum([rx.nu[nm] * ex * un.mol for rx, ex in zip(rxn_list, ext_guess)])
            Qa_list = np.array([rxn.calc_Qa_nn(nn, split=True) for rxn in rxn_list])
            Qap_list = Qa_list[:,0]
            Qar_list = Qa_list[:,1]
            eq_list = [(Ka/Qar - Qap).magnitude for Ka, Qap, Qar in zip(Ka_list, Qap_list, Qar_list)]
            return eq_list
        ext_list = fsolve(sol_ee, ext_guess) *un.mol
        return ext_list

class mixfit:
    def __init__(self, x1_arr, M_arr, curve_func="Not Given" ):
        """
        Takes x1_arr, M_arr corresponding to the mole fraction and total molar property.
        Will optionally take a curve_func, for example a quadratic. Not yet implemented. Currently defaults to a quadratic fit.
        """
        if curve_func == "Not Given":
            self.curve_func = self.default_curve
            self.nParams=3
        else:
            self.curve_func = curve_func
            print("Warning: system is not currently set up to handle other curve fits.")
        self.M1 = M_arr[-1]
        self.M2 = M_arr[0]
        self.fit, self.covar = curve_fit(lambda x1, p: self.M_nondim(x1, *p), x1_arr, M_arr.magnitude, p0=np.zeros(self.nParams))

    @staticmethod
    def default_curve(x1, a, b, c):
        return (a + b*x1 + c*x1*x1)

    # This function needs to have a function signature corresponding to the curve_func passed.
    # I would like to automate that better, perhaps with the help of curve_fit.
    def M_nondim(self, x1, param):
        x2 = 1-x1
        return x1*self.M1.magnitude + x2*self.M2.magnitude + x1*x2*self.curve_func(x1, *param)

    @un.wraps("mL/mol", (None, ""), strict=False)
    def calc_M(self, x1):
        return self.M_nondim(x1, *self.fit)
    def M_id_nondim(self, x1):
        return x1*self.M1.magnitude + (1-x1)*self.M2.magnitude
    @un.wraps("mL/mol", (None, ""), strict=False)
    def calc_M_id(self, x1):
        return self.M_id_nondim(x1)
    def dMdx_nondim(self, x1):
        return muc.ddx(lambda x: self.M_nondim(x, *self.fit), x1, .1)
    @un.wraps("mL/mol", (None, ""), strict=False)
    def calc_dMdx(self, x1):
        return self.dMdx_nondim(x1)
    def M1b_nondim(self, x1):
        M = self.M_nondim(x1, *self.fit)
        dMdx = self.dMdx_nondim(x1)
        x2 = 1-x1
        return M + x2*dMdx
    @un.wraps("mL/mol", (None, ""), strict=False)
    def calc_M1b(self, x1):
        return self.M1b_nondim(x1)
    def M2b_nondim(self, x1):
        M = self.M_nondim(x1, *self.fit)
        dMdx = self.dMdx_nondim(x1)
        x2 = 1-x1
        return M - x1*dMdx
    @un.wraps("mL/mol", (None, ""), strict=False)
    def calc_M2b(self, x1):
        return self.M2b_nondim(x1)

class Activity:
    def __init__(self, kind, params, phase="l", index=1):
        """
        kind: kind of activity model
            For phase='l' : 'mrg1', 'mrg2', 'Wilson', 'WilsonLL', 'vanLaar'
            For phase='s' : '1'
            For phase='g' : 'ig', 'im' . 'ig' indicates ideal gas, 'im' indicates ideal mixture but nonideal gas
        params: single parameter or tuple of parameters to unpack. 
            mrg1: A
            mrg2: A12, A21
            Wilson: a12, a21, V1, V2. #L12 and L21 are computed by calc_gamma12. (L for \Lambda)
            WilsonLL: L12, L21 (L for \Lambda)
            vanLaar: A12', A21'
            ig: P
            im: *either* fi0, fugacity, *or* fugacity_func(T,P). Checks if callable to decide which.
        phase: defaults to l, for compatibility with older code.
        index: defaults to 1. Used for binary liquid activities: should be either 1 or 2. At present assumes Poynting factor is 1.
        No returns.
        """
        self.kind = kind
        self.args = (kind, params)
        self.phase = phase
        self.index = index
        if self.phase == "l":
            if self.kind == "mrg1":
                self.A = params
                self.calc_gamma12 = self.calc_gamma12_mrg1
            elif self.kind == "mrg2":
                self.A12, self.A21 = params
                self.calc_gamma12 = self.calc_gamma12_mrg2
            elif self.kind == "Wilson":
                self.V1, self.V2, self.a12, self.a21 = params
                self.calc_gamma12 = self.calc_gamma12_Wilson
            elif self.kind == "WilsonLL":
                self.L12, self.L21 = params
                self.calc_gamma12 = self.calc_gamma12_WilsonLL
            elif self.kind == "vanLaar":
                self.A12p, self.A21p = params
                self.calc_gamma12 = self.calc_gamma12_vanLaar
            else:
                raise ValueError("Invalid type of activity model.")
            self.calc_a = self.calc_a_liq
        elif self.phase == "g":
            if self.kind == "ig":
                self.calc_a = self.calc_a_ig
                self.P = params
            elif self.kind == "im":
                self.calc_a = self.calc_a_im
                self.f_func = callable(params)
                if self.f_func:
                    self.calc_f = params
                else:
                    self.fi0 = params
            else:
                raise ValueError("Invalid type of activity model.")
        elif self.phase == "s":
            if self.kind == "1":
                self.calc_a = self.calc_a_s1
            else:
                raise ValueError("Invalid type of activity model.")
        else:
            raise ValueError("Invalid phase passed to initialize activity model.")

    def calc_a(self):
        print("No function for activity was set, but it was called.")
        raise NotImplementedError("Activity model not yet implemented: " + self.phase + self.kind)

    # --------------------------------------------------
    # Activity calculations (a, not gamma)
    # Unlike the gamma calculations below, these do not assume a binary mixture, so they are performed
    # for an individual species.
    def calc_a_s1(self, xi, T="Not Used", P="Not Used"):
        return 1
    
    def calc_a_ig(self, xi, T="Not Used", P="Not Used"):
        a = xi*self.P/(1*un.bar)
        return a

    def calc_a_im(self, xi, T="Not Given", P="Not Given"):
        if not self.f_func:
            a = xi*self.fi0/(1*un.bar)
            return a
        else:
            if T=="Not Given" or P=="Not Given":
                raise ValueError("Missing either T or P for calculating fugacity.")
            f = self.calc_f(T, P)
            a = xi*f/(1*un.bar)
            return a

    def calc_a_liq(self, xi, T="Not Given", P="Not Used"):
        if self.index ==1:
            x1 = xi
            x2 = 1-x1
        elif self.index == 2:
            x2 = xi
            x1 = 1-x2
        gamma1, gamma2 = self.calc_gamma12(x1, T=T)
        Poynting = 1
        if self.index == 1:
            a = xi*gamma1*Poynting
        elif self.index == 2:
            a = xi*gamma2*Poynting
        else:
            raise ValueError("Index other than 1 or 2 passed for a binary mixture of liquids.")
        return a
        
            
    # ---------------------------------------------------------
    # Activity coefficient models for binary mixture
    # All functions return both gamma1 and gamma2.
    def calc_gamma12_mrg1(self, x1, T="Not Given"):
        x2 = 1-x1
        gam1 = np.exp(x2*x2*self.A)
        gam2 = np.exp(x1*x1*self.A)
        return gam1, gam2

    def calc_gamma12_mrg2(self, x1, T="Not Given"):
        x2 = 1-x1
        gam1 = np.exp(x2*x2 * (self.A12 + 2*(self.A21 - self.A12)*x1))
        gam2 = np.exp(x1*x1 * (self.A21 + 2*(self.A12 - self.A21)*x2))
        return gam1, gam2

    def calc_gamma12_Wilson(self, x1, T):
        # if T == "Not Given":
        #     raise ValueError("Wilson VLE needs a temperature for gamma.")
        x2 = 1-x1
        L12, L21 = self.calc_Wilson_LL(self.V1, self.V2, self.a12, self.a21, T)
        der = (L12/(x1 + x2*L12) - L21/(x2 + x1*L21))
        ret1 = x1 + x2*L12
        ret2 = x2 + x1*L21
        gam1 = np.exp(x2*der)/ret1
        gam2 = np.exp(-x1*der)/ret2
        return gam1, gam2
    @staticmethod
    def calc_Wilson_LL(V1, V2, a12, a21, T):
        L12 = V2/V1 * np.exp(-a12/muc.R/T)
        L21 = V1/V2 * np.exp(-a21/muc.R/T)
        return L12, L21
    def calc_gamma12_WilsonLL(self, x1, T="Not Given"):
        x2 = 1-x1
        L12 = self.L12
        L21 = self.L21
        der = (L12/(x1 + x2*L12) - L21/(x2 + x1*L21))
        ret1 = x1 + x2*L12
        ret2 = x2 + x1*L21
        gam1 = np.exp(x2*der)/ret1
        gam2 = np.exp(-x1*der)/ret2
        return gam1, gam2

    def calc_gamma12_vanLaar(self, x1, T="Not Given"):
        x2 = 1-x1
        gam1 = np.exp(self.A12p *(1 + self.A12p*x1/self.A21p/x2)**-2)
        gam2 = np.exp(self.A21p *(1 + self.A21p*x2/self.A12p/x1)**-2)
        return gam1, gam2

    def calc_GERT(self, *args):
        x2 = 1-x1
        gam1, gam2 = self.calc_gamma12(*args)
        return x1*np.log(gam1) + x2*np.log(gam2)
    def calc_GmixRT(self, x1):
        x2 = 1-x1
        GERT = self.calc_GERT(x1)
        term = x1*np.log(x1) + x2*np.log(x2)
        return GERT + term

class vle:
    def __init__(self, kind, T="Not Given", P="Not Given"):
        """
        Class to handle binary VLE calculations. To generate a Txy diagram, pass P; for a Pxy diagram, pass T.
        Allowed kinds are baby, teen, and adult (pass a string).
        If using teen or adult, next call vle.set_act_model.
        For all, call vle.set_Psat .
        """
        self.kind = kind
        self.Tbool=True
        self.Pbool=True
        self.T = T
        self.P = P
        if self.T == "Not Given":
            self.Tbool = False
        if self.P == "Not Given":
            self.Pbool = False
        self.x1_arr = np.linspace(0, 1)

        if self.kind == "baby":
            self.calc_P = self.calc_P_baby
        elif self.kind == "teen":
            self.calc_P = self.calc_P_teen
        elif self.kind == "adult":
            self.calc_P = self.calc_P_adult
        else:
            raise ValueError("Incorrect VLE kind given. Should be 'baby', 'teen', or 'adult'.")

    def set_Psat(self, func1, func2):
        """
        Takes two callables. Each should return Psat for the fluid.
        Proper function signature: func(T) 
        """
        self.Psat1 = func1
        self.Psat2 = func2
            
    def set_act_model(self, kind, params):
        """
        kind: 'mrg1', 'mrg2', or 'Wilson'.
        params: tuple of parameters to unpack. 
            For mrg1: A
            For mrg2: A12, A21
            For Wilson: a12, a21, V1, V2. L12 and L21 are computed by calc_gamma12. (L for \Lambda)
        No returns.
        Wraps the Activity initialization function, and calls it self.act.
        """
        if self.kind == "baby":
            print("Warning: Setting an unused activity model for VLE.")
        self.act_kind = kind
        self.act = Activity(self.act_kind, params)

    def calc_P_baby(self, x1, T, calc_y=False):
        """
        Takes x1 and T. Optional parameter calc_y=False.
        If calc_y is True, also computes y1 and returns (P_tot, y1).
        If calc_y is False, returns P_tot.
        """
        x2 = 1-x1
        P1 = x1*self.Psat1(T)
        P2 = x2*self.Psat2(T)
        if calc_y:
            return P1 + P2, P1/(P1 + P2)
        else:
            return P1 + P2
    def calc_P_teen(self, x1, T, calc_y=False):
        """
        Takes x1 and T. Optional parameter calc_y=False.
        If calc_y is True, also computes y1 and returns (P_tot, y1).
        If calc_y is False, returns P_tot.
        """
        x2 = 1-x1
        gam1, gam2 = self.act.calc_gamma12(x1, T)
        P1 = x1*self.Psat1(T) * gam1
        P2 = x2*self.Psat2(T) * gam2
        if calc_y:
            return P1 + P2, P1/(P1 + P2)
        else:
            return P1 + P2

    def calc_Pxy(self, numPoints=101, x1=.5, T="Not Given", xspan=(0,1)):
        if T == "Not Given" and not self.Tbool:
            raise ValueError("Cannot perform Pxy calc without T. Either initialize VLE with one, or pass to calc_Pxy.")
        elif T == "Not Given":
            T = self.T
        if numPoints == 1:
            P, y1 = self.calc_P(x1, T, True)
            return P, x1, y1
        else:
            x1_arr = np.linspace(*xspan, numPoints)
            P_arr, y1_arr = self.calc_P(x1_arr, T, True)
            return P_arr, x1_arr, y1_arr

    # This function was, very mysteriously, crashing the Jupyter kernel without throwing any Python errors.
    # The end result is that I run a single fsolve across the entire x1 array, instead of individually.
    # I do not know why this works and the alternatives (commented out below) did not.
    def calc_Txy(self, numPoints=101, x1=.5, P="Not Given", xspan=(0,1)):
        if P == "Not Given" and not self.Pbool:
            raise ValueError("Cannot perform Txy calc without P. Either initialize VLE with one, or pass to calc_Txy.")
        elif P == "Not Given":
            P = self.P
        Tguess = self.T if self.Tbool else 300*un.K
        if numPoints == 1:
            T = fsolve(lambda t: (self.calc_P(x1, t*Tguess.units, False) - P).magnitude, Tguess.magnitude)[0]*Tguess.units
            P, y1 = self.calc_P(x1, T, True)
            return T, x1, y1
        else:
            x1_arr = np.linspace(*xspan, numPoints)
            T_arr = fsolve(lambda t: (self.calc_P(x1_arr, t*Tguess.units, False) - P).magnitude, x1_arr+Tguess.magnitude)*Tguess.units
            P_arr, y1_arr = self.calc_P(x1_arr, T_arr, True)
            return T_arr, x1_arr, y1_arr


            # T_list = []
            # P_list = []
            # y_list = []
            # for x1 in x1_arr:
            #     T = fsolve(lambda t: (self.calc_P(x1, t*Tguess.units, False) - P).magnitude, Tguess.magnitude)[0]*Tguess.units
            #     P, y1 = self.calc_P(x1, T, True)
            #     T_list.append(T)
            #     P_list.append(P)
            #     y_list.append(y1)
            # T_arr = muc.list_unit(T_list)
            # P_arr = muc.list_unit(P_list)
            # y1_arr = muc.list_unit(y_list)
            # return T_arr, x1_arr, y1_arr
                
            # def sol_T(t, x1):
            #     t *= Tguess.units
            #     P_RHS = self.calc_P(x1, t, False)
            #     return (P - P_RHS).magnitude
            # T_arr = [fsolve(lambda t: sol_T(t, x1), Tguess.magnitude)[0]*Tguess.units for x1 in x1_arr]
            # T_arr = muc.list_unit(T_arr)
            # P_arr, y1_arr = self.calc_P(x1_arr, T_arr, True)
            # return T_arr, x1_arr, y1_arr

    def calc_bbl_yP(self, x1, guess, T="Not Given"):
        """
        Takes x1, a guess for the solver, and optionally a value of T to use in the calculation.
        The guess has the form (y1, P), without any units attached.
        Returns y1, P .
        """
        if T == "Not Given":
            T = self.T
        P, new_x, y1 = self.calc_Pxy(1, x1, T)
        return y1, P
    def calc_bbl_yT(self, x1, guess, P="Not Given"):
        """
        Takes x1, a guess for the solver, and optionally a value of P to use in the calculation.
        The guess has the form (y1, T), without any units attached.
        Returns y1, T.
        """
        if P == "Not Given":
            P = self.P
        T, new_x, y1 = self.calc_Txy(1, x1, P)
        return y1, T

    def calc_dew_xP(self, y1, guess, T="Not Given"):
        """
        Takes x1, a guess for the solver, and optionally a value of T to use in the calculation.
        The guess has the form (y1, P), without any units attached.
        Returns x1, P.
        """
        if T == "Not Given":
            T = self.T
        def sol_xP(x1P):
            x1, P = x1P
            P_RHS, y1_RHS = self.calc_P(x1, T, True)
            eq1 = (P*P_RHS.units) - P_RHS 
            eq2 = y1 - y1_RHS
            return eq1.magnitude, eq2.magnitude
        xP_sol = fsolve(sol_xP, guess)
        x1 = xP_sol[0]
        P = self.calc_P(x1, T, False)
        err = sol_xP((x1, P.magnitude))
        P_err = (err[0]*P.units) / P
        y_err = (err[1])/y1
        if np.abs(P_err.to("").magnitude) > .001 or np.abs(y_err) > .001:
            print(f"Dew P calc error: {P_err.magnitude*100:.1f}% in P, {y_err*100:.1f}% in y")
        # P = Px_sol[1]*P_RHS.units
        return x1, P
    def calc_dew_xT(self, y1, guess, P="Not Given"):
        """
        Takes x1, a guess for the solver, and optionally a value of T to use in the calculation.
        The guess has the form (y1, P), without any units attached.
        Returns x1, T.
        """
        if P == "Not Given":
            P = self.P
        def sol_xT(x1T):
            x1, T = x1T
            T *= un.K
            P_RHS, y1_RHS = self.calc_P(x1, T, True)
            eq1 = P - P_RHS
            eq2 = y1 - y1_RHS
            return eq1.magnitude, eq2.magnitude
        xT_sol = fsolve(sol_xT, guess)
        x1 = xT_sol[0]
        T = xT_sol[1]*un.K
        err = sol_xT((x1, T.magnitude))
        T_err = (err[0]*P.units) / P
        y_err = (err[1])/y1
        if np.abs(T_err.to("").magnitude) > .001 or np.abs(y_err) > .001:
            print(f"Dew T calc error: {T_err.magnitude*100:.1f}% in P, {y_err*100:.1f}% in y")
        return x1, T

    def calc_flash(self, z1, x_guess):
        """
        Uses the object's values of T and P to perform a flash calculation.
        """
        def sol_x(x):
            P_RHS = self.calc_P(x, self.T, False)
            eq1 = self.P - P_RHS
            return eq1.magnitude
        x1 = fsolve(sol_x, x_guess)[0]
        P_RHS, y1 = self.calc_P(x1, self.T, True)
        x2 = 1-x1
        y2 = 1-y1
        if not (x1<z1 and z1<y1) and not (y1<z1 and z1<x1):
            print("Flash calculation found an unphysical x1 to get the right P: z1 is not between x1 and y1")
        if x1 < 0 or x2 < 0 or y1 < 0 or y2 < 0:
            print("Flash calculation found an unphysical x1 to get the right P.")
            print(f"Calculated P: {P_RHS:.1f}")
        l_frac = (z1-y1)/(x1-y1)
        v_frac = 1-l_frac
        ret = {
               "l_comp":(x1,x2),
               "v_comp":(y1.magnitude,y2.magnitude),
               "l_frac":l_frac,
               "v_frac":v_frac
        }
        return ret
        
class lle:
    def __init__(self, act_kind, act_params):
        self.act_kind = act_kind
        self.act_params = act_params
        self.act = Activity(self.act_kind, self.act_params)
    
    def calc_equilibrium(self, guess=(.1, .9), T="Not Given"):
        """
        Takes guess and optionally T.
        Guess has form (x1a, x1b).
        Returns x1a, x1b.
        """
        if T=="Not Given" and self.act_kind == "Wilson":
            raise ValueError("For Wilson activity, need a given temperature.")
        def sol_eqs(x1ab):
            x1a, x1b = x1ab
            gam1a, gam2a = self.act.calc_gamma12(x1a, T)
            gam1b, gam2b = self.act.calc_gamma12(x1b, T)
            x2a = 1-x1a
            x2b = 1-x1b
            eq1 = x1a*gam1a - x1b*gam1b
            eq2 = x2a*gam2a - x2b*gam2b
            return eq1, eq2
        x1a, x1b = fsolve(sol_eqs, guess)
        if (x1a < 0) or (x1a > 1) or (x1b<0) or (x1b>1):
            print("Solver found an unphysical x1a or x1a:", x1a, x1b)
        return x1a, x1b
        
class vlle:
    def __init__(self, act_kind, act_params, Psat1, Psat2, T="Not Given"):
        self.act_kind = act_kind
        self.act_params = act_params
        self.act = Activity(self.act_kind, self.act_params)
        
        self.Tbool=True
        self.T = T
        if self.T == "Not Given":
            self.Tbool = False

        self.lle = lle(*self.act.args)
        self.vlea = vle("teen", T=T)
        self.vleb = vle("teen", T=T)
        self.vlea.set_act_model(*self.act.args)
        self.vleb.set_act_model(*self.act.args)
        self.vlea.set_Psat(Psat1, Psat2)
        self.vleb.set_Psat(Psat1, Psat2)


    def calc_lle(self, guess=(.1, .9)):
        x1a, x1b = self.lle.calc_equilibrium(guess)
        self.x1a = x1a
        self.x1b = x1b
        return x1a, x1b

    def calc_Pys(self, T="Not Given"):
        """
        Assumes calc_lle has already been called.
        """
        if T == "Not Given" and not self.Tbool:
            raise ValueError("Cannot perform P* or y* calc without T. Either initialize VLLE with one, or pass to calc_Pys.")
        elif T == "Not Given":
            T = self.T
        x1a = self.x1a
        x1b = self.x1b
        x2a = 1-x1a
        x2b = 1-x1b
        Pa, y1a = self.vlea.calc_P(x1a, T, True)
        Pb, y1b = self.vleb.calc_P(x1b, T, True)
        self.Ps = (Pb*y1b) + (Pa*(1-y1a))
        self.y1s = Pa/self.Ps
        return self.Ps, self.y1s
    def calc_Pxy(self, numPoints=101, x1=.5, T="Not Given"):
        if T == "Not Given" and not self.Tbool:
            raise ValueError("Cannot perform Pxy calc without T. Either initialize VLLE with one, or pass to calc_Pxy.")
        elif T == "Not Given":
            T = self.T
        xa_num = int(np.ceil(100*self.x1a))+1
        xb_num = 101 - xa_num
        Pa_arr, x1a_arr, y1a_arr = self.vlea.calc_Pxy(numPoints=xa_num, xspan=(1e-6,self.x1a))
        Pb_arr, x1b_arr, y1b_arr = self.vleb.calc_Pxy(numPoints=xb_num, xspan=(self.x1b, 1-1e-6))
        
        P_arr = np.concatenate((Pa_arr, Pb_arr))
        x1_arr = np.concatenate((x1a_arr, x1b_arr))
        y1_arr = np.concatenate((y1a_arr, y1b_arr))
        return P_arr, x1_arr, y1_arr


# Class loosely based on form of mixrn.Reaction
class reac_equil:
    """
    As of 2 April 2020, class is loosely based on mixrxn.Reaction class, but not actually dependent.
    To calculate an equilibrium, you need to first call:
    set_n0, set_phases, set_act_model, set_G0, [set_H0]
    set_n0 is called once, with a list.
    set_phases is called once, with a list.
    set_act_model is called for each species.
    set_G0 is called once, with G0 and Tref.
    set_H0 is called once, with H0 and Tref.

    Available calculations:
    calc_Qa(xi, TP)
    calc_Ka(T='Not Given', DCpR_func='Not Given')
    calc_ext(ext_guess) runs a solver to compute equilibrium
    calc_nn_phase(ext)
    calc_nn(ext)
    calc_nfrac(ext)
    calc_X() assumes an extent already calculated

    For multi-reaction equilibrium, use this class with thermo.calc.multi_reac_equil.
    """
    def __init__(self, names, nus):
        """
        Takes list of stoich. coefficients and list of species names. Stores values.
        """
        self.names = names
        self.nspec = len(names)
        if len(nus) != self.nspec:
            raise ValueError("Wrong number of nu values for number of species.")
        self.act = {}
        self.nu = {}
        for i, n in enumerate(names):
            self.nu[n] = nus[i]

    def set_extra_reac(self, names, nus):
        """

        """

    def set_phases(self, phases):
        if len(phases) != self.nspec:
            raise ValueError("Wrong number of phases for number of species.")
        self.phase = {}
        for n, p in zip(self.names, phases):
            self.phase[n] = p

    def set_n0(self, n_list):
        if len(n_list) != self.nspec:
            raise ValueError("Wrong number of mole values for number of species.")
        self.nn0 = {}
        for nm, nn in zip(self.names, n_list):
            self.nn0[nm] = nn
        self.n0 = sum(n_list)

    def set_act_model(self, spec, kind, params, phase="Not Given", index="Not Given"):
        """
        Arguments: species name, phase, params
        Phase: either 's', 'g', or 'l'
        For kinds and params, see thermo.Activity class
        """
        if phase == "Not Given":
            phase = self.phase[spec]
        self.act[spec] = Activity(kind, params, phase, index)
        
    def set_G0(self, G0rxn, Tref):
        self.G0rxn = G0rxn
        self.TGrf = Tref

    def set_H0(self, H0rxn, Tref):
        self.H0rxn = H0rxn
        self.THrf = Tref

    def calc_Ka(self, T="Not Given", DCpR_func="Not Given"):
        self.K0 = np.exp(-self.G0rxn/muc.R/self.TGrf)
        if T == "Not Given":
            return self.K0
        self.K1 = np.exp(self.H0rxn/muc.R/self.THrf*(1-self.THrf/T))
        if DCpR_func == "Not Given":
            return self.K0 * self.K1
        int1 = quad(lambda t: (DCpR_func(t*un.K)).magnitude, self.THrf.to("K").magnitude, T.to("K").magnitude)[0]
        int2 = quad(lambda t: (DCpR_func(t*un.K)/t).magnitude, self.THrf.to("K").magnitude, T.to("K").magnitude)[0]
        self.K2 = np.exp(-1/T*int1*T.units + int2)
        return self.K0 * self.K1 * self.K2

    def calc_nn_phase(self, ext):
        nn = {}
        n = 0
        n_phase = {"g":0, "s":0, "l":0}
        for nm in self.names:
            nn[nm] = self.nn0[nm] + self.nu[nm]*ext
            n += nn[nm]
            n_phase[self.phase[nm]] += nn[nm]
        return nn, n_phase

    def calc_nfrac(self, ext=None):
        """
        If no extent is given, uses the value from the last extent calculation
        """
        if ext is None:
            ext = self.ext
        nn, n_phase = self.calc_nn_phase(ext)
        nfrac = {}
        for nm in self.names:
            # Compute separate mole fractions for each phase
            nfrac[nm] = nn[nm] / n_phase[self.phase[nm]]
        return nfrac

    def calc_Qa(self, ext, TP=None, debug=False, split=False):
        """
        Takes ext (extent of reaction), T, P.
        Returns Qa, the product of all activities raised to stoichiometric coefficients.
        Optional argument: split. Defaults False; if True, returns Qa_prod and Qa_reac separately.
            (Helpful for making equations converge easier.)
        """
        if TP is None:
            passTP = False
        else:
            passTP = True
        a = {}
        Qa = 1
        Qa_prod = 1
        Qa_reac = 1
        nfrac = self.calc_nfrac(ext)
        for nm in self.names:
            # nfrac[nm] = nn[nm] / n
            if passTP:
                a[nm] = self.act[nm].calc_a(nfrac[nm], *TP)
            else:
                a[nm] = self.act[nm].calc_a(nfrac[nm])
            Qa *= a[nm]**self.nu[nm]
            if split:
                if self.nu[nm] < 0:
                    Qa_reac *= a[nm]**self.nu[nm]
                elif self.nu[nm] > 0:
                    Qa_prod *= a[nm]**self.nu[nm] 
                else:
                    pass

        # print(nfrac, n_phase)
        if not debug:
            if not split:
                return Qa
            else:
                return Qa_prod, Qa_reac
        else:
            return Qa, nn, nfrac, a

    def calc_Qa_nn(self, nn, TP=None, split=False, debug=False):
        """
        Takes nn (dict of molar amounts), T, P.
        Useful for manually constructing a separate solver, esp. for multiple reactions.
        Returns Qa, the product of all activities raised to stoichiometric coefficients.
        Optional argument: split. Defaults False; if True, returns Qa_prod and Qa_reac separately.
            (Helpful for making equations converge easier.)
        """
        if TP is None:
            passTP = False
        else:
            passTP = True
        a = {}
        Qa = 1
        Qa_prod = 1
        Qa_reac = 1
        nfrac = {}
        n_phase = {"g":0, "s":0, "l":0}
        for nm in self.names:
            n_phase[self.phase[nm]] += nn[nm]
        for nm in self.names:
            nfrac[nm] = nn[nm] / n_phase[self.phase[nm]]
            if passTP:
                a[nm] = self.act[nm].calc_a(nfrac[nm], *TP)
            else:
                a[nm] = self.act[nm].calc_a(nfrac[nm])
            Qa *= a[nm]**self.nu[nm]
            if split:
                if self.nu[nm] < 0:
                    Qa_reac *= a[nm]**self.nu[nm]
                elif self.nu[nm] > 0:
                    Qa_prod *= a[nm]**self.nu[nm] 
                else:
                    pass
        # print(nfrac, n_phase)
        if not debug:
            if not split:
                return Qa
            else:
                return Qa_prod, Qa_reac
        else:
            return Qa, nn, nfrac, a

    def calc_ext(self, ext_guess, TP=None, DCpR_func=None, pint_strip=True):
        if DCpR_func is not None:
            Ka = self.calc_Ka(TP[0], DCpR_func)
        elif TP is not None:
            Ka = self.calc_Ka(TP[0])
        else:
            Ka = self.calc_Ka()
        if pint_strip:
            def sol_ext(ext):
                ext *= ext_guess.units
                Qap, Qar = self.calc_Qa(ext, TP, split=True)
                # Ka = self.calc_Ka(T)
                return (Ka/Qap - Qar).magnitude
        else:
            def sol_ext(ext):
                Qap, Qar = self.calc_Qa(ext, TP, split=True)
                # Ka = self.calc_Ka(T)
                return (Ka/Qap - Qar)
        ext = fsolve(sol_ext, ext_guess.magnitude)[0]*ext_guess.units
        # ext = root(sol_ext, ext_guess.magnitude, method="anderson").x*ext_guess.units
        if ext == ext_guess:
            print("Careful, the fsolve for extent of reaction gave back the guess value.")
        self.ext = ext
        return ext
            
    def calc_X(self, spec):
        """
        Assumes calc_ext has already been called, and uses stored value for extent of reaction.
        """
        X = -self.nu[spec]*self.ext/self.nn0[spec]
        return X



