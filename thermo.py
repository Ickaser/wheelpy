import wheelpy.muc as muc
un = muc.uReg
R = muc.R
from scipy.optimize import fsolve, curve_fit
from scipy.integrate import quad
import numpy as np

class EOS:

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
    def book_Cp(T, ABCD, pint=True):
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


class mixfit:
    def __init__(self, x1_arr, M_arr, curve_func="Not Given" ):
        """
        Takes x1_arr, M_arr corresponding to the mole fraction and total molar property.
        Will optionally take a curve_func, for example a quadratic. Not yet implemented. Currently defaults to a quadratic fit.
        """
        if curve_func == "Not Given":
            self.curve_func = self.default_curve
        else:
            self.curve_func = curve_func
            print("Warning: system is not currently set up to handle other curve fits.")
        self.M1 = M_arr[-1]
        self.M2 = M_arr[0]
        self.fit, self.covar = curve_fit(self.M_nondim, x1_arr, M_arr.magnitude)

    @staticmethod
    def default_curve(x1, a, b, c):
        return (a + b*x1 + c*x1*x1)

    # This function needs to have a function signature corresponding to the curve_func passed.
    # I would like to automate that better, perhaps with the help of curve_fit.
    def M_nondim(self, x1, a, b, c):
        x2 = 1-x1
        return x1*self.M1.magnitude + x2*self.M2.magnitude + x1*x2*self.curve_func(x1, a, b, c)

    @un.wraps("mL/mol", (None, ""), strict=False)
    def calc_M(self, x1):
        return self.M_nondim(x1, *self.fit)
    def M_id_nondim(self, x1):
        return x1*self.M1.magnitude + (1-x1)*self.M2.magnitude
    @un.wraps("mL/mol", (None, ""), strict=False)
    def calc_M_id(self, x1):
        return self.M_id_nondim(x1)
    def dMdx_nondim(self, x1):
        # u0 = M1.magnitude - M2.magnitude + a
        # u1 = 2*b-2*a
        # u2 = 3*c-3*b
        # u3 = -4*c
        # polycoeff = [u3, u2, u1, u0]
        # return np.polyval(polycoeff, x1)
        return muc.ddx(lambda x: self.M_nondim(x, *self.fit), x1, .1)
    @un.wraps("mL/mol", (None, ""), strict=False)
    def calc_dMdx(self, x1):
        return self.dMdx_nondim(x1)
    def M1b_nondim(self, x1):
        M = self.M_nondim(x1, *self.fit)
        dMdx = self.dMdx_nondim(x1)
        x2 = 1-x1
        return M + x2*dMdx
        # u0 = M1.magnitude + a
        # u1 = 2*b-2*a
        # u2 = a-4*b+3*c
        # u3 = 2*b-6*c
        # u4 = 3*c
        # polycoeff = [u4, u3, u2, u1, u0]
        # return np.polyval(polycoeff, x1)
    @un.wraps("mL/mol", (None, ""), strict=False)
    def calc_M1b(self, x1):
        return self.M1b_nondim(x1)
    def M2b_nondim(self, x1):
        M = self.M_nondim(x1, *self.fit)
        dMdx = self.dMdx_nondim(x1)
        x2 = 1-x1
        return M - x1*dMdx
        # u0 = M2.magnitude
        # u1 = 0
        # u2 = a-b
        # u3 = 2*b-2*c
        # u4 = 3*c
        # polycoeff = [u4, u3, u2, u1, u0]
        # return np.polyval(polycoeff, x1)
    @un.wraps("mL/mol", (None, ""), strict=False)
    def calc_M2b(self, x1):
        return self.M2b_nondim(x1)



