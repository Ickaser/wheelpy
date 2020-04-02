# Once that exists, add nn and ny initialization options, then Fill should populate either direction
import wheelpy.muc as muc
un = muc.uReg

class Species:
    def __init__(self, name):
        self.name = name

    def set_Hf(self, Hf):
        self.Hf = Hf

    def set_Cp_const(self, Cp):
        self.Cp = Cp
        self.Cp_const = True

    def set_Cp_func(self, coeff, coeff_kind):
        self.Cp_kind = coeff_kind
        if self.Cp_kind == 1 and len(coeff) > 4:
                raise ValueError("Wrong number of coefficients")
        elif self.Cp_kind == 2 and len(coeff) != 3:
                raise ValueError("Wrong number of coefficients")
        elif self.Cp_kind == 3 and len(coeff) != 4:
                raise ValueError("Wrong number of coefficients")
        self.Cp_coeff = coeff
        self.Cp_kind = coeff_kind
        if self.Cp_kind == 1:
            self.Cp_order = len(self.Cp[self.names[0]])
        self.Cp_const = False

    def calc_DH(self, T1, T2, pint_strip = False):
        """
        Integrates C_p dT , for some common correlations.
        Arguments: T1, T2, pint_strip (bool).
        pint_strip indicates whether to strip units for temperature calc, always assumes in correct units
        kind 1: degree <=3 polynomial, coefficients in rising order. 
        kind 2: a + b*T + c*T**-2
        kind 3: a + b*T + c*T**2 + d*T**-2 (default)
        """
        if pint_strip:
            temp_unit = T1.units
            T1 = T1.magnitude
            T2 = T2.magnitude

        if self.Cp_const:
            Tarray = T2 - T1
        elif self.Cp_kind == 1:
            Tarray = [(T2**i - T1**i)/i for i in range(1, self.Cp["order"]+1)]
        elif self.Cp_kind == 2:
            Tarray = [T2-T1, (T2*T2 - T1*T1)/2, (1/T2 - 1/T1)*-1]
        elif self.Cp_kind == 3:
            Tarray = [T2-T1, (T2*T2 - T1*T1)/2, (T2**3 - T1**3)/3, (1/T2 - 1/T1)*-1]

        if self.Cp_const:
            DH = Tarray * self.Cp
        else:
            DH = sum([t*c for t,c in zip(Tarray, self.Cp_coeff)])
        if pint_strip:
            DH *= temp_unit
        return DH

class Mixture:
    
    def __init__(self, names, vals, mFlow = None, kind = "mx"):
        """
        Takes names, vals, total flow rate, and kind (either mx or mm).
        For mx, give total flow rate and fractions.
        For mm, give total flow rate and species flow rates.
        """
        
        self.names = names
        self.vals = vals
        self.initType = kind
        self.mFrac = {}
        self.mFlows = {}
        self.nFrac = {}
        self.nFlows = {}
        
        if self.initType == "mx":
            for u, i in zip(names, vals):
                self.mFrac[u] = i
            self.mFlow = mFlow
        elif self.initType == "mm":
            for u, i in zip(names, vals):
                self.mFlows[u] = i
            self.mFlow = mFlow
        elif self.initType == "nx":
            for u, i in zip(names, vals):
                self.nFrac[u] = i
            self.nFlow = mFlow
        elif self.initType == "nn":
            for u, i in zip(names, vals):
                self.nFlows[u] = i
            self.nFlow = mFlow

        self.s = {}
        for n in self.names:
            self.s[n] = Species(n)
            
                
    def fill(self, molar = False, pint=True):
        """
        Takes one optional argument, specifying whether to use MW function and compute molar flows
        If using MW feature, ensure that mass flows are in either g or kg.
        Returns nothing.
        Looks for None values, replaces them.
        ONLY DESIGNED TO HANDLE ONE UNKNOWN AT A TIME.
        """
        
        # -----------------------------
        self.molar_fill = molar
        trueType = self.initType
        changedType = False
        if self.initType == "nx" or self.initType == "nn":
            self.from_molar()
            changedType = True
        if self.initType == "mm":
            
            # compute total flow if necessary
            if self.mFlow == None:
                self.mFlow = sum(self.mFlows.values()) 
            # compute an individual mass flow if necessary
            else: 
                for i in self.mFlows.keys():
                    if self.mFlows[i] == None:
                        self.mFlows[i] = self.mFlow - sum(filter(None, self.mFlows.values()))
           # compute mass fractions 
            for u, i in zip(self.mFlows.keys(), self.mFlows.values()):
                self.mFrac[u] = i / self.mFlow
        # -----------------------------
        elif self.initType == "mx":

            # compute an individual mass fraction if necessary
            for i in self.names:
                if self.mFrac[i] == None:
                    self.mFrac[i] = 1- sum(filter(None, list(self.mFrac.values())))
                    break
                    
            # compute mass flow rates by species
            for i in self.names:
                self.mFlows[i] = self.mFrac[i] * self.mFlow
                
        # compute molar weights, then flows, then fractions
        if molar and (self.initType=="mx" or self.initType=="mm"):
            self.MW = {}
            for s in self.mFlows.keys():
                self.MW[s] = muc.per_tab.MW(s, pint=pint)
                self.nFlows[s] = self.mFlows[s] / self.MW[s]
            self.nFlow = sum(self.nFlows.values())
            self.nFrac = {}
            for s in self.mFlows.keys():
                self.nFrac[s] = self.nFlows[s] / self.nFlow           
        if changedType:
            self.initType = trueType

    #--------------------

    def from_molar(self):
        if self.initType == "nx":
            names = self.names
            for n in names:
                self.nFlows[n] = self.nFrac[n]*self.nFlow
        elif self.initType == "nn":
            names = self.names
            self.nFlow = sum(self.nFlows.values())
            for n in names:
                self.nFrac[n] = self.nFlows[n]/self.nFlow
        self.MW = {}
        for s in self.nFlows.keys():
            self.MW[s] = muc.per_tab.MW(s)
            self.mFlows[s] = self.nFlows[s] * self.MW[s]
        self.mFlow = None
        self.initType = "mm"
# -----------------------------------------
    # alternate constructor: for making alternate list
    @classmethod
    def SubMixture(cls, mixObj, ind):
        """
        Takes two argument: Mixture object, then tuple of indices to pull.
        Designed to take mFlows and calculate fractions and total, so make sure object you pass is filled.
        ex. mixture1, (0, 1, 2)
        """
        impKeys = list(mixObj.mFlows.keys())
        
        mFlows = {}
        for i in ind:
            mFlows[impKeys[i]] = mixObj.mFlows[impKeys[i]]
            
        newObj = cls(mFlows.keys(), mFlows.values(), kind = "mm")
        newObj.fill()
        return newObj
        
#         self.mFrac = {}
#         for i in self.mFlows.keys():
#             self.mFrac[i] = self.mFlows[i] / self.mFlow
         
             
    @staticmethod
    def solveMFlow(mix1, mix2, mixEnd, spec, known = "1"):
        """
        Takes two input mixes and one output mix, a species name, and a known mFlow "1" or "end".
        Returns nothing; modifies mFlow of the given mixtures.
        """
        if known == "1":
            mix2.mFlow = (mix1.mFlow*mixEnd.mFrac[spec] - mix1.mFlow*mix1.mFrac[spec])/(mix2.mFrac[spec] - mixEnd.mFrac[spec])
            mixEnd.mFlow = mix1.mFlow + mix2.mFlow
        elif known == "end":
            mix1.mFlow = (mixEnd.mFlow*mixEnd.mFrac[spec] - mixEnd.mFlow*mix2.mFrac[spec])/(mix1.mFrac[spec] - mix2.mFrac[spec])
            mix2.mFlow = mixEnd.mFlow - mix1.mFlow
                    
# ---------------------------------------------

    def Convert(self, rxn, spec, X):
        """
        Called from a mixture. Assumes that all mFlows are in moles.
        Arguments: rxn, a Reaction object; spec, a species identifier; X, conversion for the given species
        Returns a new mixture, having all the same species, after reaction with the given conversion and stoichiometric coefficients.
        """
        names = self.mFlows.keys()
        oldFlows = self.mFlows
        xi = -oldFlows[spec] * X / rxn.nu[spec]
        newFlows = []
        for s in names:
            newFlows.append(oldFlows[s] + xi * rxn.nu[s])
        newMix = Mixture(names, newFlows, kind = "mm")
        newMix.fill()
        return newMix

    def Bin_sep(self, x1, y1, spec, molar=False):
        """
        Takes two fractions and the name of the reference species.
        Assumes a binary mixture.
        Returns two new objects, which satisfy a mole balance on the reference species.
        Object order matches the fraction order passed
        (Algorithm derived from flash calculation, and therefore uses L and V notation.)
        """
        if len(self.names) != 2:
            raise ValueError("Not a binary mixture.")
        x2 = 1-x1
        y2 = 1-y1
        if molar:
            z1 = self.nFrac[spec]
            m1 = self.nFlow
        else:
            z1 = self.mFrac[spec]
            m1 = self.mFlow
        frac_l = (z1-y1)/(x1-y1)
        frac_v = (z1-x1)/(y1-x1)
        names = self.names
        ml = m1*frac_l
        mv = m1*frac_v
        
        if molar:
            l = Mixture(names, (x1, x2), ml, kind="nx")
            l.fill(self.molar_fill)
            v = Mixture(names, (y1, y2), mv, kind="nx")
            v.fill(self.molar_fill)
        else:
            l = Mixture(names, (x1, x2), ml, kind="mx")
            l.fill(self.molar_fill)
            v = Mixture(names, (y1, y2), mv, kind="mx")
            v.fill(self.molar_fill)
        return l, v

                
    def Extract(self, spec):
        """
        Takes an iterable of species names.
        Returns a new mixture, based on the mFlows, with only the species given as arguments.
        """
        newMix = {}
        for s in spec:
            newMix[s] = self.mFlows[s]
        newObj = Mixture(newMix.keys(), newMix.values(), kind = "mm")
        newObj.fill()
        return newObj
    
    def Remove(self, spec):
        """
        Takes an iterable of species names.
        Returns a new mixture, based on the mFlows, without the species given as arguments.
        """
        newMix = {}
        for s in self.mFlows.keys():
            if s not in spec:
                newMix[s] = self.mFlows[s]
        newObj = Mixture(newMix.keys(), newMix.values(), kind = "mm")
        newObj.fill()
        return newObj
    
    def __add__(self, other):
        if isinstance(other, Mixture):
            spec = list(self.mFlows.keys())
            for s in list(other.mFlows.keys()):
                if s not in spec:
                    spec.append(s)
            newFlows = [self.mFlows.get(s, 0) + other.mFlows.get(s, 0) for s in spec]
            # for s in spec:
                # newFlows[s] = self.mFlows[s] + other.mFlows[s]
            newMix = Mixture(spec, newFlows, sum(newFlows), kind = "mm")
            newMix.fill()
            return newMix
        else:
            print("You added a mixture to something else.")

    def __sub__(self, other):
        if isinstance(other, Mixture):
            spec = list(self.mFlows.keys())
            for s in list(other.mFlows.keys()):
                if s not in spec:
                    print("Warning: You subtracted a species that isn't there.")
            newFlows = [self.mFlows.get(s, 0) - other.mFlows.get(s, 0) for s in spec]
            # for s in spec:
                # newFlows[s] = self.mFlows[s] + other.mFlows[s]
            newMix = Mixture(spec, newFlows, sum(newFlows), kind = "mm")
            newMix.fill()
            return newMix
        else:
            print("You subtracted something else from a mixture.")
#---------------------------
    def set_Hf(self, Hf_list):
        for n, h in zip(self.names, Hf_list):
            self.s[n].set_Hf(h) 

    def set_Cp_const(self, Cp):
        """
        Takes list of single coefficients, corresponding to constant Cp values for each species.
        Must be in same order as when reaction was set.
        Wraps the Species function of the same name.
        """
        for n, c in zip(self.names, Cp):
            self.s[n].set_Cp_const(c)

    def set_Cp_func(self, coeff_list, coeff_kind):
        """
        Takes list of lists of coefficients, corresponding to Cp correlation for each species.
        Must be in same order as when reaction was set.
        Wraps Species function of the same name.
        kind 1: degree <=3 polynomial, coefficients in rising order. 
        kind 2: a + b*T + c*T**-2
        kind 3: a + b*T + c*T**2 + d*T**-2 (default)
        """
        names = self.names
        for n, c in zip(names, coeff_list):
            self.s[n].set_Cp_func(c, coeff_kind) 

    def calc_H(self, T, Tref, use_Hf, pint_strip = True):
        """
        Arguments: T, Tref, pint_strip
        T: mixture temperature, with proper units for type of Cp correlation
        Tref: reference temperature for heats of formation
        pint_strip: defaults to True. If True, strips pint-style units from temperatures before using Cp correlation.
        Returns total H of mixture, referenced to reference temperature of Hf values. Also stored as self.H.
        """
        self.H_by_spec = {}
        if use_Hf:
            for n in self.mFlows.keys():
                self.H_by_spec[n] = self.mFlows[n] * (self.s[n].Hf + self.s[n].calc_DH(Tref, T, pint_strip) )
        else:
            for n in self.mFlows.keys():
                self.H_by_spec[n] = self.mFlows[n] * (self.s[n].calc_DH(Tref, T, pint_strip) )
        self.H = sum(self.H_by_spec.values())						
        return self.H
# -------------------------------------------------------------
    def print(self, kind = None, dec = 3):
        """
        Prints out one type of the values stored. 
        Optional arguments: kind, dec, pint_strip
        kind: kind of values to print. 'nn', 'mm', 'nx', or 'mx'. Defaults to type initialized with mixture.
        dec: decimal points to print. Defaults to 3.
        """
        if kind == None:
            kind = self.initType
        if kind == "mx":
            for s in self.mFrac.keys():
                print(f"{s}: {self.mFrac[s]:.{dec}f}, ", end = "")
            print(f"total flow: {self.mFlow:.{dec}f}")
        elif kind == "nx":
            for s in self.mFrac.keys():
                print(f"{s}: {self.nFrac[s]:.{dec}f}, ", end = "")
            print(f"total flow: {self.nFlow:.{dec}f}")
        elif kind == "mm":
            for s in self.mFrac.keys():
                print(f"{s}: {self.mFlows[s]:.{dec}f}, ", end = "")
            print(f"total flow: {self.mFlow:.{dec}f}")
        elif kind == "nn":
            for s in self.mFrac.keys():
                print(f"{s}: {self.nFlows[s]:.{dec}f}, ", end = "")
            print(f"total flow: {self.nFlow:.{dec}f}")

class Reaction:
    def __init__(self, names, nus):
        """
        Takes list of stoich. coefficients and list of species names. Stores values.
        """
        self.names = names
        self.s = {}
        self.nu = {}
        for i, n in enumerate(names):
            self.s[n] = Species(n)
            self.nu[n] = nus[i]

    def set_Hf(self, Hf_list, Tref):
        for n, h in zip(self.names, Hf_list):
            self.s[n].set_Hf(h) 
        self.Tref = Tref

    def set_H0rxn(self, H0rxn, Tref):
        # Kludge. Set the heat of formation of one reactant equal to heat of reaction, all else to zero
        Hf_list = [0*H0rxn.units for n in self.names]
        Hf_list[0] = H0rxn/self.nu[self.names[0]]
        self.set_Hf(Hf_list, Tref)

        self.Tref = Tref

    def set_Cp_const(self, Cp):
        """
        Takes list of single coefficients, corresponding to constant Cp values for each species.
        Must be in same order as when reaction was set.
        Wraps the Species function of the same name.
        """
        for n, c in zip(self.names, Cp):
            self.s[n].set_Cp_const(c)

    def set_Cp_func(self, coeff_list, coeff_kind):
        """
        Takes list of lists of coefficients, corresponding to Cp correlation for each species.
        Must be in same order as when reaction was set.
        Wraps Species function of the same name.
        kind 1: degree <=3 polynomial, coefficients in rising order. 
        kind 2: a + b*T + c*T**-2
        kind 3: a + b*T + c*T**2 + d*T**-2 (default)
        """
        names = self.names
        for n, c in zip(names, coeff_list):
            self.s[n].set_Cp_func(c, coeff_kind) 

    def calc_H0rxn(self):
        self.H0rxn = 0
        for n in self.names:
            self.H0rxn += self.nu[n] * self.s[n].Hf
        return self.H0rxn

    def calc_Hrxn(self, T, Tref="Not Given", pint_strip = False):
        try:
            self.Hrxn = self.H0rxn
        except:
            self.Hrxn = self.calc_H0rxn()
        if Tref == "Not Given":
            Tref = self.Tref
        for n in self.names:
            self.Hrxn += self.nu[n] * self.s[n].calc_DH(Tref, T)
        return self.Hrxn


    

