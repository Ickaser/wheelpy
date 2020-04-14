from iapws import IAPWS97
import wheelpy.muc as muc
uReg = muc.uReg

class Steam:
    """
    Wrapper by Nathan Barrett for iapws.
    See https://pypi.org/project/iapws/ for documentation
    """
    def __init__(self,T="NotGiven",P="NotGiven",x="NotGiven",Print=False):
        """
        T    : Temperature
        P    : Pressure
        x    : Steam Quality
        H    : Enthalpy
        S    : Entropy
        l    : wavelength of light emmited
        V    : Specific volume
        """
        PossibleVars = [T,P,x]
        varsGiven = []
        for i in range(len(PossibleVars)):
            if PossibleVars[i] != "NotGiven":
                varsGiven.append(i)
        if len(varsGiven) != 2:
            raise Exception("Error! Steam Class definition can only take two arguments. " + str(len(varsGiven)) + " were provided.")
        
        if varsGiven == [0,1]:
            T = T.to(uReg.K).magnitude
            P = P.to(uReg.MPa).magnitude
            self.steam = IAPWS97(T=T,P=P)

        elif varsGiven == [0,2]:
            T = T.to(uReg.K).magnitude
            self.steam = IAPWS97(T=T,x=x)
        elif varsGiven == [1,2]:
            P = P.to(uReg.MPa).magnitude
            self.steam = IAPWS97(P=P,x=x)
            
        propertyFound = False
        
        try:
            self.H = (self.steam.h) * uReg.kJ / uReg.kg
            propertyFound = True
        except:
            self.H = "NotFound"
        try:
            self.S = (self.steam.s) * uReg.kJ / uReg.kg / uReg.K
            propertyFound = True
        except:
            self.S = "NotFound"
        try:
            self.P = (self.steam.P) * uReg.MPa
            propertyFound = True
        except:
            self.P = "NotFound"
        try:
            self.T = (self.steam.T) * uReg.K
            propertyFound = True
        except:
            self.T = "NotFound"
        try:
            self.x = self.steam.x
            propertyFound = True
            if Print:
                if self.x > 1:
                    print("This is Saturated Vapor.")
                elif self.x > 0:
                    print("This is in vapor liquid equilibirum.")
                elif self.x <= 0:
                    print("This is Saturated Liquid.")
        except:
            self.x = "NotFound"
            print("WARNING! Steam condition could not be determined.")
        try:
            self.l = (self.steam.l) * uReg.nm
            propertyFound = True
        except:
            self.l = "NotFound"
        try:
            self.V = (self.steam.v) * uReg.m**3 / uReg.kg
            propertyFound = True
        except:
            self.V = "NotFound"
            
        if not propertyFound:
            raise Exception("There is no data for the conditions given.")
