# Copyright (C) 2018 Thomas Allen Knotts IV - All Rights Reserved          #
# This file, benzene.py, is a python library of the                        #
# thermophysical properties of benzene.  The properties, both              #
# constant and temperature-dependent, are taken from the DIPPR(R) Sample   #
# database which can be accessed at <https://dippr.byu.edu>.               #
# The vapor density at 1 atm comes from a spline of the data found in the  #
# 7th ed. of "The Fundamentals of Heat and Mass Transfer" by Bergman et.   #
# al. Density-related properties for the vapor phase use the DIPPR(R)      #
# and this spline function.                                                #
#                                                                          #
# benzeneproperties.py is distributed in the hope that it will be useful,  #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU General Public License for more details.                             #
#                                                                          #
# All published work which utilizes this library, or other property data   #
# from the DIPPR(R) database, should include the citation below.           #
# R. L. Rowley, W. V. Wilding, J. L. Oscarson, T. A. Knotts, N. F. Giles,  #
# DIPPR® Data Compilation of Pure Chemical Properties, Design Institute    #
# for Physical Properties, AIChE, New York, NY (2017).                     #
#                                                                          #
# All published work which utilzes the data obtained from the vdnsat       #
# function should also have a copy of and cite the following.              #
# T. L. Bergman, A. S. Lavine, F. P. Incropera, D. P. Dewitt, Fundamentals #
# of Heat and Mass Transfer 7th ed., John Wiley & Sons, Hoboken, NJ,       #
# (2011).                                                                  #
#                                                                          #
# ======================================================================== #
# benzeneproperties.py                                                     #
#                                                                          #
# Thomas A. Knotts IV                                                      #
# Brigham Young University                                                 #
# Department of Chemical Engineering                                       #
# Provo, UT  84606                                                         #
# Email: thomas.knotts@byu.edu                                             #
# ======================================================================== #
# Version 1.0 - February 2018                                              #
# ======================================================================== #

# ======================================================================== #
# benzeneproperties.py                                                     #
#                                                                          #
# This library contains functions for the properties of benzene.           #
#                                                                          #
# The library can be loaded into python via the following command:         #
# import benzeneproperties as benzene                                      #
#                                                                          #
# When imported in this way, the properties can be accessed as:            #
# benzene.tc for the critical temperature and benzene.vp(t) for the vapor  #
# pressure at temperature t where t is in units of K.                      #
# A complete list of properties, and the associated units, are found       #
# below.                                                                   #
#                                                                          #
# Function    Return Value                             Input Value         #
# ---------   --------------------------------------   -----------------   #
# tc          critical temperature in K                none                #
# pc          critical pressure in Pa                  none                #
# mw          molecular weight in kg/mol               none                #
# ldn(t)      liquid density in kg/m**3                temperature in K    #
# lcp(t)      liquid heat capacity in J/mol/K          temperature in K    #
# ltc(t)      liquid thermal conductivity in W/m/K     temperature in K    #
# vp(t)       liquid vapor pressure in Pa              temperature in K    #
# hvp(t)      heat of vaporization in J/mol            temperature in K    # 
# lpr(t)      liquid Prandtl number                    temperature in K    #
# lvs(t)      liquid viscosity in Pa*s                 temperature in K    #
# nu(t)       liquid kinematic viscosity in m**2/s     temperature in K    #
# tstat(p)    temperature at saturation in K           temperature in K    #
# vvs(t)      vapor viscosity in Pa*s                  temperature in K    #
# vtc(t)      vapor therm. conductiv. in W/m/K         temperature in K    #
# vPr(t)      vapor Prandtl number                     unitless            #
# ======================================================================== #

import numpy as np
from scipy.optimize import fsolve
from scipy   import interpolate

# critical temperature
tc = 562.05 # units of K

# critical pressure
pc = 4.895e6 # units of Pa

# molecular weight
mw = 0.07811184 # units of kg/mol
  
def ldn(t): # liquid density
    A = 1.0259
    B = 0.26666
    C = 562.05
    D = 0.28394
    E = 0
   
    y = A/(B**(1+(1-t/C)**D))
    y = y * 1000 # convert from kmol/m^3 to mol/m^3
    y = y * mw # convert from mol/m^3 to kg/m^3
    return y # units of kg/m^3
  
def lcp(t): # liquid heat capacity
    A = 162940
    B = -344.94
    C = 0.85562 
    D = 0
    E = 0
    y = A + B * t + C * t**2 + D * t**3 + E * t**4
    y = y / 1000 # convert from J/kmol/K to J/mol/K
    return y # units of J/mol/K

def ltc(t): # liquid thermal conductivity
    A = 0.23444
    B = -0.00030572
    C = 0
    D = 0
    E = 0
    y = A + B * t + C * t**2 + D * t**3 + E * t**4
    return y # units of W/m/K

def vp(t): # liquid vapor pressure
    A = 83.107
    B = -6486.2
    C = -9.2194
    D = 0.0000069844
    E = 2.0000E+00
    y = np.exp(A + B / t + C * np.log(t) + D * t**E)
    return y # units of Pa
    
def hvp(t): # heat of vaporization
    A = 50007000
    B = 0.65393
    C = -0.65393
    D = 0.029569
    tr = t/tc
    y = A * (1.0-tr)**(B + C * tr + D * tr**2)
    y = y / 1000 # convert from J/kmol to J/mol
    return y # J/mol
    
def lvs(t): # liquid viscosity
    A = 7.5117
    B = 294.68
    C = -2.794
    D = 0
    E = 0
    y = np.exp(A + B / t + C * np.log(t) + D * t**E)
    return y # units of Pa*s

def nu(t): # kinematic liquid viscosity
    return lvs(t)/ldn(t) # m**2/s

def lpr(t): # liquid Prandtl number
    return lcp(t)*lvs(t)/ltc(t)/mw # unitless

def ftsat(t,p): # function to calculate tsat with fsolve
    return vp(t) - p

def tsat(p): # saturation temperature (K) at pressure P (Pa)
    x = 700 # guess in K
    y = fsolve(ftsat,x,p)
    return(y[0]) # K
    
def vvs(t): # vapor viscosity
    A = 0.00000003134
    B = 0.9676
    C = 7.9
    return (A*t**B)/(1+C/t) # Pa*s

def vtc(t): # vapor thermal conductivity
    A = 0.00001652
    B = 1.3117
    C = 491
    return (A*t**B)/(1+C/t) # W/m/K   

def icp(t): # ideal gas heat capacity
    A = 33257.8886
    B = 51444.739266
    C = -761.088083
    D = 139737.490488
    E = 1616.907907
    F = 56829.10351
    G = 4111.398275
    y = A + B*(C/t)**2*np.exp(C/t)/(np.exp(C/t)-1)**2 + D*(E/t)**2*np.exp(E/t)/(np.exp(E/t)-1)**2 + F*(G/t)**2*np.exp(G/t)/(np.exp(G/t)-1)**2
    y = y / 1000 # convert from J/kmol/K to J/mol/K
    return y # units of J/mol/K

def vpr(t): # liquid Prandtl number
    return icp(t)*vvs(t)/vtc(t)/mw # unitless



  
