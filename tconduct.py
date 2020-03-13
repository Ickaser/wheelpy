# Copyright (C) 2020 Thomas Allen Knotts IV - All Rights Reserved          #
# This file, tconduct.py, is a python library with functions that return   #
# the positive roots (eigenvalues) of the transcendental equations, and    #
# the corresponding coefficients, needed to evaluate the exact (series)    #
# solutions to 1D transient conduction in a plane wall, long cylinder,     #
# and sphere.  These equations are found in Chapter 5 of Berman et al.     #
# (see citation below)                                                     #
#                                                                          #
# tconduct.py is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU General Public License for more details.                             #
#                                                                          #
# All published work which utilzes the data obtained from this library     #
# have a copy of and cite the following.                                   #
# T. L. Bergman, A. S. Lavine, F. P. Incropera, D. P. Dewitt, Fundamentals #
# of Heat and Mass Transfer 7th ed., John Wiley & Sons, Hoboken, NJ,       #
# (2011).                                                                  #
#                                                                          #
# ======================================================================== #
# transientconduction.py                                                   #
#                                                                          #
# Thomas A. Knotts IV                                                      #
# Brigham Young University                                                 #
# Department of Chemical Engineering                                       #
# Provo, UT  84606                                                         #
# Email: thomas.knotts@byu.edu                                             #
# ======================================================================== #
# Version 1.0 - February 2020                                              #
# ======================================================================== #

# ======================================================================== #
# tconduct.py                                                              #
#                                                                          #
# The library can be loaded into python via the following command:         #
#                                                                          #
# import tconduct as tc                                                    #
#                                                                          #
# When imported in this way, the functions can be called as:               #
#                                                                          #
# tc.zwall(n,Bi)                                                           #
#                                                                          #
# which returns the nth zeta value for a plane wall at Biot number Bi.     #
# For thermal conductivty k and convective heat transfer coefficient h,    #
# the Biot number required in these functions is calculated as follows.    #
# Plane wall: Bi = hL/k where L is 0.5 the width of the wall               #
# Infinite cylinder: Bi = hr/k where r is the radius of the cylinder       #
# Sphere: Bi = hr/k where r is the radius of the sphere                    #
#                                                                          #
# A complete list of functions are found below.                            #
# Input Values                                                             #
#   n          an integer greater than zero specifying the nth term in     #
#              the series solution to the HDE                              #
#   Bi         the Biot number for the problem                             #
# Functions and Return Values                                              #
#   zwal       nth zeta value for Bi for a plane wall                      #
#   zcyl       nth zeta value for Bi for an infinite cylinder              #
#   zsph       nth zeta value for Bi for a sphere                          #
#   Cwal       nth coefficient for Bi for a plane wall                     #
#   Ccyl       nth coefficient for Bi for an infinite cylinder             #
#   Csph       nth coefficient for Bi for a sphere                         #
#                                                                          #
# Function     Explanation; Parameter Restrictions                         #
# ----------   ----------------------------------------------------------- #
# zwal(n,Bi)   zeta for plane wall; n=[1,2000], Bi=[0.0001,500]            #
# zcyl(n,Bi)   zeta for infinite cylinder; n=[1,2000]; Bi=[0.0001,194]     #
# zsph(n,Bi)   zeta for sphere; n=[1,2000], Bi=[0.0001,8500]               #
# Cwal(n,Bi)   coefficient for plane wall; same restrictions as zwall      #
# Ccyl(n,Bi)   coefficient for infinite cyliner; same restrictions as zcyl #
# Csph(n,Bi)   coefficient for sphere; same restrictions as zsph           #
# ======================================================================== #

import numpy          as np
import scipy.special  as bsl
from scipy.optimize   import ridder

# define the funtions that need to be solved
# these are placed in the f(x)=0 form

def fwal(x,Bi): # function for a plane wall
    return x*np.tan(x)-Bi
	
def fcyl(x,Bi): # function for an infinite cylinder
    return x*bsl.j1(x)/bsl.j0(x)-Bi	
	
def fsph(x,Bi): # function for a sphere
    return 1.0-x/np.tan(x)-Bi
	
# define the zeta functions for each geometry
def zwal(n,Bi):
    try:
        Bi.ito("")
    except:
        pass
    if n == 1:
        lb=0
        ub=(2*n-1)*np.pi/2
        return ridder(fwal,lb,ub,args=(Bi,))
    else:
        lb=(2*n-1)*np.pi/2-np.pi*0.9999     
        ub=(2*n-1)*np.pi/2*0.9999
        return ridder(fwal,lb,ub,args=(Bi,))

def zcyl(n,Bi):
    try:
        Bi.ito("")
    except:
        pass
    if n==1:
        lb=0
        ub=2.4
        return ridder(fcyl,lb,ub,args=(Bi,))
    else:
        lb=0.005+(n-1)*np.pi
        ub=2.35+(n-1)*np.pi
        return ridder(fcyl,lb,ub,args=(Bi,))
		
def zsph(n,Bi):
    try:
        Bi.ito("")
    except:
        pass
    lb=n*np.pi-0.9999*np.pi
    ub=n*np.pi*0.9999
    return ridder(fsph,lb,ub,args=(Bi,))
	
# define the C functions for each geometry
def Cwal(n,Bi):
    try:
        Bi.ito("")
    except:
        pass
    z=zwal(n,Bi)
    return(4*np.sin(z)/(2*z+np.sin(2*z)))

def Ccyl(n,Bi):
    try:
        Bi.ito("")
    except:
        pass
    z=zcyl(n,Bi)
    return(2*bsl.j1(z)/(z*(bsl.j0(z)**2+bsl.j1(z)**2)))

def Csph(n,Bi):
    try:
        Bi.ito("")
    except:
        pass
    z=zsph(n,Bi)
    return(4*(np.sin(z)-z*np.cos(z))/(2*z-np.sin(2*z)))
	
