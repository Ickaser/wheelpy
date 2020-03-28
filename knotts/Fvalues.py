
# coding: utf-8

""" 
## Available Functions
 * `P(tsi,tti,tto)`: Dimensionless ratio needed for calculation
 * `R(tsi,tso,tti,tto)`: Dimensionless ratio needed for calculation
 * `F1s2t(tsi,tso,tti,tto)`: Returns the F value for exchangers with 1 shell pass and 2, 4, or any other multiple of 2 tube passes.
 * `F2s4t(tsi,tso,tti,tto)`: Returns the F value for exchangers with 2 shell passes and 4, 8, or any other multiple of 4 tube passes.
 
The arguments to the functions are:
 * `tsi`: the temperature of the shell side inlet stream in K
 * `tso`: the temperature of the shell side outlet stream in K
 * `tti`: the temperature of the tube side inlet stream in K
 * `tto`: the temperature of the tube side outlet stream in K
 
## Citation
 These functions are coded from equations found in:
 
 Bowman, Mueller, and Nagle, *Mean Temperature Difference in Design*, **Transactions of the A.S.M.E**, Vol. 62, No. 4, 1940.
 
 Only functions for 1 and 2 shell pass exchangers are coded in this library.  Expressions for other types may be found in the original reference.
     
"""
# ## Load Libraries

import numpy as np

# ## Define the Dimensionless Ratios
def P(tsi,tti,tto):
    return (tto-tti)/(tsi-tti)

def R(tsi,tso,tti,tto):
    return (tsi-tso)/(tto-tti)

def F1s2t(tsi,tso,tti,tto):
    y=np.log10((2/P(tsi,tti,tto)-1-R(tsi,tso,tti,tto)+np.sqrt(R(tsi,tso,tti,tto)**2+1))/(2/P(tsi,tti,tto)-1-R(tsi,tso,tti,tto)-np.sqrt(R(tsi,tso,tti,tto)**2+1)))
    if R(tsi,tso,tti,tto) != 1:
        x=np.sqrt(R(tsi,tso,tti,tto)**2+1)/(R(tsi,tso,tti,tto)-1)*np.log10((1-P(tsi,tti,tto))/(1-P(tsi,tti,tto)*R(tsi,tso,tti,tto)))
    else:
        x=np.sqrt(R(tsi,tso,tti,tto)**2+1)*P(tsi,tti,tto)/2.3/(1-P(tsi,tti,tto))
    return x/y

def F2s4t(tsi,tso,tti,tto):
    y=np.log10((2/P(tsi,tti,tto)-1-R(tsi,tso,tti,tto)+2/P(tsi,tti,tto)*np.sqrt((1-P(tsi,tti,tto))*(1-P(tsi,tti,tto)*R(tsi,tso,tti,tto)))+np.sqrt(R(tsi,tso,tti,tto)**2+1))/(2/P(tsi,tti,tto)-1-R(tsi,tso,tti,tto)+2/P(tsi,tti,tto)*np.sqrt((1-P(tsi,tti,tto))*(1-P(tsi,tti,tto)*R(tsi,tso,tti,tto)))-np.sqrt(R(tsi,tso,tti,tto)**2+1)))
    if R(tsi,tso,tti,tto) != 1:
        x=np.sqrt(R(tsi,tso,tti,tto)**2+1)/2/(R(tsi,tso,tti,tto)-1)*np.log10((1-P(tsi,tti,tto))/(1-P(tsi,tti,tto)*R(tsi,tso,tti,tto)))
    else:
        x=np.sqrt(R(tsi,tso,tti,tto)**2+1)*P(tsi,tti,tto)/2/2.3/(1-P(tsi,tti,tto))
    return x/y


