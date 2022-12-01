#########################################################################################################
#########################################################################################################
#
# Please refer to <arXiv link> for explaination of variables and instructions for using the code
#
#########################################################################################################
#########################################################################################################

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt



#########################################################################################################
# The model of inflation is defined in this section
#########################################################################################################

# This term defines one unit of time [ T = t * m_p * S ] where t is the actual cosmic time
S = 5e-5 


# parameters used in the potential function
M = 5.9e-6 
v0 = 0.5*M**2


# dimensionless potential function and its derivatives
def f(x):
    return x**2

def dfdx(x):
    return 2*x

def d2fdx2(x):
    return 2



#########################################################################################################
# In this section we set the initial conditions for inflation and define our dynamical equations 
# We solve the equations using the function scipy.integrate.odeint 
# After getting a numerical solution for the dynamical quantities, we define the derived quantities
#########################################################################################################


### The dynamical variables are defined as follows:
#
# x : dimensionless field value [ \phi / m_p ]
# y : dimensionless field velocity [ dx/dT or \dot\phi / (m_p ^2 * S) ]
# A : dimensionless scale factor [ a * m_p * S ]
# z : dimensionless hubble parameter [ H / (S * m_p) ]


### initial conditions for the given system

xi = 17.5 # initial field value
yi = 0 # initial field velocity
zi = np.sqrt(yi**2/6 + (v0*f(xi)/(3*S**2))) # initial value of hubble parameter
Ai = 1e-3 # initial value of scale factor


# the system of differential equations to be solved
def sys(var, T):
    [x, y, z, A] = var

    # Note that all derivatives are taken wrt the scaled, dimenstionless cosmic time T

    dxdT = y
    dydT = -3*z*y - v0*dfdx(x)/S**2 
    dzdT = -0.5*y**2 # = -z**2 + (v0*f(x)/S**2 - y**2)/3 
    dAdT = A*z
    
    return [dxdT, dydT, dzdT, dAdT]


# the period of time over which the system is integrated 
# note: time is in units of S*m_p 
T = np.linspace(0, 1000, 100000)


# invoking the ODE solver
sol = odeint(sys, [xi,yi,zi,Ai], T, rtol=3e-14, atol=2e-35, mxstep=900000000)
x, y, z, A = np.transpose(sol)
phi, phi_t, H = x, y*S, z*S


### for consistency check
# We compute the evolution of z in two ways: one by solving the differential equation for z (as done already) and the other by computing x and y and then using their relation to compute z (denoted by z1, as given below)
#   If our analysis is correct, both the solutions should match exactly. If they don't match then it's an indication that there's a mistake somewhere in the script
#   This mistake could be due to an incorrect expression in the potential function, system of equations, other expression
z1 = np.sqrt(y**2 /6 + (v0*f(x)/(3*S**2))) 


### Since the expansion of the universe is expected to be near-exponential, it is convenient to express the evolution of scale factor in terms of the power of e by which it increases (e-folds)
# The number of e-folds can be used as a time-axis for studying the behaviour of various quantities
N = np.log(A/Ai) # number of e-folds of expansion elapsed
Nt = 77.4859 # number of e-folds elapsed when inflation ends (value needs to be fixed from the behaviour of epsH)
Ne = Nt - N # number of e-folds of expansion remaining before the end of inflation


### We can define other quantities in terms of our main dynamical quantites that provide us useful information about the system

## slow-roll parameters
epsH = -(-z**2 + ((v0*f(x)/S**2 - y**2))/3)/z**2 # indicated whether the system is in a state of inflation
etaH = -(-3*z*y - v0*dfdx(x)/S**2)/(y*z) # gives the rate of change of epsH

## observable quantities (under slow-roll apparoximation)
ns = 1 + 2*etaH - 4*epsH # 
r = 16*epsH
Ps = (S*z)**2 / (8 * np.pi**2 * epsH)
Pt = 2*(S*z)**2 / (np.pi**2)

## These quantities will be useful for studying the dynamics of quantum fluctuations
aH = A*z # rate of change of scale factor
meff = 2.5*y**2 + 2*y*(-3*z*y - v0*dfdx(x)/S**2 )/z + 2*z**2 + 0.5*y**4/z**2 - v0*d2fdx2(x)/S**2 # effective mass of fluctuations (time-dependent)



#########################################################################################################
# This section is for plotting the behaviour of various quantities
# Uncomment each subsection for the corresponding plot
#########################################################################################################

''' # to check if inflation is starting and ending
plt.plot(T,N,'r', lw=2)
plt.xlabel('T')
plt.ylabel('N')
plt.show()
'''
''' # consistency check
plt.plot(N,z, 'r', lw=3, label='from system')
plt.plot(N,z1, 'k--', lw=3, label='from conservation')
plt.xlabel('N (e-folds)')
plt.ylabel('z')
plt.legend()
plt.show()
'''
''' # phase behaviour
plt.plot(phi, phi_t, 'r', lw=2)
plt.show()
'''
#''' # to check if inflation is adequate
plt.plot(N,epsH, 'r')
plt.axhline(1, color='black')
plt.xlabel('N (e-folds)')
plt.ylabel('eps_H')
plt.show()
#'''
''' # to adjust the value of V_0
plt.plot(Ne,Ps, 'r')
plt.axvline(60, color='black')
plt.yscale('log')
plt.xlim(0,70)
#plt.ylim(2e-9,2.4e-9)
plt.xlabel('Ne')
plt.ylabel('P_s')
plt.show()
'''


#########################################################################################################
# This section is for saving the data in a text file
# Uncomment this after the values of all parameters have been fixed
#########################################################################################################

#np.savetxt('data/inf_bg_data.txt',np.c_[T,N,Ne,x,y,z,aH,epsH,etaH,meff,Ps,Pt])
#print('\n\t--- Data saved successfully : inf_bg_data.txt ---\n')

#########################################################################################################
#########################################################################################################
