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

# It is required to execute the script 'inf_dyn_background.py' and save the data in a text file before this script can be executed
# We input the initial conditions and horizon exit for various scales from the background data
# Please change the filename in this line if you have saved the data with a different name or at a different location
data = np.loadtxt('data/inf_bg_data.txt') # row: T,N,Ne,x,y,z,aH,epsH,etaH,meff,Ps,Pt. 12 columns

# This function returns the index of the row in the data file where the number of e-folds before the 
#   end of inflation (Ne) attains a certain value specified by the argument 
def i(Ne):
    return np.max(np.where(data[:,2]>=Ne))



#########################################################################################################
# The model of inflation is defined in this section
#########################################################################################################

# This term defines one unit of time 
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
# In this section we set the initial conditions for both background as well as fluctuations
# The background data file is used to input initial conditions for background quantities
# We solve the dynamical equations, including the dimensionless Mukhanov-Sasaki equation using the
#   function scipy.integrate.odeint 
#########################################################################################################


### The dynamical variables are defined as follows:
#
# background:
#
# x : dimensionless field value [ \phi / m_p ]
# y : dimensionless field velocity [ dx/dT or \dot\phi / (m_p ^2 * S) ]
# A : dimensionless scale factor [ a * m_p * S ]
# z : dimensionless hubble parameter [ H / (S * m_p) ]
#
#
# fluctuations:
#
# v : real part of the Mukhanov-Sasaki variable [ v_k ] for scalar fluctuations [ \zeta_k ]
# u : imaginary part of the Mukhanov-Sasaki variable for scalar fluctuations
# h : real part of the Mukhanov-Sasaki variable for tensor fluctuations [ h_k ]
# g : imaginary part of the Mukhanov-Sasaki variable for tensor fluctuations




### Select a mode that leaves the horizon a certain number of e-folds before the end of inflation 
# Note that each mode corresponds to one point in the final power spectrum
Nk = 60 

k = data[i(Nk),6] # frequency of the mode that exits horizon 'Nk' e-folds before the end of inflation
kp = data[i(60),6] # frequency of the mode corresponding to the CMB scale i.e. 60 e-folds before the end of inflation


### initial conditions for the given system

# For initial conditions of background variables, we need to select the values of x,y and A when the mode
#   in question was sub-Hubble. A safe guess for this is 5 e-folds before horizon exit
xi = data[i(Nk+5),3]
yi = data[i(Nk+5),4]
zi = np.sqrt(yi**2/6 + (v0*f(xi)/(3*S**2)))
Ai = 1e-3 * np.exp(77.4859  - (Nk+5))

# Initial conditions for the fluctuations (in Mukhanov-Sasaki variables) are given by the Bunch-Davies vacuum
vi = (1/np.sqrt(2*k)) # real part of the field that describes scalar fluctuations
ui = 0 # imaginary part of the field that describes scalar fluctuations
v_Ti = 0
u_Ti = -k*(1/np.sqrt(2*k))/Ai
# Same applies for tensor fluctuations too
hi = (1/np.sqrt(2*k)) # real part of the field that describes tensor fluctuations
gi = 0 # imaginary part of the field that describes tensor fluctuations
h_Ti = 0
g_Ti = -k*(1/np.sqrt(2*k))/Ai



### the system of differential equations to be solved
def sys(var, T):
    [x, y, z, A, v, v_T, u, u_T, h, h_T, g, g_T] = var
    
    #background
    dxdT = y
    dydT = -3*z*y - v0*dfdx(x)/S**2 
    dzdT = -0.5*y**2 #-z**2 + (v0*f(x)/S**2 - y**2)/3 # 
    dAdT = A*z

    # scalar fluctuations
    dvdT = v_T
    dv_TdT = -z*v_T + v*(2.5*y**2 + 2*y*(-3*z*y - v0*dfdx(x)/S**2 )/z + 2*z**2 + 0.5*y**4/z**2 - v0*d2fdx2(x)/S**2 - k**2/A**2)
    dudT = u_T
    du_TdT = -z*u_T + u*(2.5*y**2 + 2*y*(-3*z*y - v0*dfdx(x)/S**2 )/z + 2*z**2 + 0.5*y**4/z**2 - v0*d2fdx2(x)/S**2 - k**2/A**2)
    
    # tensor fluctuations
    dhdT = h_T
    dh_TdT = -z*h_T - h*(k**2/A**2 - 2*z**2 + 0.5*y**2)
    dgdT = g_T
    dg_TdT = -z*g_T - g*(k**2/A**2 - 2*z**2 + 0.5*y**2)

    return [dxdT, dydT, dzdT, dAdT, dvdT, dv_TdT, dudT, du_TdT, dhdT, dh_TdT, dgdT, dg_TdT]


# This term defines one unit of time 
S = 5e-5

# the period of time over which the system is integrated 
T = np.linspace(0, 200, 10000)


sol = odeint(sys, [xi,yi,zi,Ai,vi,v_Ti,ui,u_Ti,hi,h_Ti,gi,g_Ti], T, rtol=3e-14, atol=2e-35, mxstep=900000000)
x, y, z, A, v, v_T, u, u_T, h, h_T, g, g_T = np.transpose(sol)

N = np.log(A/Ai)
Nt = 77.4859 
Ne = Nt - N 

aH = A*z
aHk = aH/k
meff = 2.5*y**2 + 2*y*(-3*z*y - v0*dfdx(x)/S**2 )/z + 2*z**2 + 0.5*y**4/z**2 - v0*d2fdx2(x)/S**2

# slow-roll parameters
epsH = -(-z**2 + ((v0*f(x)/S**2 - y**2))/3)/z**2
etaH = -(-3*z*y - v0*dfdx(x)/S**2)/(y*z)

# observable quantities (under slow-roll apparoximation)
ns = 1 + 2*etaH - 4*epsH
r = 16*epsH
Ps = (S*z)**2 / (8 * np.pi**2 * epsH)
Pt = 2*(S*z)**2 / (np.pi**2)

# values of power spectra
zeta2 = (v**2 + u**2)/(2*epsH*(A/S)**2) # scalar fluctuations
P_S = (k**3 * zeta2)/(2*np.pi**2)
h2 = (h**2 + g**2)/((A/S)**2) # tensor fluctuations
P_T = 4*(k**3 * h2)/(np.pi**2)



#########################################################################################################
# This section is for plotting the behaviour of various quantities
# Uncomment each subsection for the corresponding plot
#########################################################################################################

''' # horizon exit behaviour of the power spectrum
plt.plot(aHk,P_S,'r', lw=2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('aH/k')
plt.ylabel('P_S')
plt.show()
'''



#########################################################################################################
# This section is for saving the data in a text file
# Uncomment this after the values of all parameters have been fixed
#
# For saving data of different modes, only the last part of the file name can be changed.
# Single decimal place precision is adequate for Ne of horizon exit even for highly scale-dependent
#   power spectra
#########################################################################################################

np.savetxt('data/modes/inf_MS_data_%.1f.txt'%Nk,np.c_[Ne,aHk,zeta2,P_S,h2,P_T])
print('\n\t--- Data saved successfully : inf_MS_data_%.1f.txt ---\n'%Nk)

#########################################################################################################
#########################################################################################################
