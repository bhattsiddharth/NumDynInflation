from math import *
import numpy as np
import scipy.integrate as sci
import scipy.misc as scm
import matplotlib.pyplot as plt


#data = np.loadtxt("data/Starobinsky_MS_scales.txt") # Ne,x,y,z,A,aH 6


def i(Ne):
    return np.max(np.where(data[:,0]>=Ne))

S = 5e-5 # for scaling time


# parameters used in the potential function
b = np.power(2.0/3.0, 0.5) # \sqrt(2/3)
v0 = 0.96*np.power(10.0, -10)

def f(x):
    return np.power(1 - np.exp(-b*x),2)

def dfdx(x):
    return 2*b*np.exp(-b*x)*(1 - np.exp(-b*x))

def d2fdx2(x):
    return 2*b**2*np.exp(-2*b*x)*(2 - np.exp(b*x))

for Nk in [60]: #np.arange(0,72,2):
    kp = 1#data[i(60),5] #388.711
    k = 1#data[i(Nk),5] #53.493

    # Main system of equations to be solved
    def sys(var, T):
        [x, y, z, A, v, vv, u, vu, h, vh, g, vg] = var
        
        dxdT = y
        dydT = -3*z*y - v0*dfdx(x)/S**2 
        dzdT = -0.5*y**2 #-z**2 + (v0*f(x)/S**2 - y**2)/3 # 
        dAdT = A*z
        dvdT = vv
        dvvdT = -z*vv + v*(2.5*y**2 + 2*y*(-3*z*y - v0*dfdx(x)/S**2 )/z + 2*z**2 + 0.5*y**4/z**2 - v0*d2fdx2(x)/S**2 - k**2/A**2)
        dudT = vu
        dvudT = -z*vu + u*(2.5*y**2 + 2*y*(-3*z*y - v0*dfdx(x)/S**2 )/z + 2*z**2 + 0.5*y**4/z**2 - v0*d2fdx2(x)/S**2 - k**2/A**2)
        dhdT = vh
        dvhdT = -z*vh - h*(k**2/A**2 - 2*z**2 + 0.5*y**2)
        dgdT = vg
        dvgdT = -z*vg - g*(k**2/A**2 - 2*z**2 + 0.5*y**2)

        return [dxdT, dydT, dzdT, dAdT, dvdT, dvvdT, dudT, dvudT, dhdT, dvhdT, dgdT, dvgdT]


    T = np.linspace(0, 1000, 10000)

    xi = 5.75#data[i(Nk+5),1]
    yi = 0#data[i(Nk+5),2]
    zi = np.sqrt(yi**2/6 + (v0*f(xi)/(3*S**2)))
    Ai = 1e-3 #* np.exp(79.4135 - (Nk+5))

    vi = (1/np.sqrt(2*k))
    ui = 0
    vvi = 0
    vui = -k*(1/np.sqrt(2*k))/Ai

    hi = (1/np.sqrt(2*k))
    gi = 0
    vhi = 0
    vgi = -k*(1/np.sqrt(2*k))/Ai

    sol = sci.odeint(sys, [xi,yi,zi,Ai,vi,vvi,ui,vui,hi,vhi,gi,vgi], T, rtol=3e-14, atol=2e-35, mxstep=900000000, full_output=True)
    x, y, z, A, v, vv, u, vu, h, vh, g, vg = np.transpose(sol)

    Z = S*z 
    Y = S*y
    X = x

    N = np.log(A/Ai)
    Nt = 79.4135
    Ne = Nt - N 

    aH = A*z
    aHk = aH/k
    meff = 2.5*y**2 + 2*y*(-3*z*y - v0*dfdx(x)/S**2 )/z + 2*z**2 + 0.5*y**4/z**2 - v0*d2fdx2(x)/S**2

    # slow-roll parameters
    epsH = -(-z**2 + (v0*f(x)/S**2 - y**2)/3)/z**2
    etaH = -(-3*z*y - v0*dfdx(x)/S**2)/(y*z)

    # observables
    ns = 1 + 2*etaH - 4*epsH
    r = 16*epsH
    ##As = (24*np.pi**2)**(-1) * v0*f(x)/epsH
    Ps = (S*z)**2 / (8 * np.pi**2 * epsH)
    ##At = (1.5*np.pi**2)**(-1) * v0*f(x)
    Pt = 2*(S*z)**2 / (np.pi**2)

    zeta2 = (v**2 + u**2)/(2*epsH*(A/S)**2)
    P_S = (k**3 * zeta2)/(2*np.pi**2)
    h2 = (h**2 + g**2)/((A/S)**2)
    P_T = 4*(k**3 * h2)/(np.pi**2)

    #plt.plot(T, N, 'r')
    #plt.plot(N,epsH, 'r')
    #plt.plot(Ne,Ps,'r')
    #plt.yscale('log')
    #plt.show()

    np.savetxt('data/Starobinsky_MS_background.txt',np.c_[T,N,Ne,x,y,z,aH,epsH,etaH,meff,Ps,Pt])
    print('\n\t--- Data saved successfully: Starobinsky_MS_background.txt ---\n')

    #np.savetxt('data/Starobinsky_MS_amplitude.txt',np.c_[Ne,aHk,zeta2,P_S,h2,P_T])
    #print('\n\t--- Data saved successfully: Starobinsky_MS_amplitude.txt ---\n')

    #file_name = 'data/modes/Starobinsky_MS_amplitude%.1f.txt'%Nk
    #np.savetxt(file_name,np.c_[Ne,aHk,P_S,P_T])
    #print('\t--- Data saved successfully: %s ---\n'%file_name)

    
