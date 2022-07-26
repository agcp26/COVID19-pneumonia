# -*- coding: utf-8 -*-
"""
This code plots the solutions of the coinfection model in the time interval
[0 tmax] using the initial conditions X0.
"""

from numpy import *
import pylab as pyl
from scipy import integrate

tmax=40000      # Maximal time
numt=8000       # Number of time nodes
t = linspace(0, tmax, numt)

# Initial conditions
X0 = array([8.33e7, 1e5, 1e3, 1e5, 1e3, 0, 0, 0, 0, 0, 0, 0, 0.8])

# We define fixed parameter values

Lamb = 2000     # People/day
mu = 2.4e-5     # 1/(people*day)
sigm = 1/100.   # 1/day
gamm = 1/12.    # 1/day
gamma1 = 1/20.  # 1/day
theta = 1/14.   # 1/day
theta1 = 1/24.  # 1/day
b1 = 2.0e-9     # 1/(people*day)
b2 = 0.1        # 1/(conc.*day)
alpha1 = 1.0e-8 # 1/(people*day)
delta = 0.001   # 1/day
delta0 = 0.005  # 1/day
delta1 = 0.01   # 1/day
delta2 = 0.2    # 1/day
eta = 0.12      # 1/day
eta1 = 0.1      # 1/day
phi = 1/14.     # 1/day
phi1 = 1/30.    # 1/day
phi2 = 1/40.    # 1/day
p = 1.0e-5      # (Concentration of bacteria)/(people*day)
kappa = 1.      # Concentration of bacteria
m = 0.01        # 1/day

k1 = gamm+eta+mu
k2 = mu+delta1+theta
k3 = mu+sigm
k4 = mu+delta+phi
k5 = gamma1+eta+eta1+mu+delta0+phi1
k6 = mu+delta1+delta2+phi2
k7 = mu+delta+phi+sigm

################################################
###        CASE 1: Rc>1, Rp<1, Rb<1:

alph = 3.0e-9   # 1/(people*day)
b = 1e-10       # 1/(people*day)
r = 0.004       # 1/day

RC = alph*Lamb/(mu*k1)
RP = b*Lamb/(mu*k4)
RB = r/m
print("\nCase 1:")
print("RC =",RC)
print("RP =",RP)
print("RB =",RB)

def dX_dt(X, t=0):
    """ Returns the growth rate of each subpopulation. """
    fSS = Lamb+sigm*X[3]-mu*X[0]-alph*X[0]*(X[1]+X[5]+X[9])-b*X[0]*(X[4]+X[5]+X[7]);
    fSI = alph*X[0]*(X[1]+X[5]+X[9])-k1*X[1]-b1*X[1]*(X[4]+X[5]+X[7]);
    fSH = -X[12]*X[2]*b2-X[2]*k2+X[1]*eta;
    fSR = gamm*X[1]+theta*X[2]-k3*X[3]-b*X[3]*(X[4]+X[5]+X[7]);
    fIS = sigm*X[7]+b*X[0]*(X[4]+X[5]+X[7])-alpha1*X[4]*(X[1]+X[5]+X[9])-k4*X[4];
    fII = b1*X[1]*(X[4]+X[5]+X[7])+alpha1*X[4]*(X[1]+X[5]+X[9])-k5*X[5];
    fIH = (eta+eta1)*X[5]+b2*X[2]*X[12]-theta1*X[6]-k6*X[6];
    fIR = b*X[3]*(X[4]+X[5]+X[7])+gamma1*X[5]+theta1*X[6]-k7*X[7];
    fRS = sigm*X[11]+phi*X[4]-mu*X[8]-alph*X[8]*(X[1]+X[5]+X[9]);
    fRI = phi1*X[5]+alph*X[8]*(X[1]+X[5]+X[9])-k1*X[9];
    fRH = X[6]*phi2-X[10]*k2+X[9]*eta;
    fRR = X[7]*phi+X[10]*theta+X[9]*gamm-X[11]*k3;
    fB = p*X[6]+r*X[12]*(1-X[12]/kappa)-m*X[12];
    return array([ fSS, fSI, fSH, fSR, fIS, fII,
                   fIH, fIR, fRS, fRI, fRH, fRR, fB ])

X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
infodict['message']                     # >>> 'Integration successful.'

solSS, solSI, solSH, solSR, solIS, solII, solIH, solIR, solRS, solRI, solRH, solRR, solB = X.T

f1 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solSS/10, 'b-', label='XSS/10')
pyl.plot(t, solSI, 'r--', label='XSI')
pyl.plot(t, solSH, 'm-.', label='XSH')
pyl.plot(t, solSR, 'g:', label='XSR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Susceptible to bacterial pneumonia')
pyl.xlim(0, 5000)
pyl.rc('font', size=8)
pyl.axes([.3, .25, .4, .33])
pyl.plot(t, solSI, 'r--', label='XSI')
pyl.plot(t, solSH, 'm-.', label='XSH')
pyl.plot(t, solSR, 'g:', label='XSR')
pyl.xlim(8000, tmax)
pyl.ylim(0, 85000)
pyl.grid()
f1.savefig("Rc-a.pdf", bbox_inches='tight')

f2 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solIS, 'b-', label='XIS')
pyl.plot(t, solII, 'r--', label='XII')
pyl.plot(t, solIH, 'm-.', label='XIH')
pyl.plot(t, solIR, 'g:', label='XIR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Infected by bacterial pneumonia')
pyl.xlim(0, 2000)
pyl.rc('font', size=8)
pyl.axes([.3, .25, .4, .55])
pyl.plot(t, solIS, 'b-', label='XIS')
pyl.plot(t, solII, 'r--', label='XII')
pyl.plot(t, solIH, 'm-.', label='XIH')
pyl.plot(t, solIR, 'g:', label='XIR')
pyl.xlim(4000, tmax)
pyl.ylim(0, 2000)
pyl.grid()
f2.savefig("Rc-b.pdf", bbox_inches='tight')

f3 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solRS/10, 'b-', label='XRS/10')
pyl.plot(t, solRI, 'r--', label='XRI')
pyl.plot(t, solRH, 'm-.', label='XRH')
pyl.plot(t, solRR, 'g:', label='XRR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Recovered from bacterial pneumonia')
pyl.xlim(0, 10000)
pyl.rc('font', size=8)
pyl.axes([.25, .5, .45, .35])
pyl.plot(t, solRI, 'r--', label='XRI')
pyl.plot(t, solRH, 'm-.', label='XRH')
pyl.plot(t, solRR, 'g:', label='XRR')
pyl.xlim(4000, tmax)
pyl.ylim(0, 20000)
pyl.grid()
f3.savefig("Rc-c.pdf", bbox_inches='tight')

f4 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solB  , 'm--', label='B')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Concentration of bacteria in environment')
pyl.xlim(0, 30000)
f4.savefig("Rc-d.pdf", bbox_inches='tight')

################################################
###        CASE 2: Rc>1, Rp>1, Rb<1:

alph = 3.0e-9   # 1/(people*day)
b = 9e-10       # 1/(people*day)
r = 0.004       # 1/day

RC = alph*Lamb/(mu*k1)
RP = b*Lamb/(mu*k4)
RB = r/m
print("\nCase 2:")
print("RC =",RC)
print("RP =",RP)
print("RB =",RB)

def dX_dt(X, t=0):
    """ Returns the growth rate of each subpopulation. """
    fSS = Lamb+sigm*X[3]-mu*X[0]-alph*X[0]*(X[1]+X[5]+X[9])-b*X[0]*(X[4]+X[5]+X[7]);
    fSI = alph*X[0]*(X[1]+X[5]+X[9])-k1*X[1]-b1*X[1]*(X[4]+X[5]+X[7]);
    fSH = -X[12]*X[2]*b2-X[2]*k2+X[1]*eta;
    fSR = gamm*X[1]+theta*X[2]-k3*X[3]-b*X[3]*(X[4]+X[5]+X[7]);
    fIS = sigm*X[7]+b*X[0]*(X[4]+X[5]+X[7])-alpha1*X[4]*(X[1]+X[5]+X[9])-k4*X[4];
    fII = b1*X[1]*(X[4]+X[5]+X[7])+alpha1*X[4]*(X[1]+X[5]+X[9])-k5*X[5];
    fIH = (eta+eta1)*X[5]+b2*X[2]*X[12]-theta1*X[6]-k6*X[6];
    fIR = b*X[3]*(X[4]+X[5]+X[7])+gamma1*X[5]+theta1*X[6]-k7*X[7];
    fRS = sigm*X[11]+phi*X[4]-mu*X[8]-alph*X[8]*(X[1]+X[5]+X[9]);
    fRI = phi1*X[5]+alph*X[8]*(X[1]+X[5]+X[9])-k1*X[9];
    fRH = X[6]*phi2-X[10]*k2+X[9]*eta;
    fRR = X[7]*phi+X[10]*theta+X[9]*gamm-X[11]*k3;
    fB = p*X[6]+r*X[12]*(1-X[12]/kappa)-m*X[12];
    return array([ fSS, fSI, fSH, fSR, fIS, fII,
                   fIH, fIR, fRS, fRI, fRH, fRR, fB ])

X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
infodict['message']                     # >>> 'Integration successful.'

solSS, solSI, solSH, solSR, solIS, solII, solIH, solIR, solRS, solRI, solRH, solRR, solB = X.T

f5 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solSS/10, 'b-', label='XSS/10')
pyl.plot(t, solSI, 'r--', label='XSI')
pyl.plot(t, solSH, 'm-.', label='XSH')
pyl.plot(t, solSR, 'g:', label='XSR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Susceptible to bacterial pneumonia')
pyl.xlim(0, 8000)
pyl.rc('font', size=8)
pyl.axes([.3, .25, .4, .33])
pyl.plot(t, solSI, 'r--', label='XSI')
pyl.plot(t, solSH, 'm-.', label='XSH')
pyl.plot(t, solSR, 'g:', label='XSR')
pyl.xlim(8000, tmax)
pyl.ylim(0, 75000)
pyl.grid()
f5.savefig("RcRp-a.pdf", bbox_inches='tight')

f6 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solIS/10, 'b-', label='XIS/10')
pyl.plot(t, solII, 'r--', label='XII')
pyl.plot(t, solIH, 'm-.', label='XIH')
pyl.plot(t, solIR, 'g:', label='XIR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Infected by bacterial pneumonia')
pyl.xlim(0, 2000)
pyl.rc('font', size=8)
pyl.axes([.3, .25, .4, .55])
pyl.plot(t, solIS/10, 'b-', label='XIS/10')
pyl.plot(t, solII, 'r--', label='XII')
pyl.plot(t, solIH, 'm-.', label='XIH')
pyl.plot(t, solIR, 'g:', label='XIR')
pyl.xlim(4000, 40000)
pyl.ylim(0, 2000)
pyl.grid()
f6.savefig("RcRp-b.pdf", bbox_inches='tight')

f7 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solRS/10, 'b-', label='XRS/10')
pyl.plot(t, solRI, 'r--', label='XRI')
pyl.plot(t, solRH, 'm-.', label='XRH')
pyl.plot(t, solRR, 'g:', label='XRR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Recovered from bacterial pneumonia')
pyl.xlim(0, 15000)
pyl.rc('font', size=8)
pyl.axes([.25, .5, .45, .35])
pyl.plot(t, solRI, 'r--', label='XRI')
pyl.plot(t, solRH, 'm-.', label='XRH')
pyl.plot(t, solRR, 'g:', label='XRR')
pyl.xlim(4000, tmax)
pyl.ylim(0, 35000)
pyl.grid()
f7.savefig("RcRp-c.pdf", bbox_inches='tight')

f8 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solB  , 'm--', label='B')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Concentration of bacteria in environment')
pyl.xlim(0, 30000)
f8.savefig("RcRp-d.pdf", bbox_inches='tight')

################################################
###        CASE 3: Rc>1, Rp<1, Rb>1:

alph = 3.0e-9   # 1/(people*day)
b = 1e-10       # 1/(people*day)
r = 0.08       # 1/day

RC = alph*Lamb/(mu*k1)
RP = b*Lamb/(mu*k4)
RB = r/m
print("\nCase 3:")
print("RC =",RC)
print("RP =",RP)
print("RB =",RB)

def dX_dt(X, t=0):
    """ Returns the growth rate of each subpopulation. """
    fSS = Lamb+sigm*X[3]-mu*X[0]-alph*X[0]*(X[1]+X[5]+X[9])-b*X[0]*(X[4]+X[5]+X[7]);
    fSI = alph*X[0]*(X[1]+X[5]+X[9])-k1*X[1]-b1*X[1]*(X[4]+X[5]+X[7]);
    fSH = -X[12]*X[2]*b2-X[2]*k2+X[1]*eta;
    fSR = gamm*X[1]+theta*X[2]-k3*X[3]-b*X[3]*(X[4]+X[5]+X[7]);
    fIS = sigm*X[7]+b*X[0]*(X[4]+X[5]+X[7])-alpha1*X[4]*(X[1]+X[5]+X[9])-k4*X[4];
    fII = b1*X[1]*(X[4]+X[5]+X[7])+alpha1*X[4]*(X[1]+X[5]+X[9])-k5*X[5];
    fIH = (eta+eta1)*X[5]+b2*X[2]*X[12]-theta1*X[6]-k6*X[6];
    fIR = b*X[3]*(X[4]+X[5]+X[7])+gamma1*X[5]+theta1*X[6]-k7*X[7];
    fRS = sigm*X[11]+phi*X[4]-mu*X[8]-alph*X[8]*(X[1]+X[5]+X[9]);
    fRI = phi1*X[5]+alph*X[8]*(X[1]+X[5]+X[9])-k1*X[9];
    fRH = X[6]*phi2-X[10]*k2+X[9]*eta;
    fRR = X[7]*phi+X[10]*theta+X[9]*gamm-X[11]*k3;
    fB = p*X[6]+r*X[12]*(1-X[12]/kappa)-m*X[12];
    return array([ fSS, fSI, fSH, fSR, fIS, fII,
                   fIH, fIR, fRS, fRI, fRH, fRR, fB ])

X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
infodict['message']                     # >>> 'Integration successful.'

solSS, solSI, solSH, solSR, solIS, solII, solIH, solIR, solRS, solRI, solRH, solRR, solB = X.T

f9 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solSS/10, 'b-', label='XSS/10')
pyl.plot(t, solSI, 'r--', label='XSI')
pyl.plot(t, solSH, 'm-.', label='XSH')
pyl.plot(t, solSR, 'g:', label='XSR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Susceptible to bacterial pneumonia')
pyl.xlim(0, 5000)
pyl.rc('font', size=8)
pyl.axes([.3, .25, .4, .33])
pyl.plot(t, solSI, 'r--', label='XSI')
pyl.plot(t, solSH, 'm-.', label='XSH')
pyl.plot(t, solSR, 'g:', label='XSR')
pyl.xlim(8000, tmax)
pyl.ylim(0, 85000)
pyl.grid()
f9.savefig("RcRb-a.pdf", bbox_inches='tight')

f10 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solIS, 'b-', label='XIS')
pyl.plot(t, solII, 'r--', label='XII')
pyl.plot(t, solIH, 'm-.', label='XIH')
pyl.plot(t, solIR, 'g:', label='XIR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Infected by bacterial pneumonia')
pyl.xlim(0, 2000)
pyl.rc('font', size=8)
pyl.axes([.3, .25, .4, .55])
pyl.plot(t, solIS, 'b-', label='XIS')
pyl.plot(t, solII, 'r--', label='XII')
pyl.plot(t, solIH, 'm-.', label='XIH')
pyl.plot(t, solIR, 'g:', label='XIR')
pyl.xlim(4000, tmax)
pyl.ylim(0, 2000)
pyl.grid()
f10.savefig("RcRb-b.pdf", bbox_inches='tight')

f11 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solRS/10, 'b-', label='XRS/10')
pyl.plot(t, solRI, 'r--', label='XRI')
pyl.plot(t, solRH, 'm-.', label='XRH')
pyl.plot(t, solRR, 'g:', label='XRR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Recovered from bacterial pneumonia')
pyl.xlim(0, 10000)
pyl.rc('font', size=8)
pyl.axes([.25, .5, .45, .35])
pyl.plot(t, solRI, 'r--', label='XRI')
pyl.plot(t, solRH, 'm-.', label='XRH')
pyl.plot(t, solRR, 'g:', label='XRR')
pyl.xlim(4000, tmax)
pyl.ylim(0, 20000)
pyl.grid()
f11.savefig("RcRb-c.pdf", bbox_inches='tight')

f12 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solB  , 'm--', label='B')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Concentration of bacteria in environment')
pyl.xlim(0, 10000)
f12.savefig("RcRb-d.pdf", bbox_inches='tight')

################################################
###        CASE 4: Rc>1, Rp>1, Rb>1:

alph = 3.0e-9   # 1/(people*day)
b = 9e-10       # 1/(people*day)
r = 0.08       # 1/day

RC = alph*Lamb/(mu*k1)
RP = b*Lamb/(mu*k4)
RB = r/m
print("\nCase 4:")
print("RC =",RC)
print("RP =",RP)
print("RB =",RB)

def dX_dt(X, t=0):
    """ Returns the growth rate of each subpopulation. """
    fSS = Lamb+sigm*X[3]-mu*X[0]-alph*X[0]*(X[1]+X[5]+X[9])-b*X[0]*(X[4]+X[5]+X[7]);
    fSI = alph*X[0]*(X[1]+X[5]+X[9])-k1*X[1]-b1*X[1]*(X[4]+X[5]+X[7]);
    fSH = -X[12]*X[2]*b2-X[2]*k2+X[1]*eta;
    fSR = gamm*X[1]+theta*X[2]-k3*X[3]-b*X[3]*(X[4]+X[5]+X[7]);
    fIS = sigm*X[7]+b*X[0]*(X[4]+X[5]+X[7])-alpha1*X[4]*(X[1]+X[5]+X[9])-k4*X[4];
    fII = b1*X[1]*(X[4]+X[5]+X[7])+alpha1*X[4]*(X[1]+X[5]+X[9])-k5*X[5];
    fIH = (eta+eta1)*X[5]+b2*X[2]*X[12]-theta1*X[6]-k6*X[6];
    fIR = b*X[3]*(X[4]+X[5]+X[7])+gamma1*X[5]+theta1*X[6]-k7*X[7];
    fRS = sigm*X[11]+phi*X[4]-mu*X[8]-alph*X[8]*(X[1]+X[5]+X[9]);
    fRI = phi1*X[5]+alph*X[8]*(X[1]+X[5]+X[9])-k1*X[9];
    fRH = X[6]*phi2-X[10]*k2+X[9]*eta;
    fRR = X[7]*phi+X[10]*theta+X[9]*gamm-X[11]*k3;
    fB = p*X[6]+r*X[12]*(1-X[12]/kappa)-m*X[12];
    return array([ fSS, fSI, fSH, fSR, fIS, fII,
                   fIH, fIR, fRS, fRI, fRH, fRR, fB ])

X, infodict = integrate.odeint(dX_dt, X0, t, full_output=True)
infodict['message']                     # >>> 'Integration successful.'

solSS, solSI, solSH, solSR, solIS, solII, solIH, solIR, solRS, solRI, solRH, solRR, solB = X.T

f13 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solSS/10, 'b-', label='XSS/10')
pyl.plot(t, solSI, 'r--', label='XSI')
pyl.plot(t, solSH, 'm-.', label='XSH')
pyl.plot(t, solSR, 'g:', label='XSR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Susceptible to bacterial pneumonia')
pyl.xlim(0, 5000)
pyl.rc('font', size=8)
pyl.axes([.3, .25, .4, .33])
pyl.plot(t, solSI, 'r--', label='XSI')
pyl.plot(t, solSH, 'm-.', label='XSH')
pyl.plot(t, solSR, 'g:', label='XSR')
pyl.xlim(8000, tmax)
pyl.ylim(0, 85000)
pyl.grid()
f13.savefig("RcRpRb-a.pdf", bbox_inches='tight')

f14 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solIS/10, 'b-', label='XIS/10')
pyl.plot(t, solII, 'r--', label='XII')
pyl.plot(t, solIH, 'm-.', label='XIH')
pyl.plot(t, solIR, 'g:', label='XIR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Infected by bacterial pneumonia')
pyl.xlim(0, 2000)
pyl.rc('font', size=8)
pyl.axes([.3, .25, .4, .55])
pyl.plot(t, solIS/10, 'b-', label='XIS/10')
pyl.plot(t, solII, 'r--', label='XII')
pyl.plot(t, solIH, 'm-.', label='XIH')
pyl.plot(t, solIR, 'g:', label='XIR')
pyl.xlim(4000, tmax)
pyl.ylim(0, 2000)
pyl.grid()
f14.savefig("RcRpRb-b.pdf", bbox_inches='tight')

f15 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solRS/10, 'b-', label='XRS/10')
pyl.plot(t, solRI, 'r--', label='XRI')
pyl.plot(t, solRH, 'm-.', label='XRH')
pyl.plot(t, solRR, 'g:', label='XRR')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Recovered from bacterial pneumonia')
pyl.xlim(0, 10000)
pyl.rc('font', size=8)
pyl.axes([.25, .5, .45, .35])
pyl.plot(t, solRI, 'r--', label='XRI')
pyl.plot(t, solRH, 'm-.', label='XRH')
pyl.plot(t, solRR, 'g:', label='XRR')
pyl.xlim(4000, tmax)
pyl.ylim(0, 35000)
pyl.grid()
f15.savefig("RcRpRb-c.pdf", bbox_inches='tight')

f16 = pyl.figure()
pyl.rc('font', size=10)
pyl.plot(t, solB  , 'm--', label='B')
pyl.grid()
pyl.legend(loc='best')
pyl.xlabel('time')
pyl.ylabel('population')
pyl.title('Concentration of bacteria in environment')
pyl.xlim(0, 10000)
f16.savefig("RcRpRb-d.pdf", bbox_inches='tight')