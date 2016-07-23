import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.optimize import curve_fit


r_1 = 50 # mm
r_2 = 70 # mm

L_0=(r_1+r_2)/2+np.sqrt(((r_1+r_2)/2)**2-r_1*r_2)
l_0=(r_1+r_2)/2-np.sqrt(((r_1+r_2)/2)**2-r_1*r_2)
L_1=r_1+r_2
print(
"""
--------------------------------------------------------------------------------------------
Minimale Länge für r_1={} und r_2={}:  L={} oder {}
Maximale Länge :                       L=0 oder {}
--------------------------------------------------------------------------------------------
""".format(r_1, r_2, L_0, l_0, L_1))

def f(x, r_1, r_2):
    return 1-x*(r_1+r_2)/(r_1*r_2)+x**2/(r_1*r_2)

x = np.linspace(0, 140, 1000)
plt.plot(x, f(x, r_1, r_2), 'r-', label='Verlauf von $g_1g_2$ für $r_1$={}cm, $r_2$={}cm'.format(r_1, r_2))
plt.xlim(x.min(), x.max())
plt.ylim(f(x, r_1, r_2).min()-0.4, f(x, r_1, r_2).max())
plt.axhline(y=1, xmin=0, xmax=140, linewidth=1, color='k')
plt.axhline(y=0, xmin=0, xmax=140, linewidth=1, color='k')
plt.ylabel(r"$g_1g_2$")
plt.xlabel(r"$\mathrm{Resonanzlänge}\; L/\mathrm{mm}$")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("plot_laser_{}_{}.pdf".format(r_1, r_2))
plt.close()

r_1 = 140 # cm
r_2 = 140 # cm

L_0=(r_1+r_2)/2+np.sqrt(((r_1+r_2)/2)**2-r_1*r_2)
l_0=(r_1+r_2)/2-np.sqrt(((r_1+r_2)/2)**2-r_1*r_2)
L_1=r_1+r_2
print(
"""
--------------------------------------------------------------------------------------------
Minimale Länge für r_1={} und r_2={}:  L={} oder {}
Maximale Länge :                       L=0 oder {}
--------------------------------------------------------------------------------------------
""".format(r_1, r_2, L_0, l_0, L_1))

def f(x, r_1, r_2):
    return 1-x*(r_1+r_2)/(r_1*r_2)+x**2/(r_1*r_2)


x = np.linspace(0, 280, 1000)
plt.plot(x, f(x, r_1, r_2), 'r-', label='Verlauf von $g_1g_2$ für $r_1$={}cm, $r_2$={}cm'.format(r_1, r_2))
plt.xlim(x.min(), x.max())
plt.ylim(f(x, r_1, r_2).min(), f(x, r_1, r_2).max()+0.2)
plt.axhline(y=1, xmin=0, xmax=200, linewidth=1, color='k')
plt.ylabel(r"$g_1g_2$")
plt.xlabel(r"$\mathrm{Resonanzlänge}\; L/\mathrm{mm}$")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("plot_laser_{}_{}.pdf".format(r_1, r_2))
plt.close()
#----------------------------------------------------------------------------

########################################################################################################
################################## P o l a r i s a t i o n #############################################
# Messung der Polarisation des Laserlichtes mit einem Polarisationsfilter. Geplottet wird die Stromstärke einer Diode gegen den Winkel des
# Filters in °.
w, I = np.loadtxt('Daten_polarisation.txt', unpack = True)
wrad = w*2*np.pi/360

def f(x, a, xo):
    return a*(np.cos(x+xo))**2

params, cov = curve_fit(f, wrad, I)
a  = params[0]
xo  = params[1]
a_err = np.sqrt(cov[0][0])
xo_err = np.sqrt(cov[1][1])

print(
"""
--------------------------------------------------------------------------------------------
Ausgleichsrechnung für die POLARISATION des Laserlichtes der Form a*np.cos(x+xo)**2
a = {}+-{} µA
xo = {}+-{} rad
--------------------------------------------------------------------------------------------
""".format(a, a_err, xo, xo_err))

t = np.linspace(wrad.min(), wrad.max(), 1000)
plt.plot(wrad, I, 'b.', label='Daten')
plt.plot(t, f(t, a, xo), 'r-', label='Ausgleichsrechnung')
plt.xlim(wrad.min()-0.2, wrad.max()+0.2)
plt.ylabel(r"Stromstärke $I/$µA")
plt.xlabel(r"Polarisationswinkel $/rad$")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("plot_polarisation.pdf")
plt.close()

########################################################################################################
################################## Stabilität konkav ###################################################

x, I = np.loadtxt('Daten_kk.txt', unpack = True)

def f(x, m, b):
    return m*x+b

params, cov = curve_fit(f, x, I)
m  = params[0]
b  = params[1]
m_err = np.sqrt(cov[0][0])
b_err = np.sqrt(cov[1][1])

t = np.linspace(x.min(), x.max(), 1000)
plt.plot(x, I, 'b.', label='Daten')
plt.plot(t, f(t, m, b), 'r-', label='Ausgleichsgeraden')
plt.xlim(x.min()-5, x.max()+5)
plt.ylabel(r"Stromstärke $I/$µA")
plt.xlabel(r"Resonatorlänge $/cm$")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("plot_kk.pdf")
plt.close()

########################################################################################################
################################## M O D E  0 0 ########################################################

# Messung der Moden des Laserlichtes im Resonator über einen dünnen Draht. Geplottet wird die Stromstärke einer Diode gegen die Position der Diode im Modenbild (willkürliche 0-Setzung bei Messbeginn).
x, I = np.loadtxt('Daten_Mode00.txt', unpack = True)

def f(x, a, b, c):
    return a*np.exp(-2*((x-b)**2)/(c**2))

params, cov = curve_fit(f, x, I)
a = params[0]
a_err = np.sqrt(cov[0][0])
b = params[1]
b_err = np.sqrt(cov[1][1])
c = params[2]
c_err = np.sqrt(cov[2][2])

print(
"""
--------------------------------------------------------------------------------
Parameter der Ausgleichsrechnung: (a*e**(-(x-b)**2/(2*c**2))) Mode 00
a = {}+-{}
b = {}+-{}
c = {}+-{}
--------------------------------------------------------------------------------
""".format(a, a_err, b, b_err, c, c_err))

t = np.linspace(x.min(), x.max(), 1000)
plt.plot(x, I, 'b.', label='Daten')
plt.plot(t, f(t, a, b, c), 'r-', label='Ausgleichsrechnung')
plt.xlim(x.min()-2, x.max()+2)
plt.ylabel(r"Stromstärke $I/$µA")
plt.xlabel(r"Abstand $/cm$")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("plot_Mode00.pdf")
plt.close()

########################################################################################################
################################## M O D E  1 0 ########################################################

x, I = np.loadtxt('Daten_Mode10.txt', unpack = True)
################### WHY U NO WORK?!?! #########################################
def f(x, I0, x0, w):
    return I0*(8*((x-x0)/w)**2)*np.exp(-2*((x-x0)/w)**2)

params, cov = curve_fit(f, x, I, p0=[33,20,13])
I0 = params[0]
I0_err = np.sqrt(cov[0][0])
x0 = params[1]
x0_err = np.sqrt(cov[1][1])
w = params[2]
w_err = np.sqrt(cov[2][2])

print(
"""
---------------------------------------------------------------------------------------------------
Parameter der Ausgleichsrechnung: (I0*((8*(x-x0)**2)/(w**2))*np.exp(-((x-x0)**2)/(w**2))) Mode 10
I0 = {}+-{}
x0 = {}+-{}
w = {}+-{}
---------------------------------------------------------------------------------------------------
""".format(I0, I0_err, x0, x0_err, w, w_err))

t = np.linspace(x.min(), x.max(), 1000)
plt.plot(t, f(t, I0, x0, w), 'r-', label='Ausgleichsrechnung')
plt.plot(x, I, 'b.', label='Daten')
plt.xlim(x.min()-2, x.max()+2)
plt.ylabel(r"Stromstärke $I/$µA")
plt.xlabel(r"Abstand $/cm$")
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("plot_Mode10.pdf")
plt.close()
