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
print('Minimale Länge (g_1g_2=0) für r_1={} und r_2={}, L={} oder {}'.format(r_1, r_2, L_0, l_0))
print('Maximale Länge (g_1g_2=1), L=0 oder {}'.format(L_1))

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
plt.tight_layout()
plt.savefig("plot_laser_{}_{}.pdf".format(r_1, r_2))
plt.close()

r_1 = 140 # cm
r_2 = 140 # cm

L_0=(r_1+r_2)/2+np.sqrt(((r_1+r_2)/2)**2-r_1*r_2)
l_0=(r_1+r_2)/2-np.sqrt(((r_1+r_2)/2)**2-r_1*r_2)
L_1=r_1+r_2
print('Minimale Länge (g_1g_2=0) für r_1={} und r_2={}, L={} oder {}'.format(r_1, r_2, L_0, l_0))
print('Maximale Länge (g_1g_2=1), L=0 oder {}'.format(L_1))

def f(x, r_1, r_2):
    return 1-x*(r_1+r_2)/(r_1*r_2)+x**2/(r_1*r_2)


x = np.linspace(0, 200, 1000)
plt.plot(x, f(x, r_1, r_2), 'r-', label='Verlauf von $g_1g_2$ für $r_1$={}cm, $r_2$={}cm'.format(r_1, r_2))
plt.xlim(x.min(), x.max())
plt.ylim(f(x, r_1, r_2).min(), f(x, r_1, r_2).max()+0.2)
plt.axhline(y=1, xmin=0, xmax=200, linewidth=1, color='k')
plt.ylabel(r"$g_1g_2$")
plt.xlabel(r"$\mathrm{Resonanzlänge}\; L/\mathrm{mm}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot_laser_{}_{}.pdf".format(r_1, r_2))
plt.close()
#----------------------------------------------------------------------------

# Messung der Polarisation des Laserlichtes mit einem Polarisationsfilter. Geplottet wird die Stromstärke einer Diode gegen den Winkel des Filters in °.
w, I = np.loadtxt('Daten_polarisation.txt', unpack = True)

plt.plot(w, I, 'b.', label='Daten')
plt.xlim(w.min()-5, w.max()+5)
plt.ylabel(r"Stromstärke $I/[\mu A]$")
plt.xlabel(r"Polarisationswinkel $/[^\circ]$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot_polarisation.pdf")
plt.close()

# Messung der Polarisation des Laserlichtes mit einem Polarisationsfilter. Geplottet wird die Stromstärke einer Diode gegen den Winkel des Filters in °.
x, I = np.loadtxt('Daten_kk.txt', unpack = True)

def f(x, m, b):
    return m*x+b

params, cov = curve_fit(f, x, I)
m  = params[0]
b  = params[1]
m_err = np.sqrt(cov[0][0])
b_err = np.sqrt(cov[1][1])

plt.plot(x, I, 'rx', label='Daten')
#plt.plot(x, f(x, m, b), 'rx', label='Ausgleichsgeraden')
plt.xlim(x.min()-5, x.max()+5)
plt.ylabel(r"Stromstärke $I/[\mu A]$")
plt.xlabel(r"Resonatorlänge $/[cm]$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot_kk.pdf")
plt.close()

# Messung der Moden des Laserlichtes im Resonator über einen dünnen Draht. Geplottet wird die Stromstärke einer Diode gegen die Position der Diode im Modenbild (willkürliche 0-Setzung bei Messbeginn).
x, I = np.loadtxt('Daten_Mode00.txt', unpack = True)

def f(x, a, b, c):
    return a*e**(-(x-b)**2/(2*c**2))


plt.plot(x, I, 'b.', label='Daten')
plt.xlim(x.min()-2, x.max()+2)
plt.ylabel(r"Stromstärke $I/[\mu A]$")
plt.xlabel(r"Abstand $/[mm]$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot_Mode00.pdf")
plt.close()

x, I = np.loadtxt('Daten_Mode10.txt', unpack = True)

plt.plot(x, I, 'b.', label='Daten')
plt.xlim(x.min()-2, x.max()+2)
plt.ylabel(r"Stromstärke $I/[\mu A]$")
plt.xlabel(r"Abstand $/[mm]$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot_Mode10.pdf")
plt.close()
