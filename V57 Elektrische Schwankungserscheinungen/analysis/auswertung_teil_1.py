import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from funktionen import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

#TODO: Fehler auf die Spannungen berücksichtigen

# Einfachschaltung - Eichmessung

f, VN, U = np.genfromtxt('einfachschaltung_kalibrationsmessung.txt', unpack='True')
f = 1000 * f
U = correctSelfNoise(U, VN)
U = normVoltage(U, 1, VN)

def durchlasskoeff(U, V_l = 10, U_sin = 150e-3):
    return U * np.sqrt(2)/(V_l * U_sin**2)

integral = np.trapz(1/max(durchlasskoeff(U))*durchlasskoeff(U), f)
print(integral) #Fehler? Erstmal 2%

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(f/1000, 1/max(durchlasskoeff(U))*durchlasskoeff(U), 'bx')
plt.xlabel(r"Frequenz"r'$\,\nu\,/\,\mathrm{kHz}$')
plt.ylabel(r"Durchlasskoeffizient"r'\,$\beta$')
plt.ylim(0,1.1)
plt.grid()
plt.text(0.8, 0.9,'Integral={:3f}'.format(integral), horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
plt.tight_layout()
plt.savefig('durchlasskoeffizient.pdf')
plt.cla()
plt.clf()

# Einfachschaltung - schwacher Widerstand

R, VN, U = np.genfromtxt('einfachschaltung_schwacher_widerstand.txt', unpack='True')

U = correctSelfNoise(U, VN)
U = normVoltage(U, 1000, VN)

plt.plot(R, U, 'rx', label='Messwerte')
plt.xlabel(r"Widerstand"r'$\,R\,/\,\mathrm{\Omega}$')
plt.ylabel(r"Spannung"r'\,$U\,/\,\mathrm{V}$')
plt.grid()
params, cov = curve_fit(linearFunction, R, U)
x = np.linspace(min(R), max(R), 2)
plt.plot(x, linearFunction(x, *params), 'b-', label='Lineare Regression')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('einfachschaltung_schwacher_widerstand.pdf')
plt.cla()
plt.clf()

m_w = ufloat(params[0], np.sqrt(cov[0][0]))
T = ufloat(293, 5) #K
Δν = ufloat(3533.4, 3533.4/50)
k_w = params[0]*1/(4*T*Δν)

print('''
Schwacher Widerstand:
----------------------------
m   = {}+-{} [V^2/R]
b   = {}+-{} [V^2]
k_w = {}
'''.format(params[0], np.sqrt(cov[0][0]), params[1], np.sqrt(cov[1][1]), k_w))


# Einfachschaltung - starker Widerstand

R, VN, U = np.genfromtxt('einfachschaltung_starker_widerstand.txt', unpack='True')
U = normVoltage(U, 1, VN)
plt.plot(R, U, 'rx')


params, cov = curve_fit(linearFunction, np.delete(R, [7, 20]), np.delete(U, [7, 20]))
x = np.linspace(min(R), max(R), 2)
plt.plot(x, linearFunction(x, *params), 'b-')
plt.cla()
plt.clf()

m_s = ufloat(params[0], np.sqrt(cov[0][0]))
k_s = params[0]*1/(4*T*Δν)

print('''
Starker Widerstand:
----------------------------
m  = {}+-{} [V^2/R]
b  = {}+-{} [V^2]
k_s = {}
'''.format(params[0], np.sqrt(cov[0][0]), params[1], np.sqrt(cov[1][1]), k_s))

# Korrelatorschaltung - Eichmessung
f, VN, U = np.genfromtxt('korrelatorschaltung_kalibrationsmessung.txt', unpack='True')
f = 1000*f
U = normVoltage(U, 1, VN)
plt.plot(f, U, 'rx')


params, cov = curve_fit(correlatorCalibrationCurve, f, durchlasskoeff(U), p0=[0.1, 1, 5000])
x = np.linspace(min(f), max(f), 100)
plt.plot(x, correlatorCalibrationCurve(x, *params), 'b-')
#plt.show()


# Korrelatorschaltung - schwacher Widerstand

R = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 995]
U = [0.06, 0.09, 0.068, 0.08992, 0.10912, 0.132, 0.1528, 0.179, 0.206, 0.225, 0.251, 0.270, 0.297, 0.313, 0.330, 0.364, 0.385, 0.408, 0.425, 0.450]
plt.plot(R, U, 'rx')
#plt.show()

# Korrelatorschaltung - starker Widerstand

R = [5070, 10130, 15120, 20000, 25000, 30000, 40100, 44900, 50300, 54300, 60300, 64900, 70200, 74900, 80300, 85100, 90100, 95200]
U = [2.172, 4.292, 6.376, 8.624, 11.536, 14.352, 21, 23.2, 27, 29.5, 32.9, 35, 37.9, 40.5, 44.1, 47.3, 53.2, 56.4]
plt.plot(R, U, 'rx')
#plt.show()
plt.cla()
plt.clf()
