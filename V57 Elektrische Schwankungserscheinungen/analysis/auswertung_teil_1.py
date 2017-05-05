import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from funktionen import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

#TODO: Fehler auf die Spannungen berücksichtigen

# Einfachschaltung - Eichmessung.

f, VN, U = np.genfromtxt('data/einfachschaltung_kalibrationsmessung.txt', unpack='True')
f = 1000 * f
U = correctSelfNoise(U, VN)
U = normVoltage(U, VN)
print(U)
integral = np.trapz(1/max(durchlasskoeff(U))*durchlasskoeff(U), f)
print('Integral = {}'.format(integral)) #Fehler? Erstmal 2%

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

R, VN, U = np.genfromtxt('data/einfachschaltung_schwacher_widerstand.txt', unpack='True')

U = correctSelfNoise(U, VN)
U = normVoltage(U, VN)

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
Δν = ufloat(integral, integral/50)
k_w = params[0]*1/(4*T*Δν)

print('''
Frequenzband: Δν = {}

Schwacher Widerstand:
----------------------------
m   = {}+-{} [V^2/R]
b   = {}+-{} [V^2]
k_w = {}
'''.format(Δν, params[0], np.sqrt(cov[0][0]), params[1], np.sqrt(cov[1][1]), k_w))

# Rauschzahl
# -----------------------------------------------------------------------------
i = 9
F_korr = Rauschzahl(U[i], R[i], T, Δν)
print('''
Rauschzahl der Einfachschaltung für den Widerstand R={} Ohm und eine Temperatur von T={}K: {}
'''.format(R[i], T, F_korr))

# Einfachschaltung - starker Widerstand
# -----------------------------------------------------------------------------

R, VN, U = np.genfromtxt('data/einfachschaltung_starker_widerstand.txt', unpack='True')
U = normVoltage(U, VN)

params, cov = curve_fit(linearFunction, np.delete(R, [7, 20]), np.delete(U, [7, 20]))
x = np.linspace(min(R), max(R), 2)
plt.plot(R, U, 'rx', label='Messwerte')
plt.plot(x, linearFunction(x, *params), 'b-', label='Lineare Regression')
plt.xlabel(r"Widerstand"r'$\,R\,/\,\mathrm{\Omega}$')
plt.ylabel(r"Spannung"r'\,$U\,/\,\mathrm{V}$')
plt.grid()
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('einfachschaltung_starker_widerstand.pdf')
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
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Korrelatorschaltung - Eichmessung
# -----------------------------------------------------------------------------
f, VN, U = np.genfromtxt('data/korrelatorschaltung_kalibrationsmessung.txt', unpack='True')
f = 1000*f
U = normVoltage(U, VN, VS=10)

integral = np.trapz(1/max(durchlasskoeff(U))*durchlasskoeff(U), f)
Δν = ufloat(integral, integral/50)

fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(f/1000, 1/max(durchlasskoeff(U))*durchlasskoeff(U), 'bx')
plt.xlabel(r"Frequenz"r'$\,\nu\,/\,\mathrm{kHz}$')
plt.ylabel(r"Durchlasskoeffizient"r'\,$\beta$')
plt.ylim(0,1.1)
plt.grid()
plt.text(0.8, 0.9,'Integral={:3f}'.format(integral), horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
plt.tight_layout()
plt.savefig('durchlasskoeffizient_korr.pdf')
plt.cla()
plt.clf()

print('''
Korrelatorschaltung:
--------------------------------
Integral = {}
Δν = {} '''.format(integral, Δν))


# Korrelatorschaltung - schwacher Widerstand
# -----------------------------------------------------------------------------
R, VN, U = np.genfromtxt('data/korrelatorschaltung_schwacher_widerstand.txt', unpack='True')

U = normVoltage(U, VN, VS=10)

plt.plot(R, U, 'rx', label='Messwerte')
plt.xlabel(r"Widerstand"r'$\,R\,/\,\mathrm{\Omega}$')
plt.ylabel(r"Spannung"r'\,$U\,/\,\mathrm{V}$')
plt.grid()
params, cov = curve_fit(linearFunction, R, U)
x = np.linspace(min(R), max(R), 2)
plt.plot(x, linearFunction(x, *params), 'b-', label='Lineare Regression')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('korrelatorschaltung_schwacher_widerstand.pdf')
plt.cla()
plt.clf()

m_w = ufloat(params[0], np.sqrt(cov[0][0]))
k_w = params[0]*1/(4*T*Δν)

print('''
Frequenzband: Δν = {}

Schwacher Widerstand:
----------------------------
m   = {}+-{} [V^2/R]
b   = {}+-{} [V^2]
k_w = {}
'''.format(Δν, params[0], np.sqrt(cov[0][0]), params[1], np.sqrt(cov[1][1]), k_w))

# Rauschzahl
# -----------------------------------------------------------------------------
i = 10
F_korr = Rauschzahl(U[i], R[i], T, Δν)
print('''
Rauschzahl der Korrelatorschaltung für den Widerstand R={} Ohm und eine Temperatur von T={}K: {}
'''.format(R[i], T, F_korr))


# Korrelatorschaltung - starker Widerstand
# -----------------------------------------------------------------------------
R, VN, U = np.genfromtxt('data/korrelatorschaltung_starker_widerstand.txt', unpack='True')

U = normVoltage(U, VN, VS=10)

plt.plot(R, U, 'rx', label='Messwerte')
plt.xlabel(r"Widerstand"r'$\,R\,/\,\mathrm{\Omega}$')
plt.ylabel(r"Spannung"r'\,$U\,/\,\mathrm{V}$')
plt.grid()
params, cov = curve_fit(linearFunction, R, U)
x = np.linspace(min(R), max(R), 2)
plt.plot(x, linearFunction(x, *params), 'b-', label='Lineare Regression')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('korrelatorschaltung_starker_widerstand.pdf')
plt.cla()
plt.clf()

m_w = ufloat(params[0], np.sqrt(cov[0][0]))
k_w = params[0]/(4*T*Δν)

print('''
Frequenzband: Δν = {}

Starker Widerstand:
----------------------------
m   = {}+-{} [V^2/R]
b   = {}+-{} [V^2]
k_s = {}
'''.format(Δν, params[0], np.sqrt(cov[0][0]), params[1], np.sqrt(cov[1][1]), k_w))
