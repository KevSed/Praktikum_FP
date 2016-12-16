import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

#TODO: Fehler auf die Spannungen berücksichtigen

# Einfachschaltung - Eichmessung

f, VN, U = np.genfromtxt('eichmessung.txt', unpack='True')
f = 1000 * f
U = korrEigenrauschen(U, VN)
U = normVoltage(U, 1, VN)

def durchlasskoeff(U, V_l = 10, U_sin = 150e-3):
    return U * np.sqrt(2)/(V_l * U_sin**2)

integral = np.trapz(durchlasskoeff(U), f)
print(integral)

plt.plot(f, durchlasskoeff(U), 'rx')
#plt.show()
plt.cla()
plt.clf()

# Einfachschaltung - schwacher Widerstand

R, VN, U = np.genfromtxt('einfachschaltung_r_daten.txt', unpack='True')

plt.plot(R, U, 'rx')
#plt.show()

params, cov = curve_fit(linear, R, U)
x = np.linspace(min(R), max(R), 2)
plt.plot(x, linear(x, *params), 'b-')
#plt.show()
plt.cla()
plt.clf()
print('''
Schwacher Widerstand:
----------------------------
m = {}+-{} [V^2/R]
b = {}+-{} [V^2]
'''.format(params[0], np.sqrt(cov[0][0]), params[1], np.sqrt(cov[1][1])))

# Einfachschaltung - starker Widerstand

R, VN, U = np.genfromtxt('einfachschaltung_R_daten.txt', unpack='True')
U = normVoltage(U, 1, VN)

plt.plot(R, U, 'rx')
#plt.show()

params, cov = curve_fit(linear, np.delete(R, [7, 20]), np.delete(U, [7, 20]))
x = np.linspace(min(R), max(R), 2)
plt.plot(x, linear(x, *params), 'b-')
#plt.show()
plt.cla()
plt.clf()
print('''
Schwacher Widerstand:
----------------------------
m = {}+-{} [V^2/R]
b = {}+-{} [V^2]
'''.format(params[0], np.sqrt(cov[0][0]), params[1], np.sqrt(cov[1][1])))

# Korrelatorschaltung - Eichmessung
f, VN, U = np.genfromtxt('eichmessung_korrelator.txt', unpack='True')
f = 1000*f
U = normVoltage(U, 1, VN)
plt.plot(f, U, 'rx')

#plt.show()

params, cov = curve_fit(korrEich, f, durchlasskoeff(U), p0=[0.1, 1, 5000])
x = np.linspace(min(f), max(f), 100)
print(params)
plt.plot(x, korrEich(x, *params), 'b-')
plt.show()


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
