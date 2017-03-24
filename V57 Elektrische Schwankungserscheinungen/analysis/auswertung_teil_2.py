import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties import unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from funktionen import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

################################################################################
### Oxydkathode
################################################################################

ν, VS, VN, U, ΔU = np.genfromtxt('oxydkathode.txt', unpack = True)
Δν = frequency2bandwidth(ν, 'Oxydkathode')
U = unp.uarray(U, ΔU)
U = normVoltage(U, VN, VS=VS)
R = 2200
W = U / (R ** 2 * Δν)
print(W*1e21)
plt.plot(ν, noms(W), 'rx', label = 'Messdaten')
plt.xlabel(r'$\nu\,/\,\mathrm{Hz}$')
plt.ylabel(r'$W(\nu)\,/\,\mathrm{A^2s}$')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('oxydkathode_frequenzspektrum.pdf')
plt.cla()
plt.clf()

params, cov = curve_fit(linearFunction, np.log(ν[-16:]), np.log(noms(W[-16:])), sigma = np.log(stds(W[-16:])))
    # Wähle hier und im Folgenden mit [-16:] alle Werte mit f < 1000 Hz aus
errors = np.sqrt(np.diag(cov))
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
print('Der Exponent alpha wird zu {} bestimmt'.format(m))

x = np.linspace(1, 2000, 2)
plt.plot(ν[-16:], noms(W[-16:]), 'rx', label = 'Messdaten')
plt.plot(x, x ** noms(m) * np.exp(noms(b)), 'b-', label = 'Ausgleichsgerade')
    # Genau so, weil ln(y)=m*ln(x)+b <=> y=x**m*exp(b)
plt.xlabel(r'$\nu\,/\,\mathrm{Hz}$')
plt.ylabel(r'$W(\nu)\,/\,\mathrm{A^2s}$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1, 2000)
plt.grid()
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('oxydkathode_frequenzspektrum_linearer_teil.pdf')
plt.cla()
plt.clf()

################################################################################
### Reinmetallkathode
################################################################################

U   = [10, 18, 25, 35, 43, 60, 75, 92, 110, 125]
I60 = [0.4, 0.45, 0.45, 0.49, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
I62 = [0.6, 0.7, 0.7, 0.7, 0.7, 0.72, 0.72, 0.75, 0.75, 0.76]
I64 = [0.8, 0.9, 0.9, 0.92, 0.92, 0.97, 0.98, 1.0, 1.0, 1.0]
I65 = [0.9, 1.2, 1.21, 1.24, 1.25, 1.28, 1.3, 1.31, 1.31, 1.33]

plt.plot(U, I60, 'rx', label = r'$I_{\mathrm{H}}=0.60\mathrm{A}$')
plt.plot(U, I62, 'bx', label = r'$I_{\mathrm{H}}=0.62\mathrm{A}$')
plt.plot(U, I64, 'gx', label = r'$I_{\mathrm{H}}=0.64\mathrm{A}$')
plt.plot(U, I65, 'kx', label = r'$I_{\mathrm{H}}=0.65\mathrm{A}$')
plt.xlabel(r'$U\,/\,\mathrm{V}$')
plt.ylabel(r'$I\,/\,\mathrm{mA}$')
plt.ylim(0.4, 1.6)
plt.grid()
plt.legend(loc = 'best', ncol = 2, fontsize = 16)
plt.savefig('reinmetallkathode_linien.pdf')
plt.cla()
plt.clf()

ν, VS, VN, U, ΔU = np.genfromtxt('reinmetallkathode.txt', unpack = True)
Δν = frequency2bandwidth(ν, 'Reinmetallkathode')
U = unp.uarray(U, ΔU)
print(normVoltage(U, VN, VS=VS, VZ=(50*1000*10))*10)
U = normVoltage(U, VN, VS=VS)
R = 4680
W = U / (R ** 2 * Δν)
print(Δν)
print(W*1e21)

plt.plot(ν, noms(W), 'rx', label = 'Messdaten')
plt.xlabel(r'$\nu\,/\,\mathrm{Hz}$')
plt.ylabel(r'$W(\nu)\,/\,\mathrm{A^2s}$')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('reinmetallkathode_frequenzspektrum.pdf')
plt.cla()
plt.clf()

I = U / R ** 2
I0 = 1e-3

Δν1 = np.append(Δν[Δν > 14000], Δν[Δν == 11600])
I1 = np.append(I[Δν > 14000], I[Δν == 11600])
Δν2 = np.append(Δν[Δν < 11500], Δν[Δν == 12100])
I2 = np.append(I[Δν < 11500], I[Δν == 12100])

params1, cov1 = curve_fit(linearFunction, noms(Δν1), noms(I1))
errors1 = np.sqrt(np.diag(cov1))
m1 = ufloat(params1[0], errors1[0])
b1 = ufloat(params1[1], errors1[1])
print(m1)
print(b1)
print('Die Elementarladung wird im hochfrequenten Bereich zu {} bestimmt'.format(m1 / (2 * I0)))

params2, cov2 = curve_fit(linearFunction, noms(Δν2), noms(I2))
errors2 = np.sqrt(np.diag(cov2))
m2 = ufloat(params2[0], errors2[0])
b2 = ufloat(params2[1], errors2[1])
print(m2)
print(b2)
print('Der Elementarladung wird im niederfrequenten Bereich zu {} bestimmt'.format(m2 / (2 * I0)))

x1 = np.linspace(noms(min(Δν1)), noms(max(Δν1)), 2)
x2 = np.linspace(noms(min(Δν2)), noms(max(Δν2)), 2)
plt.plot(noms(Δν1), noms(I1), 'rx')
plt.plot(noms(Δν2), noms(I2), 'bx')
plt.plot(x1, linearFunction(x1, *params1), 'r-')
plt.plot(x2, linearFunction(x2, *params2), 'b-')
plt.xlabel(r'$\Delta\nu\,/\,\mathrm{Hz}$')
plt.ylabel(r'$\bar{I}^2\,/\,\mathrm{A^2}$')
plt.grid()
plt.tight_layout()
plt.savefig('reinmetallkathode_elementarladung.pdf')
plt.cla()
plt.clf()
