import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from funktionen import *
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

################################################################################
### Oxydkathode
################################################################################

ν, VS, VN, U = np.genfromtxt('oxydkathode.txt', unpack = True)
Δν = frequency2bandwidth(ν, 'Oxydkathode')
U = normVoltage(U, VS, VN)
    # Vielleicht zwecks Übersichtlichkeit besser auf VZ = 100 normieren. Mal schauen...
R = 2200
W = U / (R ** 2 * Δν)

plt.plot(ν, noms(W), 'rx', label = 'Messdaten')
plt.xlabel(r'$\nu\,/\,\mathrm{Hz}$')
plt.ylabel(r'$W(\nu)\,/\,\mathrm{A^2s}$')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend(loc = 'best')
plt.tight_layout()
plt.title('Oxydkathode - Frequenzspektrum')
plt.savefig('oxydkathode_frequenzspektrum.pdf')
plt.cla()
plt.clf()

params, cov = curve_fit(linearFunction, np.log(ν[-16:]), np.log(noms(W[-16:])))
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
plt.title('Oxydkathode - Frequenzspektrum (linearer Teil)')
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
plt.title('Reinmetallkathode - Charakteristische Linien')
plt.savefig('reinmetallkathode_linien.pdf')
plt.cla()
plt.clf()

ν, VS, VN, U = np.genfromtxt('reinmetallkathode.txt', unpack = True)
Δν = frequency2bandwidth(ν, 'Reinmetallkathode')
U = normVoltage(U, VS, VN)
R = 4680
W = U / (R ** 2 * Δν)

plt.plot(ν, noms(W), 'rx', label = 'Messdaten')
plt.xlabel(r'$\nu\,/\,\mathrm{Hz}$')
plt.ylabel(r'$W(\nu)\,/\,\mathrm{A^2s}$')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend(loc = 'best')
plt.tight_layout()
plt.title('Reinmetallkathode - Frequenzspektrum')
plt.savefig('reinmetallkathode_frequenzspektrum.pdf')
plt.cla()
plt.clf()

sys.exit()






















# Bestimmung der Elementarladung

I = (U/1000000) / R **2
params, cov = curve_fit(linear, noms(Δν), I)
print(params)
x = np.linspace(0, 25000)
plt.plot(noms(Δν), I, 'rx')
#plt.plot(x, linear(x,*params), 'b-')
plt.xlabel(r'$\Delta\nu\,/\,\mathrm{Hz}$')
plt.ylabel(r'$\bar{I}^2\,/\,\mathrm{A^2}$')
plt.title('Reinmetallkathode - Elementarladung')
plt.savefig('reinmetallkathode_elementarladung.pdf')
plt.cla()
plt.clf()
#plt.show()
