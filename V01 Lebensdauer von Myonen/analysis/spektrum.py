import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

n   = np.genfromtxt('analysis/spektrum.Spe', skip_header=12, skip_footer=14, unpack=True)
bkg = 5.221370
ch  = range(0,512)

plt.step(ch, n, where='mid', label='Messwerte')
plt.xlabel('Kanal')
plt.ylabel('Anzahl an Impulse')
plt.xlim(0,511)
plt.grid()
plt.legend(loc='best', numpoints=1)
plt.tight_layout()
plt.savefig('analysis/spektrum1.pdf')
plt.clf()

rem = np.array([0,1,2,12,13,14,15,16,17,18])
#rem = np.append(rem, np.arange(441,511))
rem = np.append(rem, np.arange(343,512)) # Korrektur 1. Abgabe

n   = np.delete(n, rem)
ch  = np.delete(ch, rem)

def f(x,a,b,c):
    return a*np.exp(-x/b)+c

#params, cov = curve_fit(f, ch[1:332], n[1:332], sigma=np.sqrt(n[1:332]))
params, cov = curve_fit(f, ch, n, sigma=np.sqrt(n)) # Korrektur 1. Abgabe
errors = np.sqrt(np.diag(cov))
N = ufloat(params[0], errors[0])
l = ufloat(params[1], errors[1])
U = ufloat(params[2], errors[2])

temp = np.genfromtxt('analysis/kalibrierung.txt')
m = ufloat(temp[0], temp[2])
b = ufloat(temp[1], temp[3])
print('Fitparameter: \n N = {} \n l = {} \n U = {}'.format(N, 1/l, U))
print('Berechnete Myonlebensdauer: \n ({})µs'.format(m*l))

x = np.linspace(0,332,1000)
plt.errorbar(ch, n, yerr=np.sqrt(n), fmt='rx', ms=5, mew=0.6, elinewidth=0.6, capsize=2, capthick=0.6, label='Messwerte')
plt.plot([0,332], [bkg,bkg], 'b-', label='Erwarteter Untergrund')
plt.plot(x, f(x, *params), 'g-', label='Ausgleichskurve')

plt.xlabel('Kanal')
plt.ylabel('Anzahl an Impulse')
plt.xlim(0,332)
plt.grid()
plt.legend(loc='best', numpoints=1)
plt.tight_layout()
plt.savefig('analysis/spektrum2.pdf')

plt.yscale('log')
plt.savefig('analysis/spektrum3.pdf')
