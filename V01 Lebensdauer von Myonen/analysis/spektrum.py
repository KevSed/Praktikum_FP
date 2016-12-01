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
rem = np.append(rem, np.arange(343,512)) # Korrektur 1. Abgabe

n   = np.delete(n, rem)
ch  = np.delete(ch, rem)

temp = np.genfromtxt('analysis/kalibrierung.txt')
m = ufloat(temp[0], temp[2])
b = ufloat(temp[1], temp[3])

# Korrektur 1. Abgabe
def ch2time(ch):
    return noms(m)*ch+noms(b)

def f(x,a,b,c):
    return a*np.exp(-x/b)+c

params, cov = curve_fit(f, ch, n, sigma=np.sqrt(n)) # Korrektur 1. Abgabe
errors = np.sqrt(np.diag(cov))
N = ufloat(params[0], errors[0])
l = ufloat(params[1], errors[1])
U = ufloat(params[2], errors[2])

print('Fitparameter: \n N = {} \n l = {} \n U = {}'.format(N, 1/l, U))
print('Berechnete Myonlebensdauer: \n ({})µs'.format(m*l))

fig, ax1 = plt.subplots()
ax2 = ax1.twiny()
fig.subplots_adjust(bottom=0.20)

ax2.set_frame_on(True)
ax2.patch.set_visible(False)
ax2.xaxis.set_ticks_position('bottom')
ax2.xaxis.set_label_position('bottom')
ax2.spines['bottom'].set_position(('outward', 50))

x = np.linspace(0,343,1000)
ax1.errorbar(ch, n, yerr=np.sqrt(n), fmt='rx', ms=5, mew=0.6, elinewidth=0.6, capsize=2, capthick=0.6, label='Messwerte')
ax1.plot(x, f(x, *params), 'g-', label='Ausgleichskurve')
ax1.plot([], [], 'b-', label='Erwarteter Untergrund')
ax2.plot([ch2time(0),ch2time(343)], [bkg,bkg], 'b-', label='Erwarteter Untergrund')

ax1.set_xlabel('Kanal')
ax2.set_xlabel(r'$\displaystyle\frac{\Delta t}{\mathrm{\mu s}}$')
ax1.set_xlim(0,343)
ax2.set_xlim(ch2time(0),ch2time(343))
ax1.legend(loc='best', numpoints=1)
ax1.set_ylabel('Anzahl an Impulse')
plt.tight_layout()
plt.savefig('analysis/spektrum2.pdf')

plt.yscale('log')
plt.savefig('analysis/spektrum3.pdf')
