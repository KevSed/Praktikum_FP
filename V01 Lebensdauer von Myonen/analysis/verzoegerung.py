import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

t    = np.array([-10, -6, -5, -4, -3, -2, -1,  0,  1,  2,2.5,  3,3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8,8.5,  9, 10]) # ns
rate = np.array([586,593,579,599,605,599,562,620,577,577,574,604,558,646,621,616,668,627,620,639,597,582,568,550,527]) # 1/20s

plt.errorbar(t, rate/20, yerr=np.sqrt(rate)/20, fmt='rx', ms=13, mew=2, elinewidth=1, capsize=5, capthick=1, label='Messwerte')
plt.xlabel(r'$\displaystyle\frac{T_\mathrm{VZ}}{\mathrm{ns}}$')
plt.ylabel(r'$\displaystyle\frac{N}{\mathrm{s}}$')
plt.xlim(-11,11)
plt.grid()
plt.legend(loc='best', numpoints=1)
plt.tight_layout()
plt.savefig('analysis/verzoegerung.pdf')
