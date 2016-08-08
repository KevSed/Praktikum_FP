import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def f(x):
    return [x for ind, x in enumerate(x) if True]

r, Δr, Θ_deg, ΔΘ_deg, Θ_rad, ΔΘ_rad, d, Δd, d1di, Δd1di, N, diamond, diffdiamond, Δdiffdiamond, a, Δa, c2, Δc2 = np.genfromtxt('data_probe8.txt', unpack = True)
x = np.linspace(0,12,13)

plt.errorbar(diffdiamond[1:], x[1:], xerr=Δdiffdiamond[1:], fmt='rx', label='Differenzen')
plt.grid()
plt.xlabel(r"Ein x-Label")
plt.ylabel(r"Ein y-Label")
plt.legend(loc='best')
plt.xlim(-0.3,0.3)
plt.tight_layout()
plt.savefig("../build/plot_probe8_1.pdf")
plt.close()

def g(x, m, b):
    return m*x+b

params, cov = curve_fit(g, f(c2), f(a))
m = params[0]
Δm = np.sqrt(cov[0][0])
b = params[1]
Δb = np.sqrt(cov[1][1])

x = np.linspace(0,1,2)
plt.errorbar(f(c2), f(a), xerr=f(Δc2), yerr=f(Δa), fmt='kx', label='Ein Label')
plt.plot(x, g(x, m, b), 'b-')
plt.grid()
plt.xlabel(r"Ein x-Label")
plt.ylabel(r"Ein y-Label")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("../build/plot_probe8_2.pdf")
plt.close()
