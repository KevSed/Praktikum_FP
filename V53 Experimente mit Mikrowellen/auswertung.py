import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#================#
#     Teil A     #
#================#

# Daten
mod1U = [190,200,205];
mod1A = [0,1.08,0];

mod2U = [90,96,105];
mod2A = [0,0.92,0];

mod3U = [110,120,124];
mod3A = [0,1.13,0];

# Ausgleichsrechnung
def f(x, a, b, c):
    return a * x**2 + b * x + c
    
params_mod1, cov_mod1 = curve_fit(f, mod1U, mod1A)
params_mod2, cov_mod2 = curve_fit(f, mod2U, mod2A)
params_mod3, cov_mod3 = curve_fit(f, mod3U, mod3A)

print(params_mod1)
print(params_mod2)
print(params_mod3)

# Plotten
x = np.linspace(90, 205, 2000)
plt.plot(mod1U, mod1A, 'rx', label = 'Messdaten')
plt.plot(mod2U, mod2A, 'rx') 
plt.plot(mod3U, mod3A, 'rx') 
plt.plot(x, f(x, *params_mod1), 'k-', label = 'Regressionskurven')
plt.plot(x, f(x, *params_mod2), 'k-')
plt.plot(x, f(x, *params_mod3), 'k-')
plt.ylim(0,1.5)
plt.grid()
plt.xlabel(r'$U\,[\mathrm{V}]$')
plt.ylabel(r'$A\,[\mathrm{V}]$')
plt.legend(loc = "best")
plt.tight_layout()
plt.savefig("build\plot_modenkurve.pdf")
plt.close()

#================#
#     Teil B     #
#================#

# Daten
Mic = [0,1.05,1.5,1.77,2.02,2.29];
SWR = [0,2,4,6,8,10];
The = [0,3,4.5,6,8,10];

# Plotten
plt.plot(Mic, SWR, 'rx', label = 'SWR-Meter')
plt.plot(Mic, The, 'b+', label = 'Eichkurve') 
plt.grid()
plt.xlabel(r'Mikrometereinstellung$\,[\mathrm{mm}]$')
plt.ylabel(r'Dämpfung$\,[\mathrm{dB}]$')
plt.legend(loc = "best")
plt.tight_layout()
plt.savefig("build\plot_daempfungskurve.pdf")
plt.close()