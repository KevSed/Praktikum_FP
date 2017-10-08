import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties
from uncertainties import ufloat
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from uncertainties.umath import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

# KONSTANTEN
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
p_0 = 1013.25 # mbar
T_0 = 15 # °C
A = 1e2
R = ufloat(8.3144598, 0.0000048) #J mol-1 K-1
T = 24.8 #°C CO2
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

######################################################################
# calculate contrast from voltage measurments for different
# angles of the polarisation filter. fit of a sin(2x)².
######################################################################

winkel, U_min, U_max = np.genfromtxt('kontrast.txt', unpack='True')
def contrast(max, min):
    return (max-min)/(max+min)

def fit(a, x):
    return a*(np.sin(2*x))**2

def fit2(phi, A, delta):
    return A * abs(np.sin(2 * phi + delta))

params, cov = curve_fit(fit2, np.radians(winkel), contrast(U_max, U_min))
A = params[0]
A_err = np.sqrt(cov[0][0])
delta = params[1]
delta_err = np.sqrt(cov[1][1])

x = np.linspace(0, 180, 10000)
plt.plot(winkel, contrast(U_max, U_min), 'rx', label='Messwerte')
plt.plot(x, fit2(np.radians(x), *params), label='Ausgleichskurve')
plt.xlabel(r' Winkel $\varphi$ in Grad')
plt.ylabel('Kontrast')
plt.ylim(0, 1)
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('kontrast.pdf')
plt.close()


######################################################################
# calculate refraction index of the glass slobe from interference
# counts for different angles of the glass slobe (10° steps).
######################################################################
#
#number, counts = np.genfromtxt('glas.txt', unpack='True')
#phi_1 = number*10
#phi_1 = np.radians(phi_1)
#phi_2 = phi_1 + np.radians(10)
#
vac_wavelen = 632.99e-9
#thickness = 0.5e-3
## ref_index = 1/(1 - counts * vac_wavelen/(thickness * (phi_2**2 - phi_1**2)))
#alpha = counts*vac_wavelen/(2*thickness)
#ref_index = (alpha**2+2*(1-np.cos(phi_1))*(1-alpha))/(2*(1-np.cos(phi_1)-alpha))
#
#print('''
#Refraction index of glass slabs:   {:.4f} ± {:.4f}'''
#      .format(np.mean(ref_index), ref_index.std(ddof = 1)))
#
################################################################################
# Calculate refraction index of CO2 chamber from interference
# counts for different pressures of mentioned gas. (50mbar steps)
################################################################################

p, M1, M2, M3 = np.genfromtxt('CO2.txt', unpack='True')
L = ufloat(100e-3, 0.1e-3)

n1 = M1 * vac_wavelen / L + 1
n2 = M2 * vac_wavelen / L + 1
n3 = M3 * vac_wavelen / L + 1

n1_err = stds(n1)
n2_err = stds(n2)
n3_err = stds(n3)
n1 = noms(n1)
n2 = noms(n2)
n3 = noms(n3)

print(
'''
Refraction index of CO2:   {:.5f} for (24.7°, 18mbar)
Refraction index of CO2:   {:.5f} for (24.8°,  4mbar)
Refraction index of CO2:   {:.5f} for (24.8°,  3mbar)
'''.format(np.mean(n1), np.mean(n2), np.mean(n3))
)

fit1_params, fit1_cov = np.polyfit(p, n1, 1, cov = True)
errors1 = np.sqrt(np.diag(fit1_cov))
fit2_params, fit2_cov = np.polyfit(p, n2, 1, cov = True)
errors2 = np.sqrt(np.diag(fit2_cov))
fit3_params, fit3_cov = np.polyfit(p, n3, 1, cov = True)
errors3 = np.sqrt(np.diag(fit3_cov))

m = np.array([fit1_params[0], fit2_params[0], fit3_params[0]])
b = np.array([fit1_params[1], fit2_params[1], fit3_params[1]])

print(
'''
Results of the regression for CO2:
Measurement 1: a = {:.5g} ± {:.7g}   1 - b = {:.5g} ± {:.7g}
Measurement 2: a = {:.5g} ± {:.7g}   1 - b = {:.5g} ± {:.7g}
Measurement 3: a = {:.5g} ± {:.7g}   1 - b = {:.5g} ± {:.7g}
Mean: a = {:.5g} ± {:.7g}   1 - b = {:.5g} ± {:.7g}
'''.format(fit1_params[0], errors1[0], 1 - fit1_params[1], errors1[1],
           fit2_params[0], errors2[0], 1 - fit2_params[1], errors2[1],
           fit3_params[0], errors3[0], 1 - fit3_params[1], errors3[1],
           np.mean(m), np.std(m), 1 - np.mean(b), np.std(b))
)

m = ufloat(np.mean(m), np.std(m))
b = ufloat(np.mean(b), np.std(b))

n_norm = m * T / T_0 * p_0 + b

print(
'''
Calculated refraction index at 1018hPa:
CO2: n = {}
Calculated refraction index at normal conditions:
CO2: n = {}
'''.format(unp.sqrt(m * 1018 + b), unp.sqrt(n_norm))
)

def linear(x, m, b):
    return m * x + b

x = np.linspace(p[0], p[-1], 2)

plt.errorbar(p, n1, 100 * n1_err, fmt = 'x',
             label = 'Messung 1')
plt.plot(x, linear(x, *fit1_params), label = 'Ausgleichsgerade')
plt.xlabel(r'Druck $p$ / mbar')
plt.ylabel(r'Brechungsindex $n$')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('CO2_1.pdf')
plt.close()

plt.errorbar(p, n2, 100 * n2_err, fmt = 'x',
             label = 'Messung 2')
plt.plot(x, linear(x, *fit2_params), label = 'Ausgleichsgerade')
plt.xlabel(r'Druck $p$ / mbar')
plt.ylabel(r'Brechungsindex $n$')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('CO2_2.pdf')
plt.close()

plt.errorbar(p, n3, 100 * n3_err, fmt = 'x',
             label = 'Messung 3')
plt.plot(x, linear(x, *fit3_params), label = 'Ausgleichsgerade')
plt.xlabel(r'Druck $p$ / mbar')
plt.ylabel(r'Brechungsindex $n$')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('CO2_3.pdf')
plt.close()

################################################################################
# Calculate refraction index of air chamber from interference
# counts for different pressures of mentioned gas. (50mbar steps)
################################################################################

p, M1, M2, M3 = np.genfromtxt('luft.txt', unpack='True')
L = ufloat(100e-3, 0.1e-3)

n1 = M1 * vac_wavelen / L + 1
n2 = M2 * vac_wavelen / L + 1
n3 = M3 * vac_wavelen / L + 1

n1_err = stds(n1)
n2_err = stds(n2)
n3_err = stds(n3)
n1 = noms(n1)
n2 = noms(n2)
n3 = noms(n3)

print(
'''
Refraction index of air:   {:.5f} for (24.7°, 18mbar)
Refraction index of air:   {:.5f} for (24.8°,  4mbar)
Refraction index of air:   {:.5f} for (24.8°,  3mbar)
'''.format(np.mean(n1), np.mean(n2), np.mean(n3))
)

fit1_params, fit1_cov = np.polyfit(p, n1, 1, cov = True)
errors1 = np.sqrt(np.diag(fit1_cov))
fit2_params, fit2_cov = np.polyfit(p, n2, 1, cov = True)
errors2 = np.sqrt(np.diag(fit2_cov))
fit3_params, fit3_cov = np.polyfit(p, n3, 1, cov = True)
errors3 = np.sqrt(np.diag(fit3_cov))

m = np.array([fit1_params[0], fit2_params[0], fit3_params[0]])
b = np.array([fit1_params[1], fit2_params[1], fit3_params[1]])

print(
'''
Results of the regression for air:
Measurement 1: a = {:.5g} ± {:.7g}   1 - b = {:.5g} ± {:.7g}
Measurement 2: a = {:.5g} ± {:.7g}   1 - b = {:.5g} ± {:.7g}
Measurement 3: a = {:.5g} ± {:.7g}   1 - b = {:.5g} ± {:.7g}
Mean: a = {:.5g} ± {:.7g}   1 - b = {:.5g} ± {:.7g}
'''.format(fit1_params[0], errors1[0], 1 - fit1_params[1], errors1[1],
           fit2_params[0], errors2[0], 1 - fit2_params[1], errors2[1],
           fit3_params[0], errors3[0], 1 - fit3_params[1], errors3[1],
           np.mean(m), np.std(m), 1 - np.mean(b), np.std(b))
)

m = ufloat(np.mean(m), np.std(m))
b = ufloat(np.mean(b), np.std(b))

n_norm = m * T / T_0 * p_0 + b

print(
'''
Calculated refraction index at 1018hPa:
air: n = {}
Calculated refraction index at normal conditions:
air: n = {}
'''.format(unp.sqrt(m * 1018 + b), unp.sqrt(n_norm))
)

def linear(x, m, b):
    return m * x + b

x = np.linspace(p[0], p[-1], 2)

plt.errorbar(p, n1, 100 * n1_err, fmt = 'x',
             label = 'Messung 1')
plt.plot(x, linear(x, *fit1_params), label = 'Ausgleichsgerade')
plt.xlabel(r'Druck $p$ / mbar')
plt.ylabel(r'Brechungsindex $n$')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('Luft_1.pdf')
plt.close()

plt.errorbar(p, n2, 100 * n2_err, fmt = 'x',
             label = 'Messung 2')
plt.plot(x, linear(x, *fit2_params), label = 'Ausgleichsgerade')
plt.xlabel(r'Druck $p$ / mbar')
plt.ylabel(r'Brechungsindex $n$')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('Luft_2.pdf')
plt.close()

plt.errorbar(p, n3, 100 * n3_err, fmt = 'x',
             label = 'Messung 3')
plt.plot(x, linear(x, *fit3_params), label = 'Ausgleichsgerade')
plt.xlabel(r'Druck $p$ / mbar')
plt.ylabel(r'Brechungsindex $n$')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig('Luft_3.pdf')
plt.close()
