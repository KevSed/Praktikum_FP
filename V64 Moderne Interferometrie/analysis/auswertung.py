import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties
from uncertainties import ufloat
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

print(contrast(U_max, U_min))

params, cov = curve_fit(fit2, np.radians(winkel), contrast(U_max, U_min))
A = params[0]
A_err = np.sqrt(cov[0][0])
delta = params[1]
delta_err = np.sqrt(cov[1][1])

print(A)
print(A_err)
print(delta)
print(delta_err)

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

number, counts = np.genfromtxt('glas.txt', unpack='True')
phi_1 = number*10
phi_1 = np.radians(phi_1)
phi_2 = phi_1 + np.radians(10)

vac_wavelen = 632.99e-9
thickness = 0.5e-3
# ref_index = 1/(1 - counts * vac_wavelen/(thickness * (phi_2**2 - phi_1**2)))
alpha = counts*vac_wavelen/(2*thickness)
ref_index = (alpha**2+2*(1-np.cos(phi_1))*(1-alpha))/(2*(1-np.cos(phi_1)-alpha))

print('''
Refraction index of glass slabs:   {:.4f} ± {:.4f}'''
      .format(np.mean(ref_index), ref_index.std(ddof = 1)))

######################################################################
# calculate refraction index of CO2 chamber from interference
# counts for different pressures of mentioned gas. (50mbar steps)
######################################################################

pressure, counts1, counts2, counts3 = np.genfromtxt('CO2.txt', unpack='True')
xplotlin = np.linspace(pressure[0], pressure[-1], 2)
length = ufloat(100e-3, 0.1e-3)

ref_index1 = (counts1 * vac_wavelen / length + 1)**(.5)
ref_index2 = (counts2 * vac_wavelen / length + 1)**(.5)
ref_index3 = (counts3 * vac_wavelen / length + 1)**(.5)

ref_index1_err = np.zeros(len(ref_index1))
ref_index2_err = np.zeros(len(ref_index1))
ref_index3_err = np.zeros(len(ref_index1))
refindex1 = np.zeros(len(ref_index1))
refindex2 = np.zeros(len(ref_index1))
refindex3 = np.zeros(len(ref_index1))

for i in range(0,len(ref_index1)):
    ref_index1_err[i] = stds(ref_index1[i])
    refindex1[i] = noms(ref_index1[i])

for i in range(0,len(ref_index2)):
    ref_index2_err[i] = stds(ref_index2[i])
    refindex2[i] = noms(ref_index2[i])

for i in range(0,len(ref_index3)):
    ref_index3_err[i] = stds(ref_index3[i])
    refindex3[i] = noms(ref_index3[i])

print('''
Refraction index of CO2:   {:.5f} for (24.7°, 18mbar)
Refraction index of CO2:   {:.5f} for (24.8°, 4mbar)
Refraction index of CO2:   {:.5f} for (24.8, 3mbar)'''.format(np.mean(ref_index1),np.mean(ref_index2),np.mean(ref_index3)))


def lin(a, b, x):
    return a*x+b

fit1_params, fit1_cov = np.polyfit(pressure, refindex1, 1, cov=True)
errors1 = np.sqrt(np.diag(fit1_cov))

fit2_params, fit2_cov = np.polyfit(pressure, refindex2, 1, cov=True)
errors2 = np.sqrt(np.diag(fit2_cov))

fit3_params, fit3_cov = np.polyfit(pressure, refindex3, 1, cov=True)
errors3 = np.sqrt(np.diag(fit3_cov))

steig = np.array([fit1_params[0], fit2_params[0], fit3_params[0]])
achse = np.array([fit1_params[1], fit2_params[1], fit3_params[1]])

print('''
Results of the regression for CO2:
Messung 1: a = {:.5g} ± {:.7g}   b = {:.5g} ± {:.7g}
Messung 2: a = {:.5g} ± {:.7g}   b = {:.5g} ± {:.7g}
Messung 3: a = {:.5g} ± {:.7g}   b = {:.5g} ± {:.7g}
Mean: a = {:.5g} ± {:.7g}   b = {:.5g} ± {:.7g}
'''.format(fit1_params[0], errors1[0], 1-fit1_params[1], errors1[1], fit2_params[0], errors2[0], 1-fit2_params[1], errors2[1], fit3_params[0], errors3[0], 1-fit3_params[1], errors3[1], np.mean(steig), np.std(steig), 1-np.mean(achse), np.std(achse)))

m = ufloat(np.mean(steig), np.std(steig))
b = ufloat(np.mean(achse), np.std(achse))

n_norm = m*T/T_0*p_0+b
#n_norm = m*p_0*((3*A)/(R*T_0))+b

print('''
Calculated refraction index at 1018hPa:
CO2: n = {:.8g} ± {:.8g}
Calculated refraction index at normalconditions:
CO2: n = {:.8g}
'''.format(np.mean(steig)*1018+np.mean(achse), np.std(achse)+np.std(steig),n_norm))

x = np.linspace(0, 1000, 1000)
plt.errorbar(pressure, refindex1, 100*ref_index1_err, fmt='x', label=r'Messung 1 (24.7°, 18mbar) Fehler 100fach vergrößert.')
plt.plot(x, lin(*fit1_params, x), label='Ausgleichsgeraden')
plt.title(r'Druckverteilung des Brechungsindex von $CO_2$')
plt.xlabel('Druck p/mbar')
plt.ylabel('Brechungsindex n')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('CO2_1.pdf')

plt.close()
plt.errorbar(pressure, refindex2, 100*ref_index2_err, fmt='x', label=r'Messung 2 (24.8°, 4mbar) Fehler 100fach vergrößert.')
plt.plot(x, lin(*fit2_params, x), label='Ausgleichsgeraden')
plt.title(r'Druckverteilung des Brechungsindex von $CO_2$')
plt.xlabel('Druck p/mbar')
plt.ylabel('Brechungsindex n')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('CO2_2.pdf')
plt.close()
plt.errorbar(pressure, refindex3, 100*ref_index3_err, fmt='x', label=r'Messung 3 (24.8, 3mbar) Fehler 100fach vergrößert.')
plt.plot(x, lin(*fit3_params, x), label='Ausgleichsgeraden')
plt.title(r'Druckverteilung des Brechungsindex von $CO_2$')
plt.xlabel('Druck p/mbar')
plt.ylabel('Brechungsindex n')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('CO2_3.pdf')
plt.close()

######################################################################
# calculate refraction index of air chamber from interference
# counts for different pressures of mentioned gas. (50mbar steps)
######################################################################

pressure, counts1, counts2, counts3 = np.genfromtxt('luft.txt', unpack='True')
length = ufloat(100e-3, 0.1e-3)

ref_index1 = (counts1 * vac_wavelen / length + 1)**(.5)
ref_index2 = (counts2 * vac_wavelen / length + 1)**(.5)
ref_index3 = (counts3 * vac_wavelen / length + 1)**(.5)

ref_index1_err = np.zeros(len(ref_index1))
ref_index2_err = np.zeros(len(ref_index1))
ref_index3_err = np.zeros(len(ref_index1))
refindex1 = np.zeros(len(ref_index1))
refindex2 = np.zeros(len(ref_index1))
refindex3 = np.zeros(len(ref_index1))

for i in range(0,len(ref_index1)):
    ref_index1_err[i] = stds(ref_index1[i])
    refindex1[i] = noms(ref_index1[i])

for i in range(0,len(ref_index2)):
    ref_index2_err[i] = stds(ref_index2[i])
    refindex2[i] = noms(ref_index2[i])

for i in range(0,len(ref_index3)):
    ref_index3_err[i] = stds(ref_index3[i])
    refindex3[i] = noms(ref_index3[i])

fit1_params, fit1_cov = np.polyfit(pressure, refindex1, 1, cov=True)
errors1 = np.sqrt(np.diag(fit1_cov))

fit2_params, fit2_cov = np.polyfit(pressure, refindex2, 1, cov=True)
errors2 = np.sqrt(np.diag(fit2_cov))

fit3_params, fit3_cov = np.polyfit(pressure, refindex3, 1, cov=True)
errors3 = np.sqrt(np.diag(fit3_cov))

steig = np.array([fit1_params[0], fit2_params[0], fit3_params[0]])
achse = np.array([fit1_params[1], fit2_params[1], fit3_params[1]])

ref_index1 = (counts1 * vac_wavelen / length + 1)**(.5)
ref_index2 = (counts2 * vac_wavelen / length + 1)**(.5)
ref_index3 = (counts3 * vac_wavelen / length + 1)**(.5)

print('''
Refraction index of air:   {:.8f}
Refraction index of air:   {:.8f}
Refraction index of air:   {:.8f} '''.format(np.mean(ref_index1),np.mean(ref_index2),np.mean(ref_index3)))

print('''
Results of regression for air:
Measurement 1: a = {:.5g} ± {:.7g}   b = {:.5g} ± {:.7g}
Measurement 2: a = {:.5g} ± {:.7g}   b = {:.5g} ± {:.7g}
Measurement 3: a = {:.5g} ± {:.7g}   b = {:.5g} ± {:.7g}
Mean: a = {:.5g} ± {:.7g}   b = {:.5g} ± {:.7g}
'''.format(fit1_params[0], errors1[0], 1-fit1_params[1], errors1[1], fit2_params[0], errors2[0], 1-fit2_params[1], errors2[1], fit3_params[0], errors3[0], 1-fit3_params[1], errors3[1], np.mean(steig), np.std(steig), 1-np.mean(achse), np.std(achse)))

m = ufloat(np.mean(steig), np.std(steig))
b = ufloat(np.mean(achse), np.std(achse))

n_norm = m*T/T_0*p_0+b

print('''
Calculated refraction index at 1018hPa:
Air: n = {:.8g} ± {:.8g}
Calculated refraction index at normalconditions:
Air: n = {:.8g}
'''.format(np.mean(steig)*1018+np.mean(achse), np.std(achse)+np.std(steig),n_norm))

plt.errorbar(pressure, refindex1, 100*ref_index1_err, fmt='x', markersize=4, elinewidth=1, label=r'Messung 1 Fehler 100fach vergrößert.')
plt.plot(x, lin(*fit1_params, x), label='Ausgleichsgeraden')
plt.title(r'Druckverteilung des Brechungsindex von Luft')
plt.xlabel('Druck p/mbar')
plt.ylabel('Brechungsindex n')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Luft1.pdf')
plt.close()
plt.errorbar(pressure, refindex2, 100*ref_index2_err, fmt='x', markersize=4, elinewidth=1, label=r'Messung 2 Fehler 100fach vergrößert.')
plt.plot(x, lin(*fit2_params, x), label='Ausgleichsgeraden')
plt.title(r'Druckverteilung des Brechungsindex von Luft')
plt.xlabel('Druck p/mbar')
plt.ylabel('Brechungsindex n')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Luft2.pdf')
plt.close()
plt.errorbar(pressure, refindex3, 100*ref_index3_err, fmt='x', markersize=4, elinewidth=1, label=r'Messung 3 Fehler 100fach vergrößert.')
plt.plot(x, lin(*fit3_params, x), label='Ausgleichsgeraden')
plt.title(r'Druckverteilung des Brechungsindex von Luft')
plt.xlabel('Druck p/mbar')
plt.ylabel('Brechungsindex n')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('Luft3.pdf')
plt.close()
