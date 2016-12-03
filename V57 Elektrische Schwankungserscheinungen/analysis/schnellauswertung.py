import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

#TODO: Fehler auf die Spannungen berücksichtigen

def frequency2bandwidth(f, key):
  if (key == 'Reinmetallkathode'):
    if (f > 100e3):
      if (f == 120e3): return ufloat(11.6e3, 0.0e3)
      if (f == 140e3): return ufloat(13.0e3, 0.0e3)
      if (f == 160e3): return ufloat(14.4e3, 0.0e3)
      if (f == 180e3): return ufloat(16.1e3, 0.0e3)
      if (f == 200e3): return ufloat(18.8e3, 0.0e3)
      if (f == 220e3): return ufloat(21.0e3, 0.2e3)
      if (f == 240e3): return ufloat(20.9e3, 0.2e3)
      if (f == 260e3): return ufloat(20.6e3, 0.2e3)
      if (f == 280e3): return ufloat(21.3e3, 0.2e3)
      if (f == 300e3): return ufloat(23.4e3, 0.4e3)
      if (f == 320e3): return ufloat(24.0e3, 0.4e3)
      if (f == 340e3): return ufloat(24.4e3, 0.4e3)
      if (f == 360e3): return ufloat(24.4e3, 0.4e3)
      if (f == 380e3): return ufloat(24.4e3, 0.4e3)
      if (f == 400e3): return ufloat(24.3e3, 0.4e3)
      if (f == 440e3): return ufloat(23.0e3, 0.4e3)
      if (f == 460e3): return ufloat(21.5e3, 0.4e3)
      print('Fehler: Keine passende Bandbreite für f = {}Hz gefunden.'.format(f))
    if (f > 50e3):
      return ufloat(0.109 * f + 1200, 0.0)
    if (f > 10e3):
      return ufloat(0.135 * f + 50, 0.0)
    if (f > 100):
      return ufloat(0.140 * f + 0.7, 0.0)
    if (f > 2):
      return ufloat(0.15 * f - 0.3, 0.0)

  if (key == 'Oxydkathode'):
    if (f > 100e3):
      if (f == 120e3): return ufloat(12.2e3, 0.0e3)
      if (f == 140e3): return ufloat(14.4e3, 0.0e3)
      if (f == 160e3): return ufloat(16.3e3, 0.0e3)
      if (f == 180e3): return ufloat(18.5e3, 0.0e3)
      if (f == 200e3): return ufloat(21.5e3, 0.2e3)
      if (f == 220e3): return ufloat(23.6e3, 0.2e3)
      if (f == 240e3): return ufloat(24.8e3, 0.4e3)
      if (f == 260e3): return ufloat(24.5e3, 0.4e3)
      if (f == 280e3): return ufloat(27.5e3, 0.5e3)
      if (f == 300e3): return ufloat(29.5e3, 0.3e3)
      if (f == 320e3): return ufloat(33.0e3, 0.5e3)
      if (f == 340e3): return ufloat(33.5e3, 0.5e3)
      if (f == 360e3): return ufloat(35.2e3, 0.3e3)
      if (f == 380e3): return ufloat(35.6e3, 0.2e3)
      if (f == 400e3): return ufloat(35.8e3, 0.2e3)
      if (f == 440e3): return ufloat(35.2e3, 0.2e3)
      if (f == 460e3): return ufloat(34.3e3, 0.2e3)
      print('Fehler: Keine passende Bandbreite für f = {}Hz gefunden.'.format(f))
    if (f > 50e3):
      return ufloat(0.115 * f + 1050, 0.0)
    if (f > 10e3):
      return ufloat(0.135 * f + 50, 0.0)
    if (f > 100):
      return ufloat(0.140 * f + 0.7, 0.0)
    if (f > 2):
      return ufloat(0.15 * f - 0.3, 0.0)

def normVoltage(U, VS, VN):
  # Berechne Gesamtverstärkung (ohne Quadrierung)
  V = VS * VN
  # Normiere alle Spannungen auf Verstärkung 1
  return U * (1 / V) ** 2

def linear(x, m, b):
  return m * x + b

# Einfachschaltung - Eichmessung

f = [1.078, 1.586, 2.036, 2.530, 3.021, 3.485, 3.955, 4.444, 4.963, 5.424, 5.915, 6.432, 6.916, 7.445, 8.016, 8.412, 8.895, 9.432, 9.936, 3.267, 3.725, 4.205, 4.773, 5.263, 5.793]
U = [0.0204, 0.02225, 0.0785, 0.36, 1.265, 2.86, 4.79, 5.86, 5.48, 4.32, 2.97, 1.856, 1.169, 0.7325, 0.443, 0.31575, 0.208, 0.13975, 0.098, 1.99, 3.89, 5.51, 5.78, 4.77, 3.29]

plt.plot(f, U, 'rx')
plt.show()

# Einfachschaltung - schwacher Widerstand

R = [50, 102, 152, 207, 251, 303, 353, 404, 451, 501, 555, 605, 652, 704, 752, 807, 855, 908, 949, 995]
U = [0.049, 0.064, 0.078, 0.093, 0.104, 0.114, 0.129, 0.144, 0.152, 0.167, 0.185, 0.194, 0.208, 0.222, 0.237, 0.250, 0.261, 0.274, 0.283, 0.296]

plt.plot(R, U, 'rx')
plt.show()

# Einfachschaltung - starker Widerstand

R = [1330, 5320, 10240, 15030, 19600, 25100, 30800, 35100, 40600, 44700, 50700, 56300, 60000, 64800, 70100, 76600, 80200, 85700, 91700, 95200, 97300]
U = [0.384, 1.378, 2.816, 4.012, 5.584, 6.764, 8.425, 17.3, 11.8, 12.3, 13.675, 17.85, 19.5, 18.275, 21.3, 26.2, 24.5, 26.2, 30.1, 33.6, 36.7]
plt.plot(R,U, 'rx')
plt.show()

# Korrelatorschaltung - Eichmessung

f = [1.088, 1.609, 2.071, 2.564, 3.034, 3.523, 4.033, 4.480, 4.943, 5.484, 5.976, 6.424, 6.920, 7.410, 7.926, 8.481, 9.064, 9.559, 9.940, 4.211, 4.750, 5.246, 5.755, 5.072, 5.182]
U = [0.92, 1.68, 2.8, 5.072, 9.08, 18.4, 47.0, 145.25, 717.0, 265.0, 87.25, 43.0, 24.8, 17.45, 11.5, 9.16, 7.0, 6.0, 528, 70, 354, 600, 134.75, 863, 727]
plt.plot(f, U, 'rx')
plt.show()

# Korrelatorschaltung - schwacher Widerstand

R = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 995]
U = [0.06, 0.09, 0.068, 0.08992, 0.10912, 0.132, 0.1528, 0.179, 0.206, 0.225, 0.251, 0.270, 0.297, 0.313, 0.330, 0.364, 0.385, 0.408, 0.425, 0.450]
plt.plot(R, U, 'rx')
plt.show()

# Korrelatorschaltung - starker Widerstand

R = [5070, 10130, 15120, 20000, 25000, 30000, 40100, 44900, 50300, 54300, 60300, 64900, 70200, 74900, 80300, 85100, 90100, 95200]
U = [2.172, 4.292, 6.376, 8.624, 11.536, 14.352, 21, 23.2, 27, 29.5, 32.9, 35, 37.9, 40.5, 44.1, 47.3, 53.2, 56.4]
plt.plot(R, U, 'rx')
plt.show()

# Frequenzspektrum - Oxydkathode

ν, VS, VN, U = np.genfromtxt('oxydkathode.txt', unpack=True)
Δν = [frequency2bandwidth(x, 'Oxydkathode') for x in ν]
U = normVoltage(U, VS, VN)
R = 2200
W = [x / (R ** 2 * y) for x, y in zip(U, Δν)]

plt.plot(ν, noms(W), 'rx', label='Messdaten')
plt.xlabel(r'$\nu\,/\,\mathrm{Hz}$')
plt.ylabel(r'$W(\nu)\,/\,\mathrm{A^2s}$')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.show()

# Bestimmung des Exponenten alpha

params, cov = curve_fit(linear, np.log(ν), np.log(noms(W)))
errors = np.sqrt(np.diag(cov))
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
print(m)

# Frequenzspektrum - Reinmetallkathode

ν, VS, VN, U = np.genfromtxt('reinmetallkathode.txt', unpack=True)
Δν = [frequency2bandwidth(x, 'Reinmetallkathode') for x in ν]
U = normVoltage(U, VS, VN)
R = 4680
W = [x / (R ** 2 * y) for x, y in zip(U, Δν)]

plt.plot(ν, noms(W), 'rx')
plt.xscale('log')
plt.yscale('log')
plt.show()

# Bestimmung der Elementarladung

I = U / R **2
plt.plot(ν, I, 'rx')
plt.show()
