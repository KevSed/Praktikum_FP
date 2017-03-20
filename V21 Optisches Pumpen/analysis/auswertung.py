import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.constants as const
import uncertainties.unumpy as unp

from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms
from uncertainties.unumpy import std_devs as stds

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

# Funktionen
def B(I, N, R):
    return µ * (8 / np.sqrt(125)) * I * N / R

def linear(x, m, b):
    return m * x + b

def gF(m):
    return 4 * np.pi * me / (e0 * m) * 1000 # 1000 wegen Steigung in kHz

def gJ(J, S, L):
    return (3.0023 * J * (J + 1) + 1.0023 * (S * (S + 1) - L * (L + 1)) ) / (2 * J * (J + 1))

def nucleonSpin(gJ, gF):
    return gJ / (2 * gF) - 1 / 2

def U(g, B):
    return g * µB * B

def U2(g, B, Mf, EHy):
    return g**2 * µB**2 * B**2 * (1 - Mf) / EHy

def exp(x, a, b, c):
    return a * np.exp(b * x) + c

def hyp(x, a, b, c):
    return a + b / (x - c)

# Naturkonstanten
µ  = const.mu_0
me = ufloat(const.physical_constants['electron mass'][0], const.physical_constants['electron mass'][2])
e0 = ufloat(const.physical_constants['elementary charge'][0], const.physical_constants['elementary charge'][2])
µB = ufloat(const.physical_constants['Bohr magneton'][0], const.physical_constants['Bohr magneton'][2])

# Messwerte
f = np.array([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
Us = np.array([[4.75, 3.20, 5.50, 2.65, 2.64, 2.05, 0.87, 3.24, 1.22, 2.44],  #1.Dip
               [5.95, 5.60, 9.00, 7.40, 8.55, 9.05, 9.12, 7.16, 6.33, 7.61]]) #2.Dip
Uh = np.array([[0.00, 0.10, 0.10, 0.22, 0.27, 0.34, 0.42, 0.42, 0.53, 0.55],  #1.Dip
               [0.00, 0.10, 0.10, 0.22, 0.27, 0.34, 0.42, 0.55, 0.65, 0.70]]) #2.Dip

# Windungszahl N und Radius R der Sweep-, Horizontal- und Vertikalspule
Ns = 11
Nh = 154
Nv = 20
Rs = 0.1639  #Meter
Rh = 0.1579  #Meter
Rv = 0.11735 #Meter

# Rechnungen
Bv_Erde = B(2.45 * 0.1, Nv, Rv)
print('Erdmagnetfeld (vertikal):', Bv_Erde * 1e6, 'µT')

print('Stromstärken der Sweepspule in mA: \n', Us * 0.1 * 1000)
print('Stromstärken der Horizontalspule in mA: \n', Uh * 0.3 * 1000)
Bs = B(Us * 0.1, Ns, Rs)
Bh = B(Uh * 0.3, Nh, Rh)
Bg = Bs + Bh
print('Gesamtmagnetfeldstärken in µT: \n', np.round(Bg * 1e6, decimals = 2))

params87, cov87 = curve_fit(linear, f, Bg[0])
errors87 = np.sqrt(np.diag(cov87))
m87 = ufloat(params87[0], errors87[0])
b87 = ufloat(params87[1], errors87[1])

params85, cov85 = curve_fit(linear, f, Bg[1])
errors85 = np.sqrt(np.diag(cov85))
m85 = ufloat(params85[0], errors85[0])
b85 = ufloat(params85[1], errors85[1])

print('Lineare Regression Rb-85: \n\t Steigung:', m85, 'T/kHz \n\t Achsenabschnitt:', b85, 'T')
print('Lineare Regression Rb-87: \n\t Steigung:', m87, 'T/kHz \n\t Achsenabschnitt:', b87, 'T')

gF85 = gF(m85)
gF87 = gF(m87)
print('gF85:', gF85)
print('gF87:', gF87)
print('Verhältnis gF85/gF87:', gF85/gF87)
print('Verhältnis gF87/gF85:', gF87/gF85)

gJ = gJ(0.5, 0.5, 0)
print('gJ:', gJ)

I85 = nucleonSpin(gJ, gF85)
I87 = nucleonSpin(gJ, gF87)
print('I85:', I85)
print('I87:', I87)

EHy85 = 2.01e-24 #Joule
EHy87 = 4.53e-24 #Joule
Mf = 0
print('Linear 85:', U(gF85, Bg[1][9]), 'J, Linear 87:', U(gF87, Bg[0][9]) ,'J')
print('Quadratisch 85:', U2(gF85, Bg[1][9], Mf, EHy85), 'J, Quadratisch 87:', U2(gF87, Bg[0][9], Mf, EHy87) ,'J')

print('Abgelesenes Isotopenverhältnis Rb-85/Rb-87:', ufloat(4.4,0.2)/ufloat(2.2,0.2))

plt.plot(f, Bg[1] * 10**3, 'bx', label = r'Rb-85')
plt.plot(f, Bg[0] * 10**3, 'rx', label = r'Rb-87')
plt.plot(f, linear(f, *params87) * 10**3, 'r-')
plt.plot(f, linear(f, *params85) * 10**3, 'b-')

plt.ylabel(r'$B\,/\,\mathrm{mT}$')
plt.xlabel(r'$\nu\,/\,\mathrm{kHz}$')
plt.legend()
plt.tight_layout()
plt.savefig('magnetfeld.pdf')
plt.cla()
plt.clf()

df = pd.read_csv('daten.CSV', header = None, names = ['t', 'U'], usecols = [3, 4])
df = df.loc[(df['t'] >= 0) & (df['t'] <= 0.025)]
df['U'] = df['U'] - df['U'].max()
t = df['t']
U = df['U']

params85, cov85 = curve_fit(exp, t, U)
errors85 = np.sqrt(np.diag(cov85))
a85 = ufloat(params85[0], errors85[0])
b85 = ufloat(params85[1], errors85[1])
c85 = ufloat(params85[2], errors85[2])

x = np.linspace(min(t), max(t), 1000)
plt.plot(df['t'], df['U'], 'b.', markersize = 1)
plt.plot(x, exp(x, *params85))
#plt.savefig('anstieg.pdf')
plt.cla()
plt.clf()

U = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 9.756])
T87 = np.array([6.6, 5.6, 1.49, 1.14, 3.58, 3.38, 3.56, 3.7, 3.52, 0.788])
T85 = np.array([13.2, 6.56, 3.28, 3.3, 1.32, 1.66, 1.92, 1.69, 3.58, 3.08])
n87 = np.array([3, 5, 2, 2, 8, 9, 11, 13, 13, 3])
n85 = np.array([4, 4, 3, 4, 2, 3, 4, 4, 9, 8])
T87 /= n87
T85 /= n85
print(T87)
print(T85)
param_bounds = ([0, 0, 0], [np.inf, np.inf, np.inf])
params85, covariance85 = curve_fit(hyp, U, T85, p0 = (0.3, 2.7, 0.3), bounds = param_bounds)
errors85 = np.sqrt(np.diag(covariance85))
a85 = ufloat(params85[0], errors85[0])
b85 = ufloat(params85[1], errors85[1])
c85 = ufloat(params85[2], errors85[2])

params87, covariance87 = curve_fit(hyp, U, T87, p0 = (0.15, 1.8, 0.15), bounds = param_bounds)
errors87 = np.sqrt(np.diag(covariance87))
a87 = ufloat(params87[0], errors87[0])
b87 = ufloat(params87[1], errors87[1])
c87 = ufloat(params87[2], errors87[2])

print(a87, b87, c87)
print(a85, b85, c85)
print('b85/b87', b85/b87)

x = np.logspace(-1, 1.1)
plt.plot(x, hyp(x, *params85), 'r-')
plt.plot(x, hyp(x, *params87), 'b-')
plt.plot(U, T87, 'bx', label = 'Rb-85')
plt.plot(U, T85, 'rx', label = 'Rb-87')

plt.xlabel(r'$U\,/\,\mathrm{V}$')
plt.ylabel(r'$T\,/\,\mathrm{ms}$')
#plt.xlim(0, 10.5)
plt.ylim(0, 7)
plt.legend()
plt.tight_layout()
plt.savefig('periodendauer.pdf')
plt.cla()
plt.clf()
