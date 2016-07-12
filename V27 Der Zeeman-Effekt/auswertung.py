import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit

# Daten laden
I  = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]; # Strom in A
B1 = [6,70,142,212,288,335,422,492,560,630,698,785,873,947,1020,1104]; # B-Feld bei steigendem Strom in mT
B2 = [8,76,145,213,284,346,422,488,562,642,733,810,890,967,1040,1104]; # B-Feld bei fallendem Strom in mT

B  = [7,73,143.5,212.5,286,340.5,422,490,561,636,615.5,797.5,881.5,957,1030,1104]

# Auswertung
def f(x, m, b): return m*x+b #mT
params, cov = curve_fit(f, I, B)
m = params[0]
dm = np.sqrt(cov[0][0])
b = params[1]
db = np.sqrt(cov[1][1])
print('m = {:4f}'.format(m))
print('m_err = {:4f}'.format(dm))
print('b = {:4f}'.format(b))
print('b_err = {:4f}'.format(db))

# Graphen plotten
x = np.linspace(0, 16)
plt.plot(I, B1, 'rx', label='steigender Strom')
plt.plot(I, B2, 'b+', label='fallender Strom')
plt.plot(x, f(x, m, b), 'k-', label='Ausgleichsgerade')

plt.legend(loc = 'best')
plt.xlabel(r'Stromstärke $I$ [A]')
plt.ylabel(r'Magnetfeld $B$ [mT]')
plt.ylim(-40, 1200)
plt.grid()
plt.savefig("build/plot_magnetfeld.pdf")

# Konstanten und Daten laden für die Bestimmung der Lande-Faktoren
h = 6.626e-34 # Planck'sches Wirkungsquantum in [J/s]
c = 299.7e6 # Lichtgeschwindigkeit in [m/s]

Lr = 48.91e-12 # Dispersionsgebiet für die rote Linie in [m]
Lb = 26.95e-12 # Dispersionsgebiet für die blaue Linie in [m]
lb = 480e-9 # Wellenlänge der blauen Linie in [m]
lr = 643.8e-9 # Wellenlänge der roten Linie in [m]
m_B = 9.274e-24 # Bohr'sches Magneton in [J/T]

# Rote sigma Linie
rs, rs_b = np.loadtxt('data/rs.txt', unpack = True)

B_rs = f(10, m, b)*10**(-3) #T
print('Magnetfeldstärke der roten sigma Linie: {:3f} T'.format(B_rs))
dl_rs = 0.5*(rs_b/rs)*Lr
################################################################################

# Blaue sigma Linie
bs, bs_b = np.loadtxt('data/bs.txt', unpack = True)

B_bs = f(5, m, b)*10**(-3) #T
print('Magnetfeldstärke der blauen sigma Linie: {:3f} T'.format(B_bs))
dl_bs = 0.5*(bs_b/bs)*Lb
################################################################################

# BLaue pi Linie
bp, bp_b = np.loadtxt('data/bp.txt', unpack = True)

B_bp = f(15, m, b)*10**(-3) #T
print('Magnetfeldstärke der blauen pi Linie: {:3f} T'.format(B_bp))
dl_bp = 0.5*(bp_b/bp)*Lb
################################################################################

np.savetxt(fname='dl_rs.txt', X=dl_rs, header='.', footer='', comments='# Berechnete dlambdas fuer rot sigma')
np.savetxt(fname='dl_bs.txt', X=dl_bs, header='.', footer='', comments='# Berechnete dlambdas fuer blau sigma')
np.savetxt(fname='dl_bp.txt', X=dl_bp, header='.', footer='', comments='# Berechnete dlambdas fuer blau pi')

g_rs = h*c*dl_rs/(lr**2*m_B*B_rs)
g_bs = h*c*dl_bs/(lb**2*m_B*B_bs)
g_bp = h*c*dl_bp/(lb**2*m_B*B_bp)

np.savetxt(fname='g_rs.txt', X=g_rs, header='.', footer='', comments='Berechnete Lande-Faktoren fuer rot sigma')
np.savetxt(fname='g_bs.txt', X=g_bs, header='.', footer='', comments='Berechnete Lande-Faktoren fuer blau sigma')
np.savetxt(fname='g_bp.txt', X=g_bp, header='.', footer='', comments='Berechnete Lande-Faktoren fuer blau pi')

g_rs_mean = np.mean(g_rs)
g_rs_var = np.sqrt(np.var(g_rs))
g_bs_mean = np.mean(g_bs)
g_bs_var = np.sqrt(np.var(g_bs))
g_bp_mean = np.mean(g_bp)
g_bp_var = np.sqrt(np.var(g_bp))

print('Mittelwert für den Lande-Faktor für die rote sigma Linie = {:4f}+-{:6f}'.format(g_rs_mean, g_rs_var))
print('Mittelwert für den Lande-Faktor für die blaue sigma Linie = {:4f}+-{:6f}'.format(g_bs_mean, g_bs_var))
print('Mittelwert für den Lande-Faktor für die blaue pi Linie = {:4f}+-{:6f}'.format(g_bp_mean, g_bp_var))
