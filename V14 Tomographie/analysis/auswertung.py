import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from uncertainties import ufloat

# Leermessung
hist_1 = np.genfromtxt('Leermessung_5.Spe', skip_header=12, skip_footer=16)

plt.bar(range(0, 511), hist_1)
plt.xlim(0, 120)
plt.title('Messung des leeren Würfels')
plt.xlabel('Kanal')
plt.ylabel('Ereignisse')
plt.savefig('leermessung.pdf')

# Nullmessung
I_0 = ufloat(101383 / 449.80, np.sqrt(101383) / 449.80)

hist_0 = np.genfromtxt('Nullmessung.Spe', skip_header=12, skip_footer=16)

plt.bar(range(0, 511), hist_0)
plt.xlim(0, 120)
plt.title('Messung ohne Würfel')
plt.xlabel('Kanal')
plt.ylabel('Ereignisse')
plt.savefig('nullmessung.pdf')

# Geometriematrix
b = np.sqrt(2)

A = np.matrix([[1, 0, 0, 1, 0, 0, 1, 0, 0],
               [0, 1, 0, 0, 1, 0, 0, 1, 0],
               [0, 0, 1, 0, 0, 1, 0, 0, 1],
               [1, 1, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 1, 1, 1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 1, 1, 1],
               [0, b, 0, 0, 0, b, 0, 0, 0],
               [b, 0, 0, 0, b, 0, 0, 0, b],
               [0, 0, 0, b, 0, 0, 0, b, 0],
               [0, 0, 0, 0, 0, b, 0, b, 0],
               [0, 0, b, 0, b, 0, b, 0, 0],
               [0, b, 0, b, 0, 0, 0, 0, 0]])

B = A.T*A

#
W = np.diag([1, 1, 1, 1, 1])
