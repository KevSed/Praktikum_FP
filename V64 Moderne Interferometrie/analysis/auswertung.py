import numpy as np
import matplotlib.pyplot as plt

winkel, U_min, U_max = np.genfromtxt('kontrast.txt', unpack='True')

######################################################################
# Kontrastberechnung aus Messwerten der Spannungen für verschiedene
# Winkel des Polarisationsfilters
######################################################################

def kontrast(max, min):
    return (max-min)/(max+min)

plt.plot(winkel, kontrast(U_max, U_min), 'rx', label='Messwerte')
plt.xlabel(r'Winkel $\phi$ in °')
plt.ylabel('Kontrast')
plt.title('Winkelverteilung des Kontrastes')
plt.grid()
plt.legend(loc='best')
plt.savefig('kontrast.pdf')
