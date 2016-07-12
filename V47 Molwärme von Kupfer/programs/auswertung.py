import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

'''
T, alpha = np.loadtxt('data_alpha.txt', unpack = True)

def f(x, m, b):
    return m*1/x+b

params, cov = curve_fit(f, T, alpha)

m = params[0]
b = params[1]

Δm = np.sqrt(cov[0][0])
Δb = np.sqrt(cov[1][1])

print(m)
print(Δm)
print(b)
print(Δb)

x = np.linspace(70, 300, 1000)

plt.plot(1/T, alpha, 'rx', label='Messwerte')
plt.plot(1/x, f(x, m, b), 'k-', label='Regressionsgerade')
plt.grid()
plt.xlabel(r"$T^{-1}$ in $\mathrm{K}^{-1}$")
plt.ylabel(r"$\alpha\cdot 10^{-6}$ in $\mathrm{K}^{-1}$")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("../build/plot_alpha.pdf")
plt.close()
'''


R,U,I,t,T,dT,alpha,dalpha,diffT,ddiffT,Cp,dCp,Cv,dCv=np.loadtxt('data_C.txt',unpack=True)
plt.errorbar(T,Cp,yerr=dCp,fmt='bx',label='Werte für $C_{\mathrm{P}}$')
plt.errorbar(T,Cv,yerr=dCv,fmt='rx',label='Werte für $C_{\mathrm{V}}$')
plt.legend(loc='best')
plt.tight_layout()
plt.grid()
plt.xlabel(r"$T$ in $\mathrm{K}$")
plt.ylabel(r"$C_{\mathrm{P,V}}$ in $\mathrm{JK}^{-1}\mathrm{mol}^{-1}$")
plt.tight_layout()
plt.savefig("../build/plot_Cv.pdf")
plt.close()
