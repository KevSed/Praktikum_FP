import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

t  = [1,2,3,4,5,6,7,8,9] # µs

x  = []
x += [36164*[ 22]]                         # 1µs
x += [42143*[ 44]]                         # 2µs
x += [46389*[ 66]]                         # 3µs
x += [ 4017*[ 87]+25749*[ 88]+ 2313*[ 89]] # 4µs
x += [ 6846*[109]+17056*[110]+ 8042*[111]] # 5µs
x += [ 8245*[131]+13694*[132]+10312*[133]] # 6µs
x += [ 9585*[153]+12647*[154]+12760*[155]] # 7µs
x += [ 9991*[175]+10732*[176]+13436*[177]] # 8µs
x += [ 8855*[197]+ 8293*[198]+12519*[199]] # 9µs

ch = []
for i in range(len(x)):
  ch.append(ufloat(np.mean(x[i]), np.std(x[i])))

def f(x, m, b):
  return m*x+b

params, cov = curve_fit(f, noms(ch), t)
errors = np.sqrt(np.diag(cov))
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
print('Fitparameter: \n m = {} \n b = {}'.format(m, b))

x = np.linspace(0,220,2)
plt.plot(x, f(x, *params), 'b-', label='Lineare Regression')
plt.plot(noms(ch[0:3]), t[0:3], 'rx', ms=13, mew=1)
plt.errorbar(noms(ch[3::]), t[3::], xerr=stds(ch[3::]), fmt='rx', ms=13, mew=1, elinewidth=1, capsize=5, capthick=1, label='Messwerte')
plt.xlabel(r'Kanal')
plt.ylabel(r'$\displaystyle\frac{\Delta t}{\mathrm{\mu s}}$')
plt.xlim(0,220)
plt.grid()
plt.legend(loc='best', numpoints=1)
plt.tight_layout()
plt.savefig('analysis/kalibrierung.pdf')

np.savetxt('analysis/kalibrierung.txt', np.array([noms(m), noms(b), stds(m), stds(b)]), delimiter=' ')
