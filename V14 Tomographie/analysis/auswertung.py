import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

np.set_printoptions(precision=4)

# Dateneinlesen
counts = np.genfromtxt('../data/Leermessung_5.Spe', skip_header=12,
                       skip_footer=16)
counts_0 = np.genfromtxt('../data/Nullmessung.Spe', skip_header=12,
                         skip_footer=16)
I_l, c_l, t_l = np.genfromtxt('../data/Leermessung.txt', unpack=True)
I_2, c_2, t_2 = np.genfromtxt('../data/würfel2.txt', unpack=True)
I_3, c_3, t_3 = np.genfromtxt('../data/würfel3.txt', unpack=True)
I_5, c_5, t_5 = np.genfromtxt('../data/würfel5.txt', unpack=True)

n = {'leer': counts, 'error_leer': np.sqrt(counts), 'null': counts_0,
     'error_null': np.sqrt(counts_0)}
df = pd.DataFrame(n)

# Geometriematrix
b = np.sqrt(2)

# A wird ja mit den Absorptionskoeffizienten multipliziert, um dann die Raten
# auszugeben. Demnach entsprechen die Komponenten der einzelnen Zeilen
# den Strecken durch die Unter-Würfel, die bei der jeweiligen Projektion
# durchstrahlt werden. Daher lassen sich die Zeilen bei uniformer Verteilung
# auch aufsummieren (wie bei Würfel 2 und 3).
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

# Geoimetriematrix für Würfel 2 und 3 (vollständig aus einem Material)
A_2 = np.matrix([[3],
                [3],
                [3*b],
                [3*b]])


def kleinsteQuadrate(y, A, W):
    temp = np.dot(np.linalg.inv(np.dot(A.T, np.dot(W, A))), A.T)
    a = np.dot(temp, np.dot(W, y))
    a_err = np.linalg.inv(np.dot(A.T, np.dot(W, A)))
    return a, np.sqrt(np.diag(a_err))


# Leermessung
plt.bar(range(0, 511), df['leer'], yerr=df['error_leer'])
plt.xlim(0, 120)
plt.title('Messung des leeren Würfels')
plt.xlabel('Kanal')
plt.ylabel('Ereignisse')
plt.savefig('plots/leermessung.pdf')
plt.close()

print('Spektrum der Leermessung geplottet...')

# Nullmessung

plt.bar(range(0, 511), df['null'], yerr=df['error_null'])
plt.xlim(0, 120)
plt.title('Messung ohne Würfel')
plt.xlabel('Kanal')
plt.ylabel('Ereignisse')
plt.savefig('plots/nullmessung.pdf')
plt.close()

print('(Zusätzlich: Spektrum der Nullmessung geplottet...)')

c_no = 101383
t_no = 449.80  # s
rate_no = c_no / t_no

print('Messung ohne Würfel mit Rate: {:2f} +- {:2f} counts/s'
      .format(rate_no, np.sqrt(c_no)/t_no))

# I_0s für die Projektionen 5, 6, 11, 8, 9
rate_l = c_l / t_l
err_rate_l = np.sqrt(c_l) / t_l

# Würfel 2 (5, 2, 8, 11)
rate_2 = c_2 / t_2
err_rate_2 = np.sqrt(c_2) / t_2


# Würfel 3 (5, 2, 8, 11)
rate_3 = c_3 / t_3
err_rate_3 = np.sqrt(c_3) / t_3


# Würfel 5 (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
rate_5 = c_5 / t_5
err_rate_5 = np.sqrt(c_5) / t_5


# Mitteln für leeren Würfel
I_leer = np.array(np.zeros(len(rate_5)))
err_I_leer = np.array(np.zeros(len(rate_5)))

for i in range(len(rate_5)):
    if i < 6:
        I_leer[i] = (rate_l[0] + rate_l[1]) / 2
        err_I_leer[i] = (err_rate_l[0] + err_rate_l[1]) / 2
    if (i == 6) or (i == 8) or (i == 9) or (i == 11):
        I_leer[i] = rate_l[4]
        err_I_leer[i] = err_rate_l[4]
    if (i == 7) or (i == 10):
        I_leer[i] = (rate_l[2] + rate_l[3]) / 2
        err_I_leer[i] = (err_rate_l[2] + err_rate_l[3]) / 2

# Umrechnen in ln(I_0/N_j)
# Proj_1 steht für dir Projektionen, die bei Würfel 2 und 3 verwendet wurden
I_Proj_1 = np.array([I_leer[1], I_leer[4], I_leer[7], I_leer[10]])
err_I_Proj_1 = np.array([err_I_leer[1], err_I_leer[4], err_I_leer[7], err_I_leer[10]])

I_2 = np.log(c_2 / I_Proj_1)
I_3 = np.log(c_3 / I_Proj_1)
I_5 = np.log(c_5 / I_leer)

err_I_2 = np.sqrt((np.sqrt(c_2) / c_2)**2 + (err_I_Proj_1 / I_Proj_1)**2)
err_I_3 = np.sqrt((np.sqrt(c_3) / c_3)**2 + (err_I_Proj_1 / I_Proj_1)**2)
err_I_5 = np.sqrt((np.sqrt(c_5) / c_5)**2 + (err_I_leer / I_leer)**2)


print('''
###############################################################

      ~~~ Raten der verschiedenen Würfel ~~~

      Leermessung: (Projektionen 5, 6, 11, 8, 9)
      -----------------------------------------------
      Werte = {}
      Fehler = {}
      
      Leermessung mit Mittelung: (alle 12 Projektionen)
      -----------------------------------------------
      Werte = {}
      Fehler = {}
      
      Würfel 2: (Projektionen 5, 2, 8, 11)
      -----------------------------------------------
      Werte = {}
      Fehler = {}
      
      Würfel 3: (Projektionen 5, 2, 8, 11)
      -----------------------------------------------
      Werte = {}
      Fehler = {}
      
      Würfel 5: (alle 12 Projektionen)
      -----------------------------------------------
      Werte = {}
      Fehler = {}

###############################################################
      '''
      .format(rate_l, err_rate_l, I_leer, err_I_leer, rate_2, err_rate_2,
              rate_3, err_rate_3, rate_5, err_rate_5))

print('''
      ~~~ Logarithmen der Raten der verschiedenen Würfel ~~~

      Würfel 2:
      -----------------------------------------------
      Werte = {}
      Fehler = {}
      
      Würfel 3:
      -----------------------------------------------
      Werte = {}
      Fehler = {}
      
      Würfel 5:
      -----------------------------------------------
      Werte = {}
      Fehler = {}
 
###############################################################
'''.format(I_2, err_I_2, I_3, err_I_3, I_5, err_I_5))

# Gewichtungsmatrizen
W_2 = np.diag(1/err_I_2**2)
W_3 = np.diag(1/err_I_3**2)
W_5 = np.diag(1/err_I_5**2)

# Über kleinsteQuadrate mu berechnen
mu_2, err_mu_2 = kleinsteQuadrate(I_2, W=W_2, A=A_2)
mu_3, err_mu_3 = kleinsteQuadrate(I_3, W=W_3, A=A_2)
mu_5, err_mu_5 = kleinsteQuadrate(I_5, W=W_5, A=A)

print(
    '''
      ~~~ Absorptionskoeffizienten der verschiedenen Würfel ~~~

      Würfel 2:
      -----------------------------------------------
      Werte = {}
      Fehler = {}
      
      Würfel 3:
      -----------------------------------------------
      Werte = {}
      Fehler = {}
      
      Würfel 5:
      -----------------------------------------------
      Werte = {}
      Fehler = {}
    '''.format(mu_2, err_mu_2, mu_3, err_mu_3, mu_5, err_mu_5)
)

