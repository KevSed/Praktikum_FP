import numpy as np
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)

@np.vectorize
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

@np.vectorize
def correctSelfNoise(U, VN):
    if (VN ==    1): return U - 0.0031
    if (VN ==    2): return U - 0.0030
    if (VN ==    5): return U - 0.0031
    if (VN ==   10): return U - 0.0031
    if (VN ==   20): return U - 0.0030
    if (VN ==   50): return U - 0.0030
    if (VN ==  100): return U - 0.0030
    if (VN ==  200): return U - 0.0039
    if (VN ==  500): return U - 0.0116
    if (VN == 1000): return U - 0.0362

def normVoltage(U, VS, VN, VZ = 1):
    return U * (VZ / (VS * VN)) ** 2

def linearFunction(x, m, b):
    return m * x + b

def correlatorCalibrationCurve(f, A, Q, f0):
    return A /Q ** 2 * 1 / ((f / f0) ** 2 + (f0 / f) ** 2 + 1 / Q ** 2 - 2)
