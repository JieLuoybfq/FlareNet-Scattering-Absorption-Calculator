from math import exp
from math import pi
from math import log
from math import sin
import numpy as np


def LogNormal (D_median, Sigma, Dpi, Dpi1):  # return number between 0 to 1
    A = (1 / (log (Sigma) * (2 * pi) ** 0.5)) * exp (-1 * (((log (Dpi)-log (D_median)) ** 2) / (2 * (log (Sigma)) ** 2))) * (log (Dpi1)-log (Dpi))
    return A


def RDG_Def_Scattering (K, theta, F, N, Dp, Df, Kf):
    xp = Dp * K / 2
    q = 2 * K * sin (theta / 2)
    Rg = Dp / 2 * (N / Kf) ** (1 / Df)
    if q * Rg <= 0.2:
        sq = 1
    elif q * Rg > 0.2 and q * Rg <= 1.6:
        sq = exp (-1 * (q * Rg) ** 2 / 3)
    elif q * Rg > 1.6:
        sq = (q * Rg) ** (-Df)
    C = (N ** 2) * (F * xp ** 6) / (K ** 2) * sq
    return C


def RDG_Absorption (K, N, Dp, E):
    A = -1 * np.asarray (N) * 4 * pi * K * E * (np.asarray (Dp) / 2) ** 3
    return A, A / N


def RDG_Total_Scatteing (K, F, N, Dp, Betha):
    Df=1.8
    Rg = np.asarray (Dp) / 2 /Betha
    Sigma_m = 8 / 3 * pi * (K ** 4) * ((np.asarray (Dp) / 2) ** 6) * F
    G = (1+((4 / (3 * Df)) * (K * Rg) ** 2)) ** (-0.5 * Df)
    Sigma_agg = (np.asarray (N) ** 2) * Sigma_m * G
    return Sigma_agg, Sigma_agg / N


def Effective_Density (K, Dm, dp):
    A = K * dp ** (Dm-3)
    return A


def Particle_Mass (rho, da):
    A = rho * pi * (da ** 3) / 6
    return A


def Absorption_Eff (Abs_Cross, Abs_Cross_Prim, da, dp):
    a = Abs_Cross / ((pi / 4) * np.asarray (da) ** 2)
    b = Abs_Cross_Prim / ((pi / 4) * np.asarray (dp) ** 2)
    return a, b
