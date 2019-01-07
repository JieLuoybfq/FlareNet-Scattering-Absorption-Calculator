import logging
from math import exp
from math import pi
from math import log
from math import sin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path
import os

####### Plotting Parameters
rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'STIXGeneral'


##############

def Primary_Particle_Size_meter(da_meter, dp100_meter, Dtem):
    A = dp100_meter * ((da_meter * 1e9 / 100) ** Dtem)
    return A


def Figure_Plotter_Saver_1Lines_X_Log_Y_Linear(Address, X_Array, X_Label, Y_array, Y_Legend, Y_label1, Plot_Title, label_font_size=12, Plot_Title_Size=12, Figure_DPI=1200, alpha=0.3, Marker_Size=3):
    try:
        fig, ax1 = plt.subplots()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax1.plot(X_Array, Y_array, 'k-', label=Y_Legend, alpha=0.8, markersize=Marker_Size)
        ax1.set_xlabel(X_Label, fontsize=label_font_size)
        ax1.set_xscale("log")
        ax1.set_ylabel(Y_label1, fontsize=label_font_size)
        ax1.grid(True, which='major', axis="both", alpha=0.5)
        ax1.legend(bbox_to_anchor=(0.5, 0.67), loc='center left', fontsize='large')
        plt.title(Plot_Title, fontsize=Plot_Title_Size)
        plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
        plt.clf()
        plt.close()

    except Exception as e:
        logging.exception(e)
        raise


def Number_Primary_Particle(da_meter, dp_meter, Soot_Prefactor_k, Soot_Fractal_D, Continuum=True):
    if Continuum == True:
        A = ((da_meter / dp_meter) / Soot_Prefactor_k) ** (1 / Soot_Fractal_D)
        return round(A)
    else:
        A = ((da_meter / dp_meter) / Soot_Prefactor_k) ** (1 / Soot_Fractal_D)
        return round(A)


def Primary_LogNormal_Generator(dp_Median_meter, Sigma_g_Primary, Sigma_g_Number, Number_Points, da_meter, Soot_Prefactor_k_mc, Soot_Fractal_D_mc):

    Dmax = dp_Median_meter * (Sigma_g_Primary ** Sigma_g_Number)
    Dmin = dp_Median_meter * (Sigma_g_Primary ** (-1 * Sigma_g_Number))
    Diam_total = np.linspace(Dmin, Dmax, num=Number_Points)
    Number_Primary = []
    Probability = []
    Sum = 0

    for i in range(len(Diam_total)-1):
        Number_Primary.append(Number_Primary_Particle(da_meter, Diam_total[i], Soot_Prefactor_k_mc, Soot_Fractal_D_mc))
    for i in range(len(Diam_total)-1):
        Probability.append(LogNormal_Distribution(dp_Median_meter, Sigma_g_Primary, Diam_total[i + 1], Diam_total[i]))
        Sum += Probability[i]

    return Diam_total[:-1], Number_Primary, Probability


def LogNormal_Distribution(Median, SigmaG, Dp2, Dp1):  # return number between 0 to 1
    try:
        A = (1 / (log(SigmaG) * (2 * pi) ** 0.5)) * exp(-1 * (((log(Dp1) - log(Median)) ** 2) / (2 * (log(SigmaG)) ** 2))) * (log(Dp2) - log(Dp1))
        return A
    except Exception as e:
        logging.exception(e)
        raise


def RDG_Def_Scattering(K, theta, F, N, Dp, Df, Kf):
    xp = Dp * K / 2
    q = 2 * K * sin(theta / 2)
    Rg = Dp / 2 * (N / Kf) ** (1 / Df)
    if q * Rg <= 0.2:
        sq = 1
    elif q * Rg > 0.2 and q * Rg <= 1.6:
        sq = exp(-1 * (q * Rg) ** 2 / 3)
    elif q * Rg > 1.6:
        sq = (q * Rg) ** (-Df)
    C = (N ** 2) * (F * xp ** 6) / (K ** 2) * sq
    return C


def File_Pointer(Main, FolderName, FileName, Extension):
    try:
        Path1 = Path(Main) / Path(FolderName)
        if not os.path.exists(Path1):
            os.makedirs(Path1)
        File_Address = Path1 / Path(FileName + "." + Extension)
        return File_Address
    except Exception as e:
        logging.exception(e)
        raise


def RDG_Absorption(K, N, Dp, E):
    A = -1 * np.asarray(N) * 4 * pi * K * E * (np.asarray(Dp) / 2) ** 3
    return A, A / N


def RDG_Total_Scatteing(K, F, N, Dp, Betha):
    Df = 1.8
    Rg = np.asarray(Dp) / 2 / Betha
    Sigma_m = 8 / 3 * pi * (K ** 4) * ((np.asarray(Dp) / 2) ** 6) * F
    G = (1 + ((4 / (3 * Df)) * (K * Rg) ** 2)) ** (-0.5 * Df)
    Sigma_agg = (np.asarray(N) ** 2) * Sigma_m * G
    return Sigma_agg, Sigma_agg / N


def Effective_Density(K, Dm, dp):
    A = K * dp ** (Dm - 3)
    return A


def Particle_Mass(rho, da):
    A = rho * pi * (da ** 3) / 6
    return A


def Absorption_Eff(Abs_Cross, Abs_Cross_Prim, da, dp):
    a = Abs_Cross / ((pi / 4) * np.asarray(da) ** 2)
    b = Abs_Cross_Prim / ((pi / 4) * np.asarray(dp) ** 2)
    return a, b
