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
    try:

        A = dp100_meter * ((da_meter * 1e9 / 100) ** Dtem)
        return A

    except Exception as e:
        logging.exception(e)
        raise


def Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, X_Array, X_Label, Y_array, Y_Legend, Y_label1, Plot_Title, label_font_size=12, Plot_Title_Size=12, Figure_DPI=1200, alpha=0.3, Marker_Size=3):
    try:

        fig, ax1 = plt.subplots()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax1.plot(X_Array, Y_array, 'r*', label=Y_Legend, alpha=0.8, markersize=Marker_Size)
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


def Number_Primary_Particle(da_meter, dp_meter, Soot_Prefactor_k, Soot_Fractal_D):
    try:

        A = ((da_meter / dp_meter) / Soot_Prefactor_k) ** (1 / Soot_Fractal_D)
        return round(A)

    except Exception as e:
        logging.exception(e)
        raise


def Bins_LogN_Distributed(Median_Diameter, Sigma_G, Sigma_G_Bound, Total_Number_Bins):
    try:

        Bound_D_Max = Median_Diameter * (Sigma_G ** Sigma_G_Bound)
        Bound_D_Min = Median_Diameter * (Sigma_G ** (-1 * Sigma_G_Bound))
        Diameter_Meter = []
        Diameter_Log = []
        Diameter_Nano = []
        Keyhan = 10 ** 9

        D_Ratio = (Bound_D_Max / Bound_D_Min) ** (1 / (Total_Number_Bins - 1))

        for i in range(0, Total_Number_Bins):
            d1 = Bound_D_Min * (D_Ratio ** i)
            Diameter_Meter.append(d1)
            Diameter_Log.append(log(d1))
            Diameter_Nano.append(d1 * Keyhan)

        return sorted(Diameter_Meter, key=float), sorted(Diameter_Log, key=float), sorted(Diameter_Nano, key=float)

    except Exception as e:
        logging.exception(e)
        raise


def Primary_LogN_Generator(dp_Median_meter, Sigma_g_Primary, Sigma_g_Number, Number_Points, da_meter, Soot_Prefactor_k, Soot_Fractal_D):
    try:
        if Sigma_g_Primary != 1:
            Dmax = dp_Median_meter * (Sigma_g_Primary ** Sigma_g_Number)
            Dmin = dp_Median_meter * (Sigma_g_Primary ** (-1 * Sigma_g_Number))
            Diam_total = np.linspace(Dmin, Dmax, num=Number_Points)
        elif Sigma_g_Primary == 1:
            Dmax = dp_Median_meter * (1.5 ** 2)
            Dmin = dp_Median_meter * (1.5 ** (-1 * 2))
            Diam_total = np.linspace(Dmin, Dmax, num=Number_Points)
            k = 3
        Number_Primary = []
        Probability = []
        Sum = 0

        for i in range(len(Diam_total) - 1):
            Number_Primary.append(Number_Primary_Particle(da_meter, Diam_total[i], Soot_Prefactor_k, Soot_Fractal_D))
        for i in range(len(Diam_total) - 1):
            Probability.append(LogN_Distribution(dp_Median_meter, Sigma_g_Primary, Diam_total[i + 1], Diam_total[i]))
            Sum += Probability[i]

        return Diam_total[:-1], Number_Primary, Probability

    except Exception as e:
        logging.exception(e)
        raise


def LogN_Distribution(Median, SigmaG, Dp2, Dp1):  # return number between 0 to 1
    try:
        if SigmaG != 1:
            A = (1 / (log(SigmaG) * (2 * pi) ** 0.5)) * exp(-1 * (((log(Dp1) - log(Median)) ** 2) / (2 * (log(SigmaG)) ** 2))) * (log(Dp2) - log(Dp1))
        elif SigmaG == 1:
            if Dp2 > Median and Median > Dp1:
                return 1
            else:
                return 0

        return A

    except Exception as e:
        logging.exception(e)
        raise


def RDG_Def_Scattering(K, N, Dp, q, F, D_RDG, K_RDG, check, Number):
    try:

        Monomer_Differential = F * ((Dp / 2) ** 6) * (K ** 4)
        Rg = ((N / K_RDG) ** (1 / D_RDG)) * (Dp / 2)
        S_q = RDG_Structure_Factor(q=q, Rg=Rg, D_FD=D_RDG)
        Aggregate_Differential = (N ** 2) * Monomer_Differential * S_q
        if check == 1:
            if Number == 194:
                logging.debug(f"RDG_Def_Scattering---Number={Number},Monomer_Differential={Monomer_Differential},Rg={Rg},S_q={S_q},Aggregate_Differential={Aggregate_Differential}")
        return Aggregate_Differential

    except Exception as e:
        logging.exception(e)
        raise


def Scattering_Wave_Vector(WaveLength_meter, Theta_radian):
    try:

        A = (4 * pi / WaveLength_meter) * sin(Theta_radian / 2)
        return A

    except Exception as e:
        logging.exception(e)
        raise


def RDG_Structure_Factor(q, Rg, D_FD, C=1):
    try:

        Check = q * Rg
        if Check <= 1:
            A = 1 - ((Check ** 2) / 3)  # Guinier equation (Guinier 1939;Guinier et al. 1955; Teixeira 1986).
        elif Check > 1:
            A = C * (Check) ** (-1 * D_FD)
        return A

    except Exception as e:
        logging.exception(e)
        raise


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
    try:

        A = -1 * N * 4 * pi * K * E * ((Dp / 2) ** 3)
        return A

    except Exception as e:
        logging.exception(e)
        raise


def RDG_Total_Scattering(K, N, Dp, F, D_RDG, K_RDG, check, Number):
    try:

        Monomer_Total = (8 / 3) * pi * (K ** 4) * ((Dp / 2) ** 6) * F
        Rg = ((N / K_RDG) ** (1 / D_RDG)) * (Dp / 2)
        G = (1 + ((4 / (3 * D_RDG)) * ((K * Rg) ** 2))) ** (-0.5 * D_RDG)
        Aggregate_Total = (N ** 2) * Monomer_Total * G
        if check == 1:
            logging.debug(f"RDG_Total_Scattering---Number={Number},Monomer_Total={Monomer_Total},Rg={Rg},G={G},Aggregate_Total={Aggregate_Total}")
        return Aggregate_Total

    except Exception as e:
        logging.exception(e)
        raise


def Effective_Density(K, Dm, da):
    try:

        A = K * da ** (Dm - 3)
        return A

    except Exception as e:
        logging.exception(e)
        raise


def Differential_Solid_Angle(Phi_radian, Theta_Diff, Phi_Diff):
    try:

        A = Theta_Diff * Phi_Diff * sin(Phi_radian)
        return A

    except Exception as e:
        logging.exception(e)
        raise


def Mass_Calc(rho, da):
    try:

        A = rho * pi * (da ** 3) / 6
        return A

    except Exception as e:
        logging.exception(e)
        raise


def Absorption_Eff(Abs_Cross, da):
    try:

        a = Abs_Cross / ((pi / 4) * (da ** 2))
        return a

    except Exception as e:
        logging.exception(e)
        raise
