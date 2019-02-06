import logging
from math import exp
from math import pi
from math import log
from math import sin
from math import cos
import numpy as np
from matplotlib import rcParams
from pathlib import Path
import os
import matplotlib.pyplot as plt
import csv
from mpl_toolkits.mplot3d import Axes3D

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


def Dictionary_ToCSV(address, dictionary):
    try:

        with open(address, 'w') as csv_file:
            writer = csv.writer(csv_file)
            for key, value in dictionary.items():
                if isinstance(value, list) == True:
                    value.insert(0, key)

                    # for item in Total:
                    writer.writerow(value)
                    # for i in range(len(value)):
                    #     writer.write([value[i]])
                    # writer.write("\n")

                else:
                    writer.writerow([key, value])

    except Exception as e:
        logging.exception(e)
        raise


def Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address, Identifier, X_Array, Y_array, X_Min=None, X_Max=None, Y_Min=None, Y_Max=None, X_Label=None, Y_Legend=None, Y_label1=None, Plot_Title=None, label_font_size=12, Plot_Title_Size=12, Figure_DPI=1200, alpha_Y=0.7, Marker_Size=5):
    try:

        fig, ax1 = plt.subplots()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        Number_Situation = len(Identifier)

        Linestyles = ['-', '--', '-.', ':']
        Count_Linestyles = len(Linestyles) - 1
        Markerstyles = [',', '+', '.', 'o', '*']
        Count_Markerstyles = len(Markerstyles) - 1
        Marker_Size_Step = Marker_Size / (Number_Situation + 8)
        alpha_Y_Step = alpha_Y / (Number_Situation + 8)
        i1 = 0
        i2 = 0
        for i in range(Number_Situation):
            ax1.plot(X_Array[Identifier[i]], Y_array[Identifier[i]], label=Identifier[i], alpha=alpha_Y, marker=Markerstyles[i2], linestyle=Linestyles[i1], markersize=Marker_Size)
            alpha_Y -= alpha_Y_Step
            Marker_Size -= Marker_Size_Step
            i1 += 1
            if i1 == Count_Linestyles:
                i1 = 0
                i2 += 1
                if i2 == Count_Markerstyles:
                    i2 = 0

        if X_Label != None:
            ax1.set_xlabel(X_Label, fontsize=label_font_size)
        if X_Min != None and X_Max != None:
            ax1.set_xlim(X_Min, X_Max)
        if Y_Min != None and Y_Max != None:
            ax1.set_ylim(Y_Min, Y_Max)
        ax1.set_xscale("log")
        if Y_label1 != None:
            ax1.set_ylabel(Y_label1, fontsize=label_font_size)
        ax1.grid(True, which='major', axis="both", alpha=0.5)
        ax1.legend(bbox_to_anchor=(1.05, 0.70), loc='center left', fontsize='medium')
        if Plot_Title != None:
            plt.title(Plot_Title, fontsize=Plot_Title_Size)
        plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
        plt.clf()
        plt.close()

    except Exception as e:
        logging.exception(e)
        raise


def Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, X_Array, Y_array, X_Min=None, X_Max=None, Y_Min=None, Y_Max=None, X_Label=None, Y_Legend=None, Y_label1=None, Plot_Title=None, label_font_size=12, Plot_Title_Size=12, Figure_DPI=1200, alpha_Y=0.9, Marker_Size=3):
    try:

        fig, ax1 = plt.subplots()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax1.plot(X_Array, Y_array, 'r-', label=Y_Legend, alpha=alpha_Y, markersize=Marker_Size)
        if X_Label != None:
            ax1.set_xlabel(X_Label, fontsize=label_font_size)
        if X_Min != None and X_Max != None:
            ax1.set_xlim(X_Min, X_Max)
        if Y_Min != None and Y_Max != None:
            ax1.set_ylim(Y_Min, Y_Max)
        ax1.set_xscale("log")
        if Y_label1 != None:
            ax1.set_ylabel(Y_label1, fontsize=label_font_size)
        ax1.grid(True, which='major', axis="both", alpha=0.5)
        ax1.legend(bbox_to_anchor=(0.7, 0.9), loc='center left', fontsize='large')
        if Plot_Title != None:
            plt.title(Plot_Title, fontsize=Plot_Title_Size, y=1.08)
        plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
        plt.clf()
        plt.close()

    except Exception as e:
        logging.exception(e)
        raise


def Fig_Plot_Save_1Lines_X_Log_Y_Log(Address, X_Array, X_Label, Y_array, Ybottom, YTop, Y_Legend, Y_label1, Plot_Title, label_font_size=12, Plot_Title_Size=12, Figure_DPI=1200, alpha=0.3, Marker_Size=3):
    try:

        fig, ax1 = plt.subplots()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax1.plot(X_Array, Y_array, 'r-', label=Y_Legend, alpha=0.8, markersize=Marker_Size)
        ax1.set_xlabel(X_Label, fontsize=label_font_size)
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_ylim(Ybottom, YTop)
        ax1.set_ylabel(Y_label1, fontsize=label_font_size)
        ax1.grid(True, which='major', axis="both", alpha=0.5)
        ax1.legend(bbox_to_anchor=(0.6, 0.95), loc='center left', fontsize='large')
        plt.title(Plot_Title, fontsize=Plot_Title_Size)
        plt.savefig(Address, format='jpg', dpi=Figure_DPI, bbox_inches='tight')
        plt.clf()
        plt.close()

    except Exception as e:
        logging.exception(e)
        raise


def Fig_Plot_3D_Show_XCte(X, Y_2D, Z_2D):
    try:

        fig, ax = plt.subplots()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for i in range(len(X)):
            dummy = []
            for p in range(len(Y_2D[i])):
                dummy.append(X[i])
            ax.plot(dummy, Y_2D[i][:], Z_2D[i][:])
        plt.show()
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


def Diff_Integral_Phi(Scattered_Theta, Phi_Radian, Theta_Rad, Theta_Diff, Phi_Diff, Formula=1):
    try:

        A = 0
        for phi in range(len(Phi_Radian)):
            A += Scattered_Theta * sin(Theta_Rad) * Theta_Diff * Phi_Diff * (1 - (((cos(Phi_Radian[phi])) * sin(Theta_Rad)) ** 2))
        return A

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
            Diameter_Log.append(log(d1 * Keyhan))
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


def RDG_Def_Scattering(K, N, Dp, q, F, D_RDG, K_RDG, Formula=1):
    try:

        Monomer_Differential = F * ((Dp / 2) ** 6) * (K ** 4)
        Rg = ((N / K_RDG) ** (1 / D_RDG)) * (Dp / 2)
        S_q = RDG_Structure_Factor(q=q, Rg=Rg, D_RDG=D_RDG, K_RDG=K_RDG, Formula=Formula)
        Aggregate_Differential = (N ** 2) * Monomer_Differential * S_q

        logging.debug(f"RDG_Def_Scattering,Monomer_Differential={Monomer_Differential},Rg={Rg},S_q={S_q},Aggregate_Differential={Aggregate_Differential}")
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


def RDG_Structure_Factor(q, Rg, D_RDG, K_RDG, Formula=1, C=1):
    try:

        if Formula == 1:  # Sorensen 2001
            C = 1.35 / K_RDG
            Check = q * Rg
            if Check <= 1:
                A = 1 - ((Check ** 2) / 3)  # Guinier equation (Guinier 1939;Guinier et al. 1955; Teixeira 1986).
            elif Check > 1:
                A = C * (Check) ** (-1 * D_RDG)
        elif Formula == 2:  # Yang 2005
            A = (1 + ((8 * (q * Rg) ** 2) / (3 * D_RDG)) + (q * Rg) ** 8) ** (-1 * D_RDG / 8)

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


def RDG_Total_Scattering(K, N, Dp, F, D_RDG, K_RDG, Formula=1):
    try:

        if Formula == 1:  # Sorensen 2001
            Monomer_Total = (8 / 3) * pi * (K ** 4) * ((Dp / 2) ** 6) * F
            Rg = ((N / K_RDG) ** (1 / D_RDG)) * (Dp / 2)
            G = (1 + ((4 / (3 * D_RDG)) * ((K * Rg) ** 2))) ** (-0.5 * D_RDG)

            Aggregate_Total = (N ** 2) * Monomer_Total * G
            logging.debug(f"RDG_Total_Scattering,Monomer_Total={Monomer_Total},Rg={Rg},G={G},Aggregate_Total={Aggregate_Total}_Sorensen(2001)")
        if Formula == 2:  # Yang 2005
            Monomer_Total = (8 / 3) * pi * (K ** 4) * ((Dp / 2) ** 6) * F
            Rg = ((N / K_RDG) ** (1 / D_RDG)) * (Dp / 2)
            Betha = 3 * D_RDG / (8 * (K * Rg) ** 2)
            if Betha >= 1:
                G = 1 - (2 * ((K * Rg) ** 2) / 3)
            elif Betha < 1:
                G = ((Betha / 2) * (3 - 3 * Betha + 2 * (Betha ** 2))) - ((((K * Rg * Betha) ** 2) / 3) * (3 - 4 * Betha + 3 * (Betha ** 2))) + (
                        ((2 * K * Rg) ** (-1 * D_RDG)) * ((3 / (2 - D_RDG)) - (12 / ((6 - D_RDG) * (4 - D_RDG))) - (3 * (Betha ** (1 - D_RDG / 2)) * ((1 / (2 - D_RDG)) - (2 * Betha / (4 - D_RDG)) + (2 * (Betha ** 2) / (6 - D_RDG))))))

            Aggregate_Total = (N ** 2) * Monomer_Total * G
            logging.debug(f"RDG_Total_Scattering,Monomer_Total={Monomer_Total},Rg={Rg}, Betha={Betha},G={G},Aggregate_Total={Aggregate_Total}_Yang(2005)")
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


def Effective_Density_K_FromRho100nm(Eff_dm, Eff_rho_100nm):
    try:
        da = 100 * (10 ** -9)

        K = Eff_rho_100nm / (da ** (Eff_dm - 3))
        return K

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
