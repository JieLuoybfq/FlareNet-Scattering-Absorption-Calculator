# Programmed by Keyhan Babaee Under Prof. Steven Rogak supervision, https://github.com/KeyhanB
# Version 0.4
# Jan 2, 2018
import Functions as FN
import os
import numpy as np
from math import pi
from math import log
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import logging

####### Logging Parameters
logging.basicConfig(format='%(asctime)s, %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', filename='Log.txt', level=logging.INFO, filemode='w')
##############

if __name__ == "__main__":

    logging.info("Program Started!")
    script_dir = os.path.dirname(os.path.realpath('__file__'))
    Graph_Folder = "Graph Output"

    ####################### Aggregate diameter bounds
    # Bound_D_Min = 50e-9  # smallest diameter (m)
    # Bound_D_Max = 10e-6  # largest diameter (m)
    # Bound_D_Bins = 199  # Number of bins

    ####################### Particle Sample Setting
    Sample_Total_Number_Concentration = 100  # /cm^3
    Sample_LogN_D_Median_Min = 200  # Smallest diameter of computation (nm)
    Sample_LogN_D_Median_Max = 800  # Largest diameter of computation (nm)
    Sample_LogN_D_Median_Bins = 2  # Number of the Steps
    Sample_LogN_Sigma_Min = 1.03  # Smallest Sigma G
    Sample_LogN_Sigma_Max = 2.5  # Largest Sigma G
    Sample_LogN_Sigma_Bins = 2  # Number of the Steps
    Sample_Sigma_Bound = 2.5  # Number of Sigma G to cover
    Sample_Sigma_Bins = 300  # Number of bins

    ####################### Varation of Primary Particle Size within the sample
    Primary_Diameter_100nm = 16.48 * 1e-9  # Primary particle diameter of 100 nm aggregate
    Primary_D_TEM = 0.45
    Primary_Sigma_da_CTE = 1.02  # Sigma G around specific particle primary particle size
    Primary_Sigma_da_CTE_Bound = 2.5  # Number of Sigma G to cover
    Primary_Sigma_da_CTE_Nt = 19 + 1  # Number of bins
    Primary_Betha = 0.9

    ####################### Refractive Index
    Soot_Refractive_Index = 1.95 - 0.79j  # Complex refractive index
    Soot_Primary_Agg_dp_sigma = 1  # Geometric std of primary particles within aggregate

    ####################### Continuum regime
    Soot_Fractal_D_mc = 0.52  # Fractal Properties for continuum regime
    Soot_Prefactor_k_mc = 0.85 + (0.03 * Soot_Primary_Agg_dp_sigma ** (4.4))  # Fractal continuum regime

    ####################### Free molecular regime
    Soot_Fractal_D_mfm = 0.46  # Fractal Properties free molecular regime
    Soot_Prefactor_k_mfm = 0.94 + (0.03 * Soot_Primary_Agg_dp_sigma ** (4.8))  # Fractal free molecular regime

    ####################### Effective Density
    Eff_k = 0.418  # effective density variable for the result to be in kg/m^3
    Eff_dm = 2.56  # effective density variable for the result to be in kg/m^3

    ####################### Wave Parameter
    Wave_Length = 550 * 1e-9  # nm
    Wave_Number = 2 * pi / Wave_Length

    ####################### Scattering Parameters
    Theta_Number = 48
    Theta_Radian = np.linspace(0, 2 * pi, num=Theta_Number)

    #######################
    Soot_Complex = (((Soot_Refractive_Index ** 2) - 1) / ((Soot_Refractive_Index ** 2) + 2))
    Soot_FM = (abs(Soot_Complex)) ** 2
    Soot_EM = Soot_Complex.imag
    Sample_Sigma_List = np.linspace(Sample_LogN_Sigma_Min, Sample_LogN_Sigma_Max, num=Sample_LogN_Sigma_Bins)
    Sample_D_Median_List = np.linspace(Sample_LogN_D_Median_Min, Sample_LogN_D_Median_Max, num=Sample_LogN_D_Median_Bins)
    #######################
    #
    # Diameter_Meter = []
    # Diameter_Log = []
    # Diameter_Nano = []
    # Keyhan = 10 ** 9  # Conversion factor
    # D_ratio = (Bound_D_Max / Bound_D_Min) ** (1 / (Bound_D_Bins - 1))  # Generating particle diameters
    # for i in range(Bound_D_Bins):
    #     d = Bound_D_Min * D_ratio ** (i)
    #     Diameter_Meter.append(d)
    #     Diameter_Log.append(log(d))
    #     Diameter_Nano.append(d * Keyhan)
    #######################

    # Total_Absorption_Cross = [[0 for j in range(LogNormal_Sigma_Number_Sample)] for i in range(LogNormal_D_Median_Number_Sample)]
    # Scattering_Cross = [[[0 for k in range(Nd)] for j in range(LogNormal_Sigma_Number_Sample)] for i in range(LogNormal_D_Median_Number_Sample)]
    # Scattering_Cross_Primary = [[[0 for k in range(Nd)] for j in range(LogNormal_Sigma_Number_Sample)] for i in range(LogNormal_D_Median_Number_Sample)]
    # Total_Scattering_Cross = [[0 for j in range(LogNormal_Sigma_Number_Sample)] for i in range(LogNormal_D_Median_Number_Sample)]
    # Differential_Scattering_Cross = [[[[0 for p in range(Theta_Number)] for k in range(Nd)] for j in range(LogNormal_Sigma_Number_Sample)] for i in range(LogNormal_D_Median_Number_Sample)]
    # Aggregate_Mass = [[[0 for k in range(Nd)] for j in range(LogNormal_Sigma_Number_Sample)] for i in range(LogNormal_D_Median_Number_Sample)]
    # Total_Particle_Mass = [[0 for j in range(LogNormal_Sigma_Number_Sample)] for i in range(LogNormal_D_Median_Number_Sample)]
    # Mass_Absorption_Eff = [[0 for j in range(LogNormal_Sigma_Number_Sample)] for i in range(LogNormal_D_Median_Number_Sample)]

    X_Label1 = 'Mobility-equivalent Diameter (nm)'
    Y_Label1 = 'dN/dLogDp, (#/cm' + "$^{}$".format(3) + ')'
    Y_Label2 = "Primary Particle Diameter-dp(m)"
    Y_Label3 = "Absorption Cross Section(m^2)"
    Y_Label4 = "Absorption Efficiency"

    for i in range(Sample_LogN_D_Median_Bins):
        for j in range(Sample_LogN_Sigma_Bins):

            Diameter_Meter, Diameter_Log, Diameter_Nano = FN.Bins_LogN_Distributed(Median_Diameter=Sample_D_Median_List[i] * 10 ** (-9), Sigma_G=Sample_Sigma_List[j], Sigma_G_Bound=Sample_Sigma_Bound, Total_Number_Bins=Sample_Sigma_Bins)
            LogN_Sample_SizeDistribution_Plain = []
            LogN_Sample_SizeDistribution = []

            Situation = f"Med={round(Sample_D_Median_List[i], 2)}-Sig={round(Sample_Sigma_List[j], 2)}"
            # SUM=0

            # Number Concentration
            for k in range(Sample_Sigma_Bins - 1):
                LogN_Sample_SizeDistribution_Plain.append(FN.LogN_Distribution(Median=Sample_D_Median_List[i], SigmaG=Sample_Sigma_List[j], Dp2=Diameter_Nano[k + 1], Dp1=Diameter_Nano[k]))
                LogN_Sample_SizeDistribution.append(LogN_Sample_SizeDistribution_Plain[k] * Sample_Total_Number_Concentration)
                # SUM+=LogN_Sample_SizeDistribution_Plain[k]

            Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_NC" + "-" + Situation, Extension="jpg")
            FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, Diameter_Nano[:-1], X_Label1, LogN_Sample_SizeDistribution, "Number Concentration,\nTotal Conc.= " + str(round(Sample_Total_Number_Concentration, 1)), Y_Label1, Situation)

            # Primary Particle Size
            Soot_Primary_Diameter_Median_meter = []
            for k in range(Sample_Sigma_Bins - 1):
                Soot_Primary_Diameter_Median_meter.append(FN.Primary_Particle_Size_meter(da_meter=Diameter_Meter[k], dp100_meter=Primary_Diameter_100nm, Dtem=Primary_D_TEM))
            Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_PPS" + "-" + Situation, Extension="jpg")
            FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, Diameter_Nano[:-1], X_Label1, Soot_Primary_Diameter_Median_meter, "Primary Particle Size", Y_Label2, Situation)

            # Finding Primary Particle Number and Size
            Primary_Diameter_Bank = []
            Primary_Number_Bank = []
            Primary_Probability_Bank = []

            # fig, ax = plt.subplots()
            # fig = plt.figure()
            # ax = fig.add_subplot(111, projection='3d')

            for k in range(Sample_Sigma_Bins - 1):
                Primary_Diameter, Primary_Number, Primary_Probability = FN.Primary_LogN_Generator(dp_Median_meter=Soot_Primary_Diameter_Median_meter[k], Sigma_g_Primary=Primary_Sigma_da_CTE, Sigma_g_Number=Primary_Sigma_da_CTE_Bound, Number_Points=Primary_Sigma_da_CTE_Nt, da_meter=Diameter_Meter[k], Soot_Prefactor_k=Soot_Prefactor_k_mc,
                                                                                                  Soot_Fractal_D=Soot_Fractal_D_mc)
                Primary_Diameter_Bank.append(Primary_Diameter)
                Primary_Number_Bank.append(Primary_Number)
                Primary_Probability_Bank.append(Primary_Probability)
            #     dummy = []
            #     for p in range(Primary_Sigma_da_CTE_Nt - 1):
            #         dummy.append(Diameter_Nano[k])
            #     ax.plot(dummy, Primary_Diameter_Bank[k][:], Primary_Probability_Bank[k][:])
            # plt.show()

            # Absorption RDG
            Absorption_Cross_Section_Sample = []
            Absorption_Cross_Section = 0

            for k in range(Sample_Sigma_Bins - 1):
                Absorption_Cross_Section_Specific_Particle_Size = 0
                for p in range(Primary_Sigma_da_CTE_Nt - 1):
                    Absorption_Cross_Section_Agg = FN.RDG_Absorption(K=Wave_Number, N=Primary_Number_Bank[k][p], Dp=Primary_Diameter_Bank[k][p], E=Soot_EM)  # within  aggregate
                    Absorption_Cross_Section_Specific_Particle_Size += Absorption_Cross_Section_Agg * Primary_Probability_Bank[k][p]
                Absorption_Cross_Section_Sample.append(Absorption_Cross_Section_Specific_Particle_Size)
                Absorption_Cross_Section += Absorption_Cross_Section_Sample[k] * LogN_Sample_SizeDistribution[k]

            Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_ABSCross" + "-" + Situation, Extension="jpg")
            FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, Diameter_Nano[:-1], X_Label1, Absorption_Cross_Section_Sample, "Absorption Cross Section=\n" + str(Absorption_Cross_Section), Y_Label3, Situation)

            Absorption_Efficiency_Sample = []

            for k in range(Sample_Sigma_Bins - 1):
                Absorption_Efficiency_Sample.append(FN.Absorption_Eff(Abs_Cross=Absorption_Cross_Section_Sample[k], da=Diameter_Meter[k]))

            Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_ABSEff" + "-" + Situation, Extension="jpg")
            FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, Diameter_Nano[:-1], X_Label1, Absorption_Efficiency_Sample, "Absorption Efficiency", Y_Label4, Situation)

            # Scattering RDG
            a = 3

            # Absorption_Cross_Section_Agg, Absorption_Cross_Section_Primary = FN.RDG_Absorption(Wave_Number, Primary_Number, Primary_Diameter, Soot_EM)  # within  aggregate
            # Absorption_Cross_Eff_1, Absorption_Cross_Primary_Eff_1 = FN.Absorption_Eff(Absorption_Cross_Section_Agg, Absorption_Cross_Section_Primary, diam[k], Primary_Diameter)
            # Absorption_Cross_Section_Agg_Total = np.sum(Absorption_Cross_Section_Agg * Primary_Probability)
            # Absorption_Cross_Section_Primary_Total = np.sum(Absorption_Cross_Section_Primary * Primary_Probability)
            # Absorption_Cross_Eff_Total = np.sum(Absorption_Cross_Eff_1 * Primary_Probability)
            # Absorption_Cross_Primary_Eff_Total = np.sum(Absorption_Cross_Primary_Eff_1 * Primary_Probability)
            # Scattering_Cross_Section_Agg, Scattering_Cross_Section_Primary = FN.RDG_Total_Scatteing(Wave_Number, Soot_FM, Primary_Number, Primary_Diameter, Betha)  # within  aggregate
            # Scattering_Cross_Section_Agg_Total = np.sum(Scattering_Cross_Section_Agg * Primary_Probability)
            # Scattering_Cross_Section_Primary_Total = np.sum(Scattering_Cross_Section_Primary * Primary_Probability)
            # Aggregate_Mass[i][j][k] = FN.Particle_Mass(FN.Effective_Density(eff_k, eff_dm, diam[k]), diam[k])
            # # for p in range (Theta_Number):
            # #   Differential_Scattering_Cross[i][j][k][p] = KN.RDG_Def_Scattering (Wave_Number, Theta_Radian[p], Soot_FM, Soot_Primary_Number[i][j][k], Soot_Primary_dp_median, Soot_Fractal_Df, Soot_Prefactor_kf)
            # Total_Absorption_Cross[i][j] = Total_Absorption_Cross[i][j] + LogNormal_pdf[i][j][k] * Absorption_Cross_Section_Agg_Total * Total_Number_Concentration_Sample
            # Total_Scattering_Cross[i][j] = Total_Scattering_Cross[i][j] + LogNormal_pdf[i][j][k] * Scattering_Cross_Section_Agg_Total * Total_Number_Concentration_Sample
            # Total_Particle_Mass[i][j] = Total_Particle_Mass[i][j] + LogNormal_pdf[i][j][k] * Aggregate_Mass[i][j][k] * Total_Number_Concentration_Sample
            #
            # Absorption_Cross[i][j][k] = Absorption_Cross_Section_Agg_Total
            # Absorption_Cross_Eff[i][j][k] = Absorption_Cross_Eff_Total
            # Absorption_Cross_Primary[i][j][k] = Absorption_Cross_Section_Primary_Total
            # Absorption_Cross_Primary_Eff[i][j][k] = Absorption_Cross_Primary_Eff_Total
            # Scattering_Cross[i][j][k] = Scattering_Cross_Section_Agg_Total
            # Scattering_Cross_Primary[i][j][k] = Scattering_Cross_Section_Primary_Total
        # Mass_Absorption_Eff[i][j] = Total_Absorption_Cross[i][j] / Total_Particle_Mass[i][j]
#         plt.figure(n)
#         plt.rcParams['mathtext.fontset'] = 'stix'
#         plt.rcParams['font.family'] = 'STIXGeneral'
#         plt.title("Absorption and Scattering, " + "D_Median= " + str(LogNormal_D_Median_Sample[i] * 1e9) + "(nm), Sigma=" + str(LogNormal_Sigma_Sample[j]) + ", Mass Abs Eff= " + "{:.3E}".format(Mass_Absorption_Eff[i][j]) + " (m^2/kg)")
#
#         plt.plot(Diameter_Nano, Absorption_Cross[i][j][:], label="Tot Abs= " + "{:.3E}".format(Total_Absorption_Cross[i][j]) + "(m^2)")
#         plt.plot(Diameter_Nano, Aggregate_Mass[i][j][:], label="Total Particle Mass (Kg)= " + "{:.3E}".format(Total_Particle_Mass[i][j]))
#         plt.plot(Diameter_Nano, Absorption_Cross_Primary[i][j][:], label="Primary Particle Abs")
#         plt.plot(Diameter_Nano, Scattering_Cross[i][j][:], label="Tot Sca= " + "{:.3E}".format(Total_Scattering_Cross[i][j]) + "(m^2)")
#         plt.plot(Diameter_Nano, Scattering_Cross_Primary[i][j][:], label="Primary Particle Sca")
#         plt.xlabel('Mobility-equivalent Diameter (nm)')
#         plt.ylabel('Cross Section (m^2)')
#         plt.xscale('log')
#         plt.yscale('log')
#         plt.grid(True, which='minor', alpha=0.2)
#         plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left', fancybox=True, fontsize='small')
#         GraphName = "D_Med-" + str(i) + "Sigma" + str(j) + ".jpg"
#         plt.savefig(GraphName, format='jpg', dpi=1000, bbox_inches='tight')
#         plt.clf()
#         plt.close()
#         n = n + 1
#
#         plt.figure(n)
#         plt.rcParams['mathtext.fontset'] = 'stix'
#         plt.rcParams['font.family'] = 'STIXGeneral'
#         plt.title("Absorption Efficiency, " + "D_Median= " + str(LogNormal_D_Median_Sample[i] * 1e9) + "(nm), Sigma=" + str(LogNormal_Sigma_Sample[j]) + ", Mass Abs Eff= " + "{:.3E}".format(Mass_Absorption_Eff[i][j]) + " (m^2/kg)")
#         plt.plot(Diameter_Nano, Absorption_Cross_Eff[i][j][:], label="Aggregate")
#         plt.plot(Diameter_Nano, Absorption_Cross_Primary_Eff[i][j][:], label="Primary")
#         plt.xlabel('Mobility-equivalent Diameter (nm)')
#         plt.ylabel('Efficiency')
#         plt.xscale('log')
#         plt.yscale('log')
#         plt.grid(True, which='minor', alpha=0.2)
#         plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left', fancybox=True, fontsize='small')
#         GraphName = "D_Med-" + str(i) + "Sigma" + str(j) + "_EFF" + ".jpg"
#         plt.savefig(GraphName, format='jpg', dpi=1000, bbox_inches='tight')
#         plt.clf()
#         plt.close()
#         n = n + 1
#
# a = 3
