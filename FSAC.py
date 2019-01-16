# Programmed by Keyhan Babaee Under Prof. Steven Rogak supervision, https://github.com/KeyhanB
# Version 0.4
# Jan 2, 2018
import Functions as FN
import os
import numpy as np
from math import pi
from math import log
from math import sin
from math import cos
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

    Sample_Total_Number_Concentration = 1  # /cm^3
    Sample_LogN_D_Median_Min = 200  # Smallest diameter of computation (nm)
    Sample_LogN_D_Median_Max = 550  # Largest diameter of computation (nm)
    Sample_LogN_D_Median_Bins = 2  # Number of the Steps
    Sample_LogN_Sigma_Min = 1  # Smallest Sigma G
    Sample_LogN_Sigma_Max = 1.7  # Largest Sigma G
    Sample_LogN_Sigma_Bins = 2  # Number of the Steps
    Sample_Sigma_Bound = 2.5  # Number of Sigma G to cover
    Sample_Sigma_Bins = 300  # Number of bins

    ####################### Varation of Primary Particle Size within the sample

    # Primary_Diameter_100nm = 16.48 * 1e-9  # Primary particle diameter of 100 nm aggregate
    # Primary_D_TEM = 0.45
    Primary_D_TEM = 0
    Primary_Sigma_da_CTE = 1  # Sigma G around specific particle primary particle size
    Primary_Sigma_da_CTE_Bound = 2.5  # Number of Sigma G to cover
    Primary_Sigma_da_CTE_Nt = 19 + 1  # Number of bins
    Primary_Betha = 0.9

    ####################### Refractive Index

    Soot_Refractive_Index = 1.6 - 0.6j  # Complex refractive index
    Soot_Primary_Agg_dp_sigma = 1  # Geometric std of primary particles within aggregate

    ####################### Continuum regime

    Soot_Fractal_D_mc = 0.52  # Fractal Properties for continuum regime
    Soot_Prefactor_k_mc = 0.85 + (0.03 * Soot_Primary_Agg_dp_sigma ** (4.4))  # Fractal continuum regime
    Soot_Fractal_D_mc_RDG = 1 / Soot_Fractal_D_mc
    Soot_Prefactor_k_mc_RDG = (Primary_Betha / Soot_Prefactor_k_mc) ** (1 / Soot_Fractal_D_mc)

    ####################### Free molecular regime

    # Soot_Fractal_D_mfm = 0.46  # Fractal Properties free molecular regime
    # Soot_Prefactor_k_mfm = 0.94 + (0.03 * Soot_Primary_Agg_dp_sigma ** (4.8))  # Fractal free molecular regime
    # Soot_Fractal_D_mfm_RDG = 1 / Soot_Fractal_D_mfm
    # Soot_Prefactor_k_mfm_RDG = (Primary_Betha / Soot_Prefactor_k_mfm) ** (1 / Soot_Fractal_D_mfm)

    ####################### Effective Density

    Eff_k = 0.418  # effective density variable for the result to be in kg/m^3
    Eff_dm = 2.56  # effective density variable for the result to be in kg/m^3

    ####################### Wave Parameter

    Wave_Length = 1064 * 1e-9  # nm
    Wave_Number = 2 * pi / Wave_Length  # k
    Primary_Diameter_100nm = 0.0886 * Wave_Length / pi
    ####################### Scattering Parameters

    Theta_Number = 200
    Phi_Number = 200
    Theta_Start = 0
    Theta_Finish = pi
    Phi_Start = 0
    Phi_Finish = 2 * pi
    Theta_Radian = np.linspace(Theta_Start, Theta_Finish, num=Theta_Number)
    Phi_Radian = np.linspace(Phi_Start, Phi_Finish, Phi_Number)
    Theta_Diff = abs((Theta_Finish - Theta_Start) / Theta_Number)
    Phi_Diff = abs((Phi_Finish - Phi_Start) / Phi_Number)
    # alpha=2*pi*radius/lambda

    #######################

    Soot_Complex = (((Soot_Refractive_Index ** 2) - 1) / ((Soot_Refractive_Index ** 2) + 2))
    Soot_FM = (abs(Soot_Complex)) ** 2
    Soot_EM = Soot_Complex.imag
    Sample_Sigma_List = np.linspace(Sample_LogN_Sigma_Min, Sample_LogN_Sigma_Max, num=Sample_LogN_Sigma_Bins)
    Sample_D_Median_List = np.linspace(Sample_LogN_D_Median_Min, Sample_LogN_D_Median_Max, num=Sample_LogN_D_Median_Bins)

    #######################

    Diameter_Meter = []
    Diameter_Log = []
    Diameter_Nano = []
    Bound_D_Max = 1000 * 10 ** (-9)
    Bound_D_Min = 30 * 10 ** (-9)
    Bound_D_Bins = 300
    Keyhan = 10 ** 9  # Conversion factor
    D_ratio = (Bound_D_Max / Bound_D_Min) ** (1 / (Bound_D_Bins - 1))  # Generating particle diameters
    for i in range(Bound_D_Bins):
        d = Bound_D_Min * D_ratio ** (i)
        Diameter_Meter.append(d)
        Diameter_Log.append(log(d))
        Diameter_Nano.append(d * Keyhan)

    #######################

    # Test
    Test_Np = [5, 10, 20, 50, 100, 150, 199, 348, 546, 893]
    Test_Wave1 = 532 * 10 ** (-9)
    Test_Wave2 = 1064 * 10 ** (-9)
    Test_Wave_Number1 = 2 * pi / Test_Wave1  # k
    Test_Wave_Number2 = 2 * pi / Test_Wave2  # k
    Test_Primary_Diameter1 = 0.177 * Test_Wave1 / pi
    Test_Primary_Diameter2 = 0.0886 * Test_Wave2 / pi
    Test_Soot_Prefactor_k_mc_RDG = 2.3
    Test_Soot_Fractal_D_mc_RDG = 1.78
    Test_Soot_Refractive_Index = 1.6 - 0.6j
    Test_Soot_Complex = (((Test_Soot_Refractive_Index ** 2) - 1) / ((Test_Soot_Refractive_Index ** 2) + 2))
    Test_Soot_FM = (abs(Test_Soot_Complex)) ** 2
    Test_Soot_EM = Test_Soot_Complex.imag
    Test_Absorption_Cross_Section_Agg1 = []
    Test_Absorption_Cross_Section_Agg2 = []
    Test_Scattering_Cross_Section_Agg1 = []
    Test_Scattering_Cross_Section_Agg2 = []
    for i in range(len(Test_Np)):
        Test_Absorption_Cross_Section_Agg1.append((10 ** 18) * FN.RDG_Absorption(K=Test_Wave_Number1, N=Test_Np[i], Dp=Test_Primary_Diameter1, E=Test_Soot_EM))
        Test_Absorption_Cross_Section_Agg2.append((10 ** 18) * FN.RDG_Absorption(K=Test_Wave_Number2, N=Test_Np[i], Dp=Test_Primary_Diameter2, E=Test_Soot_EM))
        Test_Scattering_Cross_Section_Agg1.append((10 ** 18) * FN.RDG_Total_Scattering(K=Test_Wave_Number1, N=Test_Np[i], Dp=Test_Primary_Diameter1, F=Test_Soot_FM, D_RDG=Test_Soot_Fractal_D_mc_RDG, K_RDG=Test_Soot_Prefactor_k_mc_RDG, Formula=2))
        Test_Scattering_Cross_Section_Agg2.append((10 ** 18) * FN.RDG_Total_Scattering(K=Test_Wave_Number2, N=Test_Np[i], Dp=Test_Primary_Diameter2, F=Test_Soot_FM, D_RDG=Test_Soot_Fractal_D_mc_RDG, K_RDG=Test_Soot_Prefactor_k_mc_RDG, Formula=2))

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
    Y_Label5 = "Scattering Cross Section(m^2)"
    Y_Label6 = "Sample Mass(kg)"
    #
    for i in range(Sample_LogN_D_Median_Bins):
        for j in range(Sample_LogN_Sigma_Bins):

            # Diameter_Meter, Diameter_Log, Diameter_Nano = FN.Bins_LogN_Distributed(Median_Diameter=Sample_D_Median_List[i] * 10 ** (-9), Sigma_G=Sample_Sigma_List[j], Sigma_G_Bound=Sample_Sigma_Bound, Total_Number_Bins=Sample_Sigma_Bins)
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

            Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_AbsCross" + "-" + Situation, Extension="jpg")
            FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, Diameter_Nano[:-1], X_Label1, Absorption_Cross_Section_Sample, "Absorption Cross Section=\n" + str(Absorption_Cross_Section), Y_Label3, Situation)

            Absorption_Efficiency_Sample = []

            for k in range(Sample_Sigma_Bins - 1):
                Absorption_Efficiency_Sample.append(FN.Absorption_Eff(Abs_Cross=Absorption_Cross_Section_Sample[k], da=Diameter_Meter[k]))

            Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_AbsEff" + "-" + Situation, Extension="jpg")
            FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, Diameter_Nano[:-1], X_Label1, Absorption_Efficiency_Sample, "Absorption Efficiency", Y_Label4, Situation)

            # Scattering RDG
            # Differential
            Differential_Scattering_Cross_Section_Full = []
            Scattering_Cross_Section_Diff = 0
            for k in range(Sample_Sigma_Bins - 1):
                # Differential_Scattering_Cross_Section_Theta = []
                Differential_Scattering_Cross_Section_T = 0
                for t in range(Theta_Number):
                    q = FN.Scattering_Wave_Vector(WaveLength_meter=Wave_Length, Theta_radian=Theta_Radian[t])
                    Differential_Scattering_Cross_Section_Specific_Particle_Size = 0
                    for p in range(Primary_Sigma_da_CTE_Nt - 1):
                        Differential_Scattering_Cross_Section_Agg = FN.RDG_Def_Scattering(K=Wave_Number, N=Primary_Number_Bank[k][p], Dp=Primary_Diameter_Bank[k][p], q=q, F=Soot_FM, D_RDG=Soot_Fractal_D_mc_RDG, K_RDG=Soot_Prefactor_k_mc_RDG, check=Primary_Probability_Bank[k][p], Number=k)
                        Differential_Scattering_Cross_Section_Specific_Particle_Size += Differential_Scattering_Cross_Section_Agg * Primary_Probability_Bank[k][p]

                    Differential_Scattering_Cross_Section_Specific_Particle_Size = Differential_Scattering_Cross_Section_Specific_Particle_Size * sin(Theta_Radian[t]) * Theta_Diff
                    for phi in range(Phi_Number):
                        Differential_Scattering_Cross_Section_T += Differential_Scattering_Cross_Section_Specific_Particle_Size * Phi_Diff

                Differential_Scattering_Cross_Section_Full.append(Differential_Scattering_Cross_Section_T)
                Scattering_Cross_Section_Diff += Differential_Scattering_Cross_Section_Full[k] * LogN_Sample_SizeDistribution[k]

            Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_ScaDiffCross" + "-" + Situation, Extension="jpg")
            FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, Diameter_Nano[:-1], X_Label1, Differential_Scattering_Cross_Section_Full, "Scattering Cross Section (From Differential)=\n" + str(Scattering_Cross_Section_Diff), Y_Label5, Situation)

            # Total Scattering
            Scattering_Cross_Section_Total_Distribution = []
            Scattering_Cross_Section_Total = 0
            for k in range(Sample_Sigma_Bins - 1):

                Scattering_Cross_Section_Total_Specific_Particle_Size = 0
                for p in range(Primary_Sigma_da_CTE_Nt - 1):
                    Scattering_Cross_Section_Total_Agg = FN.RDG_Total_Scattering(K=Wave_Number, N=Primary_Number_Bank[k][p], Dp=Primary_Diameter_Bank[k][p], F=Soot_FM, D_RDG=Soot_Fractal_D_mc_RDG, K_RDG=Soot_Prefactor_k_mc_RDG)
                    Scattering_Cross_Section_Total_Specific_Particle_Size += Scattering_Cross_Section_Total_Agg * Primary_Probability_Bank[k][p]

                Scattering_Cross_Section_Total_Distribution.append(Scattering_Cross_Section_Total_Specific_Particle_Size)
                Scattering_Cross_Section_Total += Scattering_Cross_Section_Total_Distribution[k] * LogN_Sample_SizeDistribution[k]

            Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_ScaTotCross" + "-" + Situation, Extension="jpg")
            FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, Diameter_Nano[:-1], X_Label1, Scattering_Cross_Section_Total_Distribution, "Scattering Cross Section (Total)=\n" + str(Scattering_Cross_Section_Total), Y_Label5, Situation)

            # MAC_MSC
            Mass_Sample_Eff = []
            Density = []
            Total_Mass_Eff = 0
            for k in range(Sample_Sigma_Bins - 1):
                Density.append(FN.Effective_Density(K=Eff_k, Dm=Eff_dm, da=Diameter_Meter[k]))
                Mass_Sample_Eff.append(FN.Mass_Calc(rho=Density[k], da=Diameter_Meter[k]))
                Total_Mass_Eff += Mass_Sample_Eff[k] * LogN_Sample_SizeDistribution[k]
            Total_Mass_Eff = Total_Mass_Eff * 1000  # kg to g
            MAC = Absorption_Cross_Section / Total_Mass_Eff
            MSC = Scattering_Cross_Section_Total / Total_Mass_Eff
            Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_Mass" + "-" + Situation, Extension="jpg")
            FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address, Diameter_Nano[:-1], X_Label1, Mass_Sample_Eff, "MAC=" + str(MAC) + "\n" + "MSC=" + str(MSC), Y_Label6, Situation)
            a = 3
