# Programmed by Keyhan Babaee Under Prof. Steven Rogak supervision, https://github.com/KeyhanB
# Version 1.2
# Jan 18, 2018
import Functions as FN
import os
import numpy as np
from math import pi
from math import log

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

    Sample_Total_Number_Concentration = 5000  # /cm^3
    Sample_LogN_D_Median_Min = 200  # Smallest diameter of computation (nm)
    Sample_LogN_D_Median_Max = 350  # Largest diameter of computation (nm)
    Sample_LogN_D_Median_Bins = 1  # Number of the Steps
    Sample_LogN_Sigma_Min = 1.4  # Smallest Sigma G
    Sample_LogN_Sigma_Max = 1.8  # Largest Sigma G
    Sample_LogN_Sigma_Bins = 1  # Number of the Steps
    Sample_Sigma_Bound = 2.75  # Number of Sigma G to cover
    Sample_Sigma_Bins = 250  # Number of bins
    Sample_Diameter_Min = 20  # in nm
    Sample_Diameter_Max = 2000  # in nm
    ####################### Effective Density

    Eff_k = 0.418  # effective density variable for the result to be in kg/m^3
    Eff_dm = 2.2  # effective density variable for the result to be in kg/m^3
    Eff_rho_100nm = 502.1  # kg/m^3 for 100nm aggregate
    Eff_k = FN.Effective_Density_K_FromRho100nm(Eff_dm, Eff_rho_100nm)
    Soot_Rho_Cte = 1800

    Eff_dm_Final = 2.70
    Eff_dm_Number = 3
    Eff_dm_List = np.linspace(Eff_dm, Eff_dm_Final, num=Eff_dm_Number)
    ####################### Variation of Primary Particle Size within the sample

    Primary_Diameter_100nm = 20 * 1e-9  # Primary particle diameter of 100 nm aggregate
    # Primary_D_TEM = 0.34
    Primary_D_Alpha = 1.1
    Primary_D_TEM = (2 * Primary_D_Alpha - Eff_dm) / (2 * Primary_D_Alpha - 3)
    Primary_Sigma_da_CTE = 1  # Sigma G around specific particle primary particle size

    Primary_Sigma_da_CTE_Final = 1.4
    Primary_Sigma_da_CTE_Number = 2
    Primary_Sigma_da_CTE_List = np.linspace(Primary_Sigma_da_CTE, Primary_Sigma_da_CTE_Final, num=Primary_Sigma_da_CTE_Number)

    Primary_Sigma_da_CTE_Bound = 2.75  # Number of Sigma G to cover
    Primary_Sigma_da_CTE_Nt = 19 + 1  # Number of bins
    Primary_Betha = 0.9
    Primary_Y_Min = 5  # in nm
    Primary_Y_Max = 50  # in nm
    ####################### Refractive Index

    Soot_Refractive_Index = 1.8 - 0.6j  # Complex refractive index
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

    ####################### Wave Parameter

    Wave_Length = 500 * 1e-9  # nm
    Wave_Number = 2 * pi / Wave_Length  # k

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
    logging.info("Program Initiated!")

    #######################

    # Diameter_Meter = []
    # Diameter_Log = []
    # Diameter_Nano = []
    # Bound_D_Max = 1000 * 10 ** (-9)
    # Bound_D_Min = 30 * 10 ** (-9)
    # Bound_D_Bins = 300
    # Keyhan = 10 ** 9  # Conversion factor
    # D_ratio = (Bound_D_Max / Bound_D_Min) ** (1 / (Bound_D_Bins - 1))  # Generating particle diameters
    # for i in range(Bound_D_Bins):
    #     d = Bound_D_Min * D_ratio ** (i)
    #     Diameter_Meter.append(d)
    #     Diameter_Log.append(log(d))
    #     Diameter_Nano.append(d * Keyhan)

    #######################

    # Total_Particle_Mass = [[0 for j in range(LogNormal_Sigma_Number_Sample)] for i in range(LogNormal_D_Median_Number_Sample)]
    # Mass_Absorption_Eff = [[0 for j in range(LogNormal_Sigma_Number_Sample)] for i in range(LogNormal_D_Median_Number_Sample)]

    X_Label1 = 'Mobility-equivalent Diameter (nm)'
    Y_Label1 = 'dN/dLogDp, (#/cm' + "$^{}$".format(3) + ')'
    Y_Label2 = "Primary Particle Diameter-dp(nm)"
    Y_Label3 = "Absorption Cross Section(m^2/dLogDp)"
    Y_Label4 = "Absorption Efficiency(/dLogDp)"
    Y_Label5 = "Scattering Cross Section(m^2/dLogDp)"
    Y_Label6 = "Sample Mass(kg/dLogDp)"
    logging.info("Labels Created!")

    Save_Situation = []
    Save_Diameter = {}
    Save_LogN_Sample_SizeDistribution = {}
    Save_Soot_Primary_Diameter_Median_Nano = {}
    Save_Absorption_Cross_Section_Sample_dlnDp = {}
    Save_Absorption_Efficiency_Sample_dlnDp = {}
    Save_Differential_Scattering_Cross_Section_Full_dlnDp = {}
    Save_Scattering_Cross_Section_Total_Distribution_dlnDp = {}
    Save_Mass_Sample_Eff_dlnDp = {}
    Save_Mass_Sample_Rho_Cte_dlnDp = {}

    # Graph_Folder += f"/PD={round(Primary_Diameter_100nm * 10 ** 9, 0)}-D_TEM={round(Primary_D_TEM,2)}-PSig={round(Primary_Sigma_da_CTE,2)}-EffK={round(Eff_k,2)}-EffDm={round(Eff_dm,2)}-RI={Soot_Refractive_Index}-WL={round(Wave_Length * 10 ** 9, 0)}"

    for i in range(Sample_LogN_D_Median_Bins):
        for j in range(Sample_LogN_Sigma_Bins):
            for k in range(Eff_dm_Number):

                Eff_dm = Eff_dm_List[k]
                Eff_k = FN.Effective_Density_K_FromRho100nm(Eff_dm, Eff_rho_100nm)
                Primary_D_TEM = (2 * Primary_D_Alpha - Eff_dm) / (2 * Primary_D_Alpha - 3)

                for l in range(Primary_Sigma_da_CTE_Number):

                    Primary_Sigma_da_CTE = Primary_Sigma_da_CTE_List[l]

                    Diameter_Meter, Diameter_Log, Diameter_Nano = FN.Bins_LogN_Distributed(Median_Diameter=Sample_D_Median_List[i] * 10 ** (-9), Sigma_G=Sample_Sigma_List[j], Sigma_G_Bound=Sample_Sigma_Bound, Total_Number_Bins=Sample_Sigma_Bins)
                    LogN_Sample_SizeDistribution_Plain = []
                    LogN_Sample_SizeDistribution = []

                    Situation = f"Median Diameter={round(Sample_D_Median_List[i], 1)} (nm)- Sigma={round(Sample_Sigma_List[j], 1)}-Dm={round(Eff_dm, 2)}-k={round(Eff_k, 2)}-DTEM={round(Primary_D_TEM, 2)}-Sig={round(Primary_Sigma_da_CTE, 2)}"
                    Save_Situation.append(Situation)
                    Save_Diameter[Situation] = Diameter_Nano[:-1]

                    SUM = 0
                    SUM1 = 0
                    logging.info(f"Sample's particle bins generated:{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {Diameter_Nano}")
                    # Number Concentration
                    for k in range(Sample_Sigma_Bins - 1):
                        LogN_Sample_SizeDistribution_Plain.append(FN.LogN_Distribution(Median=Sample_D_Median_List[i], SigmaG=Sample_Sigma_List[j], Dp2=Diameter_Nano[k + 1], Dp1=Diameter_Nano[k]))
                        LogN_Sample_SizeDistribution.append(LogN_Sample_SizeDistribution_Plain[k] * Sample_Total_Number_Concentration)
                        SUM += LogN_Sample_SizeDistribution_Plain[k]
                        SUM1 += LogN_Sample_SizeDistribution[k]

                    logging.info(f"Sample's particle bins populated:{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {LogN_Sample_SizeDistribution_Plain}")
                    logging.info(f"LogNormal Check:{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {SUM}: {SUM1}")

                    Save_LogN_Sample_SizeDistribution[Situation] = LogN_Sample_SizeDistribution

                    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_NC" + "-" + Situation, Extension="jpg")
                    # FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=LogN_Sample_SizeDistribution, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_Legend="Number Concentration,\nTotal Conc.= " + str(round(Sample_Total_Number_Concentration, 1)), Y_label1=Y_Label1, Plot_Title=Situation)

                    # Primary Particle Size

                    Soot_Primary_Diameter_Median_meter = []
                    Soot_Primary_Diameter_Median_Nano = []
                    for k in range(Sample_Sigma_Bins - 1):
                        Soot_Primary_Diameter_Median_meter.append(FN.Primary_Particle_Size_meter(da_meter=Diameter_Meter[k], dp100_meter=Primary_Diameter_100nm, Dtem=Primary_D_TEM))
                        Soot_Primary_Diameter_Median_Nano.append(Soot_Primary_Diameter_Median_meter[k] * 10 ** 9)
                    logging.info(f"Primary Particle Diameter(m):{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {Soot_Primary_Diameter_Median_Nano}")

                    Save_Soot_Primary_Diameter_Median_Nano[Situation] = Soot_Primary_Diameter_Median_Nano

                    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_PPS" + "-" + Situation, Extension="jpg")
                    # FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Soot_Primary_Diameter_Median_Nano, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, Y_Min=Primary_Y_Min, Y_Max=Primary_Y_Max, X_Label=X_Label1, Y_Legend="Primary Particle Size", Y_label1=Y_Label2, Plot_Title=Situation)

                    # Finding Primary Particle Number and Size
                    Primary_Diameter_Bank = []
                    Primary_Number_Bank = []
                    Primary_Probability_Bank = []

                    #
                    SUM2 = []
                    for k in range(Sample_Sigma_Bins - 1):
                        SUM_Temp = 0
                        Primary_Diameter, Primary_Number, Primary_Probability = FN.Primary_LogN_Generator(dp_Median_meter=Soot_Primary_Diameter_Median_meter[k], Sigma_g_Primary=Primary_Sigma_da_CTE, Sigma_g_Number=Primary_Sigma_da_CTE_Bound, Number_Points=Primary_Sigma_da_CTE_Nt, da_meter=Diameter_Meter[k], Soot_Prefactor_k=Soot_Prefactor_k_mc,
                                                                                                          Soot_Fractal_D=Soot_Fractal_D_mc)
                        Primary_Diameter_Bank.append(Primary_Diameter)
                        Primary_Number_Bank.append(Primary_Number)
                        Primary_Probability_Bank.append(Primary_Probability)
                        for m in range(len(Primary_Probability)):
                            SUM_Temp += Primary_Probability[m]
                        SUM2.append(SUM_Temp)
                        logging.debug(f"Primary Particle Diameter(m):{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {Primary_Diameter}")
                        logging.debug(f"Primary Particle Number:{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {Primary_Number}")
                    # FN.Fig_Plot_3D_Show_XCte(Diameter_Nano[:-1], Primary_Diameter_Bank, Primary_Probability_Bank)
                    logging.info(f"Primary Particle Probability Check:{Situation}:Sigma Bound:{Sample_Sigma_Bound}: {SUM2}")

                    # Absorption RDG
                    Absorption_Cross_Section_Sample = []
                    Absorption_Cross_Section = 0
                    Absorption_Cross_Section_Sample_dlnDp = []

                    for k in range(Sample_Sigma_Bins - 1):
                        Absorption_Cross_Section_Specific_Particle_Size = 0
                        for p in range(Primary_Sigma_da_CTE_Nt - 1):
                            Absorption_Cross_Section_Agg = FN.RDG_Absorption(K=Wave_Number, N=Primary_Number_Bank[k][p], Dp=Primary_Diameter_Bank[k][p], E=Soot_EM)  # within  aggregate
                            Absorption_Cross_Section_Specific_Particle_Size += Absorption_Cross_Section_Agg * Primary_Probability_Bank[k][p]

                        ABS = Absorption_Cross_Section_Specific_Particle_Size * LogN_Sample_SizeDistribution[k]
                        Absorption_Cross_Section_Sample_dlnDp.append(ABS / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))
                        Absorption_Cross_Section_Sample.append(ABS)
                        Absorption_Cross_Section += ABS

                    Save_Absorption_Cross_Section_Sample_dlnDp[Situation] = Absorption_Cross_Section_Sample_dlnDp

                    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_AbsCross" + "-" + Situation, Extension="jpg")
                    # FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Absorption_Cross_Section_Sample_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_Legend="Absorption Cross Section=\n" + str(round(Absorption_Cross_Section * 10 ** 12,2)) + " um"+ "$^{}$".format(2), Y_label1=Y_Label3,
                    # Plot_Title=Situation)

                    Absorption_Efficiency_Sample = []
                    Absorption_Efficiency_Sample_dlnDp = []

                    for k in range(Sample_Sigma_Bins - 1):
                        Absorption_Efficiency_Sample.append(FN.Absorption_Eff(Abs_Cross=Absorption_Cross_Section_Sample[k], da=Diameter_Meter[k]))
                        Absorption_Efficiency_Sample_dlnDp.append(Absorption_Efficiency_Sample[k] / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))

                    Save_Absorption_Efficiency_Sample_dlnDp[Situation] = Absorption_Efficiency_Sample_dlnDp

                    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_AbsEff" + "-" + Situation, Extension="jpg")
                    # FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Absorption_Efficiency_Sample_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_Legend="Absorption Efficiency", Y_label1=Y_Label4, Plot_Title=Situation)

                    # Scattering RDG
                    # Differential
                    Differential_Scattering_Cross_Section_Full = []
                    Differential_Scattering_Cross_Section_Full_dlnDp = []
                    Scattering_Cross_Section_Diff_Tot = 0
                    qDp_Full = []
                    for k in range(Sample_Sigma_Bins - 1):
                        qDp_Temp = []
                        Differential_Scattering_Cross_Section_T = 0
                        for t in range(Theta_Number):
                            q = FN.Scattering_Wave_Vector(WaveLength_meter=Wave_Length, Theta_radian=Theta_Radian[t])
                            qDp_Temp.append(q * Soot_Primary_Diameter_Median_meter[k])
                            Differential_Scattering_Cross_Section_Specific_Particle_Size = 0
                            for p in range(Primary_Sigma_da_CTE_Nt - 1):
                                Differential_Scattering_Cross_Section_Agg = FN.RDG_Def_Scattering(K=Wave_Number, N=Primary_Number_Bank[k][p], Dp=Primary_Diameter_Bank[k][p], q=q, F=Soot_FM, D_RDG=Soot_Fractal_D_mc_RDG, K_RDG=Soot_Prefactor_k_mc_RDG, Formula=2)
                                Differential_Scattering_Cross_Section_Specific_Particle_Size += Differential_Scattering_Cross_Section_Agg * Primary_Probability_Bank[k][p]

                            Differential_Scattering_Cross_Section_T += FN.Diff_Integral_Phi(Differential_Scattering_Cross_Section_Specific_Particle_Size, Phi_Radian, Theta_Radian[t], Theta_Diff, Phi_Diff)

                        qDp_Full.append(qDp_Temp)
                        Diff_Scatter = Differential_Scattering_Cross_Section_T * LogN_Sample_SizeDistribution[k]
                        Differential_Scattering_Cross_Section_Full.append(Diff_Scatter)
                        Differential_Scattering_Cross_Section_Full_dlnDp.append(Diff_Scatter / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))
                        Scattering_Cross_Section_Diff_Tot += Diff_Scatter

                    Save_Differential_Scattering_Cross_Section_Full_dlnDp[Situation] = Differential_Scattering_Cross_Section_Full_dlnDp

                    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_ScaDiffCross" + "-" + Situation, Extension="jpg")
                    # FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Differential_Scattering_Cross_Section_Full_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1,
                    # Y_Legend="Total Scattering Cross Section (Differential)=\n" + str(round(Scattering_Cross_Section_Diff_Tot * 10 ** 12,2))+ " um"+ "$^{}$".format(2),
                    # Y_label1=Y_Label5, Plot_Title=Situation)

                    # Total Scattering
                    Scattering_Cross_Section_Total_Distribution = []
                    Scattering_Cross_Section_Total_Distribution_dlnDp = []
                    Scattering_Cross_Section_Total = 0
                    for k in range(Sample_Sigma_Bins - 1):

                        Scattering_Cross_Section_Total_Specific_Particle_Size = 0
                        for p in range(Primary_Sigma_da_CTE_Nt - 1):
                            Scattering_Cross_Section_Total_Agg = FN.RDG_Total_Scattering(K=Wave_Number, N=Primary_Number_Bank[k][p], Dp=Primary_Diameter_Bank[k][p], F=Soot_FM, D_RDG=Soot_Fractal_D_mc_RDG, K_RDG=Soot_Prefactor_k_mc_RDG, Formula=2)
                            Scattering_Cross_Section_Total_Specific_Particle_Size += Scattering_Cross_Section_Total_Agg * Primary_Probability_Bank[k][p]

                        Total_Scatter = Scattering_Cross_Section_Total_Specific_Particle_Size * LogN_Sample_SizeDistribution[k]
                        Scattering_Cross_Section_Total_Distribution.append(Total_Scatter)
                        Scattering_Cross_Section_Total_Distribution_dlnDp.append(Total_Scatter / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))
                        Scattering_Cross_Section_Total += Total_Scatter

                    Save_Scattering_Cross_Section_Total_Distribution_dlnDp[Situation] = Scattering_Cross_Section_Total_Distribution_dlnDp

                    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_ScaTotCross" + "-" + Situation, Extension="jpg")
                    # FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Scattering_Cross_Section_Total_Distribution_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_Legend="Scattering Cross Section (Total)=\n" + str(round(Scattering_Cross_Section_Total * 10 ** 12,2)) + " um"+ "$^{}$".format(2),
                    # Y_label1=Y_Label5, Plot_Title=Situation)

                    # MAC_MSC
                    Mass_Sample_Eff = []
                    Mass_Sample_Rho_Cte = []
                    Mass_Sample_Eff_dlnDp = []
                    Mass_Sample_Rho_Cte_dlnDp = []
                    Density = []
                    Total_Mass_Eff = 0
                    Total_Mass_Rho_Cte = 0

                    for k in range(Sample_Sigma_Bins - 1):
                        Density.append(FN.Effective_Density(K=Eff_k, Dm=Eff_dm, da=Diameter_Meter[k]))
                        Mass_Sample_Eff.append(FN.Mass_Calc(rho=Density[k], da=Diameter_Meter[k]))
                        Mass_Sample_Rho_Cte.append(FN.Mass_Calc(rho=Soot_Rho_Cte, da=Diameter_Meter[k]))
                        Mass_Effective = Mass_Sample_Eff[k] * LogN_Sample_SizeDistribution[k]
                        Mass_Cte = Mass_Sample_Rho_Cte[k] * LogN_Sample_SizeDistribution[k]
                        Mass_Sample_Eff_dlnDp.append(Mass_Effective / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))
                        Mass_Sample_Rho_Cte_dlnDp.append(Mass_Cte / (log(Diameter_Nano[k + 1]) - log(Diameter_Nano[k])))
                        Total_Mass_Eff += Mass_Effective
                        Total_Mass_Rho_Cte += Mass_Cte

                    Total_Mass_Eff = Total_Mass_Eff * 1000  # kg to g
                    Total_Mass_Rho_Cte = Total_Mass_Rho_Cte * 1000
                    MAC = Absorption_Cross_Section / Total_Mass_Eff
                    MSC = Scattering_Cross_Section_Total / Total_Mass_Eff
                    MSC_Diff = Scattering_Cross_Section_Diff_Tot / Total_Mass_Eff

                    MAC_Rho_Cte = Absorption_Cross_Section / Total_Mass_Rho_Cte
                    MSC_Rho_Cte = Scattering_Cross_Section_Total / Total_Mass_Rho_Cte
                    MSC_Diff_Rho_Cte = Scattering_Cross_Section_Diff_Tot / Total_Mass_Rho_Cte

                    Save_Mass_Sample_Eff_dlnDp[Situation] = Mass_Sample_Eff_dlnDp
                    Save_Mass_Sample_Rho_Cte_dlnDp[Situation] = Mass_Sample_Rho_Cte_dlnDp

                    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_Mass" + "-" + Situation, Extension="jpg")
                    # FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Mass_Sample_Eff_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_Legend="MAC=" + str(round(MAC, 2)) + "\n" + "MSC=" + str(round(MSC, 2)) + "\n" + "MSC-Diff=" + str(round(MSC_Diff, 2)), Y_label1=Y_Label6,
                    # Plot_Title=Situation)

                    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_Mass_Rho_Cte" + "-" + Situation, Extension="jpg")
                    # FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear(Address=Address, X_Array=Diameter_Nano[:-1], Y_array=Mass_Sample_Rho_Cte_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_Legend="MAC=" + str(round(MAC_Rho_Cte, 2)) + "\n" + "MSC=" + str(round(MSC_Rho_Cte, 2)) + "\n" + "MSC-Diff=" + str(round(MSC_Diff_Rho_Cte, 2)),
                    # Y_label1=Y_Label6, Plot_Title=Situation)

    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="LogN_Sample_SizeDistribution", Extension="jpg")
    FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_LogN_Sample_SizeDistribution, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_label1=Y_Label1, Plot_Title="Number Concentration")

    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Soot_Primary_Diameter_Median", Extension="jpg")
    FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Soot_Primary_Diameter_Median_Nano, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_label1=Y_Label2, Plot_Title="Primary Particle Diameter")

    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Absorption_Cross_Section_Sample", Extension="jpg")
    FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Absorption_Cross_Section_Sample_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_label1=Y_Label3, Plot_Title="Absorption Cross Section_ Total")

    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Absorption_Efficiency_Sample", Extension="jpg")
    FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Absorption_Efficiency_Sample_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_label1=Y_Label4, Plot_Title="Absorption Efficiency")

    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Differential_Scattering_Cross_Section", Extension="jpg")
    FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Differential_Scattering_Cross_Section_Full_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_label1=Y_Label5, Plot_Title="Scattering Cross Section_ Diff")

    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Scattering_Cross_Section_Total", Extension="jpg")
    FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Scattering_Cross_Section_Total_Distribution_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_label1=Y_Label5, Plot_Title="Scattering Cross Section_ Total")

    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="Mass_Sample_Eff", Extension="jpg")
    FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Mass_Sample_Eff_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_label1=Y_Label6, Plot_Title="Number Concentration")

    Address = FN.File_Pointer(Main=script_dir, FolderName=Graph_Folder, FileName="_Mass_Sample_Rho_Cte", Extension="jpg")
    FN.Fig_Plot_Save_1Lines_X_Log_Y_Dictionary_Linear(Address=Address, Identifier=Save_Situation, X_Array=Save_Diameter, Y_array=Save_Mass_Sample_Rho_Cte_dlnDp, X_Min=Sample_Diameter_Min, X_Max=Sample_Diameter_Max, X_Label=X_Label1, Y_label1=Y_Label6, Plot_Title="Number Concentration")
