from math import pi
import Functions as FN
import numpy as np

# Test

if __name__ == "__main__":

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

    Theta_Number = 300
    Phi_Number = 300
    Theta_Start = 0
    Theta_Finish = pi
    Phi_Start = 0
    Phi_Finish = 2 * pi

    Theta1 = 30 / 180 * pi
    Theta2 = 150 / 180 * pi

    Test_Soot_Complex = (((Test_Soot_Refractive_Index ** 2) - 1) / ((Test_Soot_Refractive_Index ** 2) + 2))
    Test_Soot_FM = (abs(Test_Soot_Complex)) ** 2
    Test_Soot_EM = Test_Soot_Complex.imag
    Test_Absorption_Cross_Section_Agg1 = []
    Test_Absorption_Cross_Section_Agg2 = []
    Test_Scattering_Cross_Section_Agg1 = []
    Test_Scattering_Cross_Section_Agg2 = []
    Def_ScatteringW1_T1 = []
    Def_ScatteringW1_T2 = []
    Def_ScatteringW2_T1 = []
    Def_ScatteringW2_T2 = []
    Def_Scattering_RatioW1 = []
    Def_Scattering_RatioW2 = []

    Theta_Radian = np.linspace(Theta_Start, Theta_Finish, num=Theta_Number)
    Phi_Radian = np.linspace(Phi_Start, Phi_Finish, Phi_Number)
    Theta_Diff = abs((Theta_Finish - Theta_Start) / Theta_Number)
    Phi_Diff = abs((Phi_Finish - Phi_Start) / Phi_Number)

    for i in range(len(Test_Np)):
        Test_Absorption_Cross_Section_Agg1.append((10 ** 18) * FN.RDG_Absorption(K=Test_Wave_Number1, N=Test_Np[i], Dp=Test_Primary_Diameter1, E=Test_Soot_EM))
        Test_Absorption_Cross_Section_Agg2.append((10 ** 18) * FN.RDG_Absorption(K=Test_Wave_Number2, N=Test_Np[i], Dp=Test_Primary_Diameter2, E=Test_Soot_EM))
        Test_Scattering_Cross_Section_Agg1.append((10 ** 18) * FN.RDG_Total_Scattering(K=Test_Wave_Number1, N=Test_Np[i], Dp=Test_Primary_Diameter1, F=Test_Soot_FM, D_RDG=Test_Soot_Fractal_D_mc_RDG, K_RDG=Test_Soot_Prefactor_k_mc_RDG, Formula=2))
        Test_Scattering_Cross_Section_Agg2.append((10 ** 18) * FN.RDG_Total_Scattering(K=Test_Wave_Number2, N=Test_Np[i], Dp=Test_Primary_Diameter2, F=Test_Soot_FM, D_RDG=Test_Soot_Fractal_D_mc_RDG, K_RDG=Test_Soot_Prefactor_k_mc_RDG, Formula=2))
        qW1_T1 = FN.Scattering_Wave_Vector(WaveLength_meter=Test_Wave1, Theta_radian=Theta1)
        qW1_T2 = FN.Scattering_Wave_Vector(WaveLength_meter=Test_Wave1, Theta_radian=Theta2)
        qW2_T1 = FN.Scattering_Wave_Vector(WaveLength_meter=Test_Wave2, Theta_radian=Theta1)
        qW2_T2 = FN.Scattering_Wave_Vector(WaveLength_meter=Test_Wave2, Theta_radian=Theta2)
        Def_ScatteringW1_T1.append(FN.RDG_Def_Scattering(K=Test_Wave_Number1, N=Test_Np[i], Dp=Test_Primary_Diameter1, q=qW1_T1, F=Test_Soot_FM, D_RDG=Test_Soot_Fractal_D_mc_RDG, K_RDG=Test_Soot_Prefactor_k_mc_RDG, Formula=2))
        Def_ScatteringW1_T2.append(FN.RDG_Def_Scattering(K=Test_Wave_Number1, N=Test_Np[i], Dp=Test_Primary_Diameter1, q=qW1_T2, F=Test_Soot_FM, D_RDG=Test_Soot_Fractal_D_mc_RDG, K_RDG=Test_Soot_Prefactor_k_mc_RDG, Formula=2))
        Def_Scattering_RatioW1.append(Def_ScatteringW1_T1[i] / Def_ScatteringW1_T2[i])
        Def_ScatteringW2_T1.append(FN.RDG_Def_Scattering(K=Test_Wave_Number2, N=Test_Np[i], Dp=Test_Primary_Diameter2, q=qW2_T1, F=Test_Soot_FM, D_RDG=Test_Soot_Fractal_D_mc_RDG, K_RDG=Test_Soot_Prefactor_k_mc_RDG, Formula=2))
        Def_ScatteringW2_T2.append(FN.RDG_Def_Scattering(K=Test_Wave_Number2, N=Test_Np[i], Dp=Test_Primary_Diameter2, q=qW2_T2, F=Test_Soot_FM, D_RDG=Test_Soot_Fractal_D_mc_RDG, K_RDG=Test_Soot_Prefactor_k_mc_RDG, Formula=2))
        Def_Scattering_RatioW2.append(Def_ScatteringW2_T1[i] / Def_ScatteringW2_T2[i])
    FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear("Test/Test_Wave1.jpg", Test_Np, "Number of Primary Particle", Def_Scattering_RatioW1, "Wave Length=" + str(round(Test_Wave1 * 1e9)) + " nm", "Ratio 150/30", "Scattering Cross Section Ratio")
    FN.Fig_Plot_Save_1Lines_X_Log_Y_Linear("Test/Test_Wave2.jpg", Test_Np, "Number of Primary Particle", Def_Scattering_RatioW2, "Wave Length=" + str(round(Test_Wave2 * 1e9)) + " nm", "Ratio 150/30", "Scattering Cross Section Ratio")

    qDp1 = []
    qDp2 = []
    for t in range(Theta_Number):
        qW1_T1 = FN.Scattering_Wave_Vector(WaveLength_meter=Test_Wave1, Theta_radian=Theta_Radian[t])
        qW2_T1 = FN.Scattering_Wave_Vector(WaveLength_meter=Test_Wave2, Theta_radian=Theta_Radian[t])
        qDp1.append(qW1_T1 * Test_Primary_Diameter1)
        qDp2.append(qW2_T1 * Test_Primary_Diameter2)

    Differential_Scattering_Total1 = []
    Differential_Scattering_Total2 = []
    Ratio_Diff_Total1 = []
    Ratio_Diff_Total2 = []

    for i in range(len(Test_Np)):

        Monomer_Differential1 = Test_Soot_FM * ((Test_Primary_Diameter1 / 2) ** 6) * (Test_Wave_Number1 ** 4) * Test_Np[i]
        Monomer_Differential2 = Test_Soot_FM * ((Test_Primary_Diameter2 / 2) ** 6) * (Test_Wave_Number2 ** 4) * Test_Np[i]
        Differential_Scattering_Phi1 = 0
        Differential_Scattering_Phi2 = 0
        Ratio_Diff1 = []
        Ratio_Diff2 = []

        for t in range(Theta_Number):
            qW1_T1 = FN.Scattering_Wave_Vector(WaveLength_meter=Test_Wave1, Theta_radian=Theta_Radian[t])
            qW2_T1 = FN.Scattering_Wave_Vector(WaveLength_meter=Test_Wave2, Theta_radian=Theta_Radian[t])
            Def_ScatteringW1 = FN.RDG_Def_Scattering(K=Test_Wave_Number1, N=Test_Np[i], Dp=Test_Primary_Diameter1, q=qW1_T1, F=Test_Soot_FM, D_RDG=Test_Soot_Fractal_D_mc_RDG, K_RDG=Test_Soot_Prefactor_k_mc_RDG, Formula=2)
            Def_ScatteringW2 = FN.RDG_Def_Scattering(K=Test_Wave_Number2, N=Test_Np[i], Dp=Test_Primary_Diameter2, q=qW2_T1, F=Test_Soot_FM, D_RDG=Test_Soot_Fractal_D_mc_RDG, K_RDG=Test_Soot_Prefactor_k_mc_RDG, Formula=2)
            Ratio_Diff1.append(Def_ScatteringW1 / Monomer_Differential1)
            Ratio_Diff2.append(Def_ScatteringW2 / Monomer_Differential2)
            Differential_Scattering_Phi1 += FN.Diff_Integral_Phi(Def_ScatteringW1, Phi_Radian, Theta_Radian[t], Theta_Diff, Phi_Diff)
            Differential_Scattering_Phi2 += FN.Diff_Integral_Phi(Def_ScatteringW2, Phi_Radian, Theta_Radian[t], Theta_Diff, Phi_Diff)

            Ratio_Diff_Total1.append(Ratio_Diff1)
            Ratio_Diff_Total2.append(Ratio_Diff2)
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Log(f"Test/W_1Ratio_Diff_Scattering{Test_Np[i]}.jpg", qDp1, "q Dp", Ratio_Diff1, 1, 1000, "Wave Length=" + str(round(Test_Wave1 * 1e9)) + " nm" + "\n" + "N= " + str(Test_Np[i]), "Normalized", "Scattering Cross Section Ratio")
        FN.Fig_Plot_Save_1Lines_X_Log_Y_Log(f"Test/W_2Ratio_Diff_Scattering{Test_Np[i]}.jpg", qDp2, "q Dp", Ratio_Diff2, 10, 1000, "Wave Length=" + str(round(Test_Wave2 * 1e9)) + " nm" + "\n" + "N= " + str(Test_Np[i]), "Normalized", "Scattering Cross Section Ratio")

        Differential_Scattering_Total1.append((10 ** 18) * Differential_Scattering_Phi1)
        Differential_Scattering_Total2.append((10 ** 18) * Differential_Scattering_Phi2)
    a = 3
