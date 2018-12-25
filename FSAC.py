# Programmed by Keyhan Babaee Under Prof. Steven Rogak supervision, https://github.com/KeyhanB
# Version 0.3
# Nov 12, 2018
import Functions as KN
import matplotlib.pyplot as plt
import numpy as np
import os
from math import pi
from math import log
from math import exp

Soot_Refractive_Index = 1.95-0.79j# Complex refractive index

Soot_Primary_Agg_dp_sigma = 1  # Geometric std of primary particles within agg
Soot_Fractal_D_mc = 0.52  # Fractal Properties for continuum regime
Soot_Prefactor_k_mc = 0.85+(0.03 * Soot_Primary_Agg_dp_sigma ** (4.4))  # Fractal free molecular regime
Soot_Fractal_D_mfm = 0.46  # Fractal Properties free molecular regime
Soot_Prefactor_k_mfm = 0.94+(0.03 * Soot_Primary_Agg_dp_sigma ** (4.8))  # Fractal free molecular regime
eff_k = 0.418  # effective density variable for the result to be in kg/m^3
eff_dm = 2.56  # effective density variable for the result to be in kg/m^3
Diameter_Primary_100nm = 16.48 * 1e-9  # primary diameter of 100 nm aggregate
D_TEM = 0.45
Primary_Particle_Sigma_da_cte = 1.33
Primary_Particle_Sigma_da_cte_How_Many = 2
Primary_Particle_Sigma_da_cte_number_of_points = 19
Betha = 0.9


def Primary_Particle_Size_meter (da_meter, dp100_meter, Dtem):
    A = dp100_meter * ((da_meter * 1e9 / 100) ** Dtem)
    return A


def Number_Primary_Particle (da_meter, dp_meter, Continuum=True):
    if Continuum == True:
        A = ((da_meter / dp_meter) / Soot_Prefactor_k_mc) ** (1 / Soot_Fractal_D_mc)
        return round (A)
    else:
        A = ((da_meter / dp_meter) / Soot_Prefactor_k_mfm) ** (1 / Soot_Fractal_D_mfm)
        return round (A)


def LogNormal_Generator (dp_median_meter, sigma_g, Howmany_Sigma_g, da_meter, number):
    B = Howmany_Sigma_g
    dmax = dp_median_meter * (sigma_g ** B)
    dmin = dp_median_meter * (sigma_g ** (-1 * B))
    diam1 = []
    diam2 = []
    D_ratio = (dmax / dp_median_meter) ** (1 / (number-1))
    for i in range (number):
        d1 = dp_median_meter * D_ratio ** (i)
        diam1.append (d1)
    D_ratio = (dp_median_meter / dmin) ** (1 / (number-1))
    for i in range (number-1):
        d2 = dmin * D_ratio ** (i)
        diam2.append (d2)
    Diam_total = diam2+diam1
    prob = []
    number_primary = []
    sum = 0

    for i in range (len (Diam_total)-1):
        number_primary.append (Number_Primary_Particle (da_meter, Diam_total[i]))
    for i in range (len (Diam_total)-1):
        prob.append ((1 / (log (sigma_g) * (2 * pi) ** 0.5)) * exp (-1 * (((log (Diam_total[i])-log (dp_median_meter)) ** 2) / (2 * (log (sigma_g)) ** 2))) * (log (Diam_total[i+1])-log (Diam_total[i])))
        sum += prob[i]

    return Diam_total[:-1], number_primary, prob, sum


# a,b,c,d=LogNormal_Generator(20e-9,Primary_Particle_Sigma_da_cte,2,100e-9,19)
# Soot Primary Particle Distribution Properties (Log Normal)
# Soot_Primary_Sample_dp_median = 17.4 * 1e-9  # nm median diameter
# Soot_Primary_Sample_dp_sigma = 1  # Geometric std
# Soot_Primary_Agg_dp_median = 17.4 * 1e-9  # nm median diameter
# Soot_Primary_Rho = 1770
Soot_Complex = (((Soot_Refractive_Index ** 2)-1) / ((Soot_Refractive_Index ** 2)+2))
Soot_FM = (abs (Soot_Complex)) ** 2
Soot_EM = Soot_Complex.imag
Wave_Length = 550 * 1e-9  # nm
Wave_Number = 2 * pi / Wave_Length

Total_Number_Concentration_Sample = 1e4  # /cm^3
LogNormal_D_Median_Small_Sample = 100e-9  # smallest diameter of computation (100nm)
LogNormal_D_Median_Large_Sample = 800e-9  # largest diameter of computation (800nm)
LogNormal_D_Median_Number_Sample = 3  # number of the computation
LogNormal_D_Median_Sample = np.linspace (LogNormal_D_Median_Small_Sample, LogNormal_D_Median_Large_Sample, num=LogNormal_D_Median_Number_Sample)

LogNormal_Sigma_min_Sample = 1.05
LogNormal_Sigma_max_Sample = 2.5
LogNormal_Sigma_Number_Sample = 3
LogNormal_Sigma_Sample = np.linspace (LogNormal_Sigma_min_Sample, LogNormal_Sigma_max_Sample, num=LogNormal_Sigma_Number_Sample)

Theta_Number = 48
Theta_Radian = np.linspace (0, 2 * pi, num=Theta_Number)

diam = []
Ln_Diam1 = []
Diameter_Nano = []
D_Small = 10e-9  # smallest diameter of computation (10nm)
D_Large = 1e-6  # largest diameter of computation (10um)
Nd = 199
Keyhan = 10 ** 9  # Conversion factor
D_ratio = (D_Large / D_Small) ** (1 / (Nd-1))  # Generating particle diameters
for i in range (Nd):
    d = D_Small * D_ratio ** (i)
    diam.append (d)
    Ln_Diam1.append (log (d))
    Diameter_Nano.append (d * Keyhan)

LogNormal_pdf = [[[0 for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Soot_Primary_Diameter_Median_meter = [[[0 for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Soot_Primary_Diameter_Number_Probability = [[[0 for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Absorption_Cross = [[[0 for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Absorption_Cross_Eff = [[[0 for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Absorption_Cross_Primary = [[[0 for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Absorption_Cross_Primary_Eff = [[[0 for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Total_Absorption_Cross = [[0 for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Scattering_Cross = [[[0 for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Scattering_Cross_Primary = [[[0 for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Total_Scattering_Cross = [[0 for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Differential_Scattering_Cross = [[[[0 for p in range (Theta_Number)] for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Aggregate_Mass = [[[0 for k in range (Nd)] for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Total_Particle_Mass = [[0 for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]
Mass_Absorption_Eff = [[0 for j in range (LogNormal_Sigma_Number_Sample)] for i in range (LogNormal_D_Median_Number_Sample)]

n = 1

for i in range (LogNormal_D_Median_Number_Sample):
    for j in range (LogNormal_Sigma_Number_Sample):
        for k in range (Nd-1):
            LogNormal_pdf[i][j][k] = KN.LogNormal (LogNormal_D_Median_Sample[i], LogNormal_Sigma_Sample[j], diam[k], diam[k+1])
            Soot_Primary_Diameter_Median_meter[i][j][k] = Primary_Particle_Size_meter (diam[k], Diameter_Primary_100nm, D_TEM)
            Primary_Diameter, Primary_Number, Primary_Probability, Okay = LogNormal_Generator (Soot_Primary_Diameter_Median_meter[i][j][k], Primary_Particle_Sigma_da_cte, Primary_Particle_Sigma_da_cte_How_Many, diam[k],
                                                                                               Primary_Particle_Sigma_da_cte_number_of_points)
            Absorption_Cross_Section_Agg, Absorption_Cross_Section_Primary = KN.RDG_Absorption (Wave_Number, Primary_Number, Primary_Diameter, Soot_EM)  # within  aggregate
            Absorption_Cross_Eff_1, Absorption_Cross_Primary_Eff_1 = KN.Absorption_Eff (Absorption_Cross_Section_Agg, Absorption_Cross_Section_Primary, diam[k], Primary_Diameter)
            Absorption_Cross_Section_Agg_Total = np.sum (Absorption_Cross_Section_Agg * Primary_Probability)
            Absorption_Cross_Section_Primary_Total = np.sum (Absorption_Cross_Section_Primary * Primary_Probability)
            Absorption_Cross_Eff_Total = np.sum (Absorption_Cross_Eff_1 * Primary_Probability)
            Absorption_Cross_Primary_Eff_Total = np.sum (Absorption_Cross_Primary_Eff_1 * Primary_Probability)
            Scattering_Cross_Section_Agg, Scattering_Cross_Section_Primary = KN.RDG_Total_Scatteing (Wave_Number, Soot_FM, Primary_Number, Primary_Diameter, Betha)  # within  aggregate
            Scattering_Cross_Section_Agg_Total = np.sum (Scattering_Cross_Section_Agg * Primary_Probability)
            Scattering_Cross_Section_Primary_Total = np.sum (Scattering_Cross_Section_Primary * Primary_Probability)
            Aggregate_Mass[i][j][k] = KN.Particle_Mass (KN.Effective_Density (eff_k, eff_dm, diam[k]), diam[k])
            # for p in range (Theta_Number):
            #   Differential_Scattering_Cross[i][j][k][p] = KN.RDG_Def_Scattering (Wave_Number, Theta_Radian[p], Soot_FM, Soot_Primary_Number[i][j][k], Soot_Primary_dp_median, Soot_Fractal_Df, Soot_Prefactor_kf)
            Total_Absorption_Cross[i][j] = Total_Absorption_Cross[i][j]+LogNormal_pdf[i][j][k] * Absorption_Cross_Section_Agg_Total*Total_Number_Concentration_Sample
            Total_Scattering_Cross[i][j] = Total_Scattering_Cross[i][j]+LogNormal_pdf[i][j][k] * Scattering_Cross_Section_Agg_Total*Total_Number_Concentration_Sample
            Total_Particle_Mass[i][j] = Total_Particle_Mass[i][j]+LogNormal_pdf[i][j][k] * Aggregate_Mass[i][j][k]*Total_Number_Concentration_Sample

            Absorption_Cross[i][j][k]=Absorption_Cross_Section_Agg_Total
            Absorption_Cross_Eff[i][j][k]=Absorption_Cross_Eff_Total
            Absorption_Cross_Primary[i][j][k]=Absorption_Cross_Section_Primary_Total
            Absorption_Cross_Primary_Eff[i][j][k]=Absorption_Cross_Primary_Eff_Total
            Scattering_Cross[i][j][k]=Scattering_Cross_Section_Agg_Total
            Scattering_Cross_Primary[i][j][k]=Scattering_Cross_Section_Primary_Total
        Mass_Absorption_Eff[i][j] = Total_Absorption_Cross[i][j] / Total_Particle_Mass[i][j]
        plt.figure (n)
        plt.rcParams['mathtext.fontset'] = 'stix'
        plt.rcParams['font.family'] = 'STIXGeneral'
        plt.title ("Absorption and Scattering, "+"D_Median= "+str (LogNormal_D_Median_Sample[i] * 1e9)+"(nm), Sigma="+str (LogNormal_Sigma_Sample[j])+", Mass Abs Eff= "+"{:.3E}".format (Mass_Absorption_Eff[i][j])+" (m^2/kg)")

        plt.plot (Diameter_Nano, Absorption_Cross[i][j][:], label="Tot Abs= "+"{:.3E}".format (Total_Absorption_Cross[i][j])+"(m^2)")
        plt.plot (Diameter_Nano, Aggregate_Mass[i][j][:], label="Total Particle Mass (Kg)= "+"{:.3E}".format (Total_Particle_Mass[i][j]))
        plt.plot (Diameter_Nano, Absorption_Cross_Primary[i][j][:], label="Primary Particle Abs")
        plt.plot (Diameter_Nano, Scattering_Cross[i][j][:], label="Tot Sca= "+"{:.3E}".format (Total_Scattering_Cross[i][j])+"(m^2)")
        plt.plot (Diameter_Nano, Scattering_Cross_Primary[i][j][:], label="Primary Particle Sca")
        plt.xlabel ('Mobility-equivalent Diameter (nm)')
        plt.ylabel ('Cross Section (m^2)')
        plt.xscale ('log')
        plt.yscale ('log')
        plt.grid (True, which='minor', alpha=0.2)
        plt.legend (bbox_to_anchor=(1.04, 0.5), loc='center left', fancybox=True, fontsize='small')
        GraphName = "D_Med-"+str (i)+"Sigma"+str (j)+".jpg"
        plt.savefig (GraphName, format='jpg', dpi=1000, bbox_inches='tight')
        plt.clf ()
        plt.close ()
        n = n+1

        plt.figure (n)
        plt.rcParams['mathtext.fontset'] = 'stix'
        plt.rcParams['font.family'] = 'STIXGeneral'
        plt.title ("Absorption Efficiency, "+"D_Median= "+str (LogNormal_D_Median_Sample[i] * 1e9)+"(nm), Sigma="+str (LogNormal_Sigma_Sample[j])+", Mass Abs Eff= "+"{:.3E}".format (Mass_Absorption_Eff[i][j])+" (m^2/kg)")
        plt.plot (Diameter_Nano, Absorption_Cross_Eff[i][j][:], label="Aggregate")
        plt.plot (Diameter_Nano, Absorption_Cross_Primary_Eff[i][j][:], label="Primary")
        plt.xlabel ('Mobility-equivalent Diameter (nm)')
        plt.ylabel ('Efficiency')
        plt.xscale ('log')
        plt.yscale ('log')
        plt.grid (True, which='minor', alpha=0.2)
        plt.legend (bbox_to_anchor=(1.04, 0.5), loc='center left', fancybox=True, fontsize='small')
        GraphName = "D_Med-"+str (i)+"Sigma"+str (j)+"_EFF"+".jpg"
        plt.savefig (GraphName, format='jpg', dpi=1000, bbox_inches='tight')
        plt.clf ()
        plt.close ()
        n = n+1

a = 3
