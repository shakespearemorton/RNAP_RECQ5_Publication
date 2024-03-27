import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy import integrate
from scipy.interpolate import UnivariateSpline

file_paths = ['one_fe.txt','two_fe.txt','two_p_twoseven_fe.txt']#'og_fe.txt','one_fe.txt',,'three_fe.txt']
#file_paths = ['hist_fe_275.txt']
labels = ['+1.00 Charge : 1.2 Mol', '+2.00 Charge : 2.6e-3 Mol','+2.75 Charge : 4.2e-5 Mol']
use = 0
plt.figure(figsize=(8, 6), dpi=1200)  # Adjust size and resolution
plt.rcParams['font.family'] = 'serif'  # Set global font family to serif
plt.rcParams['font.serif'] = 'DejaVu Serif'  # Specific serif font
plt.rcParams['axes.labelsize'] = 15  # Axis label size
plt.rcParams['xtick.labelsize'] = 12  # X tick label size
plt.rcParams['ytick.labelsize'] = 12
# Creating KDE plots
for file_path in file_paths:

    # Read the file and display the first few lines to understand its structure
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    
    # Resetting the lists to avoid appending the same data again
    coor, free_energy, free_energy_error, prob, prob_error = [], [], [], [], []
    
    # Processing the file and skipping lines not containing the main data
    for line in lines[1:]:  
        if line.startswith("#"):  # Skipping any lines starting with '#'
            continue
    
        entries = line.split()
        coor.append(float(entries[0]))
        free_energy.append(float('inf') if entries[1] == 'inf' else float(entries[1]))
        free_energy_error.append(float('nan') if entries[2] == 'nan' else float(entries[2]))
    
    # Convert lists to numpy arrays for easier handling
    coor = np.array(coor)
    free_energy = np.array(free_energy) #- free_energy[-1]
    free_energy[free_energy == np.inf] = np.nan
    valid_indices = ~np.isnan(free_energy) & ~np.isinf(free_energy)
    filtered_coor = coor[valid_indices]
    filtered_free_energy = free_energy[valid_indices]
    plt.scatter(filtered_coor, filtered_free_energy - filtered_free_energy[-1],label = labels[use])
    plt.xlabel('CoM Separation Distance (nm)')
    plt.ylabel('Free Energy (kcal/mol)')
    plt.grid(True)
    use +=1
    
    # Create the plot for Free Energy
    
    
    bulk_index = np.abs(filtered_coor - 4).argmin()
    bulk_free_energy = filtered_free_energy[bulk_index]
    favorable_indices = (filtered_coor <= 4) & (filtered_free_energy < bulk_free_energy)
    favorable_coor = filtered_coor[favorable_indices]
    favorable_free_energy = filtered_free_energy[favorable_indices]
    
    x = favorable_coor
    y = favorable_free_energy

    
    x_range = np.linspace(np.min(x),np.max(x),1025)
    spline = UnivariateSpline(x, y, s=5)
    
    # Calculate the spline over the same range as before
    spline_free_energy = spline(x_range)
    #plt.plot(x_range, spline_free_energy, '--', label='Univariate Spline', color='purple')
    
    int_spline, int_spline_error = integrate.quad(spline, np.min(x), np.max(x))
    
    T = 300
    R = 1.98720425864083e-3
    x = x_range
    y = spline_free_energy
    def integrand(x):
        #We define ΔG(x) = G(x) − G(0), where x = 0 (that is, the center of the binding pocket) is defined as the grid point associated with the lowest grid PMF.
        gx = spline(x) - np.min(y)
        return 4 * np.pi * x**2 * np.exp(-gx / (R * T))
    # Show the section that will be integrated
    #plt.scatter(x_range,linear_interp(x_range))
    
    ### BINDING AFFINITY CALCULATION ###
    
    #ΔG_P = −RT ln(V_P / V°) = −RT ln ∫_(pocket) e^(-ΔG(x) / RT) dV / A^3
    gp, gp_error = -R*T*np.log(integrate.quad(integrand, np.min(x), np.max(x)))
    # Continue with your calculations as before
    g_xb = y[-1] - np.min(y)
    gb = -7.42*R*T #ΔG_B = −RT ln(V_B / V°) ≈ −7.42RT (J/mol)
    gv = gp - gb # ΔG_V = ΔG_P - ΔG_B
    dg = -g_xb + gv
    kd = np.exp(dg/(R*T))
    print(kd, ' Mol')
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('fe.pdf',format='pdf')