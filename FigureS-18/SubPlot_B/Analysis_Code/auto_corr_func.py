import numpy as np
import matplotlib.pyplot as plt

#data1 = np.loadtxt('iri_five')
data2 = np.loadtxt('sixty_sri_valency.txt')
data3 = np.loadtxt('com.txt')
legends = ['IDR Valency', 'CTD Valency', 'Com Separation Distance (nm)']

def calculate_delta_A(A):
    mean_A = np.mean(A)
    delta_A = A - mean_A
    return delta_A

def autocorrelation_function(A):
    N = len(A)
    delta_A = calculate_delta_A(A)
    C_A = np.array([np.sum(delta_A[:N-i] * delta_A[i:]) for i in range(N)]) / np.sum(delta_A**2)
    return C_A

# Calculate autocorrelation function

plt.rcParams['font.family'] = 'serif'  # Set global font family to serif
plt.rcParams['font.serif'] = 'DejaVu Serif'  # Specific serif font
plt.rcParams['axes.labelsize'] = 15  # Axis label size
plt.rcParams['xtick.labelsize'] = 12  # X tick label size
plt.rcParams['ytick.labelsize'] = 12
plt.figure(dpi = 1200)
plt.grid(True)

datum = [data2,data3]
count = 0
for i in datum:
    x = np.asarray(range(len(i)))/100
    C_A = autocorrelation_function(i)
    if count == 0:
        plt.plot(x,C_A, label=legends[count],linewidth = 4)
    elif count == 1:
        plt.plot(x,C_A, label=legends[count],linewidth = 4,color = 'red')
    else:
        plt.plot(x,C_A, label=legends[count],linewidth = 4,color='purple')
    count +=1

# Plotting
plt.fill_between(x,0.37, -0.37,alpha=0.2,color='darkgoldenrod')
plt.xlabel('t ($\\mu s$)')
plt.ylabel('$C_A(t)$')
plt.xlim([0,10])
plt.legend()
plt. tight_layout()
plt.savefig('Equilibration.pdf',format='pdf')
