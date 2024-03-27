import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Assuming your data is in 'data.txt'
data = pd.read_csv('og.txt', sep='\t', header=None, names=['Step', 'r0', 'Value'])

# Grouping the data by 'r0'
grouped = data.groupby('r0')

# Setting up the style of the plot
plt.figure(figsize=(8, 6), dpi=1200)  # Adjust size and resolution
plt.rcParams['font.family'] = 'serif'  # Set global font family to serif
plt.rcParams['font.serif'] = 'DejaVu Serif'  # Specific serif font
plt.rcParams['axes.labelsize'] = 15  # Axis label size
plt.rcParams['xtick.labelsize'] = 12  # X tick label size
plt.rcParams['ytick.labelsize'] = 12
# Creating KDE plots
for r0, group in grouped:
    sns.kdeplot(group['Value'], shade=True)
plt.grid(False)
plt.xlabel('CoM Separation Distance (nm)')
plt.ylabel('Density')
plt.tight_layout()
plt.savefig('umbrellas.pdf',format='pdf')
