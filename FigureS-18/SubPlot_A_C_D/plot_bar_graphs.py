import numpy as np
import matplotlib.pyplot as plt

#names = ['1:3','1:5','1:10','50','60','70']
#idr = [4.75,6.39,8.10,7.73,7.74,7.51]
#idr_std = [2.57,2.77,3.19,3.26,3.19,3.33]
#sri = [4.27,6.68,9.1,6.99,6.91,6.99]
#sri_std = [2.46,2,3.12,2.34,2.53,2.6]
names = ['1:5','1:7','1:8','1:9','30','50','60']
com = [11.33,11.92,12.93,13.24,12.34,12.37,12.15]
com_std = [0.01,0.02,0.02,0.03,0.01,0.001,0.01]
idr = [10.9,10.7,12.7,10.64,13.08,12.22,12.02]
idr_std = [3.84,4.39,5.98,4.69,4.64,3.83,4.84]
sri = [6.2,10.1,10.15,11.9,9.3,8.6,9.0]
sri_std = [2.8,3.86,2.18,2.80,2.78,3.90,2.99]
r5 = 7.8
r5_std = 2.7

fig, ax = plt.subplots(dpi = 1200)
plt.rcParams['font.family'] = 'serif'  # Set global font family to serif
plt.rcParams['font.serif'] = 'DejaVu Serif'  # Specific serif font
plt.rcParams['axes.labelsize'] = 15  # Axis label size
plt.rcParams['xtick.labelsize'] = 12  # X tick label size
plt.rcParams['ytick.labelsize'] = 12
plt.grid(False)
# Adding bars for idr
idr_bars = ax.bar(np.arange(len(names)), idr, yerr=idr_std, capsize=5, label='IDR', alpha=0.6)


# Highlight the background for r5 std deviation
#plt.axhspan(r5 - r5_std, r5 + r5_std, facecolor='orange', alpha=0.3)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('IDR Valency')
ax.set_xticks(np.arange(len(names)))
ax.set_xticklabels(names)
for i, label in enumerate(names):

    if label in [ '50']:

        ax.text(i, -1.5, 'Copies', ha='center', va='top')
        
    if label in [ '1:7']:

        ax.text(i+0.5, -1.5, 'Ratio', ha='center', va='top')
plt. tight_layout()
plt.savefig('idr_valency.pdf',format='pdf')

fig, ax = plt.subplots(dpi = 1200)
plt.grid(False)
# Adding bars for idr
idr_bars = ax.bar(np.arange(len(names)), sri, yerr=sri_std, capsize=5, color='red', alpha=0.6)


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('CTD Valency')
ax.set_xticks(np.arange(len(names)))
ax.set_xticklabels(names)
for i, label in enumerate(names):

    if label in [ '50']:

        ax.text(i, -1.5, 'Copies', ha='center', va='top')
        
    if label in [ '1:7']:

        ax.text(i+0.5, -1.5, 'Ratio', ha='center', va='top')
plt. tight_layout()
plt.savefig('sri_valency.pdf',format='pdf')

fig, ax = plt.subplots(dpi = 1200)
plt.grid(False)
idr_bars = ax.bar(np.arange(len(names)), com, yerr=com_std, capsize=5, color='purple', alpha=0.6)


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('CoM Separation Distance (nm)')
ax.set_xticks(np.arange(len(names)))
ax.set_xticklabels(names)
plt.axhline(12.545, c='k',linestyle = '--',label='Experimental (12.5)')
for i, label in enumerate(names):

    if label in [ '50']:

        ax.text(i, -1.5+11, 'Copies', ha='center', va='top')
        
    if label in [ '1:7']:

        ax.text(i+0.5, -1.5+11, 'Ratio', ha='center', va='top')
plt.legend()
plt.ylim([10,14])
plt. tight_layout()
plt.savefig('com_sep.pdf',format='pdf')
