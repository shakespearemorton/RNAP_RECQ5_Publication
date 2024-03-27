import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import openmm.app as app
import analyticalrdf
import glob
import pandas as pd
from progressbar import progressbar
from lxml import etree


def load(fname):
    tree = etree.parse(fname)
    root = tree.getroot()
    data = []

    for marker in root.iter('marker'):
        # Extract the relevant attributes
        id_ = int(marker.attrib['id'])
        x = float(marker.attrib['x'])
        y = float(marker.attrib['y'])
        z = float(marker.attrib['z'])
        coordX = float(marker.attrib['coordX'].rstrip('px'))
        coordY = float(marker.attrib['coordY'].rstrip('px'))
        coordZ = float(marker.attrib['coordZ'].rstrip('px'))
        data.append([id_, x, y, z, coordX, coordY, coordZ])

    df = pd.DataFrame(data, columns=['id', 'x', 'y', 'z', 'coordX', 'coordY', 'coordZ'])
    df = center(df)
    return df[['x', 'y', 'z']].values

def center(df):
    mean_x = df['x'].mean()
    mean_y = df['y'].mean()
    mean_z = df['z'].mean()
    df_centered = df.copy()
    df_centered['x'] = df['x'] - mean_x
    df_centered['y'] = df['y'] - mean_y
    df_centered['z'] = df['z'] - mean_z

    return df_centered

plt.figure(dpi=1200)

files = glob.glob('../cryo-ET_data/*.cmm')
exp = []
for file in files:
    points = load(file)/10
    r = np.arange(0,100,0.1)
    gr = analyticalrdf.RDF_AnalyticNorm(points,r,0.1)
    for i in range(len(r)):
        if gr[-1-i] > 0:
            cut = r[-1-i]/2
            break
    V = (4/3) * np.pi * cut**3
    rho = len(points) / V
    rsp = (4*np.pi*r**2)
    Fsphere = rsp/(rsp*(1-(r/(2*cut))**2)*(1+r/(4*cut)))
    for i in range(len(r)):
        if r[i] > cut:
            imax = i
            break
    sgr = gr*Fsphere
    exp.append(sgr[1:-1])
stds = np.std(exp,axis=0)
rdfs = np.mean(exp,axis=0)
plt.plot(r[:imax],rdfs[:imax], linewidth=3,label='RNAPII - exp',c='rosybrown')
plt.fill_between(r[:imax],rdfs[:imax]+stds[:imax],rdfs[:imax]-stds[:imax],alpha=0.3,color='rosybrown')
x = [r[np.argmax(i[:imax])] for i in exp]
print(np.min(x[x!=0]))

