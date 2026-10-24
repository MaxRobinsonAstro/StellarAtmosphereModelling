import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use('dark_background')
plt.rcParams.update({
    'font.family' : 'serif',
    'font.serif' : 'cmr10',
    'mathtext.fontset' : 'cm',
    'axes.formatter.use_mathtext' : True,
    'axes.unicode_minus' : False,
    'font.size' : 12 })


"""
opening the fits file, combining tables, and plotting hists of disks and spectral types
"""
a = fits.open('OneDrive - Maynooth University/Project/sample.fit')

dfs = []
for i in range(1,6):
    data = a[i].data
    
    columns = {}
    for name in data.names:
        column = np.array(data[name])
        
        if column.dtype.byteorder not in ('=' , '|'):
            column = column.astype(column.dtype.newbyteorder('='))
            
        columns[name] = column
    
    dfs.append(pd.DataFrame(columns))

df = dfs[0].copy()

df_34 = pd.concat([dfs[3], dfs[4]], axis = 0)

for d in dfs[1:3]:
    d = d.loc[:, ~d.columns.duplicated()]
    df = df.merge(d, on = 'Name', how = 'left')
    
df0 = df.merge(df_34, on = 'Name', how = 'left')

disk = ['Debris', 'Hybrid', 'Protoplanetary', '--']
amount = df0['DiscClass'].value_counts()
plt.figure(figsize=(9,5))
plt.bar(disk, amount, width = 0.65, color = 'tomato', edgecolor = 'white')
plt.xlabel('Disk Classification')
plt.ylabel('Number of Disks')

df = df0[df0['DiscClass'] != 'Protoplanetary'].copy()
df = df0[df0['DiscClass'] !=  'Hybrid'].copy()
df = df0[df0['DiscClass'] !=  '--'].copy()

df['SpType'] = df['SpType'].str[0]
order = ['O', 'B', 'A', 'F', 'G', 'K', 'M']
counts = df['SpType'].value_counts().reindex(order)

plt.figure(figsize=(9,5))
plt.bar(order, counts, width = 0.75, color = 'tomato', edgecolor = 'white')
plt.xlabel('Spectral Type')
plt.ylabel('Number of Stars')

"""
3d plot with distance and orientation
"""

fig = plt.figure(figsize=(9,8))
ax = fig.add_subplot(111, projection = '3d')

ra = np.radians(df['RAJ2000'])
dec = np.radians(df['DEJ2000'])
d = df['Dist']

x = d * np.cos(dec) * np.cos(ra)
y = d * np.cos(dec) * np.sin(ra)
z = d * np.sin(dec)

ax.scatter(x, y, z)
for x_i, y_i, z_i in zip(x, y, z):
    ax.plot([0, x_i], [0, y_i], [0, z_i], alpha = 0.6, linewidth = 0.8)
ax.set_xlabel('$x$ (pc)')
ax.set_ylabel('$y$ (pc)')
ax.set_zlabel('$z$ (pc)')

"""
polar plot with RA and dec with OHP limits and colorbar of m
"""

df['M'] = 4.83 - 2.5 * np.log10(df['L']) + df['Av']
df['m'] = df['M'] -5 + 5 * np.log10(df['Dist'])

fig = plt.figure(figsize=(9,8))
ax = fig.add_subplot(111, projection = 'polar')

j = ax.scatter(df['RAJ2000'] * np.pi / 180, df['DEJ2000'], c = df['m'], cmap = 'nipy_spectral_r')
ax.set_xlabel('RA')
ax.set_ylabel('DEC')
f = 43.9317 + (25 - 90)
ax.axhline(f, c = 'tomato', linestyle = '--')

cbar = plt.colorbar(j)
cbar.set_label('Apparent Magnitude')

angs = np.arange(0, 360, 7.5)
ax.set_xticks(angs * np.pi /180)
xlabels = ['0h', '', '1h', '', '2h', '', '3h', '', '4h', '', '5h', '', '6h', '', '7h', '', '8h', '', 
           '9h', '', '10h', '', '11h', '', '12h', '', '13h', '', '14h', '', '15h', '', '16h', '', 
           '17h', '', '18h', '', '19h', '', '20h', '', '21h', '', '22h', '', '23h', '']
ax.set_xticklabels(xlabels)
ohp_min = 22.5 * 15
ohp_max = 8.5  * 15

ax.axvline(ohp_min * np.pi / 180, c = 'tomato', linestyle = '--')
ax.axvline(ohp_max * np.pi / 180, c = 'tomato', linestyle = '--')
ax.set_rlabel_position(0)
ax.set_ylim(-90, 45)
ax.grid(alpha = 0.7)

"""
OHP stars only
""" 

df = df[df['DEJ2000'] >= f].copy()
df = df[(df['RAJ2000'] >= ohp_min) | (df['RAJ2000'] <= ohp_max)].copy()

fig = plt.figure(figsize=(9,8))
ax = fig.add_subplot(111, projection = 'polar')

j = ax.scatter(df['RAJ2000'] * np.pi / 180, df['DEJ2000'], c = df['m'], cmap = 'nipy_spectral_r')
ax.set_xlabel('RA')
ax.set_ylabel('DEC')

cbar = plt.colorbar(j)
cbar.set_label('Apparent Magnitude')


ax.set_xticks(angs * np.pi /180)
ax.set_xticklabels(xlabels)
ax.set_rlabel_position(0)
ax.set_ylim(-90, 45)
ax.axvspan(ohp_min * np.pi / 180, ohp_max * np.pi / 180, color = 'maroon')
ax.axhspan(-90,f, color = 'maroon')
ax.grid(alpha = 0.7)

counts = df['SpType'].value_counts().reindex(order)

plt.figure(figsize=(9,5))
plt.bar(order, counts, width = 0.75, color = 'tomato', edgecolor = 'white')
plt.xlabel('Spectral Type')
plt.ylabel('Number of Stars')


plt.show()