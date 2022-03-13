import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np 
from matplotlib.backends.backend_pdf import PdfPages

mpl.rcParams['font.size'] = 10
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'large'
mpl.rcParams["font.family"] = "serif"
csfont = {'fontname':'Times New Roman'}
plt.rcParams['mathtext.fontset']='dejavuserif'

input_data_path = "/Users/jeager/Documents/PhD_work/ekman_ocean/input_data/ProCb"
output_data_path = "/Users/jeager/Documents/PhD_work/ekman_ocean/output_data/ProCb"
plot_path = "/Users/jeager/Documents/PhD_work/ekman_ocean/plots"

lats_file = "lats.dat"
lons_file = "lons.dat"
os.chdir(input_data_path)
lats_index = np.loadtxt(lats_file, dtype=int, delimiter="\t",usecols=[0])
lats_data = np.loadtxt(lats_file, delimiter="\t",usecols=[1])
lons_index = np.loadtxt(lons_file, dtype=int, delimiter="\t",usecols=[0])
lons_data = np.loadtxt(lons_file, delimiter="\t",usecols=[1])

T_init_file = "surface_temperature.dat"
T_init_data = np.loadtxt(T_init_file, delimiter="\t",usecols=lons_index)

os.chdir(output_data_path)
T_surf_file = "T_surf_10000_days.dat"
T_deep_file = "T_deep_10000_days.dat"
T_surf_data = np.loadtxt(T_surf_file, delimiter="\t",usecols=lons_index)
T_deep_data = np.loadtxt(T_deep_file, delimiter="\t",usecols=lons_index)
surf_diff = T_surf_data - T_init_data

fig = plt.figure(figsize=(12,10))
ax = plt.subplot(2,2,1)

c_levels = np.arange(120.,301.,10)
cntr = ax.contourf(lons_data,lats_data,T_init_data,
    levels=c_levels,
	)
ax.set_ylabel(r'Latitude ($\degree$)')
ax.set_xlabel(r'Longitude ($\degree$)')
ax.set_title("No Transport Slab Ocean")
cbar = fig.colorbar(cntr)
cbar.set_label('Temperature (K)')

ax = plt.subplot(2,2,2)

c_levels = np.arange(120.,301.,10)
cntr = ax.contourf(lons_data,lats_data,T_surf_data,
    levels=c_levels,
    )
ax.set_ylabel(r'Latitude ($\degree$)')
ax.set_xlabel(r'Longitude ($\degree$)')
ax.set_title("Surface Ocean with transport")
cbar = fig.colorbar(cntr)
cbar.set_label('Temperature (K)')

ax = plt.subplot(2,2,3)

c_levels = np.arange(120.,301.,10)
cntr = ax.contourf(lons_data,lats_data,T_deep_data,
    levels=c_levels,
    )
ax.set_ylabel(r'Latitude ($\degree$)')
ax.set_xlabel(r'Longitude ($\degree$)')
ax.set_title("Deep Ocean with transport")
cbar = fig.colorbar(cntr)


ax = plt.subplot(2,2,4)

diff = T_surf_data - T_init_data
c_levels = np.arange(-17.5,17.6,2.5)
# c_levels = np.arange(190.,331.,5)
cntr = ax.contourf(lons_data,lats_data,diff,
    levels=c_levels,
    cmap=plt.cm.RdBu_r
    )
ax.set_ylabel(r'Latitude ($\degree$)')
ax.set_xlabel(r'Longitude ($\degree$)')
ax.set_title("Surface Ocean T change")
cbar = fig.colorbar(cntr)
cbar.set_label('Temperature (K)')

plt.tight_layout()
plt.show()
os.chdir(plot_path)
fig.savefig('temperature_maps.pdf')
