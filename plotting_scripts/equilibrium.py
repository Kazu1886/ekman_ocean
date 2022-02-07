import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np 
from matplotlib.backends.backend_pdf import PdfPages

mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'large'
mpl.rcParams["font.family"] = "serif"
csfont = {'fontname':'Times New Roman'}
plt.rcParams['mathtext.fontset']='dejavuserif'

input_data_path = "/Users/jeager/Documents/PhD_work/codron_ocean/input_data/ProCb"
output_data_path = "/Users/jeager/Documents/PhD_work/codron_ocean/output_data/ProCb"
plot_path = "/Users/jeager/Documents/PhD_work/codron_ocean/plots"
lats_file = "lats.dat"
lons_file = "lons.dat"
T_surf_data = []
T_deep_data = []

t_num=151
t_step=100
t = []
times = []
for i in range(0,t_num):
    time=i*t_step
    t+=[time]
    times+=[str(time)]

os.chdir(input_data_path)
lats_index = np.loadtxt(lats_file, dtype=int, delimiter="\t",usecols=[0])
lats_data = np.loadtxt(lats_file, delimiter="\t",usecols=[1])
lons_index = np.loadtxt(lons_file, dtype=int, delimiter="\t",usecols=[0])
lons_data = np.loadtxt(lons_file, delimiter="\t",usecols=[1])

for time in times:
    if time=='0':
        T_surf_file = "surface_temperature.dat"
        T_deep_file = "surface_temperature.dat"
        os.chdir(input_data_path)
    else:
        T_surf_file = "T_surf_"+time+"_days.dat"
        T_deep_file = "T_deep_"+time+"_days.dat"
        os.chdir(output_data_path)
    T_surf = np.loadtxt(T_surf_file, delimiter="\t",usecols=lons_index)
    T_deep = np.loadtxt(T_deep_file, delimiter="\t",usecols=lons_index)
    T_surf = np.mean(T_surf,axis=1)
    T_surf = np.average(T_surf,axis=0,weights=np.cos(np.pi/180.*lats_data))
    T_deep = np.mean(T_deep,axis=1)
    T_deep = np.average(T_deep,axis=0,weights=np.cos(np.pi/180.*lats_data))
    T_surf_data+=[T_surf]
    T_deep_data+=[T_deep]

fig = plt.figure(1, figsize=(10,5))
ax = fig.add_subplot(1,1,1)

plt.axhline(y=T_surf_data[0], color='gray', linestyle=':')
ax.plot(t,T_surf_data,'-',label='Surface')
ax.plot(t,T_deep_data,'--',label='Deep')
ax.grid(linestyle=':', alpha=0.5)
ax.set_xlim(0.0,15000.)
ax.set_ylabel('Global Average Temperature (K)')
ax.set_xlabel('Time (days)')
ax.legend(loc='center right', fontsize=10)

plt.tight_layout()
plt.show()
os.chdir(plot_path)
fig.savefig('equilibrium.pdf')
