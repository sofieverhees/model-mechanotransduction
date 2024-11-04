import matplotlib.pyplot as plt
import numpy as np
from extra_functions import plot
from matplotlib import rcParams
import numpy as np
import pyvista as pv

myblue = (0,119/235,187/235)
myred=(187/235,85/235,102/235)
myyellow=(221/235,170/235,51/235)
mygrey=(187/235,187/235,187/235)
mygreen="#66BB55"
mymagenta="#7733DD"
myblack="#48494B"

axes_labels = 25
axes_ticks = 20
legend_labels = 20
#rcParams.update({'figure.autolayout': True})
rcParams['font.size'] = 20
#rcParams['xtick.labelbottom'] = False
#rcParams['ytick.labelleft'] = False

dt = 0.05
T = 10
D1 = 100
Es = np.array([0.001, 0.01, 0.1, 0.3, 1.0, 5.7, 10.0, 50.0, 100.0, 7000000.0])
n = len(Es)

for meshnr in [1]:

        ca = np.zeros(n)
        ca_min = np.zeros(n)
        ca_max = np.zeros(n)
        pa = np.zeros(n)
        pa_min = np.zeros(n)
        pa_max = np.zeros(n)
        ca_2D = np.zeros(n)
        ca_min_2D = np.zeros(n)
        ca_max_2D = np.zeros(n)
        pa_2D = np.zeros(n)
        pa_min_2D = np.zeros(n)
        pa_max_2D = np.zeros(n)
        pa_2xD = np.zeros(n)
        pa_min_2xD = np.zeros(n)
        pa_max_2xD = np.zeros(n)
        i = 0

        for E in Es:
            file_ca = pv.read("../../../new_results/simulations/rhomodel3D_ca_reduced_dt="+str(dt)+"_T="+str(T)+"_D="+str(D1)+"_"+str(E)+"E000199.vtu")
            cell_ca = file_ca.point_data_to_cell_data()
            ca_values = cell_ca.get_array(cell_ca.array_names[0])

            sized = file_ca.compute_cell_sizes(length=False, area=False)
            cell_volumes = np.abs(sized.cell_data["Volume"])
            tot_vol = np.sum(cell_volumes)
            ca[i] = np.sum(ca_values*cell_volumes) / tot_vol
            #print('mean ca', np.sum(ca_values*cell_volumes) / tot_vol)

            file_ca_2D = pv.read("../../../new_results/simulations/rhomodel2D_ca_reduced_dt="+str(dt)+"_T="+str(T)+"_D="+str(D1)+"_"+str(E)+"E000199.vtu")
            cell_ca_2D = file_ca_2D.point_data_to_cell_data()
            ca_2D_values = cell_ca_2D.get_array(cell_ca_2D.array_names[0])

            sized = file_ca_2D.compute_cell_sizes(length=False, area=False)
            cell_volumes = np.abs(sized.cell_data["Volume"])
            tot_vol = np.sum(cell_volumes)
            ca_2D[i] = np.sum(ca_2D_values*cell_volumes) / tot_vol
            #print('mean ca', np.sum(ca_values*cell_volumes) / tot_vol)
            
            file_pa = pv.read("../../../new_results/simulations/rhomodel3D_p_reduced_dt="+str(dt)+"_T="+str(T)+"_D="+str(D1)+"_"+str(E)+"E000199.vtu")
            cell_pa = file_pa.point_data_to_cell_data()
            pa_values = cell_pa.get_array(cell_pa.array_names[0])

            sized = file_pa.compute_cell_sizes(length=False, volume=False)
            cell_areas = np.abs(sized.cell_data["Area"])
            tot_vol = np.sum(cell_areas)
            pa[i] = (np.sum(pa_values*cell_areas) / tot_vol)*10**(-15)
            #print('mean pa', np.sum(pa_values*cell_areas) / tot_vol)

            file_pa_2D = pv.read("../../../new_results/simulations/rhomodel2D_p_reduced_dt="+str(dt)+"_T="+str(T)+"_D="+str(D1)+"_"+str(E)+"E000199.vtu")
            cell_pa_2D = file_pa_2D.point_data_to_cell_data()
            pa_2D_values = cell_pa_2D.get_array(cell_pa_2D.array_names[0])

            sized = file_pa_2D.compute_cell_sizes(length=False, volume=False)
            cell_areas = np.abs(sized.cell_data["Area"])
            tot_vol = np.sum(cell_areas)
            pa_2D[i] = (np.sum(pa_2D_values*cell_areas) / tot_vol)*10**(-15)
            #print('mean pa', np.sum(pa_values*cell_areas) / tot_vol)
            
            file_pa_2xD = pv.read("../../../new_results/simulations/rhomodel2xD_p_reduced_dt="+str(dt)+"_T="+str(T)+"_D="+str(D1)+"_"+str(E)+"E000199.vtu")
            cell_pa_2xD = file_pa_2xD.point_data_to_cell_data()
            pa_2xD_values = cell_pa_2xD.get_array(cell_pa_2xD.array_names[0])

            sized = file_pa_2xD.compute_cell_sizes(length=False, volume=False)
            cell_areas = np.abs(sized.cell_data["Area"])
            tot_vol = np.sum(cell_areas)
            pa_2xD[i] = (np.sum(pa_2xD_values*cell_areas) / tot_vol)*10**(-15)
            #print('mean pa', np.sum(pa_values*cell_areas) / tot_vol)

            tempstats = np.load('../../../new_results/temp/tempstats0_3D-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_D='+str(D1)+'_'+str(E)+'E.npy')
            tempstats_mm = np.load('../../../new_results/temp/tempstats0_3D-minmax-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_D='+str(D1)+'_'+str(E)+'E.npy')
            tempstats_2D = np.load('../../../new_results/temp/tempstats0_2D-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_D='+str(D1)+'_'+str(E)+'E.npy')
            tempstats_mm_2D = np.load('../../../new_results/temp/tempstats0_2D-minmax-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_D='+str(D1)+'_'+str(E)+'E.npy')
            tempstats_2xD = np.load('../../../new_results/temp/tempstats0_2xD-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_D='+str(D1)+'_'+str(E)+'E.npy')
            tempstats_mm_2xD = np.load('../../../new_results/temp/tempstats0_2xD-minmax-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_D='+str(D1)+'_'+str(E)+'E.npy')

            #ca[i] = tempstats[0,-1]
            ca_min[i] = tempstats_mm[0,-1]
            ca_max[i] = tempstats_mm[1,-1]
            #ca_2D[i] = tempstats_2D[0,-1]
            ca_min_2D[i] = tempstats_mm_2D[0,-1]
            ca_max_2D[i] = tempstats_mm_2D[1,-1]

            #pa[i] = tempstats[2,-1]*10**(-15)
            pa_min[i] = tempstats_mm[4,-1]*10**(-15)
            pa_max[i] = tempstats_mm[5,-1]*10**(-15)
            #pa_2D[i] = tempstats_2D[2,-1]*10**(-15)
            pa_min_2D[i] = tempstats_mm_2D[4,-1]*10**(-15)
            pa_max_2D[i] = tempstats_mm_2D[5,-1]*10**(-15)
            #pa_2xD[i] = tempstats_2xD[2,-1]*10**(-15)
            pa_min_2xD[i] = tempstats_mm_2xD[4,-1]*10**(-15)
            pa_max_2xD[i] = tempstats_mm_2xD[5,-1]*10**(-15)
            
            i += 1

        ca_err = np.array([ca - ca_min, ca_max - ca])
        pa_err = np.array([pa - pa_min, pa_max - pa])
        ca_err_2D = np.array([ca_2D - ca_min_2D, ca_max_2D - ca_2D])
        pa_err_2D = np.array([pa_2D - pa_min_2D, pa_max_2D - pa_2D])
        pa_err_2xD = np.array([pa_2xD - pa_min_2xD, pa_max_2xD - pa_2xD])


        fig = plt.figure(figsize=(12,8))
        plt.xscale('log')
        plt.xlabel('E (kPa)', fontsize=axes_labels)
        plt.ylabel('$\phi_a$ ($\mu$M)', fontsize=axes_labels)
        plt.ylim(0.2,1.2)
        plt.xticks(fontsize=axes_ticks)
        plt.yticks(fontsize=axes_ticks)
        plt.grid(alpha=0.3, which='both', linestyle=':')
        plt.errorbar(Es, ca, yerr=ca_err, marker='o', markerfacecolor=myblue, markeredgecolor=myblue, linestyle='-', color=myblue, capsize=10, label='3D stimulus')
        plt.errorbar(Es, ca_2D, yerr=ca_err_2D, marker='o', markerfacecolor=myred, markeredgecolor=myred, linestyle='-', color=myred, capsize=10, label='2xD/2D stimulus')
        plt.legend(loc="upper left", fontsize = 'large')
        plt.savefig('../../../new_results/figures/PhiTotvsEs0_allD-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_D='+str(D1)+'.png')
        plt.show()
        plt.close()

        fig = plt.figure(figsize=(12,8))
        plt.xscale('log')
        plt.xlabel('E (kPa)', fontsize=axes_labels)
        plt.ylabel(r'$\rho_a$ ($\mu$mol/$\mu$m$^2$)', fontsize=axes_labels)
        plt.ylim(0, 0.9*10**(-15))
        plt.xticks(fontsize=axes_ticks)
        plt.yticks(fontsize=axes_ticks)
        plt.grid(alpha=0.3, which='both', linestyle=':')
        plt.errorbar(Es, pa, yerr=pa_err, marker='o', markerfacecolor=myblue, markeredgecolor=myblue, linestyle='-', color=myblue, capsize=10, label='3D stimulus')
        plt.errorbar(Es, pa_2xD, yerr=pa_err_2xD, marker='o', markerfacecolor=myred, markeredgecolor=myred, linestyle='-', color=myred, capsize=10, label='2xD stimulus')
        plt.errorbar(Es, pa_2D, yerr=pa_err_2D, marker='o', markerfacecolor=mygreen, markeredgecolor=mygreen, linestyle='-', color=mygreen, capsize=10, label='2D stimulus')
        plt.legend(loc="upper left", fontsize = 'large')
        plt.savefig('../../../new_results/figures/RhoAvsEs0_allD-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_D='+str(D1)+'.png')
        plt.show()
        plt.close()
        
