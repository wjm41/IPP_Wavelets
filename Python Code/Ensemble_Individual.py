# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 10:17:47 2017

@author: user
"""
#import matplotlib.pyplot as plt
import pywt
import time as sys_time
import sys
import os
import numpy as np
import h5py
import Read_Data
import Process_Data
import Diagnostics

#from Single_Particle import Statistics
from mpi4py import MPI

tstart = sys_time.time()


#MPI communicator object:
mpi_comm = MPI.COMM_WORLD
#ID number of the particular MPI process (one out of mpi_size):
mpi_rank = mpi_comm.Get_rank()
#Total number of MPI processes:
mpi_size = mpi_comm.Get_size()

def MPI_mean_array(local_array, root=0):
    """Obtain the global mean from all the MPI tasks and collect the info on the root node."""
    global_array = np.zeros_like(local_array)
    mpi_comm.Reduce(local_array, global_array, op=MPI.SUM, root=root)
    return global_array/float(mpi_size)


def Diffusion(v_input,ref,N):
    '''Takes in an input vector array, selects time data from t=ref to end of simulation.
    Subtracts from each component their value at t=ref, and takes the rms of the vector magnitude.'''
    v_input = v_input[:,ref:N]
    
    v_diff = np.zeros_like(v_input)
    
    v_ref= v_input.T[0]
    
    v_diff= v_input.T-v_ref
     
    v_diff=np.linalg.norm(v_diff.T,axis=0)

    return MPI_mean_array(v_diff)

#Read input arguments
wavelet = str(sys.argv[1])
data_dir = str(sys.argv[2])  # root directory with the trajectory files
species = str(sys.argv[3])   # species to analyse
fraction = float(sys.argv[4])  # number between 0 and 1 to represent a fraction of N_part_total

"""
For the 'beta0p5' simulation data, use Omega_c = 0.3536 for electrons and Omega_c = 0.00553 for ions.
For the 'beta0p2' simulation data, use Omega_c = 0.528 for electrons and Omega_c = 0.00825 for ions.
"""

if data_dir=='beta0p5_eps0p4':
    if species=='electrons':
        #q=-0.015625
        q=-1
        Omega_c=0.3536
        m=1
    else:
        q=1
        Omega_c=0.00533
        m=1*64
else:
    if species=='electrons':
        q=-1
        Omega_c=0.528
        m=1
    else:
        q=1
        Omega_c=0.00825 
        m=1*64

#Choose a given wavelet and edge mode:
wavelet = pywt.Wavelet(wavelet)
mode = 'periodization'

#Get total number of particles:
cwd = os.getcwd()
dset = 't'
path = cwd + '/' + data_dir +'/'
base_name = species+ '-tracks-'    
f=h5py.File(path + base_name + dset + '.h5', 'r')
N_part_total = f.attrs['NTRACKS']
f.close()

N_part = int(fraction*N_part_total)

#Print some stuff from the process with rank 0:
if mpi_rank==0:
    print("Running on " + str(mpi_size) + " MPI processes.")
    print("Input data = "+species+" from " + data_dir)
    print("Total no. of particles =" +str(N_part_total))
    print("No. of particles used =" +str(N_part))

#Partition the total input file with N particles into mpi_size parts:
mpi_borders = np.linspace(0, N_part, mpi_size + 1).astype('int')

my_border_low = mpi_borders[mpi_rank]
my_border_high = mpi_borders[mpi_rank+1]

#Redefining N_part as number of particles analysed per processor
N_part = my_border_high - my_border_low


#Read Trajectory Data
time, U, E, B, Energy = Read_Data.Read_Trajectory_Data(my_border_low, my_border_high, data_dir, species)
#Comment: U is used to denote the proper relativistic velocity, which is U = gamma*V

N = len(time)

#Process Trajectory Data using chosen wavelet and edge mode
if mpi_rank==0:
    print('Processing Data for particles...')

#Transfrom the data and compute new quantities:
time, V, E, b, V_par, V_perp, vpar, vperp, Acc_mag, Acc_par, Acc_perp1, Acc_perp2 = Process_Data.Process_raw(time, U, E, B, Omega_c, q, m)
N_pad = len(vpar[0,:])

ref=int(np.ceil(N/4))

#vpar_Global=Diffusion(vpar,ref,N)
#vperp_Global=Diffusion(vperp,ref,N)
#
#if mpi_rank==0:
#    np.savetxt(path+'/'+species+'/Diffusion/v_par_rms_raw.dat', np.c_[time[ref:N], vpar_Global])
#    np.savetxt(path+'/'+species+'/Diffusion/v_perp_rms_raw.dat', np.c_[time[ref:N], vperp_Global])

#Compute wavelet coefficeints:
coefs_V = Process_Data.Process_coef(V, wavelet, mode)
coefs_E = Process_Data.Process_coef(E, wavelet, mode)
coefs_b = Process_Data.Process_coef(b, wavelet, mode)
coefs_V_par = Process_Data.Process_coef(V_par, wavelet, mode)
coefs_V_perp = Process_Data.Process_coef(V_perp, wavelet, mode)
coefs_vpar = Process_Data.Process_coef(vpar, wavelet, mode)
coefs_vperp = Process_Data.Process_coef(vperp, wavelet, mode)
#coefs_energy = Process_Data.Process_coef(Energy, wavelet, mode)
#coefs_Acc_mag = Process_Data.Process_coef(Acc_mag, wavelet, mode)
#coefs_Acc_par = Process_Data.Process_coef(Acc_par, wavelet, mode)
#coefs_Acc_perp1 = Process_Data.Process_coef(Acc_perp1, wavelet, mode)
#coefs_Acc_perp2 = Process_Data.Process_coef(Acc_perp2, wavelet, mode)

J = len(coefs_vpar)
 
if mpi_rank==0:
    print('Processing done!')
    print('Number of time steps in the original signal:' + str(N))
    print('Number of time steps in the padded array:' + str(N_pad))
    print('Number of scales for the orthogonal transform:' + str(J-1))

#Prepare the arrays to hold all the moments:
V_moments_vs_j = np.zeros((3,J, 3, 5)) #at 3 times, for all 3 components, versus J for 5 different moments (2nd, 3rd, 4th, 5th, 6th)
vpar_moments_vs_j = np.zeros((3,J, 5))
vperp_moments_vs_j = np.zeros((3,J, 5))
scale_j = np.zeros(J) #approximate time scale corresponding to some given J in the real physical units
    
for j in range(J):
    """Extract that selected scale for all the particles and all the vector components over 
    a given time window specified by t_ind_low and t_ind_high:"""
    
#    Mu_j_full= Diagnostics.new_extract_scale(coefs_mu, j, wavelet, mode, 0 ,N)
    V_j_full = Diagnostics.new_extract_scale(coefs_V, j, wavelet, mode, 0, N)
    E_j_full = Diagnostics.new_extract_scale(coefs_E, j, wavelet, mode, 0, N)
    b_j_full = Diagnostics.new_extract_scale(coefs_b, j, wavelet, mode, 0, N)
    
    V_par_j_full = Diagnostics.new_extract_scale(coefs_V_par, j, wavelet, mode, 0, N)
    V_perp_j_full = Diagnostics.new_extract_scale(coefs_V_perp, j, wavelet, mode, 0, N)
    
    vpar_j_full = Diagnostics.new_extract_scale(coefs_vpar, j, wavelet, mode, 0, N)
    vperp_j_full = Diagnostics.new_extract_scale(coefs_vperp, j, wavelet, mode, 0, N)
    
#    Acc_mag_j_full = Diagnostics.new_extract_scale(coefs_Acc_mag, j, wavelet, mode ,0 ,N)
#    Acc_par_j_full = Diagnostics.new_extract_scale(coefs_Acc_par, j, wavelet, mode ,0 ,N)
#    Acc_perp1_j_full = Diagnostics.new_extract_scale(coefs_Acc_perp1, j, wavelet, mode ,0 ,N)
#    Acc_perp2_j_full = Diagnostics.new_extract_scale(coefs_Acc_perp2, j, wavelet, mode ,0 ,N)
#    
    scale_j[j] = time[int(N_pad//len(coefs_vpar[j][0]))] - time[0]    
#    
    #Specify the time window over which to average:
    for i in range(3):
        t_ind_low = i*N//4
        t_ind_high = (i+2)*N//4
        
        V_j = V_j_full[:, :, t_ind_low : t_ind_high]
        E_j = E_j_full[:, :, t_ind_low : t_ind_high]
        b_j = b_j_full[:, :, t_ind_low : t_ind_high]
        V_par_j = V_par_j_full[:, :, t_ind_low : t_ind_high]
        V_perp_j = V_perp_j_full[:, :, t_ind_low : t_ind_high]
        vpar_j = vpar_j_full[:, t_ind_low : t_ind_high]
        vperp_j = vperp_j_full[:, t_ind_low : t_ind_high]    
#        
#        Acc_mag_j = Acc_mag_j_full[:, t_ind_low : t_ind_high]
#        Acc_par_j = Acc_par_j_full[:, t_ind_low : t_ind_high]
#        Acc_perp1_j = Acc_perp1_j_full[ :, t_ind_low : t_ind_high]    
#        Acc_perp2_j = Acc_perp2_j_full[ :, t_ind_low : t_ind_high]
#
#        # Evaluate the net,parallel, and perpendicular instantaneous power ( P = time deriv. of energy) at scale j 
#        P_j_net = q*(V_j[0]*E_j[0] + V_j[1]*E_j[1] + V_j[2]*E_j[2])
#        P_j_par = q*(V_par_j[0]*E_j[0] + V_par_j[1]*E_j[1] + V_par_j[2]*E_j[2])
#        P_j_perp = q*(V_perp_j[0]*E_j[0] + V_perp_j[1]*E_j[1] + V_perp_j[2]*E_j[2])
#        
#        #Calculate the PDFs
#        norm_t = 1./float(t_ind_high - t_ind_low)
#        P_net_PDF_tavg, PDF_bins = Diagnostics.PDF_nomean(P_j_net[:,0])
#        P_par_PDF_tavg, PDF_bins = Diagnostics.PDF_nomean(P_j_par[:,0])
#        P_perp_PDF_tavg, PDF_bins = Diagnostics.PDF_nomean(P_j_perp[:,0])
#        P_net_PDF_tavg *= norm_t
#        P_par_PDF_tavg *= norm_t
#        P_perp_PDF_tavg *= norm_t
#        
#        acc_mag_PDF_tavg, PDF_bins = Diagnostics.PDF_nomean(Acc_mag_j[:,0])
#        acc_par_PDF_tavg, PDF_bins = Diagnostics.PDF_nomean(Acc_par_j[:,0])
#        acc_perp1_PDF_tavg, PDF_bins = Diagnostics.PDF_nomean(Acc_perp1_j[:,0])
#        acc_perp2_PDF_tavg, PDF_bins = Diagnostics.PDF_nomean(Acc_perp2_j[:,0])
#
#        acc_mag_PDF_tavg *= norm_t
#        acc_par_PDF_tavg *= norm_t
#        acc_perp1_PDF_tavg *= norm_t
#        acc_perp2_PDF_tavg *= norm_t
#        
#        for t_ind in range(1, t_ind_high-t_ind_low):
#            P_net_PDF, PDF_bins = Diagnostics.PDF_nomean(P_j_net[:,t_ind])
#            P_par_PDF, PDF_bins = Diagnostics.PDF_nomean(P_j_par[:,t_ind])
#            P_perp_PDF, PDF_bins = Diagnostics.PDF_nomean(P_j_perp[:,t_ind])
#            P_net_PDF_tavg += P_net_PDF*norm_t 
#            P_par_PDF_tavg += P_par_PDF*norm_t
#            P_perp_PDF_tavg += P_perp_PDF*norm_t 
#            
#            acc_mag_PDF, PDF_bins = Diagnostics.PDF_nomean(Acc_mag_j[:,t_ind])
#            acc_par_PDF, PDF_bins = Diagnostics.PDF_nomean(Acc_par_j[:,t_ind])
#            acc_perp1_PDF, PDF_bins = Diagnostics.PDF_nomean(Acc_perp1_j[:,t_ind])
#            acc_perp2_PDF, PDF_bins = Diagnostics.PDF_nomean(Acc_perp2_j[:,t_ind])
#            acc_mag_PDF_tavg += acc_mag_PDF*norm_t 
#            acc_par_PDF_tavg += acc_par_PDF*norm_t
#            acc_perp1_PDF_tavg += acc_perp1_PDF*norm_t 
#            acc_perp2_PDF_tavg += acc_perp2_PDF*norm_t            
#
#        P_net_PDF_global = MPI_mean_array(P_net_PDF_tavg)
#        P_par_PDF_global = MPI_mean_array(P_par_PDF_tavg)
#        P_perp_PDF_global = MPI_mean_array(P_perp_PDF_tavg)
#        
#        acc_mag_PDF_global = MPI_mean_array(acc_mag_PDF_tavg)
#        acc_par_PDF_global = MPI_mean_array(acc_par_PDF_tavg)
#        acc_perp1_PDF_global = MPI_mean_array(acc_perp1_PDF_tavg)
#        acc_perp2_PDF_global = MPI_mean_array(acc_perp2_PDF_tavg)
#        
#        Num_bins=len(PDF_bins)/2
#        First_half_net=P_net_PDF_global[:Num_bins]
#        Second_half_net=P_net_PDF_global[Num_bins:]
#        Skewness_net=np.log10(First_half_net[::-1]/Second_half_net)
#        
#        First_half_par=P_par_PDF_global[:Num_bins]
#        Second_half_par=P_par_PDF_global[Num_bins:]
#        Skewness_par=np.log10(First_half_par[::-1]/Second_half_par)
#        
#        First_half_perp=P_perp_PDF_global[:Num_bins]
#        Second_half_perp=P_perp_PDF_global[Num_bins:]
#        Skewness_perp=np.log10(First_half_perp[::-1]/Second_half_perp)
#                          
#        if mpi_rank==0:
#            np.savetxt(path+'/'+species+'/PDFs/pdf_P_j' + str(j) + '_net_t'+str(i+1)+'.dat', np.c_[PDF_bins, P_net_PDF_global])
#            np.savetxt(path+'/'+species+'/PDFs/pdf_P_j' + str(j) + '_par_t'+str(i+1)+'.dat', np.c_[PDF_bins, P_par_PDF_global])
#            np.savetxt(path+'/'+species+'/PDFs/pdf_P_j' + str(j) + '_perp_t'+str(i+1)+'.dat', np.c_[PDF_bins, P_perp_PDF_global]) 
#            
#            np.savetxt(path+'/'+species+'/PDFs/pdf_Acc_j' + str(j) + '_mag_t'+str(i+1)+'.dat', np.c_[PDF_bins, acc_mag_PDF_global])
#            np.savetxt(path+'/'+species+'/PDFs/pdf_Acc_j' + str(j) + '_par_t'+str(i+1)+'.dat', np.c_[PDF_bins, acc_par_PDF_global])
#            np.savetxt(path+'/'+species+'/PDFs/pdf_Acc_j' + str(j) + '_perp1_t'+str(i+1)+'.dat', np.c_[PDF_bins, acc_perp1_PDF_global]) 
#            np.savetxt(path+'/'+species+'/PDFs/pdf_Acc_j' + str(j) + '_perp2_t'+str(i+1)+'.dat', np.c_[PDF_bins, acc_perp2_PDF_global])
#            
#            np.savetxt(path+'/'+species+'/Skewness/P_net_j'+str(j)+'_t'+str(i+1)+'.dat', np.c_[PDF_bins[Num_bins:], Skewness_net])
#            np.savetxt(path+'/'+species+'/Skewness/P_par_j'+str(j)+'_t'+str(i+1)+'.dat', np.c_[PDF_bins[Num_bins:], Skewness_par])
#            np.savetxt(path+'/'+species+'/Skewness/P_perp_j'+str(j)+'_t'+str(i+1)+'.dat', np.c_[PDF_bins[Num_bins:], Skewness_perp])

        #Calculate all the moments at scale j
        if mpi_rank==0:
            print("Calculating moments at scale "+str(j))
        V_moments = Diagnostics.Averaged_Moments(V_j)    
        V_moments_global = MPI_mean_array(V_moments)    
        V_moments_vs_j[i][j] = V_moments_global
        
        vpar_moments = Diagnostics.Averaged_Moments(vpar_j)    
        vpar_moments_global = MPI_mean_array(vpar_moments)    
        vpar_moments_vs_j[i][j] = vpar_moments_global
        
        vperp_moments = Diagnostics.Averaged_Moments(vperp_j)    
        vperp_moments_global = MPI_mean_array(vperp_moments)    
        vperp_moments_vs_j[i][j] = vperp_moments_global
        
    ##Calculate Ensemble Average of Mu
    #Mu_global=MPI_mean_array(np.mean(Mu_j_full,axis=0))
    
    #Output Full trace of Moments for calculating Flatness
#    if mpi_rank == 0:
#        np.savetxt(path+'/'+species+'/Moments/Mu_'+str(j)+'_trace.dat', np.c_[time, Mu_global])
#        moments_par=Diagnostics.Moment(vpar_j_full,N)
#        moments_perp=Diagnostics.Moment(vperp_j_full,N)
#        np.savetxt(path+'/'+species+'/Moments/v_par_j'+str(j)+'_trace.dat', np.c_[time, moments_par])
#        np.savetxt(path+'/'+species+'/Moments/v_perp_j'+str(j)+'_trace.dat', np.c_[time, moments_perp])
        
#    #Calculate velocity diffusion
#    vpar_j_Global=Diffusion(vpar_j_full, ref, N)
#    vperp_j_Global=Diffusion(vperp_j_full,ref, N)
#    
#    if mpi_rank==0:
#        np.savetxt(path+'/'+species+'/Diffusion/v_par_rms_j'+str(j)+'.dat', np.c_[time[ref:N], vpar_j_Global])
#        np.savetxt(path+'/'+species+'/Diffusion/v_perp_rms_j'+str(j)+'.dat', np.c_[time[ref:N], vperp_j_Global])
#        
for i in range(3):
    np.savetxt(path+'/'+species+'/Moments/moments_vpar_t'+str(i+1)+'.dat', np.c_[scale_j, vpar_moments_vs_j[i]])
    np.savetxt(path+'/'+species+'/Moments/moments_vperp_t'+str(i+1)+'.dat', np.c_[scale_j, vperp_moments_vs_j[i]])
           
#Specify the time window over which to average:
    
#    Energy_j_full = Diagnostics.new_extract_scale(coefs_energy, j, wavelet, mode, 0, N)
#    for i in range(3):
#        t_ind_low = i*N//4
#        t_ind_high = (i+2)*N//4
#    
#        Energy_j = Energy_j_full[:, t_ind_low : t_ind_high]           
#    
#        norm_t = 1./float(t_ind_high - t_ind_low)
#        E_j_PDF_tavg, PDF_bins = Diagnostics.PDF(Energy_j[:,0])
#        
#        E_j_PDF_tavg *= norm_t
#        
#        for t_ind in range(1, t_ind_high-t_ind_low):
#            E_j_PDF, PDF_bins = Diagnostics.PDF(Energy_j[:,t_ind])
#            
#            E_j_PDF_tavg += E_j_PDF*norm_t 
#            
#        E_j_PDF_global = MPI_mean_array(E_j_PDF_tavg)
#        
#        if mpi_rank == 0:
#            np.savetxt(path+'/'+species+'/PDFs/pdf_energy_j'+str(j)+'_t'+str(i+1)+'.dat', np.c_[PDF_bins, E_j_PDF_global])
#            
#        #Calculating the skewness    
#        Num_bins=len(PDF_bins)/2
#        First_half=E_j_PDF_global[:Num_bins]
#        Second_half=E_j_PDF_global[Num_bins:]
#        Skewness=np.log10(First_half[::-1]/Second_half)
#        if mpi_rank == 0:
#            np.savetxt(path+'/'+species+'/Skewness/energy_j'+str(j)+'_t'+str(i+1)+'.dat', np.c_[PDF_bins[Num_bins:], Skewness])
#            
#for i in range(3):
#    t_ind_low = i*N//4
#    t_ind_high = (i+2)*N//4
#
#    Energy_fraction = Energy[:, t_ind_low : t_ind_high]        
#
#    norm_t = 1./float(t_ind_high - t_ind_low)
#    E_PDF_tavg, PDF_bins = Diagnostics.PDF_nomean(Energy_fraction[:,0])
#    
#    E_PDF_tavg *= norm_t
#    
#    for t_ind in range(1, t_ind_high-t_ind_low):
#        E_PDF, PDF_bins = Diagnostics.PDF_nomean(Energy_fraction[:,t_ind])
#        
#        E_PDF_tavg += E_PDF*norm_t 
#        
#    E_PDF_global = MPI_mean_array(E_PDF_tavg)
#    
#    if mpi_rank == 0:
#        np.savetxt(path+'/'+species+'/PDFs/pdf_energy_raw_t'+str(i+1)+'.dat', np.c_[PDF_bins, E_PDF_global])    
#        
#    #Calculating the skewness    
#    Num_bins=len(PDF_bins)/2
#    First_half=E_j_PDF_global[:Num_bins]
#    Second_half=E_j_PDF_global[Num_bins:]
#    Skewness=np.log10(First_half[::-1]/Second_half)
#    if mpi_rank == 0:
#        np.savetxt(path+'/'+species+'/Skewness/energy_raw_t'+str(i+1)+'.dat', np.c_[PDF_bins[Num_bins:], Skewness])

if mpi_rank==0:
    print('Moments and PDFs exported')
    t2 = (sys_time.time()-tstart)
    print('Diagnostics complete; time taken = ' + str(np.round(t2/60,decimals=1)) + ' minutes')
sys.exit(0)