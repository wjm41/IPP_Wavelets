# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 14:15:15 2017

@author: user
"""
import numpy as np
from scipy.signal import tukey
import pywt
np.set_printoptions(threshold=np.nan)

      
def mod_signal(sig, prec='float32'):
  """Reduce edge effects for a periodized signal. First the mean is subtracted and a window function is 
  applied on the zero-mean signal. Then, the signal is padded with zeros at the end. Finally, 
  the mean value is added back to the windowed signal."""
  sig_avg = np.mean(sig, axis=-1)
  N = sig.shape[-1]
  ndim = sig.ndim
  # Pad the signal with as many points needed to reach a power of 2, but not with less than N//2 points for the sake of 
  # reducing edge effects:
  nPad = 2**(int(np.log2(N)) + 1) - N
  if nPad < N//2:
      nPad = 2**(int(np.log2(N)) + 2) - N
  #For scalar data (1st dimension is particle index, 2nd index is for time):
  if ndim == 2:
      return np.array(np.pad((sig - sig_avg[:, np.newaxis])*tukey(N, 0.2), ((0,0), (0, nPad)), 'constant')  
      + sig_avg[:, np.newaxis], dtype=prec)
  #For vector data (1st index is vector component, 2nd is particle index, 3rd is time):
  if ndim == 3:
      return np.array(np.pad((sig - sig_avg[:, :, np.newaxis])*tukey(N, 0.2), ((0,0), (0,0), (0, nPad)), 'constant')  
      + sig_avg[:, :, np.newaxis], dtype=prec)
      

def Process_raw(time, U, E, B, Omega_c, q, m):

    #Initialize the low-pass filter for the ExB drift:
    N = len(time)
    omega = np.fft.rfftfreq(N, (time[-1] - time[0])/(2*np.pi*N))
    omega0 = 0.9*Omega_c
    lp_filter = np.where(omega > omega0, np.zeros_like(omega), np.ones_like(omega))

    #Calculate b and E x B:
    Bsq = B[0]**2 + B[1]**2 + B[2]**2
    v_exb = np.cross(E, B, axis=0)/Bsq
    b = B/np.sqrt(Bsq)    
    
    #Low-pass filter just slightly below Omega_c (similar to a gyroaverage operation):
    v_exb = np.fft.irfft(np.fft.rfft(v_exb)*lp_filter, n=N)

    gamma_exb = 1./np.sqrt(1. - v_exb[0]**2 - v_exb[1]**2 - v_exb[2]**2)
    u_exb = gamma_exb*v_exb
    gamma = np.sqrt(1. + U[0]**2 + U[1]**2 + U[2]**2)
    boost = gamma - (U[0]*u_exb[0] + U[1]*u_exb[1] + U[2]*u_exb[2]) / (1. + gamma_exb)
    #Transfrom the proper velocity to the ExB frame using the relativistic transform for the velocity:
    U = U - boost*u_exb
    gamma = np.sqrt(1. + U[0]**2 + U[1]**2 + U[2]**2) 
    #Now, calculate the regular velocity without contributions from the ExB drift:
    V = U / gamma
    
    #Transform the E & B Fields
    n_exb = v_exb/np.linalg.norm(v_exb, axis=0)
    
    E = gamma_exb*(E + np.cross(v_exb, B, axis=0)) - (gamma_exb - 1)*(E[0]*n_exb[0] + E[1]*n_exb[1] + E[2]*n_exb[2])*n_exb

    B = gamma_exb*(B - np.cross(v_exb, E, axis=0)) - (gamma_exb - 1)*(B[0]*n_exb[0] + B[1]*n_exb[1] + B[2]*n_exb[2])*n_exb
    
    #Updating B field
    Bsq = B[0]**2 + B[1]**2 + B[2]**2
    b = B/np.sqrt(Bsq)    
    
    # Components of velocity parallel and perpendicular to the local magnetic field:
    vpar = V[0]*b[0] + V[1]*b[1] + V[2]*b[2] #vpar as a scalar (can have positive and negative values)
    V_par = vpar*b #vpar as a vector
    V_perp = V - V_par #vperp as a vector
    vperp = np.sqrt(V_perp[0]**2 + V_perp[1]**2 + V_perp[2]**2) #vperp as a (positive) scalar

    #Calculate the acceleration
    acc = q*(E + np.cross(V, B, axis=0))/m
    acc_mag=np.linalg.norm(acc, axis=0) #magnitude of acceleration as a scalar
    acc_par=acc[0]*b[0] + acc[1]*b[1] + acc[2]*b[2]
    
    n_vperp = V_perp/np.linalg.norm(V_perp, axis=0)
    acc_perp1 = acc[0]*n_vperp[0] + acc[1]*n_vperp[1] + acc[2]*n_vperp[2]
    
    n_vxb=np.cross(V,B,axis=0)/np.linalg.norm(np.cross(V,B,axis=0),axis=0)
    acc_perp2 = acc[0]*n_vxb[0] + acc[1]*n_vxb[1] + acc[2]*n_vxb[2]
            
    #Reduce edge effects for the signals:
    V = mod_signal(V)
    E = mod_signal(E)
    b = mod_signal(b)
    V_par = mod_signal(V_par)
    V_perp = mod_signal(V_perp)
    vpar = mod_signal(vpar)
    vperp = mod_signal(vperp)
    acc_mag=mod_signal(acc_mag)
    acc_par=mod_signal(acc_par)
    acc_perp1=mod_signal(acc_perp1)
    acc_perp2=mod_signal(acc_perp2)
    
#    time*=Omega_c
    
    return time, V, E, b, V_par, V_perp, vpar, vperp, acc_mag, acc_par, acc_perp1, acc_perp2
    
    
def Process_coefs(V, E, B, V_par, V_perp, vpar, vperp, ene, wavelet, mode):
    
    #Perform the wavelet transform for all (!) the particle and all (!) vector components at once: 
    coefs_V = (pywt.wavedec(V, wavelet, mode=mode))
    coefs_E = (pywt.wavedec(E, wavelet, mode=mode))
    coefs_B = (pywt.wavedec(B, wavelet, mode=mode))
    coefs_V_par = (pywt.wavedec(V_par, wavelet, mode=mode))
    coefs_V_perp = (pywt.wavedec(V_perp, wavelet, mode=mode))
    coefs_vpar = (pywt.wavedec(vpar, wavelet, mode=mode))
    coefs_vperp = (pywt.wavedec(vperp, wavelet, mode=mode))
    coefs_energy=(pywt.wavedec(ene, wavelet, mode=mode))
    
    """  Note: 

    For vector data, you select the wavelet coefficients for scale "j", for vector component "c", for particle "p"
    with:
    
    coefs[j][c,p]
    
    For scalar data, you select wavelet coefficients for scale "j", for particle "p" with:
    
    coefs[j][p]
    
    """
    
    J = len(coefs_V)
    
    return coefs_V, coefs_E, coefs_B, coefs_V_par, coefs_V_perp, coefs_vpar, coefs_vperp, coefs_energy, J
    
    
def Process_coef(input_data, wavelet, mode):
    """Obtain just one set of coefficeints."""
    
    #Perform the wavelet transform for all (!) the particle and all (!) vector components at once: 
    return pywt.wavedec(input_data, wavelet, mode=mode)


