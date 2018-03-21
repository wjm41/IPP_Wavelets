# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 10:46:38 2017

@author: user
"""
import numpy as np
import pywt


def new_extract_scale(coefs, level, wavelet, mode, t_ind_low, t_ind_high):
  """Extract a given scale corresponding to a given level of the wavelet 
  coefficients. The input are the wavelet coefficents 'coefs', the level 'level,
  the wavelet type 'wavelet' and wavelet transform edge mode 'mode'.
  Note that level=0 corresponds to the so-called approximation coefficients that
  represents the local mean of the signal. The so-called detail coefficients start 
  with level=1."""
  c = []
  max_levels = len(coefs) - 1
  for l in range(max_levels+1):
    if l==level:
      c = c + [coefs[l]]
    else:
      c = c + [np.zeros_like(coefs[l])]
  scale = pywt.waverec(c, wavelet, mode=mode)
  ndim = scale.ndim
  #Return the extracted scale only over a given time window:
  if ndim == 2:
      return scale[:, t_ind_low : t_ind_high]
  if ndim == 3:
      return scale[:, :, t_ind_low : t_ind_high]
  
def Averaged_Moments(extracted):
    '''Given a matrix of wavelet coefficients for N particles'''
    ndim = extracted.ndim
    if ndim == 3:
        moments = np.zeros((5,3))
    if ndim == 2:
        moments = np.zeros(5)
    for q in range(2,7):
        moments[q-2] = np.mean(np.mean(np.abs(extracted)**q, axis=(ndim-2)), axis=-1)
    return np.transpose(moments)

def Moment(extracted,N):
    ndim = extracted.ndim
    if ndim == 3:
        moments = np.zeros((5,3,N))
    if ndim == 2:
        moments = np.zeros((5,N))
    for q in range(2,7):
        moments[q-2] = np.mean(np.abs(extracted)**q, axis=(ndim-2))
    return np.transpose(moments)

def Flatness(Second,Fourth):
    '''Takes the 2nd and 4th moments as matrix inputs, and returns the flatness as a matrix.'''
    flatness=Fourth/Second**2
    return flatness

def Anisotropy_Magnitude(Top_2nd,Bottom_2nd):
    '''Takes two second moments as matrix inputs, and returns the magnitude of 
    the component-wise anisotropy magnitude as a matrix.'''
    Mag=Top_2nd/Bottom_2nd
    return Mag

def Anisotropy_Flucuation(Top_sigma,Bottom_sigma):
    '''Takes two spatial variabilities as matrix inputs, and returns the 
    magnitude of the component-wise anisotropy flucuation as a matrix.'''
    Fluc=Top_sigma/Bottom_sigma
    return Fluc

def Intermittency(Top_Flat,Bottom_Flat):
    '''Takes two flatnesses as matrix inputs, and returns the magnitude of 
    the component-wise anisotropic intermittency as a matrix.'''
    Lamb=(Top_Flat-1)/(Bottom_Flat-1)
    return Lamb

def Pitch_Angle(v_perp,vpar):
    '''Takes the scale-dependent perpendicular and parallel velocities as matrix inputs, 
    and returns the magnitude of the pitch angle as a matrix'''
    pitch=np.arctan(v_perp/vpar)
    return pitch

def Mu(v_perp,B):
    '''Takes the scale-dependent perpendicular velocity and magnetic field magnitude as matrix inputs,
    and returns the magnitude of the magnetic moment as a matrix'''
    mu=(v_perp**2)/B
    return mu


def PDF(data, val_range=8., dbin=0.1):
    '''Constructs a PDF distribution from an input matrix, using histogram bins of size dbin from -val_range to val_range.'''
    #Normalize with the standard deviation
    data = (data - np.mean(data)) / np.std(data)
    
    #Calculate the pdfs:
    pdf, edges = np.histogram(data, bins=np.arange(-val_range,
    val_range+dbin, dbin), density=True, range=(-val_range, val_range))
    bins = edges[:-1] + 0.5*dbin
    return pdf, bins

def PDF_nomean(data, val_range=8., dbin=0.1):
    '''Constructs a PDF distribution from an input matrix, using histogram bins of size dbin from -val_range to val_range.'''
    #Normalize with the standard deviation
    data = data / np.std(data)
    
    #Calculate the pdfs:
    pdf, edges = np.histogram(data, bins=np.arange(-val_range,
    val_range+dbin, dbin), density=True, range=(-val_range, val_range))
    bins = edges[:-1] + 0.5*dbin
    return pdf, bins
    
    
    
