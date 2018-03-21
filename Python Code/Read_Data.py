# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 13:56:16 2017

@author: user
"""

import numpy as np
import h5py
import os

def Read_Trajectory_Data(part_ind_low, part_ind_high, data_dir, species, prec='float32'):
    #Set path to the particle trajectory data:
    cwd = os.getcwd()
    path = cwd + '/' + data_dir +'/'
    base_name = species+ '-tracks-'
    
    dset = 't'
    f=h5py.File(path + base_name + dset + '.h5', 'r')

#    print('Reading time data...')
    time = np.array(f['/' + dset][0,:])
    f.close()

#    print('Reading momenta...')
    dset = 'p1'
    f=h5py.File(path + base_name + dset + '.h5', 'r')
    p1 = np.array(f['/' + dset][part_ind_low:part_ind_high, : ], dtype=prec, copy=False)
    f.close()
    dset = 'p2'
    f=h5py.File(path + base_name + dset + '.h5', 'r')
    p2 = np.array(f['/' + dset][part_ind_low:part_ind_high, : ], dtype=prec, copy=False)
    f.close()
    dset = 'p3'
    f=h5py.File(path + base_name + dset + '.h5', 'r')
    p3 = np.array(f['/' + dset][part_ind_low:part_ind_high, : ], dtype=prec, copy=False)
    f.close() 
#    print('Reading E field...')
    dset = 'E1'
    f=h5py.File(path + base_name + dset + '.h5', 'r')
    E1 = np.array(f['/' + dset][part_ind_low:part_ind_high, : ], dtype=prec, copy=False)
    f.close()
    dset = 'E2'
    f=h5py.File(path + base_name + dset + '.h5', 'r')
    E2 = np.array(f['/' + dset][part_ind_low:part_ind_high, : ], dtype=prec, copy=False)
    f.close()
    dset = 'E3'
    f=h5py.File(path + base_name + dset + '.h5', 'r')
    E3 = np.array(f['/' + dset][part_ind_low:part_ind_high, : ], dtype=prec, copy=False)
    f.close()
    
#    print('Reading B field...')
    dset = 'B1'
    f=h5py.File(path + base_name + dset + '.h5', 'r')
    B1 = np.array(f['/' + dset][part_ind_low:part_ind_high, : ], dtype=prec, copy=False)
    f.close()
    dset = 'B2'
    f=h5py.File(path + base_name + dset + '.h5', 'r')
    B2 = np.array(f['/' + dset][part_ind_low:part_ind_high, : ], dtype=prec, copy=False)
    f.close()
    dset = 'B3'
    f=h5py.File(path + base_name + dset + '.h5', 'r')
    B3 = np.array(f['/' + dset][part_ind_low:part_ind_high, : ], dtype=prec, copy=False)
    f.close()
    
    dset = 'ene'
    f=h5py.File(path + base_name + dset + '.h5', 'r')
    Energy = np.array(f['/' + dset][part_ind_low:part_ind_high, : ], dtype=prec, copy=False)
    f.close()


    #Work with vectors (!):
    U = np.array([p1, p2, p3])
    B = np.array([B1, B2, B3])
    E = np.array([E1, E2, E3])
    
    return time, U, E, B, Energy
