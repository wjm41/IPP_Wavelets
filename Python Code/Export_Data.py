# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:04:27 2017

@author: user
"""
import os

def Export(Dataset, species, quantity,scale,data_name,Data):
    '''Exports the rows of Data into a .data file that is named based on Dataset,
    species,quantity,scale level, and dataname'''
    cwd = os.getcwd()
    File_name=cwd+'/'+Dataset+'/'+species+'/'+quantity+'/'+data_name+'_'+scale+'.dat'
    File=open(File_name,'w')
    for i in Data:
        File.write(str(i)+'\n')
    File.close()
    
    
def Mass_Export(Data_File,species,V,E,B,vpar,vperp,vparsq,vperpsq,j,Thing):
    
    Export(Data_File,species,'vx',str(j),Thing,V[0])
    Export(Data_File,species,'vy',str(j),Thing,V[1])
    Export(Data_File,species,'vz',str(j),Thing,V[2])
    
    Export(Data_File,species,'ex',str(j),Thing,E[0])
    Export(Data_File,species,'ey',str(j),Thing,E[1])
    Export(Data_File,species,'ez',str(j),Thing,E[2])
    
    Export(Data_File,species,'bx',str(j),Thing,B[0])
    Export(Data_File,species,'by',str(j),Thing,B[1])
    Export(Data_File,species,'bz',str(j),Thing,B[2])
    
    Export(Data_File,species,'vparx',str(j),Thing,vpar[0])
    Export(Data_File,species,'vpary',str(j),Thing,vpar[1])
    Export(Data_File,species,'vparz',str(j),Thing,vpar[2])
    
    Export(Data_File,species,'vperpx',str(j),Thing,vperp[0])
    Export(Data_File,species,'vperpy',str(j),Thing,vperp[1])
    Export(Data_File,species,'vperpz',str(j),Thing,vperp[2])
    
    Export(Data_File,species,'vparsq',str(j),Thing,vparsq)
    Export(Data_File,species,'vperpsq',str(j),Thing,vperpsq)
    
    
def Export_PDF(Dataset, species, quantity,scale,data_name,time,PDF_Data):
    '''Exports the rows of PDF_Data into a .data file that is named based on Dataset,
    species,quantity,scale level, dataname, and the time fraction.'''
    cwd = os.getcwd()
    File_name=cwd+'/'+Dataset+'/'+species+'/'+quantity+'/'+data_name+'_'+str(time)+'_'+scale+'.dat'
    File=open(File_name,'w')
    for i in range(len(PDF_Data[0])):
        File.write(str(PDF_Data[0][i])+'\t'+str(PDF_Data[1][i])+'\n')
    File.close()
    
def Mass_Export_PDF(Data_File,species,vx,vy,vz,j,time,Thing,name):
    
    Export_PDF(Data_File,species,name+'x',str(j),Thing,time,vx)
    Export_PDF(Data_File,species,name+'y',str(j),Thing,time,vy)
    Export_PDF(Data_File,species,name+'z',str(j),Thing,time,vz)
    
