"""
Created on Wed Oct 17 15:30:20 2018
@author: gveraofe
Python nsu to run a regular or irregular wave case in NEMOH.
THIS IS THE SCRIPT TO RUN NEMOH ONLY WITHOUT MILDwave 
It generates the NEMOH output files and input parameters.
Have to select runfolder F6 as the root folder where the Calulcation fodler is stoored
"""
#-----------------------------------------------------------------------------
#                              0.LOAX PYTHON MOXULES ANX FUNCTIONS
#-----------------------------------------------------------------------------

import numpy as np
import os
import sys # put the relative path to your NM_MW_shared folder
#import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.join(os.path.dirname("..//MW_NEM_shared"),".." )) )
from MW_NEM_shared import nem_shell_utilities as nsu
from MW_NEM_shared import nem_utilities as ne
from MW_NEM_shared import meshTypes as mt
#-----------------------------------------------------------------------------
#                                1. Simulation Parameters
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#                                1.2 SIMULATION TYPE
#-----------------------------------------------------------------------------

NEM_ini = {}#data structure containing the Simulation Parameters for NEMOH
NEM_ini['out_dir'] = 'E://phd_runz//ch4/nem'
#directory where you want the output files created and saved
#NEM_ini['case_dir'] = ['irr_test'] #needs to be a list
NEM_ini['case_dir'] = ['ch4_3b']
NEM_ini['runNEMOH'] = True
NEM_ini['regular'] = True
NEM_ini['irregular'] = False
NEM_ini['directional'] = False
NEM_ini['interpolation'] = False
#-----------------------------------------------------------------------------
#                                1.2 WAVE CONDITIONS
#-----------------------------------------------------------------------------
NEM_ini['depth'] = [30.0]
if NEM_ini['regular']:
    NEM_ini['T'] = [6., 8.]#define each T as list
    NEM_ini['H'] = [2. ,2.]#define each H as list 
    NEM_ini['deg'] = [0.0]
    NEM_ini['Nf'] = 1  # GAEL why the 'NF' here?
elif NEM_ini['irregular']:
    NEM_ini['Tp'] = [6., 8.]
    NEM_ini['T'] = []
    NEM_ini['Hs'] = [2, 2]      
    NEM_ini['spectra'] = 'PM' #Define the type of Spectrum you want to use 'JS' for JONSWAP or 'PM' for Pierson-Moskovitz
    NEM_ini['fini'] = [1/(Tp+Tp*2/3) for Tp in NEM_ini['Tp']]
    NEM_ini['fend'] = [1/(Tp-Tp*1/3) for Tp in NEM_ini['Tp']]
    NEM_ini['Nf'] = 20
    NEM_ini['deg'] =  [[0.0],[0.0]]#[[0.0],[0.0]] denest list
    if NEM_ini['directional']:
        NEM_ini['deg_dir'] = "E:\IRREGULAR_WAVES_CLEAN_SCRIPTS\SHORT\EB\Tp_01.26_dx_0.08_dout_0.16_nf_300_ndir_50_t_700_depth_0.7"
        NEM_ini['deg'] = []
        NEM_ini['deg_main'] = [0.0 , 30.0]
        NEM_ini['s1'] = [15.8,15.8]

nsu.create_input(NEM_ini)
#-----------------------------------------------------------------------------
#                              1.3 GRID
#-----------------------------------------------------------------------------  
NEM_GRID = {}
NEM_GRID['Lg'] = 400.0 #lenght and width of the basin 
NEM_GRID['Wg'] = 400.0 #width of the basin
NEM_GRID['dx'] = 2.0#introduce the target deltax for NEMOH simulation
NEM_GRID['dy'] = 2.0#introduce the target deltay for NEMOH simulation
NEM_GRID['Nx'] = int(np.round(NEM_GRID['Lg']/NEM_GRID['dx']) + 1) #redifine number of gridpoints as Nxg = np.round(Lg/dx) + 1 Maximum number of GRID CELLS IN NEMOH 200
NEM_GRID['Ny'] = int(np.round(NEM_GRID['Wg']/NEM_GRID['dy']) + 1) #redifine number of gridpoints as Nyg = np.round(Wg/dy) + 1
if NEM_ini['interpolation']:
    NEM_GRID['dxN'] = 1
    NEM_GRID['dyN'] = 1
#-----------------------------------------------------------------------------
#                              1.3 BODY
#-----------------------------------------------------------------------------  
NEM_BODY = {}
NEM_BODY['diameter']=10.
NEM_BODY['draft']=2.
NEM_BODY['nbody'] = 3 #number of bodies
NEM_BODY['dof'] = [0,0,1,0,0,0]#number of degrees of freedom for each simulation
NEM_BODY['ndof'] = np.sum(NEM_BODY['dof'])
##
#NEM_BODY['xBody'] = [0 ]# #x-coordinates of each body (0,0) center of the domain
#NEM_BODY['yBody'] = [ 0] #y-coordinates of each body (0,0) center of the domain

#NEM_BODY['xBody'] = [0 , 0]# #x-coordinates of each body (0,0) center of the domain
#NEM_BODY['yBody'] = [-25., 25] #y-coordinates of each body (0,0) center of the domain

NEM_BODY['xBody']= [25,25. , -25.] #x-coordinates of each body (0,0) center of the domain
NEM_BODY['yBody'] = [ 25 ,-25., 0.]#y-coordinates of each body (0,0) center of the domain
cylmesh = [0]*NEM_BODY['nbody']

#ob1 = [0]*NEM_BODY['nbody']
#ob2 = [0]*NEM_BODY['nbody']
NEM_BODY['cG'] = -1.0 #z-coordinate of the gravity center of each body 
nPanels = 200   #number of panels for the NEMOH mesh
nsym = 0 #MESH THE BODY WITH SIMETRY AXIS
NEM_BODY['rho'] = 1025.0
NEM_BODY['PTOtype'] = 'Bpto_L'#Theory_LINEAR = Bpto_L WECSIM_LINEAR = Bpto_wsL WECSIM_HYDRA = Bpto_wsH
NEM_BODY['Bpto'] = 360000 #28.5
NEM_BODY['ptoProp'] = [0.0,NEM_BODY['Bpto'],0.0]#[Mextra,Bpto,Kextra]
#-----------------------------------------------------------------------------
#                              1.3.1 MESHING
#-----------------------------------------------------------------------------
NEMOHDir ='E:/NEM/nemoh_py_complete/NEMOH/'#DIRECTORY WHERE THE Calculation folder containting NEMOH FORTRAN codes is located
os.chdir(os.path.join(NEMOHDir,'Calculation'))#In the directory Calculation all the NEMOH runs for each frequency will be calculated
                       #Regular wave calculations folder is saved with corresponding H and T
                       #The results folder for irregular wave is changed for a folder with name Tp and Hs
if not(os.path.exists('mesh')) :    
    os.mkdir('mesh')
# and f'inresults forlder                
if not(os.path.exists('results')) :    
    os.mkdir('results')
for iB in range(NEM_BODY['nbody']):
    #ob1[iB] = mt.cylinder(0.315,0.1575+0.0082,[NEM_BODY['xBody'][iB],NEM_BODY['yBody'][iB],0.0])#diameter and draft
    #ob2[iB] = mt.hemisphere (0.315,[NEM_BODY['xBody'][iB],NEM_BODY['yBody'][iB],-0.1575-0.0082])#diameter 
    #cylmesh[iB] = mt.Mesh()
    #cylmesh[iB].combineMesh(ob1[iB],ob2[iB])
    cylmesh[iB] = mt.cylinder(NEM_BODY['diameter'],NEM_BODY['draft'],[ NEM_BODY['xBody'][iB],NEM_BODY['yBody'][iB],0.0])#diameter and draft
    if NEM_BODY['nbody'] > 1:
        mt.writeMesh(cylmesh[iB],'./mesh/axisym{0:d}'.format(iB+1))
    else:
        mt.writeMesh(cylmesh[iB],'./mesh/axisym')
        
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
ne.createMeshOpt(NEM_BODY['cG'],nPanels,nsym,rho=NEM_BODY['rho'],g=9.81,nbody=NEM_BODY['nbody'],xG=NEM_BODY['xBody'],yG=NEM_BODY['yBody'])#CALLS THE MESH FUNCION INCLUDED IN NEMOH
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
(NEM_BODY['Mass'],NEM_BODY['Kh']) = ne.calcM(rho=1026.0)
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

#del ob1,ob2,iB,nPanels,nsym,cylmesh
#TODO PLOT THE MESH GEOMETRY BEFORE RUNNING NEMOH
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

#ne.plotmesh(NEM_BODY['nbody'],NEMOHDir,NEM_BODY['diameter']/2.,NEM_BODY['xBody'],NEM_BODY['yBody'])
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
 
#-----------------------------------------------------------------------------
#                              1.4. NEMOH SIMULATION PARAMETERS (ADVANCED OPTIONS)
#-----------------------------------------------------------------------------
# Basic Options (RAO calculation)
nrFreq = 1         #each NEMOH simulation is run with one frequency TODO run NEMOH for several frequencies
# Advanced Options
NEM_advOps = {}
NEM_advOps['rhoW'] = 1025.0      #water density
NEM_advOps ['dirCheck'] =  True#Activate to change wave direction
NEM_advOps['dirStep'] = 1
NEM_advOps['irfCheck']  = False#Activate for IRF calculations
NEM_advOps['irfDur'] = 40.0
NEM_advOps['irfStep'] = 0.01
NEM_advOps['kochCheck'] = False#Activate to calculate Kochin Function
NEM_advOps['kochStart'] = 0.0
NEM_advOps['kochStop'] = 360.0
NEM_advOps['kochStep'] = 24
NEM_advOps['fsCheck'] = True #Activate to calculate free surface elevation
NEM_advOps['fsNx'] = NEM_GRID['Nx']
NEM_advOps['fsNy'] = NEM_GRID['Ny']
NEM_advOps['fsLengthX'] = NEM_GRID['Lg']
NEM_advOps['fsLengthY'] = NEM_GRID['Wg']  
NEM_advOps['Show_Console'] = False #Toggle ON or OFF the WINDOWS CONSOLE
#-----------------------------------------------------------------------------
#                              3. RUNNING NEMOH
#-----------------------------------------------------------------------------
(NEM_OUT) = nsu.runNEM(nrFreq,NEM_ini,NEM_BODY,NEM_advOps,NEM_GRID,NEMOHDir)
#-----------------------------------------------------------------------------
#                              4. GENERATING COUPLING OUTPUTS FROM NEMOH
#-----------------------------------------------------------------------------
NEM = {'NEM_ini' : NEM_ini, 'NEM_GRID' :NEM_GRID, 'NEM_BODY' : NEM_BODY, 'NEM_advOps' : NEM_advOps, 'NEM_OUT' : NEM_OUT}
# set the colormap and centre the colorbar
#-----------------------------------------------------------------------------
#                              5. SAVING THE STRUCTUREZ IN .NPZ FILE
#-----------------------------------------------------------------------------
np.savez(os.path.join(NEM_ini['out_dir'],'nm_'+ NEM_ini['case_dir'][0]),**NEM)
#-----------------------------------------------------------------------------
#                              6. GENERATING PLOTS OF THE ETAs
#-----------------------------------------------------------------------------
# TODO Call plot_nemoh_reg or irreg funcions

