import numpy as np
import pandas as pd
import os
from model1_functions import *
import matplotlib.pyplot as plt




pper = ['setup','peak','breakdown']
#initial profiles csv for set up peak and breakdown based on glider profiles
#these intial profiles are 25m
init_profs = pd.read_csv('./carp_pea_init_profs.csv')

res_df = None #initialize none to concat result dataframe to for each period
extend = False
dmax_extend = 40 #depth you want profile extended to
dz = 0.25


fig,ax=plt.subplots(1,3,figsize=(12,8),sharey=True)


for ii,per in enumerate(pper):
    print(per)
    init_df = init_profs[init_profs.per==per]
    # if we want this extended to 40m
    if extend == True:
        dmin_extend = init_df['z'].iloc[-1].astype(int)+dz
        z_extend = np.arange(dmin_extend,dmax_extend+dz,dz)
        dens_extend = [init_df['dens'].iloc[-1]]*len(z_extend)
        mldL_extend = [init_df['mldL'].iloc[-1]]*len(z_extend)
        mldU_extend = [init_df['mldU'].iloc[-1]]*len(z_extend)
        p_thick_extend = [init_df['p_thick'].iloc[-1]]*len(z_extend)
        pea_extend = [init_df['pea'].iloc[-1]]*len(z_extend)
        per_extend = [init_df['per'].iloc[-1]]*len(z_extend)
        extend_df = pd.DataFrame({'z':z_extend,'dens':dens_extend,'mldU':mldU_extend,'mldL':mldL_extend,'p_thick':p_thick_extend,'pea':pea_extend,'per':per_extend})

        init_df = pd.concat([init_df,extend_df],ignore_index=True)
    

    

    ## model

    datasave_Pstr_H=[]
    datasave_Pstr_L=[]
    datasave_L_D = []
    datasave_H_D = []

    #params from initial compisite profiles
    H = init_df['z'].iloc[-1].astype(int)


    dens = init_df['dens'].values
    mldU_g = init_df['mldU'].iloc[0].astype(int)
    mldL_g = init_df['mldL'].iloc[0].astype(int)

    #change to z up positive like carpenter
    z = np.flipud(init_df['z'].values)
    mldL_c = H-mldL_g
    mldU_c = H-mldU_g

    b=mldU_c-mldL_c #pycnocline thickness
    h = (mldU_c+mldL_c)/2 #height off bottom to center of pycnocline
    
    

    rho_s = init_df['dens'].iloc[0].astype(int) #surface dens
    rho_b = init_df['dens'].iloc[-1].astype(int) #bottom dens

    
    #params replicated from Carpenter et al., 2016
    CDL = 0.35          # Low Drag Coefficient (dimensionless)
    CDH = 1.0           # High Drage Coefficient (dimensionless)

    D = 11.28           # Diameter of Monopile Turbine foundation (m)
    l = 1000            # Turbine Spaceing (m)

    rho0 =1026 #np.trapz(dens,dx=dz,axis=0)/H #1026 #Reference Ocean Density (kg/m^3)
    g = 9.81            # Acceleration Due to Gravity (m/s^2)
    Rif = 0.17  #flux richardson number

    # Calculate Phi (kj/m2)
    phi = phi_carpenter(dens,z,dz,H)
    phi = phi*1000 #bring back to SI units for Pstr calculation

    
    #range of current velocities covering tidal to storm driven magnitudes
    u = np.arange(0.1,0.9,0.1)
    for uu in u:
        uu = np.round(uu,1)
        print(uu)
        Pstr_L, Pstr_H = Pstr (CDL,CDH,D,l,rho0,uu,H)
        datasave_Pstr_L.append(Pstr_L)
        datasave_Pstr_H.append(Pstr_H)
        tau_mix_L = Tau_mix_L (phi,H,Rif,Pstr_L,b)
        tau_mix_H = Tau_mix_H (phi,H,Rif,Pstr_H,b)
        datasave_L_D.append(tau_mix_L)
        datasave_H_D.append(tau_mix_H)

    data = {
        'Current Velocity': u,
        'Pstr CD = 0.35': datasave_Pstr_L,
        'Pstr CD = 1.0': datasave_Pstr_H,
        'Tmix CD = 0.35': datasave_L_D,
        'Tmix CD 1.0': datasave_H_D,
        'h':[h]*len(u),
        'per':[per]*len(u)
    }

    df = pd.DataFrame(data)

    if res_df is None:
        res_df = df
    else:
        res_df = pd.concat([res_df,df])
    res_df.to_csv('./all_results_H'+str(H)+'.csv',index=False)
    

plot_model_results(res_df,H)