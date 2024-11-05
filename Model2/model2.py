import numpy as np
import pandas as pd
import os
from model2_functions import *
import matplotlib.pyplot as plt



pper = ['setup','peak','breakdown']
#initial profiles csv for set up peak and breakdown based on glider profiles
#these intial profiles are 25m
init_profs = pd.read_csv('./carp_pea_init_profs.csv')
res_df =None #initialize none to concat result dataframe to for each period


for ii,per in enumerate(pper):
    print(per)
    init_df = init_profs[init_profs.per==per]
    dz=0.25
    extend=False
    dmax_extend=40 #depth you want profile extended to
    if extend == True:
        # if we want this extended to 40m
        
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

    lidx=np.where(z==mldL_c)
    uidx=np.where(z==mldU_c)
    
    delta_rho = dens[-1]-dens[0]



    CDL = 0.35          # Low Drag Coefficient (dimensionless)
    CDH = 1.0           # High Drage Coefficient (dimensionless)

    D = 11.28           # Diameter of Monopile Turbine foundation (m)
    l = 1000            # Turbine Spaceing (m)

    
    rho0 =1026 #np.trapz(dens,dx=dz,axis=0)/H #1026 #Reference Ocean Density (kg/m^3)
    g = 9.81            # Acceleration Due to Gravity (m/s^2)
    Rif = 0.17  #flux richardson number

    u = np.arange(0.01,0.81,0.01)
    
    for uu in u:
        uu = np.round(uu,3)
        print(uu)
        Pstr_L, Pstr_H = Pstr (CDL,CDH,D,l,rho0,uu,H)
        datasave_Pstr_L.append(Pstr_L)
        datasave_Pstr_H.append(Pstr_H)
        # Rate of change of pycnocline thickness over time (m/s)
        # dbdt = (2 * np.pi * Rif * P_str) / (g * delta_rho * H)  # [m/s]
        dbdt_L,dbdt_H = dbdt(Rif, Pstr_L, Pstr_H, g, delta_rho, H)
        # Mixing time scale based on unsteady model (s)
        t_mix_days_H, t_mix_days_L = t_mix_days(H,dbdt_L,dbdt_H)
        datasave_L_D.append(t_mix_days_L)
        datasave_H_D.append(t_mix_days_H)

    data = {
    'Current Velocity': u,
    'Pstr CD = 0.35': datasave_Pstr_L,
    'Pstr CD = 1.0': datasave_Pstr_H,
    'Tmix CD = 0.35': datasave_L_D,
    'Tmix CD 1.0': datasave_H_D,
    'per':[per]*len(u)}
    
    df = pd.DataFrame(data)
    
    if res_df is None:
        res_df = df
    else:
        res_df = pd.concat([res_df,df])
    res_df.to_csv('./all_results_H'+str(H)+'model2.csv',index=False)
plot_model_results(res_df,H)
