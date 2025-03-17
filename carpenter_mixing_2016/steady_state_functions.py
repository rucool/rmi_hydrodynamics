import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

def Pstr (CDL,CDH,D,l,rho0,u,H):

    A= D*H

    Pstr_L = (rho0*CDL*A*u**3)/(2*l**2)          
    Pstr_H = (rho0*CDH*A*u**3)/(2*l**2)             

    return Pstr_L, Pstr_H



def phi_carpenter(rho,z,dz,max_depth):
    """
    Based on Carpenter et al 2016
    rho - rho profile
    z - depth profile (z oriented upwards +)
    dz - step for interpolation (meters)

    returns pea kj/m^2
    """

    g = 9.8 #m/s

    
    #calculate depth averaged rho
    rho_mix=np.trapz(rho,dx=dz,axis=0)*(1/max_depth)

    drho = (rho_mix-rho)*(z)*g
    phi=np.trapz(drho[::-1],z[::-1],dx=dz,axis=0)/1000

    return phi



def Tau_mix_L (phi,H,Rif,Pstr_L,b):

    tau_mix_L = (phi*H)/(Rif*(Pstr_L)*b)
    tau_mix_L = tau_mix_L/86400


    return tau_mix_L



def Tau_mix_H (phi,H,Rif,Pstr_H,b):

    tau_mix_H = (phi*H)/(Rif*(Pstr_H)*b)
    tau_mix_H = tau_mix_H/86400


    return tau_mix_H

    
    
def plot_model_results(res_df,H):
    ss = res_df[res_df.per=='setup']
    pp = res_df[res_df.per=='peak']
    bb = res_df[res_df.per=='breakdown']
    
    
    fig,ax=plt.subplots(1,2,figsize=(10,4),sharey=True)
    ss_x=ss['Current Velocity'].values
    ss_yL=ss['Tmix CD = 0.35'].values
    ss_yH=ss['Tmix CD 1.0'].values
    
    pp_x=pp['Current Velocity'].values
    pp_yL=pp['Tmix CD = 0.35'].values
    pp_yH=pp['Tmix CD 1.0'].values
    
    bb_x=bb['Current Velocity'].values
    bb_yL=bb['Tmix CD = 0.35'].values
    bb_yH=bb['Tmix CD 1.0'].values
    
    ax[0].plot(ss_x,ss_yL ,linewidth=2.5, marker=".", markersize=20, linestyle="--",label='setup')
    ax[0].plot(pp_x,pp_yL ,linewidth=2.5, marker=".", markersize=20, linestyle="--",label='peak')
    ax[0].plot(bb_x,bb_yL ,linewidth=2.5, marker=".", markersize=20, linestyle="--",label='breakdown')
    # ax[0].legend()
    ax[0].set_title('Low Drag Case')
    ax[0].grid(True, which='both')
    ax[0].set_yscale('log')
    ax[0].set_xlabel('Current Velocity (m/s)')
    ax[0].set_ylabel('Mixing Time Period (Days)')
    
    ax[1].plot(ss_x,ss_yH ,linewidth=2.5, marker=".", markersize=20, linestyle="--",label='setup')
    ax[1].plot(pp_x,pp_yH ,linewidth=2.5, marker=".", markersize=20, linestyle="--",label='peak')
    ax[1].plot(bb_x,bb_yH ,linewidth=2.5, marker=".", markersize=20, linestyle="--",label='breakdown')
    ax[1].legend(ncol=1)
    ax[1].set_title('High Drag Case')
    ax[1].grid(True, which='both')
    ax[1].set_yscale('log')
    ax[1].set_xlabel('Current Velocity (m/s)')
    fig.suptitle('H: '+str(H))
    
    plt.tight_layout(pad=0.5)
    
    plt.savefig('./result_all_H'+str(H)+'.png',bbox_inches='tight')



