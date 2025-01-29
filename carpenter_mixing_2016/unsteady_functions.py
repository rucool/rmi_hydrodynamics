import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

def Pstr (CDL,CDH,D,l,rho0,u,H):
    A=D*H

    Pstr_L = (rho0*CDL*A*u**3)/(2*l**2)          
    Pstr_H = (rho0*CDH*A*u**3)/(2*l**2)             

    return Pstr_L, Pstr_H


def dbdt(Rif, Pstr_L, Pstr_H, g, delta_rho, H):
    # dbdt = (2 * np.pi * Rif * P_str) / (g * delta_rho * H)  # [m/s]
    dbdt_L =  (2 * np.pi * Rif * Pstr_L) / (g * delta_rho * H)  # [m/s]
    dbdt_H = (2 * np.pi * Rif * Pstr_H) / (g * delta_rho * H)  # [m/s]
    return dbdt_L, dbdt_H


def t_mix_days(H,dbdt_L,dbdt_H):
    t_mix_L = H/dbdt_L
    t_mix_days_L = t_mix_L / (3600*24)
    t_mix_H = H/dbdt_H
    t_mix_days_H = t_mix_H / (3600*24)
    return t_mix_days_H, t_mix_days_L


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
    
    ax[0].plot(ss_x,ss_yL ,linewidth=2.5, marker=".", markersize=5, linestyle="--",label='setup')
    ax[0].plot(pp_x,pp_yL ,linewidth=2.5, marker=".", markersize=5, linestyle="--",label='peak')
    ax[0].plot(bb_x,bb_yL ,linewidth=2.5, marker=".", markersize=5, linestyle="--",label='breakdown')
    # ax[0].legend()
    ax[0].set_title('Low Drag Case')
    ax[0].grid(True, which='both')
    ax[0].set_yscale('log')
    ax[0].set_xlabel('Current Velocity (m/s)')
    ax[0].set_ylabel('Mixing Time Period (Days)')
    
    ax[1].plot(ss_x,ss_yH ,linewidth=2.5, marker=".", markersize=5, linestyle="--",label='setup')
    ax[1].plot(pp_x,pp_yH ,linewidth=2.5, marker=".", markersize=5, linestyle="--",label='peak')
    ax[1].plot(bb_x,bb_yH ,linewidth=2.5, marker=".", markersize=5, linestyle="--",label='breakdown')
    ax[1].legend(ncol=1)
    ax[1].set_title('High Drag Case')
    ax[1].grid(True, which='both')
    ax[1].set_yscale('log')
    ax[1].set_xlabel('Current Velocity (m/s)')
    fig.suptitle('H: '+str(H))
    
    plt.tight_layout(pad=0.5)
    
    plt.savefig('./result_all_H'+str(H)+'.png',bbox_inches='tight')
