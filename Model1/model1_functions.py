import numpy as np
import math


def Pstr (CDL,CDH,D,l,rho0,u,h):

    A= D*h

    Pstr_L = (rho0*CDL*A*u**3)/(2*l**2)          
    Pstr_H = (rho0*CDH*A*u**3)/(2*l**2)             

    return Pstr_L, Pstr_H



def Phi (rho_s,rho_b,g,H3,H2,h):

    z = np.linspace(0,h,h)
    dz = h/h

    rho_s = np.repeat(rho_s,H3)
    rho_p = np.linspace(1022,1026,H2)

    if H3 == 20:
        if h == 20:
            rho = rho_s
            rho_mix = np.repeat(np.mean(rho),h)
        elif h == 25:
            rho = np.concatenate((rho_s,rho_p[0:5]))
            rho_mix = np.repeat(np.mean(rho),h)
        elif h >= 30:
            rho_b = np.repeat(rho_b,(h-(H3+H2)))
            rho = np.concatenate((rho_s,rho_p,rho_b))
            rho_mix = np.repeat(np.mean(rho),h)
    else :
        rho_b = np.repeat(rho_b,(h-(H3+H2)))
        rho = np.concatenate((rho_s,rho_p,rho_b))
        rho_mix = np.repeat(np.mean(rho),h)
        
    phi = (rho_mix-rho)*g*(-z)*dz
    phi = np.sum(phi)/1000
    
    return phi   



def Tau_mix_L (phi,h,Rif,Pstr_L,H2):

    tau_mix_L = (phi*h)/(Rif*(Pstr_L)*H2)
    tau_mix_L = tau_mix_L/86400


    return tau_mix_L



def Tau_mix_H (phi,h,Rif,Pstr_H,H2):

    tau_mix_H = (phi*h)/(Rif*(Pstr_H)*H2)
    tau_mix_H = tau_mix_H/86400


    return tau_mix_H
