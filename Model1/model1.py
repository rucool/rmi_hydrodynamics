## Ryan Cassin      10/03/2023

import numpy as np
import pandas as pd
import os
from NJDPU_Functions import Pstr,Phi,Tau_mix_L,Tau_mix_H

filename = 'Test_Model_1_H25_H315.csv'       # Name of Output .csv File

CDL = 0.35          # Low Drag Coefficient (demensionless)
CDH = 1.0           # High Drage Coefficient (demensionless)

D = 11.28           # Diameter of Monopile Turbine foundation (m)
l = 1000            # Turbine Spaceing (m)

rho0 = 1026         # Reference Ocean Density (kg/m^3)
rho_s = 1022
rho_b = 1026
g = 9.81            # Acceleration Due to Gravity (m/s^2)
Rif = 0.17          # Renolds Flux Number

H3 = 15             # Thickness of Surface Layer (m)
H2 = 5             # Thickness of Pycnocline (m)

H = 60              # Total Water Depth (m)
h = 20              # Min Depth (m)


# Lists to store data
datasave_Pstr_L = []
datasave_Pstr_H = []
datasave_phi = []
datasave_L_D = []
datasave_H_D = []

# Calculate Phi and The Stength of Stiring Pstr
while h <= H:
    phi = Phi (rho_s,rho_b,g,H3,H2,h)
    datasave_phi.append(phi)
    u = 0.1
    while u <= 0.8:
        Pstr_L, Pstr_H = Pstr (CDL,CDH,D,l,rho0,u,h)
        datasave_Pstr_L.append(Pstr_L)
        datasave_Pstr_H.append(Pstr_H)
        u += 0.1
    h += 5

h = 20      # Reset h value

# Calculate the Mixing Time Period (Tau_mix) for Low and High Drag Cases
while h <= H:

    if h > 0 and h <= 20:
        phi = datasave_phi[0]
        phi = phi*1000
        Pstr_L = datasave_Pstr_L[0:8]
        Pstr_H = datasave_Pstr_H[0:8]
        n = 0
        while n <= len(Pstr_L)-1:
            Pstr_L_value = Pstr_L[n]
            Pstr_H_value = Pstr_H[n]
            tau_mix_L = Tau_mix_L (phi,h,Rif,Pstr_L_value,H2)
            tau_mix_H = Tau_mix_H (phi,h,Rif,Pstr_H_value,H2)
            datasave_L_D.append(tau_mix_L)
            datasave_H_D.append(tau_mix_H)
            n += 1

    if h > 20 and h <= 25:
        phi = datasave_phi[1]
        phi = phi*1000
        Pstr_L = datasave_Pstr_L[8:16]
        Pstr_H = datasave_Pstr_H[8:16]
        n = 0
        while n <= len(Pstr_L)-1:
            Pstr_L_value = Pstr_L[n]
            Pstr_H_value = Pstr_H[n]
            tau_mix_L = Tau_mix_L (phi,h,Rif,Pstr_L_value,H2)
            tau_mix_H = Tau_mix_H (phi,h,Rif,Pstr_H_value,H2)
            datasave_L_D.append(tau_mix_L)
            datasave_H_D.append(tau_mix_H)
            n += 1

    if h > 25 and h <= 30:
        phi = datasave_phi[2]
        phi = phi*1000
        Pstr_L = datasave_Pstr_L[16:24]
        Pstr_H = datasave_Pstr_H[16:24]
        n = 0
        while n <= len(Pstr_L)-1:
            Pstr_L_value = Pstr_L[n]
            Pstr_H_value = Pstr_H[n]
            tau_mix_L = Tau_mix_L (phi,h,Rif,Pstr_L_value,H2)
            tau_mix_H = Tau_mix_H (phi,h,Rif,Pstr_H_value,H2)
            datasave_L_D.append(tau_mix_L)
            datasave_H_D.append(tau_mix_H)
            n += 1

    if h > 30 and h <= 35:
        phi = datasave_phi[3]
        phi = phi*1000
        Pstr_L = datasave_Pstr_L[24:32]
        Pstr_H = datasave_Pstr_H[24:32]
        n = 0
        while n <= len(Pstr_L)-1:
            Pstr_L_value = Pstr_L[n]
            Pstr_H_value = Pstr_H[n]
            tau_mix_L = Tau_mix_L (phi,h,Rif,Pstr_L_value,H2)
            tau_mix_H = Tau_mix_H (phi,h,Rif,Pstr_H_value,H2)
            datasave_L_D.append(tau_mix_L)
            datasave_H_D.append(tau_mix_H)
            n += 1

    if h > 35 and h <= 40:
        phi = datasave_phi[4]
        phi = phi*1000
        Pstr_L = datasave_Pstr_L[32:40]
        Pstr_H = datasave_Pstr_H[32:40]
        n = 0
        while n <= len(Pstr_L)-1:
            Pstr_L_value = Pstr_L[n]
            Pstr_H_value = Pstr_H[n]
            tau_mix_L = Tau_mix_L (phi,h,Rif,Pstr_L_value,H2)
            tau_mix_H = Tau_mix_H (phi,h,Rif,Pstr_H_value,H2)
            datasave_L_D.append(tau_mix_L)
            datasave_H_D.append(tau_mix_H)
            n += 1

    if h > 40 and h <= 45:
        phi = datasave_phi[5]
        phi = phi*1000
        Pstr_L = datasave_Pstr_L[40:48]
        Pstr_H = datasave_Pstr_H[40:48]
        n = 0
        while n <= len(Pstr_L)-1:
            Pstr_L_value = Pstr_L[n]
            Pstr_H_value = Pstr_H[n]
            tau_mix_L = Tau_mix_L (phi,h,Rif,Pstr_L_value,H2)
            tau_mix_H = Tau_mix_H (phi,h,Rif,Pstr_H_value,H2)
            datasave_L_D.append(tau_mix_L)
            datasave_H_D.append(tau_mix_H)
            n += 1

    if h > 45 and h <= 50:
        phi = datasave_phi[6]
        phi = phi*1000
        Pstr_L = datasave_Pstr_L[48:56]
        Pstr_H = datasave_Pstr_H[48:56]
        n = 0
        while n <= len(Pstr_L)-1:
            Pstr_L_value = Pstr_L[n]
            Pstr_H_value = Pstr_H[n]
            tau_mix_L = Tau_mix_L (phi,h,Rif,Pstr_L_value,H2)
            tau_mix_H = Tau_mix_H (phi,h,Rif,Pstr_H_value,H2)
            datasave_L_D.append(tau_mix_L)
            datasave_H_D.append(tau_mix_H)
            n += 1

    if h > 50 and h <= 55:
        phi = datasave_phi[7]
        phi = phi*1000
        Pstr_L = datasave_Pstr_L[56:64]
        Pstr_H = datasave_Pstr_H[56:64]
        n = 0
        while n <= len(Pstr_L)-1:
            Pstr_L_value = Pstr_L[n]
            Pstr_H_value = Pstr_H[n]
            tau_mix_L = Tau_mix_L (phi,h,Rif,Pstr_L_value,H2)
            tau_mix_H = Tau_mix_H (phi,h,Rif,Pstr_H_value,H2)
            datasave_L_D.append(tau_mix_L)
            datasave_H_D.append(tau_mix_H)
            n += 1

    if h > 55 and h <= 60:
        phi = datasave_phi[8]
        phi = phi*1000
        Pstr_L = datasave_Pstr_L[64:72]
        Pstr_H = datasave_Pstr_H[64:72]
        n = 0
        while n <= len(Pstr_L)-1:
            Pstr_L_value = Pstr_L[n]
            Pstr_H_value = Pstr_H[n]
            tau_mix_L = Tau_mix_L (phi,h,Rif,Pstr_L_value,H2)
            tau_mix_H = Tau_mix_H (phi,h,Rif,Pstr_H_value,H2)
            datasave_L_D.append(tau_mix_L)
            datasave_H_D.append(tau_mix_H)
            n += 1
            
    h = h+5

# Save Data to .csv File
u = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
u = np.tile(u,(9,1))

data = {
    'Current Velocity': u.flatten(),
    'Pstr CD = 0.35': datasave_Pstr_L,
    'Pstr CD = 1.0': datasave_Pstr_H,
    'Tmix CD = 0.35': datasave_L_D,
    'Tmix CD 1.0': datasave_H_D
}

filepath = os.path.join("./", filename)

df = pd.DataFrame(data)
df.to_csv(filepath,index=False)

print("Model 1 |", filename, "| Run Complete" )

