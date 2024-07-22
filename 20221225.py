#!/usr/bin/env python
# coding: utf-8

# In[15]:


#分层模拟，模拟北滘2.6-13.9m 2022.12.25 化学反应和基因都存在于孔隙水，TOC存在于沉积物(为wt%)数据来源2021.11.9于胜超
import numpy as np
import xlrd, xlwt
import pandas as pd
import math
import os
import matplotlib.pyplot as plt
import pdb


# define porosity dist.==========
def por(depth):                   #BJ
    beta = 0.008  # cm-1
    por_0 = 0.55
    por_00 = 0.3
    # porosity= 0.30 + 0.25 * 2.71828**(-depth/100)
    # porosity= 0.29 + 0.36 * 2.71828**(-depth/7)
    porosity = por_0 + (por_00 - por_0) * exp((-depth * beta))
    return porosity

# define diffusion coefficient==============【m2 yr-1】
def De(depth, Ions): 
    #Ds=[[61.92,96.24],[183.12,284.4],[176.64,274.32],[175.92,272.88],[99.12,153.6],[194.4,301.68]] #cm2/yr
    #if 17.6<depth<=20.7 or 21.2<depth<=24.3:
        #eff_D=Ds[Ions-1][0]        
    #else:
        #eff_D=Ds[Ions-1][1] 
    
    # Dsw=np.array([6.7,22.9,19.8,19.1,19.0,10.7,21.0,500])*1e-6 * 60 * 60 * 24 * 365 #【cm2/year】    
    Dsw = diff()
    eff_D= Dsw[Ions] / (1-2* np.log(por(depth)))
    return eff_D* 10**(-4)

# Dsw f(Temp, Sality, Presure) 
def diff(Temp = 25, Sality = 35, Presure = 1.013253):
    Dmol_CO2 = 292.2801   # 【cm2/year】
    Dmol_HCO3 = 152.2161
    Dmol_O2 = 371.2064
    Dmol_Mn = 95.66153
    Dmol_NO2 = 309.8471
    Dmol_NH4 = 285.7813
    Dmol_SO4 = 146.8013
    Dmol_H2S = 255.6641
    Dmol_PO4 = 78.81548
    Dsw=np.array([Dmol_c6, Dmol_o2, Dmol_nh4, Dmol_no2, Dmol_no3, Dmol_so4, Dmol_h2s])
    return Dsw






# define 沉积物埋藏速率==========【m/s】
def Solid_rate(depth):                   #BJ 
    Vs_0 = 1.96 * 10**(-3)   
    Vs_00 = 1.96 * 10**(-3)   
    Vs = Vs_00 * (1-por(1e100)) / (1-por(depth))  # Chuang P C et al.(2019); Rui Zhao et al. (2016)
    return Vs * 1e-2/86400/365


# define 孔隙水垂向流速==========【m/s】
def Water_rate(depth):                   #BJ
    Vp_0 = 1.72   
    Vp_00 = 1.72
    Vs_0 = 1.96 * 10**(-3)   
    Vs_00 = 1.96 * 10**(-3)    
    
    Vp = Vp_00 * por(1e100) / por(depth)                # Rui Zhao et al. (2016)

    # Vp = (por(1e10)*Vs_00 - por(0)*Vp_0) /por(depth)  # Chuang P C et al.（2019）
    return Vp * 1e-2/86400/365


# define Gibbs free energy =======【KJ/mol】
def Gibbs_G( Concentration ):
    D_G = np.zeros((10, 1), dtype=np.float64)
    # 1：C_co2, 2：C_hco3, 3：C_c6, 4：C_H, 5：C_no2, 6：C_no3, 7：C_n2, 8：C_nh4, 9：C_h2s, 10：C_so4, 11:C_DIC,  12:C_o2
    C_co2 = Concentration[0]; C_hco3 = Concentration[1]; C_c6 = Concentration[2]; C_H  = Concentration[3]; C_no2 = Concentration[4]
    C_no3 = Concentration[5];  C_n2  = Concentration[6]; C_nh4= Concentration[7]; C_h2s= Concentration[8]; C_so4 = Concentration[9]
    C_DIC = Concentration[10]; C_o2  = Concentration[11]
    D_G[0] = -479.0 + 8.31 * 298.15 * 1e-3  * np.log( C_co2 / C_c6**(1/6) / C_o2 )                                    #DeltaG_cox 
    D_G[1] = -321.6 + 8.31 * 298.15 * 1e-3  * np.log( C_co2 * C_no2**2    / C_c6**(1/6) / C_no3**2)                  #=DeltaG_narg
    D_G[2] = -594.1 + 8.31 * 298.15 * 1e-3  * np.log( C_co2 * C_n2**(2/3) / C_c6**(1/6) / C_no2**(4/3) / C_H**(4/3))#=DeltaG_nirk
    D_G[3] = -352.4 + 8.31 * 298.15 * 1e-3  * np.log( C_co2 * C_nh4**(2/3)/ C_c6**(1/6) / C_no2**(2/3) / C_H**(4/3))#=DeltaG_nrf
    D_G[4] = -76.1  + 8.31 * 298.15 * 1e-3  * np.log( C_hco3* C_h2s**(1/2)/ C_c6**(1/6) / C_so4**(1/2))             #=DeltaG_dsr
    D_G[5] = -189.9 + 8.31 * 298.15 * 1e-3  * np.log( C_no2 * C_H**2      / C_nh4       / C_o2**(3/2))              #=DeltaG_amoA
    D_G[6] = -362.7 + 8.31 * 298.15 * 1e-3  * np.log( C_n2  / C_nh4       / C_no2)                                  #=DeltaG_hzo
    D_G[7] = -201.5 + 8.31 * 298.15 * 1e-3  * np.log( C_so4**(1/4) * C_no2  * C_H**(1/2)  / C_h2s**(1/4) / C_no3)             #=DeltaG_nap
    D_G[8] = -157.4 + 8.31 * 298.15 * 1e-3  * np.log( C_no3**2 / C_no2**2       / C_o2)                           #=DeltaG_nor
    D_G[9] = -401.8 + 8.31 * 298.15 * 1e-3  * np.log( C_so4 / C_co2**2   / C_o2**2 / C_h2s / C_hco3**2)           #=DeltaG_sox
    return D_G


# define A function of free energy yield grams of biomass Y.==========【g/mol/donor】
def Biomass_Y(D_G):
    Y_G = np.zeros((10, 1), dtype=np.float64)
    # 1：cox,      2：narG,   3：nirk,   4：nrf,    5：dsr,     6：amoA,     7：hzo,    8：nap        9：nor      10：sox
    Y_G[0] = 2.08 - 0.0211 * D_G[0] * 6
    Y_G[1] = 2.08 - 0.0211 * D_G[1] * 6
    Y_G[2] = 2.08 - 0.0211 * D_G[2] * 6
    Y_G[3] = 2.08 - 0.0211 * D_G[3] * 6
    Y_G[4] = 2.08 - 0.0211 * D_G[4] * 6
    Y_G[5] = 2.08 - 0.0211 * D_G[5]
    Y_G[6] = 2.08 - 0.0211 * D_G[6]
    Y_G[7] = 2.08 - 0.0211 * D_G[7]
    Y_G[8] = 2.08 - 0.0211 * D_G[8]
    Y_G[9] = 2.08 - 0.0211 * D_G[9]
    return Y_G

# define thermodynamic potential factor F_T==========【-】
def Thermo_T(D_G):
    F_t = np.zeros((10, 1), dtype=np.float64)
    F_t = 1.0 / ( np.exp((D_G + 96.485*0.12) /(8.31 * 1e-3 * 298.15) ) + 1)
    return F_t



# define the rate of gene production Ri==========【Copies/L/s】
def Rate_R(taf, F_t, Concentration):
    G_R = np.zeros((10, 1), dtype=np.float64)
    # Concentration 1：C_co2, 2：C_hco3, 3：C_c6, 4：C_H, 5：C_no2, 6：C_no3, 7：C_n2, 8：C_nh4, 9：C_h2s, 10：C_so4, 11:C_DIC,  12:C_o2
    C_co2 = Concentration[0]; C_hco3 = Concentration[1]; C_c6 = Concentration[2]; C_H  = Concentration[3]; C_no2 = Concentration[4]
    C_no3 = Concentration[5];  C_n2  = Concentration[6]; C_nh4= Concentration[7]; C_h2s= Concentration[8]; C_so4 = Concentration[9]
    C_DIC = Concentration[10]; C_o2  = Concentration[11];
    mu=np.array([0.28, 0.151, 0.247, 0.162, 0.0636, 0.432, 0.864, 0.864, 0.432, 0.864])/86400         #s-1
    kk=np.array([[0.7e-6,0.121e-6], [0.7e-6,0.3e-6], [0.7e-6,0.3e-6], [0.7e-6,0.3e-6], [0.7e-6,3e-6,15e-6], [107e-6,18.75e-6], [5e-6,5e-6,0.2e-6], [0.121e-6,0.121e-6], [64.3e-6,16.9e-6], [0.121e-6,0.121e-6]])
    G_R[0] = taf[0] * F_t[0] * mu[0] * (C_o2 / ( C_o2 + kk[0][1] )) * ( C_c6 /(C_c6 + kk[0][0]))                         #1.cox           
    G_R[1] = taf[1] * F_t[1] * mu[1] * (C_no3 / ( C_no3 + kk[1][1] )) * ( C_c6 /(C_c6 + kk[1][0]))                       #2.narG            
    G_R[2] = taf[2] * F_t[2] * mu[2] * (C_no2 / ( C_no2 + kk[2][1] )) * ( C_c6 /(C_c6 + kk[2][0]))                       #3.nirk
    G_R[3] = taf[3] * F_t[3] * mu[3] * (C_no2 / ( C_no2 + kk[3][1] )) * ( C_c6 /(C_c6 + kk[3][0]))                       #4.nrf
    G_R[4] = taf[4] * F_t[4] * mu[4] * (C_so4 / ( C_so4 + kk[4][1] )) * ( C_c6 /(C_c6 + kk[4][0]))  * (C_o2 /( C_o2+ kk[4][2]))   #5.dsr
    G_R[5] = taf[5] * F_t[5] * mu[5] * (C_nh4 / ( C_nh4 + kk[5][0] )) * ( C_o2 /(C_o2 + kk[5][1]))                                #6.amoA
    G_R[6] = taf[6] * F_t[6] * mu[6] * (C_nh4 / ( C_nh4 + kk[6][0] )) * ( C_no2 /(C_no2 + kk[6][1])) * ( C_o2 /(C_o2 + kk[6][2]))  #7.hzo
    G_R[7] = taf[7] * F_t[7] * mu[7] * (C_no3 / ( C_no3 + kk[7][0] )) * ( C_h2s /(C_h2s + kk[7][1]))                               #8.nap 
    G_R[8] = taf[8] * F_t[8] * mu[8] * (C_no2 / ( C_no2 + kk[8][0] )) * ( C_o2 / ( C_o2 + kk[8][1] ))                              #9.nor
    G_R[9] = taf[9] * F_t[9] * mu[9] * (C_h2s / ( C_h2s + kk[9][0]))  * ( C_o2 / ( C_o2 + kk[9][1] ))                              #10.sox
    return G_R

# define the Organic Matter Remineralization ==========【wt.%/s】    
def DOC_Rc(C_DIC, C_c6, TOC_age):
    Kc = 40 / 1e3   # mol/Kg
    Rc = Kc / (C_DIC + Kc) * 0.16 * TOC_age**(-0.95) * C_c6        
    return Rc/365/86400

# define sediment age.==========【Year】
def age(depth):                  
    if depth <= 31.2:
        ag = 8.8                 
    elif  31.2 < depth <= 37.6:   
        ag = 20
    elif  37.6 < depth <= 40:
        ag = 30
    else:
        ag = 30
    return ag*1e3

# define the rate of DOC production by microbial mortality ==========【mol/L/s   &&  gene/g/s】
def DOC_Rm(taf):
    lamda = 0.001 /86400           #s-1
    Ni  = 3.75 * 1e13    #[genes/g]
    Rm = 0.4276  * lamda / 12 * (taf[0]/Ni + taf[1]/Ni + taf[2]/Ni + taf[3]/Ni + taf[4]/Ni + taf[5]/Ni + taf[6]/Ni/8 + taf[7]/Ni + taf[8]/Ni + taf[9]/Ni)
    return Rm 

# define Nitrogen as a nutrient. ==========【mol/L/s】
def Nit_U(taf, GeneR, C_nh4, C_no2):
    Un = np.zeros((3, 1), dtype=np.float64)
    K_nutrient = 134/1e9  #mol
    Ni = 3.75 * 1e13    #[genes/g]
    mu=np.array([0.28, 0.151, 0.247, 0.162, 0.0636, 0.432, 0.864, 0.864, 0.432, 0.864])/86400         #s-1
    uu = 0.1127  / 14 * (taf[0] * mu[0] /Ni + taf[1] * mu[1] /Ni + taf[2] * mu[2] /Ni + taf[3] * mu[3] /Ni + taf[4] * mu[4] /Ni+ taf[5] * mu[5] /Ni + taf[6] * mu[6] /Ni/8 + taf[7] * mu[7] /Ni + taf[8] * mu[8] /Ni + taf[9] * mu[9] /Ni)
    uu = 0.00000026
    Un[1] = (uu - Un[0]) * (C_no2 / (C_no2 + K_nutrient))
    Un[2] = uu - Un[0] - Un[1]
    return Un

# 定义生物代谢产生的化学物质. ==========【mol/L/s】
def R_nyg(GR, YG, DetG):
    Rny = np.zeros((10, 1), dtype=np.float64)
    Ni = 3.75 * 1e13    #[genes/g]
    #Gene 1：cox,      2：narG,   3：nirk,   4：nrf,    5：dsr,     6：amoA,     7：hzo,    8：nap        9：nor      10：sox    
    for i in range(0,10):
        Rny[i] = GR[i]/Ni/YG[i]
    return Rny

# 计算右端项生物代谢贡献的化学物质总和. ==========【mol/L/s】
def Sum_R(Concentration, taf, depth):
    SumR = np.zeros((7, 1), dtype=np.float64)
    Ni  = 3.75 * 1e13    #[genes/g]
    # 【mol/L】1：C_co2, 2：C_hco3, 3：C_c6, 4：C_H, 5：C_no2, 6：C_no3, 7：C_n2, 8：C_nh4, 9：C_h2s, 10：C_so4, 11:C_DIC,  12:C_o2
    C_co2 = Concentration[0]; C_hco3 = Concentration[1]; C_c6 = Concentration[2]; C_H  = Concentration[3]; C_no2 = Concentration[4]
    C_no3 = Concentration[5];  C_n2  = Concentration[6]; C_nh4= Concentration[7]; C_h2s= Concentration[8]; C_so4 = Concentration[9]
    C_DIC = Concentration[10]; C_o2  = Concentration[11]; C_POC  = Concentration[12]
    Con   = [C_co2, C_hco3, C_c6, C_H, C_no2, C_no3, C_n2, C_nh4, C_h2s, C_so4, C_DIC, C_o2]   
    DetG = Gibbs_G( Con )
    YG = Biomass_Y( DetG )
    Ft = Thermo_T( DetG )
    GR = Rate_R(taf, Ft, Con)
    Rc = DOC_Rc(C_DIC, C_POC, age(depth))
    Rn = 16/106 * Rc
    Rm  = DOC_Rm(taf)
    Un  = Nit_U(taf, GR, C_nh4, C_no2)
    Rnyg= R_nyg(GR, YG, DetG)
    #Gene 0：cox,      1：narG,   2：nirk,   3：nrf,    4：dsr,     5：amoA,     6：hzo,    7：nap        8：nor      9：sox    
    SumR[0] =  Rc/6 + Rm/6 - Rnyg[0] - Rnyg[1] - Rnyg[2] - Rnyg[3] - Rnyg[4]         
    #SumR[0] = -Rc      ##########################20200117沉积物
    SumR[1] = -6*Rnyg[0] - 1.5*Rnyg[5] - Rnyg[8] - 2*Rnyg[9]
    SumR[2] = Rn + 4*Rnyg[3] - Rnyg[5] - Rnyg[6] - Un[0]
    SumR[3] = 12*Rnyg[1] + Rnyg[5] +Rnyg[7] - Rnyg[6] - 2*Rnyg[8]- 8*Rnyg[2] - 4*Rnyg[3] - Un[1]
    SumR[4] = 2*Rnyg[8] - 12*Rnyg[1] - Rnyg[7] - Un[2]
    SumR[5] = Rnyg[7]/4 + Rnyg[9] - 3*Rnyg[4]
    SumR[6] = -Rnyg[7]/4 - Rnyg[9] + 3*Rnyg[4]
    return SumR 

def Rate_R2(Concentration, taf, depth):     #基因产生速率 ==========【cpoies/L/s】
    GR = np.zeros((10, 1), dtype=np.float64)
    Ni  = 3.75 * 10**13    #[genes/g]
    # 1：C_co2, 2：C_hco3, 3：C_c6, 4：C_H, 5：C_no2, 6：C_no3, 7：C_n2, 8：C_nh4, 9：C_h2s, 10：C_so4, 11:C_DIC,  12:C_o2
    C_co2 = Concentration[0]; C_hco3 = Concentration[1]; C_c6 = Concentration[2]; C_H  = Concentration[3]; C_no2 = Concentration[4]
    C_no3 = Concentration[5];  C_n2  = Concentration[6]; C_nh4= Concentration[7]; C_h2s= Concentration[8]; C_so4 = Concentration[9]
    C_DIC = Concentration[10]; C_o2  = Concentration[11];
    Con   = [C_co2, C_hco3, C_c6, C_H, C_no2, C_no3, C_n2, C_nh4, C_h2s, C_so4, C_DIC, C_o2]   
    DetG = Gibbs_G( Con )
    YG = Biomass_Y( DetG )
    Ft = Thermo_T( DetG )
    GR = Rate_R(taf, Ft, Con)
    return GR


#求解浓度方程组========================【mol/L】
def model_sol2(z, Ions, Concentration, Sum, N, dz):   
    dt = 1
    dz = 0.01
    A = np.zeros((N + 1, N + 1), dtype=np.float64)           #控制方程矩阵和右手矩
    rhs = np.zeros((N + 1, 1), dtype=np.float64)
    C_model =  np.zeros((N + 1, 1), dtype=np.float64)
    if Ions == 0:                         #C_c6
        rhs[0] = 0;    A[0, 0] = 1;    A[0, 1] = -1;       rhs[N] = 0;       A[N, N] = 1     ;    A[N, N - 1] = -1
        #C_model[0] = Sum[0];             C_model[N] = Sum[N];
    
    elif Ions == 1:                      #o2 
        rhs[0] = 0;     A[0, 0] = 1;    A[0, 1] = -1;       rhs[N] = 0;      A[N, N] = 1;       A[N, N - 1] = -1
        #C_model[0] = Sum[0];        C_model[N] = Sum[N];
        
    elif Ions == 2:                       #NH4
        rhs[0] = 0;     A[0, 0] = 1;    A[0, 1] = -1;       rhs[N] = 2706.1;      A[N, N] = 1
        #C_model[0] = Sum[0];        C_model[N] = 0;
        
    elif Ions == 3:                      #NO2
        rhs[0] = 3.75182E-06;     A[0, 0] = 1;         rhs[N] = 1.18286E-06;   A[N, N] = 1
        #C_model[0] = 0;             C_model[N] = 0;
        
    elif Ions == 4:                      #NO3
        rhs[0] = 1.46774E-05;       A[0, 0] = 1     ;    rhs[N] = 9.05582E-06;      A[N, N] = 1;       A[N, N - 1] = -1
        #C_model[0] = 0;             C_model[N] = Sum[N];
        
    elif Ions == 5:                       #SO4
        rhs[0] = 0;     A[0, 0] = 1;    A[0, 1] = -1;       rhs[N] = 0;      A[N, N] = 1;       A[N, N - 1] = -1
        #C_model[0] = Sum[0];        C_model[N] = Sum[N]
        
    elif Ions == 6:                       #H2S 
        #rhs[0] = 0;     A[0, 0] = 1;    A[0, 1] = -1;       rhs[N] = 0      ;   A[N, N] = 1;       A[N, N - 1] = -1
        rhs[0] = 0.19*1e3/34/1e6;    rhs[N] = 0.1*1e3/34/1e6;      A[0, 0] = 1     ;    A[N, N] = 1
        #C_model[0] = Sum[0];        C_model[N] = Sum[N]
        
    else:
        print("未找到该物质！！！")
    
    for i in range(1, N):  
        #A[i, i] =  -De(z[i]-0.5* dh, Ions)/dh**2 -De(z[i]+0.5* dh, Ions)/dh**2  - 1/dt
        #A[i, i - 1] = De(z[i]-0.5* dh, Ions)/dh**2  + Water_rate(z[i]) /(2*dh)
        #A[i, i + 1] = De(z[i]+0.5* dh, Ions)/dh**2  - Water_rate(z[i]) /(2*dh)
        #rhs[i] = -Sum[i] -Concentration[i]/dt 
        
        #A[i, i] =  -(De(z[i-1], Ions)*De(z[i], Ions))**0.5 * (por(z[i-1])*por(z[i]))**0.5 /dz**2 -(De(z[i+1], Ions)*De(z[i], Ions))**0.5 * (por(z[i+1])*por(z[i]))**0.5/dz**2  
        #A[i, i - 1] = (De(z[i-1], Ions)*De(z[i], Ions))**0.5 * (por(z[i-1])*por(z[i]))**0.5 /dz**2  + Water_rate(z[i-1])*por(z[i-1]) /(2*dz)
        #A[i, i + 1] = (De(z[i+1], Ions)*De(z[i], Ions))**0.5 * (por(z[i+1])*por(z[i]))**0.5 /dz**2 - Water_rate(z[i+1])*por(z[i+1]) /(2*dz)
        #rhs[i] = Sum[i]*por(z[i])
        #C_model[i] = (A[i, i - 1]*Concentration[i-1] + A[i, i]*Concentration[i] + A[i, i+1]*Concentration[i+1] + rhs[i]) / por(z[i])
        phiD_i     = por(z[i])* De(z[i], Ions)
        phiv_i     = por(z[i])* Water_rate(z[i])
        phiv_i_next  = por(z[i+1])* Water_rate(z[i+1])
        phiv_i_last  = por(z[i-1])* Water_rate(z[i-1])
        phiD_i_low =  (dz+dz) / (dz/(De(z[i], Ions)* por(z[i])) + dz/(De(z[i-1], Ions)* por(z[i-1])))  
        phiD_i_up =   (dz+dz) / (dz/(De(z[i], Ions)* por(z[i])) + dz/(De(z[i+1], Ions)* por(z[i+1]))) 
        phiv_i_low = (dz+dz) / (dz/(Water_rate(z[i])* por(z[i])) + dz/(Water_rate(z[i-1])* por(z[i-1])))  
        phiv_i_up =  (dz+dz) / (dz/(Water_rate(z[i])* por(z[i])) + dz/(Water_rate(z[i+1])* por(z[i+1])))
#         D_i_up = phiD_i_up*(Concentration[i+1]-Concentration[i])/(dz+dz)
#         D_i_low = phiD_i_low*(Concentration[i]-Concentration[i-1])/(dz+dz)
#         C_i_up = Concentration[i] + phiv_i_up*(Concentration[i+1]-Concentration[i])/(dz+dz)*dz/(por(z[i])*Water_rate(z[i]))
#         C_i_low = Concentration[i-1] + phiv_i_low*(Concentration[i]-Concentration[i-1])/(dz+dz)*dz/(por(z[i-1])*Water_rate(z[i-1]))
#         C_model[i] = 2/(por(z[i])*dz)*(D_i_up-D_i_low)-(phiv_i_up*C_i_up-phiv_i_low*C_i_low)/(dz*por(z[i])) + Sum[i]
        A[i, i] =  2/(dz*por(z[i]))*(-phiD_i_up/(dz+dz)-phiD_i_low/(dz+dz)) - 1/(dz*por(z[i]))*(phiv_i_up-phiv_i_up**2/phiv_i*dz/(dz+dz)) + 1/(dz*por(z[i]))*(phiv_i_low**2/phiv_i_last*dz/(dz+dz)) - 1/dt
        A[i, i - 1] = 2/(dz*por(z[i]))*phiD_i_low/(dz+dz) + 1/(dz*por(z[i]))*(phiv_i_low-phiv_i_low**2/phiv_i_last*dz/(dz+dz))
        A[i, i + 1] = 2/(dz*por(z[i]))*phiD_i_up/(dz+dz) - 1/(dz*por(z[i]))*(phiv_i_up**2/phiv_i*dz/(dz+dz))
        rhs[i] = -Sum[i] -Concentration[i]/dt

        
        
         
    # y=solve(A,b)===============
    C_model = np.dot(np.linalg.inv(A), rhs)
    return C_model

#求解基因方程组========================【copies/L】
def model_sol3(Rct, gene, geneR, z, N):
    lamda = 0.001/86400     #s-1 
    Ions = 7
    dh = 0.01
    dt = 1
    A = np.zeros((N + 1, N + 1), dtype=np.float64)           #控制方程矩阵和右手矩
    rhs = np.zeros((N + 1, 1), dtype=np.float64)
    R_model =  np.zeros((N + 1, 1), dtype=np.float64)
    
    if Rct==5:                             
        # evaluate the rhs vector============
        rhs[0] = 10840000;   rhs[N] = 714000;         A[0 , 0] = 1;      A[N , N ] = 1 
        #R_model[0] = 0;        R_model[N] = 0;
        
    elif Rct==6:                          
        # evaluate the rhs vector============
        rhs[0] = 2700000;   rhs[N] = 7870000;         A[0 , 0] = 1;      A[N , N ] = 1        
        #R_model[0] = 0;        R_model[N] = 0;
        
    else:
        # evaluate the rhs vector============
        rhs[0] = 0; rhs[N] = 0
        # evaluate the global matrix A============
        A[0, 0] = 1; A[0, 1] = -1
        A[N, N] = 1; A[N, N - 1] = -1
        #R_model[0] = -lamda * por(z[0])* gene[0] + geneR[0]*por(z[0]);       
        #R_model[N] = -lamda * por(z[N])* gene[N] + geneR[N]*por(z[N]);
        
        
    for i in range(1, N): 
#         A[i, i] = -lamda * por(z[i])
#         A[i, i + 1] = - Solid_rate(z[i+1])*por(z[i+1])/ (2*dz) 
#         A[i, i - 1] = Solid_rate(z[i-1])*por(z[i-1]) / (2*dz) 
#         rhs[i]  = geneR[i]*por(z[i])
#         R_model[i] = (A[i, i - 1]*gene[i-1] + A[i, i]*gene[i] + A[i, i+1]*gene[i+1] + rhs[i]) / por(z[i])
        

        #假设所有基因分布于孔隙水中
#         A[i, i] =  -(De(z[i-1], Ions)*De(z[i], Ions))**0.5 * (por(z[i-1])*por(z[i]))**0.5 /dh**2 -(De(z[i+1], Ions)*De(z[i], Ions))**0.5 * (por(z[i+1])*por(z[i]))**0.5/dh**2 
#         A[i, i - 1] = (De(z[i-1], Ions)*De(z[i], Ions))**0.5 * (por(z[i-1])*por(z[i]))**0.5 /dh**2  + Water_rate(z[i-1])*por(z[i-1]) /(2*dh)
#         A[i, i + 1] =(De(z[i+1], Ions)*De(z[i], Ions))**0.5 * (por(z[i+1])*por(z[i]))**0.5/dh**2 - Water_rate(z[i+1])*por(z[i+1]) /(2*dh)
#         rhs[i] = (geneR[i]-lamda * gene[i])*por(z[i])
#         R_model[i] = (A[i, i - 1]*gene[i-1] + A[i, i]*gene[i] + A[i, i+1]*gene[i+1] + rhs[i]) / por(z[i])
        phiD_i     = por(z[i])* De(z[i], Ions)
        phiv_i     = por(z[i])* Water_rate(z[i])
        phiv_i_next  = por(z[i+1])* Water_rate(z[i+1])
        phiv_i_last  = por(z[i-1])* Water_rate(z[i-1])
        phiD_i_low =  (dz+dz) / (dz/(De(z[i], Ions)* por(z[i])) + dz/(De(z[i-1], Ions)* por(z[i-1])))  
        phiD_i_up =   (dz+dz) / (dz/(De(z[i], Ions)* por(z[i])) + dz/(De(z[i+1], Ions)* por(z[i+1]))) 
        phiv_i_low = (dz+dz) / (dz/(Water_rate(z[i])* por(z[i])) + dz/(Water_rate(z[i-1])* por(z[i-1])))  
        phiv_i_up =  (dz+dz) / (dz/(Water_rate(z[i])* por(z[i])) + dz/(Water_rate(z[i+1])* por(z[i+1])))       
#         D_i_up = phiD_i_up*(gene[i+1]-gene[i])/(dz+dz)
#         D_i_low = phiD_i_low*(gene[i]-gene[i-1])/(dz+dz)
#         gene_i_up = gene[i] + phiv_i_up*(gene[i+1]-gene[i])/(dz+dz)*dz/(por(z[i])*Water_rate(z[i]))
#         gene_i_low = gene[i-1] + phiv_i_low*(gene[i]-gene[i-1])/(dz+dz)*dz/(por(z[i-1])*Water_rate(z[i-1]))
#         R_model[i] = 2/(por(z[i])*dz)*(D_i_up-D_i_low)-(phiv_i_up*gene_i_up-phiv_i_low*gene_i_low)/(dz*por(z[i])) + geneR[i]-lamda * gene[i]
        A[i, i] =  2/(dz*por(z[i]))*(-phiD_i_up/(dz+dz)-phiD_i_low/(dz+dz)) - 1/(dz*por(z[i]))*(phiv_i_up-phiv_i_up**2/phiv_i*dz/(dz+dz)) + 1/(dz*por(z[i]))*(phiv_i_low**2/phiv_i_last*dz/(dz+dz)) - 1/dt
        A[i, i - 1] = 2/(dz*por(z[i]))*phiD_i_low/(dz+dz) + 1/(dz*por(z[i]))*(phiv_i_low-phiv_i_low**2/phiv_i_last*dz/(dz+dz))
        A[i, i + 1] = 2/(dz*por(z[i]))*phiD_i_up/(dz+dz) - 1/(dz*por(z[i]))*(phiv_i_up**2/phiv_i*dz/(dz+dz))
        rhs[i] = -(geneR[i]-lamda*gene[i]) -gene[i]/dt        
        
        
    # y=solve(A,b)===============
    R_model = np.dot(np.linalg.inv(A), rhs)    
    return R_model

#求解基因最大误差========================【copies/L】
def Max_error(gene,Tao):
    det_G = np.zeros((10, 1), dtype=np.float64)
    for i in range(0, N+1):
        det_G[0] = max( abs(Gene[i][0] - Tao[i][0]),  det_G[0])
        det_G[1] = max( abs(Gene[i][1] - Tao[i][1]),  det_G[1])
        det_G[2] = max( abs(Gene[i][2] - Tao[i][2]),  det_G[2])
        det_G[3] = max( abs(Gene[i][3] - Tao[i][3]),  det_G[3])
        det_G[4] = max( abs(Gene[i][4] - Tao[i][4]),  det_G[4])
        det_G[5] = max( abs(Gene[i][5] - Tao[i][5]),  det_G[5])
        det_G[6] = max( abs(Gene[i][6] - Tao[i][6]),  det_G[6])
        det_G[7] = max( abs(Gene[i][7] - Tao[i][7]),  det_G[7])
        det_G[8] = max( abs(Gene[i][8] - Tao[i][8]),  det_G[8])
        det_G[9] = max( abs(Gene[i][9] - Tao[i][9]),  det_G[9])    
#         det_G[0] = max( abs(Gene[i][0]),  det_G[0])
#         det_G[1] = max( abs(Gene[i][1]),  det_G[1])
#         det_G[2] = max( abs(Gene[i][2]),  det_G[2])
#         det_G[3] = max( abs(Gene[i][3]),  det_G[3])
#         det_G[4] = max( abs(Gene[i][4]),  det_G[4])
#         det_G[5] = max( abs(Gene[i][5]),  det_G[5])
#         det_G[6] = max( abs(Gene[i][6]),  det_G[6])
#         det_G[7] = max( abs(Gene[i][7]),  det_G[7])
#         det_G[8] = max( abs(Gene[i][8]),  det_G[8])
#         det_G[9] = max( abs(Gene[i][9]),  det_G[9])    
    Maxi = max(det_G[0], det_G[1], det_G[2], det_G[3], det_G[4],det_G[5], det_G[6], det_G[7], det_G[8], det_G[9])
    return Maxi


#############################0.打开文件，读取数据#########################
#f=open("Gene.txt","w")


#############################1.网格剖分#########################
dz = 0.01
z = np.arange(2.25, 27.07, dz)    
N = 2481  
dt = 1
############################2.初始条件【浓度 & 基因】###########    
Concentration = np.zeros((N+1, 13), dtype=np.float64)
Concentration0 = np.zeros((N+1, 13), dtype=np.float64)
Gene = np.zeros((N+1, 10), dtype=np.float64)
Gene0 = np.zeros((N+1, 10), dtype=np.float64)
C_Name=('1: C_co2', '2: C_hco3', '3：C_c6', '4：C_H', '5：C_no2', '6：C_no3', '7：C_n2', '8：C_nh4', '9：C_h2s', '10：C_so4', '11:C_DIC',  '12:C_o2', '13:C_POC')
#民众浓度【mol/Kg】 0：C_co2, 1：C_hco3, 2：C_c6, 3：C_H, 4：C_no2, 5：C_no3, 6：C_n2, 7：C_nh4, 8：C_h2s, 9：C_so4, 10:C_DIC,  11:C_o2 ,    12:C_POC
for i in range(0, N+1):
    Concentration[i][0] = 1e-7
    Concentration[i][1] = 1e-7
    Concentration[i][2] = 1.7406692090160792e-13 *z[i]**10+( -5.13761779665674e-12 *z[i]**9)+( -3.299492097242631e-09 *z[i]**8)+( 4.1389054394304156e-07 *z[i]**7)+( -2.264804351160562e-05 *z[i]**6)+( 0.0006872252874341565 *z[i]**5)+( -0.012198587738491208 *z[i]**4)+( 0.12474859247885309 *z[i]**3)+( -0.686904587275624 *z[i]**2)+( 1.8257871630136067 *z[i])+( -1.2278103388739463 )    
    Concentration[i][3] = 1e-8
    Concentration[i][4] = -2.104057760848451e-18 *z[i]**10+( 6.934420218771268e-16 *z[i]**9)+( -9.767489583058706e-14 *z[i]**8)+( 7.681350630378094e-12 *z[i]**7)+( -3.6985136147738343e-10 *z[i]**6)+( 1.1256255934651154e-08 *z[i]**5)+( -2.1608277180044606e-07 *z[i]**4)+( 2.53764465218647e-06 *z[i]**3)+( -1.7076336967873157e-05 *z[i]**2)+( 5.799423359163791e-05 *z[i])+( -7.258477628876808e-05 )    
    Concentration[i][5] = -4.114737151006369e-18 *z[i]**10+( 1.3972857441617904e-15 *z[i]**9)+( -2.0339629216879802e-13 *z[i]**8)+( 1.6584375993822927e-11 *z[i]**7)+( -8.308714617214002e-10 *z[i]**6)+( 2.6413882877429733e-08 *z[i]**5)+( -5.319320032065825e-07 *z[i]**4)+( 6.585245157764569e-06 *z[i]**3)+( -4.698474269050189e-05 *z[i]**2)+( 0.00017058023658993992 *z[i])+( -0.00022164061478607358 )    
    Concentration[i][6] = 1e-9
    Concentration[i][7] = 4.0079606354710043e-16 *z[i]**10+( -1.3450856104693778e-13 *z[i]**9)+( 1.9426796453294307e-11 *z[i]**8)+( -1.5796148910509377e-09 *z[i]**7)+( 7.943604358946093e-08 *z[i]**6)+( -2.5567281834295015e-06 *z[i]**5)+( 5.274498022060886e-05 *z[i]**4)+( -0.0006798313738710033 *z[i]**3)+( 0.005145444523999974 *z[i]**2)+( -0.020026351877634956 *z[i])+( 0.03693768489388249 )
    Concentration[i][8] = 1e-6
    Concentration[i][9] = -2.5307806879068483e-16 *z[i]**10+( 1.1581302656815528e-13 *z[i]**9)+( -2.144527370431241e-11 *z[i]**8)+( 2.1492115480998486e-09 *z[i]**7)+( -1.29496356885863e-07 *z[i]**6)+( 4.875921300396796e-06 *z[i]**5)+( -0.00011467985708067698 *z[i]**4)+( 0.001631167643256955 *z[i]**3)+( -0.013254034192667165 *z[i]**2)+( 0.0572204118309532 *z[i])+( -0.056397032990941696 )
    Concentration[i][10] = 6.044766538676444e-13 *z[i]**10+( -1.3942290130980215e-10 *z[i]**9)+( 1.3877921778824684e-08 *z[i]**8)+( -7.823565478694807e-07 *z[i]**7)+( 2.7560825478028093e-05 *z[i]**6)+( -0.0006311396704769519 *z[i]**5)+( 0.009453810666922498 *z[i]**4)+( -0.09050594960248813 *z[i]**3)+( 0.5198921497727369 *z[i]**2)+( -1.5549435999957122 *z[i])+( 2.1299886068091514 )    
    Concentration[i][11] = 1e-9
    Concentration[i][12] = (-9e-10*z[i]**6 - 1e-7*z[i]**5 + 4e-6*z[i]**4 + 0.0005*z[i]**3 - 0.0213*z[i]**2 + 0.2553*z[i] + 0.6715)*1e6
    #基因【genes/L】 0：cox,      1：narG,   2：nirk,   3：nrf,    4：dsr,     5：amoA,     6：hzo,    7：nap        8：nor      9：sox
    Gene[i][0] = 1e5
    Gene[i][1] = 1e5
    Gene[i][2] = 1e5
    Gene[i][3] = 1e5
    Gene[i][4] = 1e5
    Gene[i][5] = 1.8149115050146092e-05 *z[i]**10+( -0.0025970258840197304 *z[i]**9)+( 0.12524748342028677 *z[i]**8)+( -0.9848821749028847 *z[i]**7)+( -123.61588309591274 *z[i]**6)+( 4925.4942433373935 *z[i]**5)+( -75163.8547815459 *z[i]**4)+( 321348.1047953144 *z[i]**3)+( 4250441.997581997 *z[i]**2)+( -53197485.408834904 *z[i])+( 169686702.3348549 )
    Gene[i][6] = -3.1496922177991296e-06 *z[i]**10+( 0.00031832652419954723 *z[i]**9)+( -0.004374674229161838 *z[i]**8)+( -0.4849318876896918 *z[i]**7)+( 11.03207967439595 *z[i]**6)+( 729.1492911791356 *z[i]**5)+( -40833.72222230589 *z[i]**4)+( 864073.0902662921 *z[i]**3)+( -9246811.681851482 *z[i]**2)+( 49431498.17036085 *z[i])+( -99993069.60356417 )
    Gene[i][8] = 1e5
    Gene[i][9] = 1e5
Concentration = np.divide(Concentration, 1e6)                        #浓度【mol/L】
Gene = Gene
Concentration0 = Concentration
Gene0 = Gene

Ttime = 0
while 1:
    Ttime = Ttime + dt
    #f.write("当前计算时间"+str(Ttime)+"\n")
    print("当前计算时间", Ttime)
############################3.picard内层迭代##########################     
    iter = 0
    while 1:
        iter = iter + 1
    ######################3.1 重置避免物质浓度 避免<0#################
        for i in range(0, N+1):                                                    
            for j in range(0,13):
                if Concentration[i][j]<0:
                    #print(C_Name[j],z[i],Concentration[i][j])
                    Concentration[i][j] = 1e-10
            for j in range(0,10):
                if Gene[i][j]<0:
                    #print(C_Name[j],z[i],Concentration[i][j])
                    Gene[i][j] = 1e-7

    ######################3.2 浓度求解##########################            
        SumR = np.zeros((7, N+1), dtype=np.float64)
        Sum0=0
        for i in range(0, N+1):
            Sum0 = Sum_R(Concentration[i,:], Gene[i,:], z[i])                           #计算细菌代谢产生的不同化学物质 Ri
            for j in range(0,7):
                SumR[j][i] = Sum0[j]                                                 #【mol/L/s】
        Cc = np.zeros((7, N+1), dtype=np.float64)                                    #计算iter次化学物质浓度,并更新
        Sum1=0;  Sum2=0;  Sum3=0;  Sum4=0;  Sum5=0;  Sum6=0;  Sum7=0
        Sum1= model_sol2(z, 0, Concentration[:,2],   SumR[0,:], N, dz)                           #C6
        Sum2= model_sol2(z, 1, Concentration[:,11],  SumR[1,:], N, dz)                           #O2
        Sum3= model_sol2(z, 2, Concentration[:,7],   SumR[2,:], N, dz)                           #NH4
        Sum4= model_sol2(z, 3, Concentration[:,4],   SumR[3,:], N, dz)                           #NO2
        Sum5= model_sol2(z, 4, Concentration[:,5],   SumR[4,:], N, dz)                           #NO3
        Sum6= model_sol2(z, 5, Concentration[:,9],   SumR[5,:], N, dz)                           #SO4
        Sum7= model_sol2(z, 6, Concentration[:,8],   SumR[6,:], N, dz)                          #H2S      

    ###########################3.3求解基因丰度######################
        GeneR = np.zeros((10, N+1), dtype=np.float64)
        Sum8=0
        for i in range(0, N+1):
            Sum8 = Rate_R2(Concentration[i,:], Gene[i,:], z[i])                    #计算不同层位的Ri，基因生长速率；【Coies/L/s】
            for j in range(0, 10):
                GeneR[j][i] = Sum8[j]
        #print(GeneR[9,:])
        Tao = np.zeros((N+1, 10), dtype=np.float64)
        Sum9=0
        for i in range(0, 10):
            Sum9 = model_sol3(i,  Gene[:,i],  GeneR[i,:], z,  N)
            for j in range(0, N+1):
                Tao[j][i] = Sum9[j]

    ###########################3.3物质浓度更新######################
        for i in range(0, N+1):         
    #         Concentration[i][2]=Concentration[i][2]+Sum1[i]                                    
    #         Concentration[i][11]=Concentration[i][11]+Sum2[i]
    #         Concentration[i][7]=Concentration[i][7]+Sum3[i]
    #         Concentration[i][4]=Concentration[i][4]+Sum4[i]
    #         Concentration[i][5]=Concentration[i][5]+Sum5[i]
    #         Concentration[i][9]=Concentration[i][9]+Sum6[i]
            Concentration[i][2]=Sum1[i]                                    
            Concentration[i][11]=Sum2[i]
            Concentration[i][7]=Sum3[i]
            Concentration[i][4]=Sum4[i]
            Concentration[i][5]=Sum5[i]
            Concentration[i][9]=Sum6[i]        
           # Concentration[i][8]=Sum7[i]   

    ###########################3.4 迭代收敛判断并更新基因丰度#########################
        Maxi = Max_error(Gene,Tao)
    #   Maxi = Max_error(Tao) 
    #   Gene = gene +Tao
        Gene = Tao

        if Maxi<1e-5  or  iter>1000:
            #f.write("The iter ="+str(iter)+",   Maxi error ="+str(Maxi)+"\n")
            print("The iter =", iter,",   Maxi error =", Maxi)   
            break
        else: 
            #f.write("The iter ="+str(iter)+",   Maxi error ="+str(Maxi)+"\n")
            print("The iter =", iter, ",   Maxi error =", Maxi) 

############################4.稳定收敛判断并输出结果##########################            
    Max = Max_error(Gene, Gene0) 
    Gene0 = Gene
    Concentration0 = Concentration
    print(Max)
    if Ttime%(100*dt)==0:
        np.savetxt('Concentration'+str(Ttime)+'.txt',(Concentration))
        np.savetxt('Gene'+str(Ttime)+'.txt',(Gene))
    if Ttime>86400 or Max<1e-2:
        np.savetxt('Concentration.txt',(Concentration))
        np.savetxt('Gene.txt',(Gene))
        break                         


# In[ ]:




