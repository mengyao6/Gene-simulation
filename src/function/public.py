import numpy as np

from src.model.chemicals import Chemicals
from src.model.gene import Genes, Mus
from src.model.DeltaG import DeltaGs

Fara = 96.485 # [KJ/V] Faraday constant
dPsi = 0.12   # [V] electric potential

R = 8.31 * 1e-3      # [KJ/k/mole] the gas constant
T = 298.15    # [-] the thermodynamic temperature

#Gene 0：cox,      1：narG,   2：nirk,   3：nrf,    4：dsr,     5：amoA,     6：hzo,    7：nap        8：nor      9：sox    
gene_names = ['cox', 'narG', 'nirk', 'nrf', 'dsr', 'amoA', 'hzo', 'nap', 'nor', 'sox']
# 0: C_co2, 1: C_hco3, 2: C_c6, 3: C_H, 4: C_no2, 5: C_no3, 6: C_n2, 7: C_nh4, 8: C_h2s, 9: C_so4, 10:C_DIC,  11:C_o2 ,    12:C_POC
Chemical_names = ['C_co2', 'C_hco3', 'C_c6', 'C_H', 'C_no2', 'C_no3', 'C_n2', 'C_nh4', 'C_h2s', 'C_so4', 'C_DIC', 'C_o2', 'C_POC']

simualtion = ['C_c6', 'C_o2', 'C_nh4', 'C_no2', 'C_no3', 'C_so4', 'C_h2s']

# define porosity dist.==========
def por(depth):                   #BJ
    beta = 0.008  # cm-1
    por_sta = 0.55     #porosity at start
    por_end = 0.3    #porosity at infinity
    # porosity= 0.30 + 0.25 * 2.71828**(-depth/100)
    # porosity= 0.29 + 0.36 * 2.71828**(-depth/7)
    porosity = por_sta + (por_end - por_sta) * np.exp((-depth * beta))
    return porosity


# Dsw f(Temp, Sality, Presure) 
def diff(Temp = 25, Sality = 35, Presure = 1.013253):
    sinyr = 60*60*24*365 # number of seconds in one year
    Dmol = Chemicals()
    Dsw = Chemicals()

    Dmol.C_co2 = 292.2801   # [m2/s]
    Dmol.C_hco3 = 152.2161
    Dmol.C_o2 = 371.2064
    Dmol.C_Mn = 95.66153
    Dmol.C_no2 = 309.8471
    Dmol.C_nh4 = 285.7813
    Dmol.C_so4 = 146.8013
    Dmol.C_h2s = 255.6641
    Dmol.C_po4 = 78.81548
    Dmol.C_c6 = 300.00  # TODO
    Dmol.C_no3 = 300.00  # TODO

    for matter in Chemical_names:
        matter_value = getattr(Dmol, matter)
        # Apply the transformation and store it in F_ts
        setattr(Dsw, matter, matter_value*sinyr*10**4)

    return Dsw


# define diffusion coefficient==============[m2/yr]
def De(depth, Ions): 
    
    # Ref Meiqing
    # Ds=[[61.92,96.24],[183.12,284.4],[176.64,274.32],[175.92,272.88],[99.12,153.6],[194.4,301.68]] #cm2/yr
    # if 17.6<depth<=20.7 or 21.2<depth<=24.3:
    #     eff_D=Ds[Ions-1][0]        
    # else:
    #     eff_D=Ds[Ions-1][1] 
    # Dsw=np.array([6.7,22.9,19.8,19.1,19.0,10.7,21.0,500])*1e-6 * 60 * 60 * 24 * 365 #[cm2/year]   
   
   
   
    # #Ref the R script from Zhaorui 
    Dsw = diff()
    eff_D=  getattr(Dsw, simualtion[Ions]) / ( 1- 2 * np.log(por(depth)))
    
    return eff_D


# define 沉积物埋藏速率==========[m/s]
def Solid_rate(depth):                   #BJ 
    Vs_start = 1.96 * 10**(-3)   
    Vs_end = 1.96 * 10**(-3)   
    Vs = Vs_end * (1-por(1e100)) / (1-por(depth))  # Chuang P C et al.(2019); Rui Zhao et al. (2016)
    return Vs * 1e-2/86400/365


# define 孔隙水垂向流速==========[m/s]
def Water_rate(depth):                   #BJ
    Vp_start = 1.72   
    Vp_end = 1.72
    Vs_start = 1.96 * 10**(-3)   
    Vs_end = 1.96 * 10**(-3)    
    
    Vp = Vp_end * por(1e100) / por(depth)                # Rui Zhao et al. (2016)

    # Vp = (por(1e10)*Vs_00 - por(0)*Vp_0) /por(depth)  # Chuang P C et al.2019
    return Vp * 1e-2/86400/365


# define Gibbs free energy =======[KJ/mol]
def Gibbs_G( Concentration: Chemicals):
    # D_G = np.zeros((10, 1), dtype=np.float64)
    # 1：C_co2, 2：C_hco3, 3：C_c6, 4：C_H, 5：C_no2, 6：C_no3, 7：C_n2, 8：C_nh4, 9：C_h2s, 10：C_so4, 11:C_DIC,  12:C_o2
    C_co2 = Concentration.C_co2; C_hco3 = Concentration.C_hco3; C_c6 = Concentration.C_c6; C_H  = Concentration.C_H; C_no2 = Concentration.C_no2
    C_no3 = Concentration.C_no3;  C_n2  = Concentration.C_n2; C_nh4= Concentration.C_nh4; C_h2s= Concentration.C_h2s; C_so4 = Concentration.C_so4
    C_DIC = Concentration.C_DIC; C_o2  = Concentration.C_o2
    D_G = Genes()
    D_G.cox = -479.0 + R * T * np.log( C_co2 / C_c6**(1/6) / C_o2 )                                    #DeltaG_cox 
    D_G.narG = -321.6 + R * T * np.log( C_co2 * C_no2**2    / C_c6**(1/6) / C_no3**2)                  #=DeltaG_narG
    D_G.nirk = -594.1 + R * T * np.log( C_co2 * C_n2**(2/3) / C_c6**(1/6) / C_no2**(4/3) / C_H**(4/3))#=DeltaG_nirk
    D_G.nrf = -352.4 + R * T * np.log( C_co2 * C_nh4**(2/3)/ C_c6**(1/6) / C_no2**(2/3) / C_H**(4/3))#=DeltaG_nrf
    D_G.dsr = -76.1  + R * T * np.log( C_hco3* C_h2s**(1/2)/ C_c6**(1/6) / C_so4**(1/2))             #=DeltaG_dsr
    D_G.amoA = -189.9 + R * T * np.log( C_no2 * C_H**2      / C_nh4       / C_o2**(3/2))              #=DeltaG_amoA
    D_G.hzo = -362.7 + R * T * np.log( C_n2  / C_nh4       / C_no2)                                  #=DeltaG_hzo
    D_G.nap = -201.5 + R * T * np.log( C_so4**(1/4) * C_no2  * C_H**(1/2)  / C_h2s**(1/4) / C_no3)    #=DeltaG_nap
    D_G.nor = -157.4 + R * T * np.log( C_no3**2 / C_no2**2       / C_o2)                           #=DeltaG_nor
    D_G.sox = -401.8 + R * T * np.log( C_so4 * C_co2**2   / C_o2**2 / C_h2s / C_hco3**2)           #=DeltaG_sox
    return D_G


# define A function of free energy yield grams of biomass Y.==========[g/mol/donor]
def Biomass_Y(D_G: Genes):
    Y_G = np.zeros((10, 1), dtype=np.float64)
    # 1：cox,      2：narG,   3：nirk,   4：nrf,    5：dsr,     6：amoA,     7：hzo,    8：nap        9：nor      10：sox
    Y_G[0] = 2.08 - 0.0211 * D_G.cox * 6
    Y_G[1] = 2.08 - 0.0211 * D_G.narG * 6
    Y_G[2] = 2.08 - 0.0211 * D_G.nirk * 6
    Y_G[3] = 2.08 - 0.0211 * D_G.nrf * 6
    Y_G[4] = 2.08 - 0.0211 * D_G.dsr * 6
    Y_G[5] = 2.08 - 0.0211 * D_G.amoA
    Y_G[6] = 2.08 - 0.0211 * D_G.hzo
    Y_G[7] = 2.08 - 0.0211 * D_G.nap
    Y_G[8] = 2.08 - 0.0211 * D_G.nor
    Y_G[9] = 2.08 - 0.0211 * D_G.sox
    return Y_G

# define thermodynamic potential factor F_T==========
def Thermo_T(D_G: Genes):
    F_ts = DeltaGs()
    # Loop through each gene and compute its transformed value
    for gene in gene_names:
        # Use getattr to access the gene value dynamically
        gene_value = getattr(D_G, gene)
        # Apply the transformation and store it in F_ts
        setattr(F_ts, gene, 1.0 / (np.exp((gene_value + Fara * dPsi) / (R * T)) + 1))

    return F_ts



# define the rate of gene production Ri==========[Copies/g/s]
def Rate_R(taf: Genes, F_t: Genes, Concentration: Chemicals):
    G_R = np.zeros((10, 1), dtype=np.float64)
    # Concentration 1：C_co2, 2：C_hco3, 3：C_c6, 4：C_H, 5：C_no2, 6：C_no3, 7：C_n2, 8：C_nh4, 9：C_h2s, 10：C_so4, 11:C_DIC,  12:C_o2
    C_co2 = Concentration.C_co2; C_hco3 = Concentration.C_hco3; C_c6 = Concentration.C_c6; C_H  = Concentration.C_H; C_no2 = Concentration.C_no2
    C_no3 = Concentration.C_no3;  C_n2  = Concentration.C_n2; C_nh4 = Concentration.C_nh4; C_h2s= Concentration.C_h2s; C_so4 = Concentration.C_so4
    C_DIC = Concentration.C_DIC; C_o2  = Concentration.C_o2; C_POC= Concentration.C_POC

    # 1：cox,      2：narG,   3：nirk,   4：nrf,    5：dsr,     6：amoA,     7：hzo,    8：nap        9：nor      10：sox
    # mu=np.array([0.28, 0.151, 0.247, 0.162, 0.0636, 0.432, 0.864, 0.864, 0.432, 0.864])/86400         #s-1
    mu = Mus()
    kk=np.asarray([[0.7e-6,0.121e-6], [0.7e-6,0.3e-6], [0.7e-6,0.3e-6], [0.7e-6,0.3e-6], [0.7e-6,3e-6,15e-6], [107e-6,18.75e-6], [5e-6,5e-6,0.2e-6], [0.121e-6,0.121e-6], [64.3e-6,16.9e-6], [0.121e-6,0.121e-6]], dtype=object)
    G_R[0] = taf.cox * F_t.cox * mu.cox * (C_o2 / ( C_o2 + kk[0][1] )) * ( C_c6 /(C_c6 + kk[0][0]))                         #1.cox           
    G_R[1] = taf.narG * F_t.narG * mu.narG * (C_no3 / ( C_no3 + kk[1][1] )) * ( C_c6 /(C_c6 + kk[1][0]))                       #2.narG            
    G_R[2] = taf.nirk * F_t.nirk * mu.nirk * (C_no2 / ( C_no2 + kk[2][1] )) * ( C_c6 /(C_c6 + kk[2][0]))                       #3.nirk
    G_R[3] = taf.nrf * F_t.nrf * mu.nrf * (C_no2 / ( C_no2 + kk[3][1] )) * ( C_c6 /(C_c6 + kk[3][0]))                       #4.nrf
    G_R[4] = taf.dsr * F_t.dsr * mu.dsr * (C_so4 / ( C_so4 + kk[4][1] )) * ( C_c6 /(C_c6 + kk[4][0]))  * (C_o2 /( C_o2+ kk[4][2]))   #5.dsr
    G_R[5] = taf.amoA * F_t.amoA * mu.amoA * (C_nh4 / ( C_nh4 + kk[5][0] )) * ( C_o2 /(C_o2 + kk[5][1]))                                #6.amoA
    G_R[6] = taf.hzo * F_t.hzo * mu.hzo * (C_nh4 / ( C_nh4 + kk[6][0] )) * ( C_no2 /(C_no2 + kk[6][1])) * ( C_o2 /(C_o2 + kk[6][2]))  #7.hzo
    G_R[7] = taf.nap * F_t.nap * mu.nap * (C_no3 / ( C_no3 + kk[7][0] )) * ( C_h2s /(C_h2s + kk[7][1]))                               #8.nap 
    G_R[8] = taf.nor * F_t.nor * mu.nor * (C_no2 / ( C_no2 + kk[8][0] )) * ( C_o2 / ( C_o2 + kk[8][1] ))                              #9.nor
    G_R[9] = taf.sox * F_t.sox * mu.sox * (C_h2s / ( C_h2s + kk[9][0]))  * ( C_o2 / ( C_o2 + kk[9][1] ))                              #10.sox
    return G_R

# define the Organic Matter Remineralization ==========[mol/g/s] 
def DOC_Rc(C_DIC, C_c6, TOC_age):
    Kc = 40 / 1e3   # mol/L
    Rc = Kc / (C_DIC + Kc) * 0.16 * TOC_age**(-0.95) * C_c6        
    return Rc/365/86400

# define sediment age.==========[Year]
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

# define the rate of DOC production by microbial mortality ==========[mol/g/s   &&  gene/g/s]
def DOC_Rm(taf: Genes):
    lamda = 0.001 /86400           #s-1
    Ni  = 3.75 * 1e13    #[genes/g]
    #Gene[genes/L 0: cox,    1: narG,   2: nirk,   3: nrf,    4: dsr,     5: amoA,     6: hzo,    7: nap        8: nor      9: sox]
    Rm = 0.4276  * lamda / 12 * (taf.cox/Ni + taf.narG/Ni + taf.nirk/Ni + taf.nrf/Ni + taf.dsr/Ni + taf.amoA/Ni + taf.hzo/Ni/8 + taf.nap/Ni + taf.nor/Ni + taf.sox/Ni)
    return Rm 

# define Nitrogen as a nutrient. ==========[mol/g/s]
def Nit_U(taf: Genes, GeneR, C_nh4, C_no2):
    Un = np.zeros((3, 1), dtype=np.float64)
    K_nutrient = 134/1e9  #mol
    Ni = 3.75 * 1e13    #[genes/g]
    # mu=np.array([0.28, 0.151, 0.247, 0.162, 0.0636, 0.432, 0.864, 0.864, 0.432, 0.864])/86400         #s-1
    mu = Mus()
    un = 0.1127  / 13 * (taf.cox * mu.cox /Ni + taf.narG * mu.narG /Ni + taf.nirk * mu.nirk /Ni + taf.nrf * mu.nrf /Ni + taf.dsr * mu.dsr /Ni+ taf.amoA * mu.amoA /Ni + taf.hzo * mu.hzo /Ni/8 + taf.nap * mu.nap /Ni + taf.nor * mu.nor /Ni + taf.sox * mu.sox /Ni)
    # un = 0.00000026
    Un[0] = un * (C_nh4 / (C_nh4 + K_nutrient))
    Un[1] = (un - Un[0]) * (C_no2 / (C_no2 + K_nutrient))
    Un[2] = un - Un[0] - Un[1]
    return Un

# 定义生物代谢产生的化学物质. ==========[mol/L/s] 
def R_nyg(GR, YG):
    Rny = DeltaGs()
    Ni = 3.75 * 1e13    #[genes/g]
    for i , gene in enumerate(gene_names):
        setattr(Rny, gene, GR[i][0] / Ni / YG[i][0])
    return Rny

# 计算右端项生物代谢贡献的化学物质总和. ==========[mol/L/s]
def Sum_R(Concentration: Chemicals, taf: Genes, depth: float):
    SumR = np.zeros((7, 1), dtype=np.float64)
    Ni  = 3.75 * 1e13    #[genes/g]
    # [mol/L】：C_co2, 2：C_hco3, 3：C_c6, 4：C_H, 5：C_no2, 6：C_no3, 7：C_n2, 8：C_nh4, 9：C_h2s, 10：C_so4, 11:C_DIC,  12:C_o2
    # C_co2 = Concentration[0]; C_hco3 = Concentration[1]; C_c6 = Concentration[2]; C_H  = Concentration[3]; C_no2 = Concentration[4]
    # C_no3 = Concentration[5];  C_n2  = Concentration[6]; C_nh4= Concentration[7]; C_h2s= Concentration[8]; C_so4 = Concentration[9]
    # C_DIC = Concentration[10]; C_o2  = Concentration[11]; C_POC  = Concentration[12]
    # Con   = [C_co2, C_hco3, C_c6, C_H, C_no2, C_no3, C_n2, C_nh4, C_h2s, C_so4, C_DIC, C_o2]   
    DetG = Gibbs_G( Concentration )
    YG = Biomass_Y( DetG )
    Ft = Thermo_T(DetG)
    GR = Rate_R(taf, Ft, Concentration )
    Rc = DOC_Rc(Concentration.C_DIC, Concentration.C_POC, age(depth))
    Rn = 16/106 * Rc
    Rm  = DOC_Rm(taf)
    Un  = Nit_U(taf, GR, Concentration.C_nh4, Concentration.C_no2)
    Rnyg= R_nyg(GR, YG)
    #Gene 0：cox,      1：narG,   2：nirk,   3：nrf,    4：dsr,     5：amoA,     6：hzo,    7：nap        8：nor      9：sox    
    SumR[0] =  Rc/6 + Rm/6 - Rnyg.cox - Rnyg.narG - Rnyg.nirk - Rnyg.nrf - Rnyg.dsr         
    #SumR[0] = -Rc      ##########################20200117沉积物
    SumR[1] = -6*Rnyg.cox - 1.5*Rnyg.amoA - 0.5*Rnyg.nor - 2*Rnyg.sox
    SumR[2] = Rn + 4*Rnyg.nrf - Rnyg.amoA - Rnyg.hzo - Un[0]
    SumR[3] = 12*Rnyg.narG + Rnyg.amoA + 4*Rnyg.nap - Rnyg.hzo - Rnyg.nor- 8*Rnyg.nirk - 4*Rnyg.nrf - Un[1]
    SumR[4] = Rnyg.nor - 12*Rnyg.narG - 4*Rnyg.nap - Un[2]
    SumR[5] = Rnyg.nap + Rnyg.sox - 3*Rnyg.dsr
    SumR[6] = -Rnyg.nap - Rnyg.sox + 3*Rnyg.dsr
    return SumR 

# 基因产生速率 ===========[cpoies/L/s]
def Rate_R2(Concentration: Chemicals, taf: Genes, depth: float):    
    GR = np.zeros((10, 1), dtype=np.float64)
    Ni  = 3.75 * 10**13    #[genes/g]
    # 1：C_co2, 2：C_hco3, 3：C_c6, 4：C_H, 5：C_no2, 6：C_no3, 7：C_n2, 8：C_nh4, 9：C_h2s, 10：C_so4, 11:C_DIC,  12:C_o2
    # C_co2 = Concentration[0]; C_hco3 = Concentration[1]; C_c6 = Concentration[2]; C_H  = Concentration[3]; C_no2 = Concentration[4]
    # C_no3 = Concentration[5];  C_n2  = Concentration[6]; C_nh4= Concentration[7]; C_h2s= Concentration[8]; C_so4 = Concentration[9]
    # C_DIC = Concentration[10]; C_o2  = Concentration[11]
    # Con   = [C_co2, C_hco3, C_c6, C_H, C_no2, C_no3, C_n2, C_nh4, C_h2s, C_so4, C_DIC, C_o2]   
    DetG = Gibbs_G( Concentration )
    Ft = Thermo_T( DetG )
    GR = Rate_R(taf, Ft, Concentration)
    return GR

#求解基因最大误差========================[copies/L]
def Max_error(Gene,Tao):
    det_g = np.max(np.abs(Gene - Tao), axis=0)
    return np.max(det_g)