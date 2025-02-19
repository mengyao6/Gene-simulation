import numpy as np
import xlrd, xlwt
import pandas as pd
import math
import os
import matplotlib.pyplot as plt
import pdb
from datetime import date

from src.function.public import Rate_R2, Sum_R, Max_error
from src.function.solving import model_sol2, model_sol3

from src.model.chemicals import Chemicals
from src.model.gene import Genes
from src.model.DeltaG import DeltaGs


#############################0. Open file, read data about gene#########################
#f=open("Gene.txt","w")

# ��������·��
save_dir = str(date.today())
os.makedirs(save_dir, exist_ok=True)


#############################1. Meshing  #########################
dz = 0.01
z = np.arange(2.25, 27.07, dz)    
N = 2481  
dt = 1
############################2.Initialisation of state variables (the concentration dimenson  mol/L; the gene dimension genes/L) ###########    
Concentration = np.zeros((N+1, 13), dtype=np.float64)
Concentration0 = np.zeros((N+1, 13), dtype=np.float64)
Gene = np.zeros((N+1, 10), dtype=np.float64)
Gene0 = np.zeros((N+1, 10), dtype=np.float64)
C_Name=('0: C_co2', '1: C_hco3', '2: C_c6', '3: C_H', '4: C_no2', '5: C_no3', '6: C_n2', '7: C_nh4', '8: C_h2s', '9: C_so4', '10:C_DIC',  '11:C_o2', '12:C_POC')
# The dta from Minzong zone [wt.%��mol/L 0: C_co2, 1: C_hco3, 2: C_c6, 3: C_H, 4: C_no2, 5: C_no3, 
# 6: C_n2, 7: C_nh4, 8: C_h2s, 9: C_so4, 10:C_DIC,  11:C_o2,    12:C_POC]
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
    #Gene[genes/L 0: cox,    1: narG,   2: nirk,   3: nrf,    4: dsr,     5: amoA,     6: hzo,    7: nap        8: nor      9: sox]
    Gene[i][0] = 1e5
    Gene[i][1] = 1e5
    Gene[i][2] = 1e5
    Gene[i][3] = 1e5
    Gene[i][4] = 1e5
    Gene[i][5] = 1.8149115050146092e-05 *z[i]**10+( -0.0025970258840197304 *z[i]**9)+( 0.12524748342028677 *z[i]**8)+( -0.9848821749028847 *z[i]**7)+( -123.61588309591274 *z[i]**6)+( 4925.4942433373935 *z[i]**5)+( -75163.8547815459 *z[i]**4)+( 321348.1047953144 *z[i]**3)+( 4250441.997581997 *z[i]**2)+( -53197485.408834904 *z[i])+( 169686702.3348549 )
    Gene[i][6] = -3.1496922177991296e-06 *z[i]**10+( 0.00031832652419954723 *z[i]**9)+( -0.004374674229161838 *z[i]**8)+( -0.4849318876896918 *z[i]**7)+( 11.03207967439595 *z[i]**6)+( 729.1492911791356 *z[i]**5)+( -40833.72222230589 *z[i]**4)+( 864073.0902662921 *z[i]**3)+( -9246811.681851482 *z[i]**2)+( 49431498.17036085 *z[i])+( -99993069.60356417 )
    Gene[i][7] = 1e5  # TODO ������ why set 0 before
    Gene[i][8] = 1e5
    Gene[i][9] = 1e5
Concentration = np.divide(Concentration, 1e6)                        #concentration [mol/L]
Gene = Gene
Concentration0 = Concentration
Gene0 = Gene


Ttime = 0
while 1:
    Ttime = Ttime + dt
    #f.write("��ǰ����ʱ��"+str(Ttime)+"\n")
    print("new time: ", Ttime)
############################3. picard�ڲ����##########################     
    iter = 0
    while 1:
        iter = iter + 1
    ######################3.1 ���ñ�������Ũ�� Ϊ��ֵ#################
        # print("Concentration �е���Сֵ:", np.min(Concentration))
        Concentration[Concentration < 0] = 1e-20

        # print("Gene �е���Сֵ:", np.min(Gene))
        Gene[Gene < 0] = 1e-20

    ###########################3.2 ��������######################
        GeneR = np.zeros((10, N+1), dtype=np.float64)
        for i in range(0, N+1):   #���㲻ͬ��λ��Ri�������������ʣ�[Coies/L/s]
            #print(GeneR[9,:])
            GeneR[:, i] = Rate_R2(Chemicals(Concentration[i,:]), Genes(Gene[i,:]), z[i])[:, 0]                    
        
        Tao = np.zeros((N+1, 10), dtype=np.float64)
        for i in range(0, 10):
            Tao[:, i] = model_sol3(i,  Gene[:,i],  GeneR[i,:], z,  N, "sedimentary")[:, 0]                     

    ######################3.3 Ũ�����##########################            
        SumR = np.zeros((7, N+1), dtype=np.float64)
        for i in range(0, N+1):  #����ϸ����л�����Ĳ�ͬ��ѧ���� Ri   #[mol/L/s]
            SumR[:, i] = Sum_R(Chemicals(Concentration[i,:]), Genes(Gene[i,:]), z[i])[:, 0]                           

        Sum1=0;  Sum2=0;  Sum3=0;  Sum4=0;  Sum5=0;  Sum6=0;  Sum7=0
        Sum1= model_sol2(z, 0, Concentration[:,2],   SumR[0,:], N, dz)                           #C6
        Sum2= model_sol2(z, 1, Concentration[:,11],  SumR[1,:], N, dz)                           #O2
        Sum3= model_sol2(z, 2, Concentration[:,7],   SumR[2,:], N, dz)                           #NH4
        Sum4= model_sol2(z, 3, Concentration[:,4],   SumR[3,:], N, dz)                           #NO2
        Sum5= model_sol2(z, 4, Concentration[:,5],   SumR[4,:], N, dz)                           #NO3
        Sum6= model_sol2(z, 5, Concentration[:,9],   SumR[5,:], N, dz)                           #SO4
        Sum7= model_sol2(z, 6, Concentration[:,8],   SumR[6,:], N, dz)                          #H2S        

    ###########################3.3����Ũ�ȸ���######################
        for i in range(0, N+1):         
            Concentration[i][2]=Sum1[i]                                    
            Concentration[i][11]=Sum2[i]
            Concentration[i][7]=Sum3[i]
            Concentration[i][4]=Sum4[i]
            Concentration[i][5]=Sum5[i]
            Concentration[i][9]=Sum6[i]        
           # Concentration[i][8]=Sum7[i]   


    ###########################3.4 ���������жϲ����»�����#########################
        Maxi = Max_error(Gene, Tao)
        Gene = Tao

        if Maxi<1e-5  or  iter>1000:
            #f.write("The iter ="+str(iter)+",   Maxi error ="+str(Maxi)+"\n")
            print("The iter =", iter,",   Maxi error =", Maxi)   
            break
        else: 
            #f.write("The iter ="+str(iter)+",   Maxi error ="+str(Maxi)+"\n")
            print("The iter =", iter, ",   Maxi error =", Maxi) 

############################4.�ȶ������жϲ�������##########################            
    Max = Max_error(Gene, Gene0) 
    Gene0 = Gene
    Concentration0 = Concentration
    print(Max)
    if Ttime%(100*dt)==0:
        np.savetxt(os.path.join(save_dir, f'Concentration-{Ttime}.txt'), (Concentration))
        np.savetxt(os.path.join(save_dir, f'Gene-{Ttime}.txt'), (Gene))
    if Ttime>86400 or Max<1e-2:
        np.savetxt(os.path.join(save_dir, f'Concentration.txt'), (Concentration))
        np.savetxt(os.path.join(save_dir, f'Gene.txt'), (Gene))        
        break 