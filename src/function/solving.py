import numpy as np


from src.function.public import por, De, Solid_rate, Water_rate

# 求解浓度方程组========================[mol/L]
def model_sol2(z, Ions, Concentration, Sum, N, dh, Chemicals_present_in = "water", Discretization_methods = 1):   
    dt = 1
    dh = 0.01
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
        print(f"cann't find this substance, {Ions}")
    
    for i in range(1, N):  

        # 1. 假设所有物质分布于沉积物中，
        if Chemicals_present_in == "sedimentary":
            # 离散方式： 空间离散方式为中心二阶差分格式，
            A[i, i] =  -De(z[i]-0.5* dh, Ions)/dh**2 -De(z[i]+0.5* dh, Ions)/dh**2  - 1/dt
            A[i, i - 1] = De(z[i]-0.5* dh, Ions)/dh**2  + Water_rate(z[i]) /(2*dh)
            A[i, i + 1] = De(z[i]+0.5* dh, Ions)/dh**2  - Water_rate(z[i]) /(2*dh)
            rhs[i] = -Sum[i] -Concentration[i]/dt 
        
        # 2. 假设所有物质分布于孔隙水中，
        if Chemicals_present_in == "water":
            # 离散方式： 参考 毛荣 在控制点上离散，同时界面通量参考Li，保证连续并修改为调和平均值，
            phiD_i     = por(z[i])* De(z[i], Ions)
            phiv_i     = por(z[i])* Water_rate(z[i])
            phiv_i_next  = por(z[i+1])* Water_rate(z[i+1])
            phiv_i_last  = por(z[i-1])* Water_rate(z[i-1])
            phiD_i_low =  (dh+dh) / (dh/(De(z[i], Ions)* por(z[i])) + dh/(De(z[i-1], Ions)* por(z[i-1])))  
            phiD_i_up =   (dh+dh) / (dh/(De(z[i], Ions)* por(z[i])) + dh/(De(z[i+1], Ions)* por(z[i+1]))) 
            phiv_i_low = (dh+dh) / (dh/(Water_rate(z[i])* por(z[i])) + dh/(Water_rate(z[i-1])* por(z[i-1])))  
            phiv_i_up =  (dh+dh) / (dh/(Water_rate(z[i])* por(z[i])) + dh/(Water_rate(z[i+1])* por(z[i+1])))
            A[i, i] =  2/(dh*por(z[i]))*(-phiD_i_up/(dh+dh)-phiD_i_low/(dh+dh)) - 1/(dh*por(z[i]))*(phiv_i_up-phiv_i_up**2/phiv_i*dh/(dh+dh)) + 1/(dh*por(z[i]))*(phiv_i_low**2/phiv_i_last*dh/(dh+dh)) - 1/dt
            A[i, i - 1] = 2/(dh*por(z[i]))*phiD_i_low/(dh+dh) + 1/(dh*por(z[i]))*(phiv_i_low-phiv_i_low**2/phiv_i_last*dh/(dh+dh))
            A[i, i + 1] = 2/(dh*por(z[i]))*phiD_i_up/(dh+dh) - 1/(dh*por(z[i]))*(phiv_i_up**2/phiv_i*dh/(dh+dh))
            rhs[i] = -Sum[i] -Concentration[i]/dt
        
        # 修改控制方程和迭代格式, 由于模型最终处于稳定态, 直接计算值dC/dt 达到稳定
        if False:
            phiD_i     = por(z[i])* De(z[i], Ions)
            phiv_i     = por(z[i])* Water_rate(z[i])
            phiv_i_next  = por(z[i+1])* Water_rate(z[i+1])
            phiv_i_last  = por(z[i-1])* Water_rate(z[i-1])
            phiD_i_low =  (dh+dh) / (dh/(De(z[i], Ions)* por(z[i])) + dh/(De(z[i-1], Ions)* por(z[i-1])))  
            phiD_i_up =   (dh+dh) / (dh/(De(z[i], Ions)* por(z[i])) + dh/(De(z[i+1], Ions)* por(z[i+1]))) 
            phiv_i_low = (dh+dh) / (dh/(Water_rate(z[i])* por(z[i])) + dh/(Water_rate(z[i-1])* por(z[i-1])))  
            phiv_i_up =  (dh+dh) / (dh/(Water_rate(z[i])* por(z[i])) + dh/(Water_rate(z[i+1])* por(z[i+1])))
            D_i_up = phiD_i_up*(Concentration[i+1]-Concentration[i])/(dh+dh)
            D_i_low = phiD_i_low*(Concentration[i]-Concentration[i-1])/(dh+dh)
            C_i_up = Concentration[i] + phiv_i_up*(Concentration[i+1]-Concentration[i])/(dh+dh)*dh/(por(z[i])*Water_rate(z[i]))
            C_i_low = Concentration[i-1] + phiv_i_low*(Concentration[i]-Concentration[i-1])/(dh+dh)*dh/(por(z[i-1])*Water_rate(z[i-1]))
            C_model[i] = 2/(por(z[i])*dh)*(D_i_up-D_i_low)-(phiv_i_up*C_i_up-phiv_i_low*C_i_low)/(dh*por(z[i])) + Sum[i]

            # 调和平均值
            #A[i, i] =  -(De(z[i-1], Ions)*De(z[i], Ions))**0.5 * (por(z[i-1])*por(z[i]))**0.5 /dh**2 -(De(z[i+1], Ions)*De(z[i], Ions))**0.5 * (por(z[i+1])*por(z[i]))**0.5/dh**2  
            #A[i, i - 1] = (De(z[i-1], Ions)*De(z[i], Ions))**0.5 * (por(z[i-1])*por(z[i]))**0.5 /dh**2  + Water_rate(z[i-1])*por(z[i-1]) /(2*dh)
            #A[i, i + 1] = (De(z[i+1], Ions)*De(z[i], Ions))**0.5 * (por(z[i+1])*por(z[i]))**0.5 /dh**2 - Water_rate(z[i+1])*por(z[i+1]) /(2*dh)
            #rhs[i] = Sum[i]*por(z[i])
            #C_model[i] = (A[i, i - 1]*Concentration[i-1] + A[i, i]*Concentration[i] + A[i, i+1]*Concentration[i+1] + rhs[i]) / por(z[i])
         
    # y=solve(A,b)===============
    C_model = np.dot(np.linalg.inv(A), rhs)
    return C_model

#求解基因方程组========================[copies/L]
def model_sol3(Rct, gene, geneR, z, N, Gene_present_in = "sedimentary", Discretization_methods = 1):
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

        # 1. 假设所有基因分布于沉积物中，
        if Gene_present_in == "sedimentary":
            # 离散方式1： 空间离散方式为中心二阶差分格式，
            A[i, i] = - lamda * (1 - por(z[i]))
            A[i, i + 1] = - Solid_rate(z[i+1]) * (1 - por(z[i+1])) / (2*dh) 
            A[i, i - 1] = Solid_rate(z[i-1]) * (1 - por(z[i-1]))  / (2*dh) 
            rhs[i]  = - geneR[i] * (1 - por(z[i]))  - gene[i] / dt * (1 - por(z[i]))


        # 2. 假设所有基因分布于孔隙水中，
        # ??? 基因的水动力弥散系数
        # 离散方式： 参考 毛荣 在控制点上离散，同时界面通量参考Li，保证连续并修改为调和平均值，
        if Gene_present_in == "water":
            phiD_i     = por(z[i])* De(z[i], Ions)
            phiv_i     = por(z[i])* Water_rate(z[i])
            phiv_i_next  = por(z[i+1])* Water_rate(z[i+1])
            phiv_i_last  = por(z[i-1])* Water_rate(z[i-1])
            phiD_i_low =  (dh+dh) / (dh/(De(z[i], Ions)* por(z[i])) + dh/(De(z[i-1], Ions)* por(z[i-1])))  
            phiD_i_up =   (dh+dh) / (dh/(De(z[i], Ions)* por(z[i])) + dh/(De(z[i+1], Ions)* por(z[i+1]))) 
            phiv_i_low = (dh+dh) / (dh/(Water_rate(z[i])* por(z[i])) + dh/(Water_rate(z[i-1])* por(z[i-1])))  
            phiv_i_up =  (dh+dh) / (dh/(Water_rate(z[i])* por(z[i])) + dh/(Water_rate(z[i+1])* por(z[i+1])))       
            A[i, i] =  2/(dh*por(z[i]))*(-phiD_i_up/(dh+dh)-phiD_i_low/(dh+dh)) - 1/(dh*por(z[i]))*(phiv_i_up-phiv_i_up**2/phiv_i*dh/(dh+dh)) + 1/(dh*por(z[i]))*(phiv_i_low**2/phiv_i_last*dh/(dh+dh)) - 1/dt
            A[i, i - 1] = 2/(dh*por(z[i]))*phiD_i_low/(dh+dh) + 1/(dh*por(z[i]))*(phiv_i_low-phiv_i_low**2/phiv_i_last*dh/(dh+dh))
            A[i, i + 1] = 2/(dh*por(z[i]))*phiD_i_up/(dh+dh) - 1/(dh*por(z[i]))*(phiv_i_up**2/phiv_i*dh/(dh+dh))
            rhs[i] = -(geneR[i]-lamda*gene[i]) -gene[i]/dt 

            # 修改控制方程和迭代格式, 由于模型最终处于稳定态, 直接计算值dTao/dt 达到稳定
            if False:
                phiD_i     = por(z[i])* De(z[i], Ions)
                phiv_i     = por(z[i])* Water_rate(z[i])
                phiv_i_next  = por(z[i+1])* Water_rate(z[i+1])
                phiv_i_last  = por(z[i-1])* Water_rate(z[i-1])
                phiD_i_low =  (dh+dh) / (dh/(De(z[i], Ions)* por(z[i])) + dh/(De(z[i-1], Ions)* por(z[i-1])))  
                phiD_i_up =   (dh+dh) / (dh/(De(z[i], Ions)* por(z[i])) + dh/(De(z[i+1], Ions)* por(z[i+1]))) 
                phiv_i_low = (dh+dh) / (dh/(Water_rate(z[i])* por(z[i])) + dh/(Water_rate(z[i-1])* por(z[i-1])))  
                phiv_i_up =  (dh+dh) / (dh/(Water_rate(z[i])* por(z[i])) + dh/(Water_rate(z[i+1])* por(z[i+1])))  

                D_i_up = phiD_i_up*(gene[i+1]-gene[i])/(dh+dh)
                D_i_low = phiD_i_low*(gene[i]-gene[i-1])/(dh+dh)
                gene_i_up = gene[i] + phiv_i_up*(gene[i+1]-gene[i])/(dh+dh)*dh/(por(z[i])*Water_rate(z[i]))
                gene_i_low = gene[i-1] + phiv_i_low*(gene[i]-gene[i-1])/(dh+dh)*dh/(por(z[i-1])*Water_rate(z[i-1]))
                R_model[i] = 2/(por(z[i])*dh)*(D_i_up-D_i_low)-(phiv_i_up*gene_i_up-phiv_i_low*gene_i_low)/(dh*por(z[i])) + geneR[i]-lamda * gene[i]


    # y=solve(A,b)===============
    R_model = np.dot(np.linalg.inv(A), rhs)    
    return R_model