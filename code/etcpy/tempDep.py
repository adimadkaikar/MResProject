import numpy as np
from scipy.optimize import fsolve
import pandas as pd
import time


def sample_data_uncertainty_with_constraint_increasing_topt(inpt,columns=None):
    if type(inpt)==tuple:
        params,seed = inpt
        np.random.seed(seed+int(time.time()))
    else: params = inpt
    '''
    # params is a dataframe with following columns:
    # Tm,Tm_std:  melting temperature. Given in K
    #
    # T90, temperature at which 90% of an enzyme is denatured. This is not mandentory. If missing, protein length will be used
    #      to calculate denaturation curve. Given in K
    # 
    # dCpt,dcpt_std: the heat capacity difference between transition state and ground state in enzyme catalyzed reaction. Given
    #                in J/mol/K
    # 
    # Topt,Topt_std: the optimal temprature at which the specific activity is maximized. Given in K
    #
    # Length, protein length. This will not be used in this function, but it will be used later in the calculation of thermal
    #         parameters.
    #
    # xx_std, corresponding uncertainty given by standard deviation.
    # 
    # columns: a list of columns to be sampled, could be any combination of ['Tm','dCpt','Topt]. 
    #          If it is None, then sample all three columns
    # 
    # The script will return an new dataframe with the same columns but with randomly sampled data
    '''
    
    sampled_params = params.copy()
    if columns is None: columns = ['Topt', 'dCpt', 'Tm']
    
    for col in columns:
        if col == 'Topt':
            sampled_params.loc[sampled_params['topt_source']=='BullShit', 'Topt'] = sampled_params.loc[sampled_params['topt_source']=='BullShit', 'Topt'] +(sampled_params.loc[sampled_params['topt_source']=='BullShit', 'Tm'] - sampled_params.loc[sampled_params['topt_source']=='BullShit', 'Topt'])*0.3
            
            sampled_params.loc[sampled_params['topt_source']!='BullShit', 'Topt'] = sampled_params.loc[sampled_params['topt_source']!='BullShit', 'Topt'] +(sampled_params.loc[sampled_params['topt_source']!='BullShit', 'Tm'] - sampled_params.loc[sampled_params['topt_source']!='BullShit', 'Topt'])*0.2
            
        if col == 'dCpt':
            lst = [np.random.normal(-4000,params.loc[ind,col+'_std']) for ind in sampled_params.index]
            sampled_params[col] = lst
        if col == 'Tm':
            for ind in sampled_params.index:
                if params.loc[ind, 'TmTag']=='Mean':
                    Tm = np.random.normal(params.loc[ind,col],params.loc[ind, 'Tm_std'])
                else:
                    if params.loc[ind, col] < 320:
                        Tm = params.loc[ind, col]
                    else:
                        Tm = np.random.normal(params.loc[ind,col],2)
                if Tm <= params.loc[ind, 'Topt']:
                    count = 0
                    while sampled_params.loc[ind, 'Tm'] and count <10:
                        Tm = np.random.normal(params.loc[ind,col],1)
                        count = count + 1
                sampled_params.loc[ind, 'Tm'] = Tm
    return sampled_params


def sample_data_uncertainty_with_constraint_small_increasing_topt(inpt,columns=None):
    if type(inpt)==tuple:
        params,seed = inpt
        np.random.seed(seed+int(time.time()))
    else: params = inpt
    '''
    # params is a dataframe with following columns:
    # Tm,Tm_std:  melting temperature. Given in K
    #
    # T90, temperature at which 90% of an enzyme is denatured. This is not mandentory. If missing, protein length will be used
    #      to calculate denaturation curve. Given in K
    # 
    # dCpt,dcpt_std: the heat capacity difference between transition state and ground state in enzyme catalyzed reaction. Given
    #                in J/mol/K
    # 
    # Topt,Topt_std: the optimal temprature at which the specific activity is maximized. Given in K
    #
    # Length, protein length. This will not be used in this function, but it will be used later in the calculation of thermal
    #         parameters.
    #
    # xx_std, corresponding uncertainty given by standard deviation.
    # 
    # columns: a list of columns to be sampled, could be any combination of ['Tm','dCpt','Topt]. 
    #          If it is None, then sample all three columns
    # 
    # The script will return an new dataframe with the same columns but with randomly sampled data
    '''
    
    sampled_params = params.copy()
    if columns is None: columns = ['Topt', 'dCpt', 'Tm']
    
    for col in columns:
        if col == 'Topt':
            for ind in sampled_params.index:
                if params.loc[ind, 'topt_source']=='BullShit':
                    topt = params.loc[ind, col] + (params.loc[ind, 'Tm']-params.loc[ind, col])*0.08
                else:
                    topt = params.loc[ind, col] + (params.loc[ind, 'Tm']-params.loc[ind, col])*0.05
                Tm = params.loc[ind, 'Tm']
                count = 0
                if topt>= Tm:
                    topt = topt - (topt - Tm)*2
                sampled_params.loc[ind, col] = topt
            
        if col == 'dCpt':
            lst = [np.random.normal(-4000,params.loc[ind,col+'_std']) for ind in sampled_params.index]
            sampled_params[col] = lst
        if col == 'Tm':
            for ind in sampled_params.index:
                if params.loc[ind, 'TmTag']=='Mean':
                    sampled_params.loc[ind, 'Tm'] = np.random.normal(params.loc[ind,col],1)
                else:
                    sampled_params.loc[ind, 'Tm'] = np.random.normal(params.loc[ind,col],0.5)
    return sampled_params


def simulate_growth_with_resampling(model,dfae,sigma,oriParams,Tadj=0):
    '''
    # model, cobra model
    # dfae, a dataframe of experimental data
    # sigma, enzyme saturation factor
    # params, a dataframe containing thermal parameters of enzymes: dHTH, dSTS, dCpu, Topt
    # Ensure that Topt is in K. Other parameters are in standard units.
    # Tadj, as descrbed in map_fNT
    #
    '''
    rs = list()
    rg = list()
    ro = list()
    Ts = dfae.index+273.15
    
    for T in Ts:
        i = 1
        while i<= 2:
            params = sample_data_uncertainty_with_constraint(oriParams, ['Topt', 'dCpt'])
            df = calculate_thermal_params(params)
            with model:
                # map temperature constraints
                map_fNT(model,T,df)
                map_kcatT(model,T,df)
                set_NGAMT(model,T)
                set_sigma(model,sigma)

                try: 
                    r = model.optimize().objective_value
                except:
                    print('Failed to solve the problem')
                    r = 0
                    g = 0
                    o = 0
            if abs(r - dfae['r'].values) < 0.1:
                break
            else:
                i = i +1
            print('Growth at ', T-273.15, 'is: ',r)
            g = model.summary().uptake_flux.loc['EX_glc__D_e','flux']
            o = model.summary().uptake_flux.loc['EX_o2_e','flux']
            rs.append(r)
            rg.append(g)
            ro.append(o)
    return rs, rg, ro


def sample_data_uncertainty_with_constraint_random_topt(inpt,columns=None):
    if type(inpt)==tuple:
        params,seed = inpt
        np.random.seed(seed+int(time.time()))
    else: params = inpt
    '''
    # params is a dataframe with following columns:
    # Tm,Tm_std:  melting temperature. Given in K
    #
    # T90, temperature at which 90% of an enzyme is denatured. This is not mandentory. If missing, protein length will be used
    #      to calculate denaturation curve. Given in K
    # 
    # dCpt,dcpt_std: the heat capacity difference between transition state and ground state in enzyme catalyzed reaction. Given
    #                in J/mol/K
    # 
    # Topt,Topt_std: the optimal temprature at which the specific activity is maximized. Given in K
    #
    # Length, protein length. This will not be used in this function, but it will be used later in the calculation of thermal
    #         parameters.
    #
    # xx_std, corresponding uncertainty given by standard deviation.
    # 
    # columns: a list of columns to be sampled, could be any combination of ['Tm','dCpt','Topt]. 
    #          If it is None, then sample all three columns
    # 
    # The script will return an new dataframe with the same columns but with randomly sampled data
    '''
    
    sampled_params = params.copy()
    if columns is None: columns = ['Topt', 'dCpt', 'Tm']
    
    for col in columns:
        if col == 'Topt':
            for ind in sampled_params.index:
                if params.loc[ind, 'topt_source']=='BullShit':
                    M = (params.loc[ind, 'Tm']-params.loc[ind,col])*0.15 + params.loc[ind,col]
                    SD = abs((params.loc[ind, 'Tm']-params.loc[ind,col])*0.05)
                    sampled_params.loc[ind, 'Topt'] = np.random.normal(M,SD)
                else: 
                    M = (params.loc[ind, 'Tm']-params.loc[ind,col])*0.09 + params.loc[ind,col]
                    SD = abs((params.loc[ind, 'Tm']-params.loc[ind,col])*0.02)
                    #print(M, SD)
                    sampled_params.loc[ind, 'Topt'] = np.random.normal(M,SD)
            
        if col == 'dCpt':
            lst = [np.random.normal(-4000,params.loc[ind,col+'_std']) for ind in sampled_params.index]
            sampled_params[col] = lst
        if col == 'Tm':
            for ind in sampled_params.index:
                if params.loc[ind, 'TmTag']=='Mean':
                    sampled_params.loc[ind, 'Tm'] = np.random.normal(params.loc[ind,col],1)
                else:
                    sampled_params.loc[ind, 'Tm'] = np.random.normal(params.loc[ind,col],0.5)
    return sampled_params
