# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import pandas as pd
import numpy as np
import pickle
from etcpy import etc
import os
from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import r2_score


# ### Load data

def load_exp_batch_data(infile):
    # return two dataframes for aerobic and anaerobic conditions
    # with temperature as index
    #
    data = dict()
    for line in open(infile):
        cont = line.split()
        data[cont[0]] = [float(item) for item in cont[1:]]
    
    def build_df(ind,col):
        df = pd.DataFrame()
        df[col] = data[col]
        df.index = data[ind]
        return df
    
    dfae_exp = build_df('Ts_ae','r_ae')
    dfan_exp = build_df('Ts_an','r_an')
    #print(dfae_exp,dfan_exp)
    return dfae_exp,dfan_exp


def load_exp_data(filename):
    # Load the csv file
    df = pd.read_csv(filename, index_col=0)
    
    # Get unique values in the 'ID' column
    unique_ids = df['ID'].unique()
    
    # Create an empty dictionary to store the split DataFrames
    split_dfs = {}
    
    # Split the DataFrame based on unique 'ID' values
    for unique_id in unique_ids:
        split_dfs[unique_id] = df[df['ID'] == unique_id].drop('ID', axis=1)
    
    # Return the dictionary of split DataFrames
    return split_dfs


def aerobic_exp_data():
    dfae_batch = exp_batch_data(os.path.join(path,'data/ExpGrowth.tsv'))
    return  {'data':dfae_batch['r_ae'].values}


def anaerobic_exp_data():
    dfae_batch,dfan_batch = load_exp_batch_data(os.path.join(path,'data/ExpGrowth.tsv'))
    return  {'data':dfan_batch['r_an'].values}


def chemostat_exp_data():
    dfchemo = pd.read_csv(os.path.join(path,'data/Chemostat_exp_data.txt'),sep='\t')
    rxn_lst = [
        'r_1714_REV',#Glucose
        'r_1672', #CO2
        'r_1761', # Ethanol
    ]
    columns = ['Glucose','CO2','Ethanol']
    exp_flux = []
    for i,rxn_id in enumerate(rxn_lst):
        exp_flux += list(dfchemo[columns[i]])
    return  {'data':np.array(exp_flux)}


# +
#path = os.path.dirname(os.path.realpath(__file__)).replace('code','')
#params = pd.read_csv(os.path.join(path,'data/model_enzyme_params.csv'),index_col=0)
#dfae_batch,dfan_batch = load_exp_batch_data(os.path.join(path,'data/ExpGrowth.tsv'))
#dfchemo = pd.read_csv(os.path.join(path,'data/Chemostat_exp_data.txt'),sep='\t',index_col=0)
# -

# ### Models 

def format_input(thermalParams):
    new_params = params.copy()
    for key,val in thermalParams.items():
        [ind,col] = key.split('_')
        new_params.loc[ind,col] = val
    
    # Update T90
    new_params['T90'] = params['T90']-params['Tm'] + new_params['Tm']

    df = etc.calculate_thermal_params(new_params)
    return df,new_params


def aerobic(thermalParams):
    # thermalParams: a dictionary with ids like uniprotid_Topt 
    df,new_params = format_input(thermalParams)
    mae = pickle.load(open(os.path.join(path,'models/aerobic.pkl'),'rb'))
    
    try: rae = etc.simulate_growth(mae,dfae_batch.index+273.15,df=df,sigma=0.5)
    except: rae = np.zeros(dfae_batch.shape[0])
    
    rae = [0 if x is None else x for x in rae]
    rae = [0 if x<1e-3 else x for x in rae]
    print(rae)
    rexp = aerobic_exp_data()['data']
    
    print('r2_batch:',r2_score(rexp,rae))
    print('MSE_ae',MSE(rexp,rae))
    return {'data':np.array(rae)}


def anaerobic(thermalParams):
    df,new_params = format_input(thermalParams)
    man = pickle.load(open(os.path.join(path,'models/anaerobic.pkl'),'rb'))

    try: ran = etc.simulate_growth(man,dfan_batch.index+273.15,df=df,sigma=0.5)
    except: ran = np.zeros(dfan_batch.shape[0])
    ran = [0 if x is None else x for x in ran]
    rexp = anaerobic_exp_data()['data']
    
    print('r2_batch_an:',r2_score(rexp,ran))
    print('MSE_an',MSE(rexp,ran))

    return  {'data':np.array(ran)}


def anaerobic_reduced(thermalParams):
    df,new_params = format_input(thermalParams)
    man = pickle.load(open(os.path.join(path,'models/anaerobic.pkl'),'rb'))
    sel_temp = [5.0,15.0,26.3,30.0,33.0,35.0,37.5,40.0]
    try: ran = etc.simulate_growth(man,np.array(sel_temp)+273.15,df=df,sigma=0.5)
    except: ran = np.zeros(len(sel_temp))
    ran = [0 if x is None else x for x in ran]
    rexp = dfan_batch.loc[sel_temp,'r_an'].values
    #anaerobic_exp_data()['data']
    
    try:
        print('ra2_batch_an:',r2_score(rexp,ran))
        print('MSE_an',MSE(rexp,ran))
    except: print('Model error:',len(rexp),len(ran))

    return  {'data':np.array(ran)}


def chemostat(thermalParams):
    df,new_params = format_input(thermalParams)
    mae = pickle.load(open(os.path.join(path,'models/aerobic.pkl'),'rb'))
    exp_flux = chemostat_exp_data()['data']
    
    growth_id = 'r_2111'
    glc_up_id = 'r_1714_REV'
    prot_pool_id = 'prot_pool_exchange'
    dilut = 0.1
    sigma = 0.5
    
    try:
        solution = etc.simulate_chomostat(mae,dilut,new_params,dfchemo.index+273.15,
                                              sigma,growth_id,glc_up_id,prot_pool_id)

        # Extract fluxes
        rxn_lst = [
                'r_1714_REV',#Glucose
                'r_1672', #CO2
                'r_1761', # Ethanol
            ]
        columns = ['Glucose','CO2','Ethanol']

        pred_flux = []
        for i,rxn_id in enumerate(rxn_lst):
            x = [s.fluxes[rxn_id] for s in solution]
            x.extend([0]*(len(dfchemo.index)-len(x)))
            pred_flux += x
        print(pred_flux)
    
    except: pred_flux = [0 for item in exp_flux]
    
    print('r2_flux:',r2_score(exp_flux,pred_flux))
    print('MSE_chemo',MSE(exp_flux,pred_flux))

    return  {'data':np.array(pred_flux)}


def simulate_at_three_conditions(args,distance_function,Yobs,min_epsilon):
    
    data_batch = aerobic(args)['data']
    d_ae = distance_function(Yobs['rae'],data_batch)
    if d_ae < min_epsilon['rae']: return False, {}, {}
    
    data_batch_an = anaerobic_reduced(args)['data']
    d_an = distance_function(Yobs['ran'],data_batch_an)
    if d_an < min_epsilon['ran']: return False, {}, {}
    
    data_chemo = chemostat(args)['data']
    d_c = distance_function(Yobs['chemostat'],data_chemo)
    if d_c < min_epsilon['chemostat']: return False, {}, {}
    return True,{'rae':data_batch,'chemostat':data_chemo,'ran':data_batch_an},{'rae':d_ae,'chemostat':d_c,'ran':d_an}


def simulate_at_three_conditions_2(args):
    data_batch = aerobic(args)['data']
    data_batch_an = anaerobic_reduced(args)['data']
    data_chemo = chemostat(args)['data']
    
    return {'rae':data_batch,'chemostat':data_chemo,'ran':data_batch_an}


