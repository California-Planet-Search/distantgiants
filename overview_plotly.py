"""
Module to make overview plot on Jump
"""


import numpy as np
import pandas as pd
import plotly.express as px
pd.options.mode.chained_assignment = None


sql_df = pd.read_csv('csv/sql_df.csv')

def obs_type(dataframe, hires_template_counts = 6e4, apf_template_counts = 9e7, radius = 5):
    """
    Takes a dataframe with observations of a set of stars and assigns each observation a type based on iodine, counts, and instrument.
    The output is the same dataframe with an extra column for observation type.
    
    Arguments:
    
        dataframe (pd.DataFrame): dataframe with all observations of a set of stars
        hires_template_counts (float): min counts to be considered a template on HIRES
        apf_template_counts (float): min counts to be considered a template on APF
        radius (int): Used to identify jitter tests. Search radius (in hours) for other jitter exposures
    
    Returns:
        dataframe (pd.DataFrame): updated dataframe with observation types
    """
    
    dataframe = dataframe.sort_values(by=['star_id', 'bjd']).reset_index(drop=True)
    dataframe['obs_type'] = np.zeros((len(dataframe),1))
    jitter = False

    for i in range(len(dataframe)):
        row = dataframe.iloc[i]
        if row['iodine_in'] == 't':
            dataframe['obs_type'][i] = 'rv'
        elif row['iodine_in'] == 'f':
            if row['instrument'] == 'hires_j':
                if row['counts'] >= hires_template_counts:
                    dataframe['obs_type'][i] = 'template'
                else:
                    dataframe['obs_type'][i] = 'recon'
            elif row['instrument'] == 'apf':
                if row['counts'] >= apf_template_counts:
                    dataframe['obs_type'][i] = 'template'
                else:
                    dataframe['obs_type'][i] = 'recon'

    
    # HIRES-only df for getting jitter
    for name in dataframe.drop_duplicates(subset='star_id')['star_id']:
        
        single_star_df_hires = dataframe.query('star_id == "{}" and instrument == "hires_j"'.format(name))
        single_star_df_hires['temp_index'] = [i for i in range(len(single_star_df_hires))]
        
        for k in range(len(single_star_df_hires))[1:-1]:
            # Look through each hires obs and see if it has an observation on either side within time radius
            diff_list = [abs(single_star_df_hires.query('temp_index == {}'.format(k))['bjd'].values-single_star_df_hires.query('temp_index == {}'.format(k-1+j))['bjd'].values) for j in range(3)]
            truth_list = [diff_list[j] < radius and single_star_df_hires.query('temp_index == {}'.format(k-1+j))['iodine_in'].values == 't' for j in range(3)]
            
            if truth_list.count(True) == 3:
                jitter = True
                
                # These are the (consecutive) temporary indices for the jitter sequence
                jit_temp_ind = [k-1, k, k+1]
                
                # These are the (not necessarily consecutive) original indices of the jitter sequence
                jit_ind = single_star_df_hires[single_star_df_hires['temp_index'].isin(jit_temp_ind)].index
    
                for i in jit_ind:
                    single_star_df_hires.loc[i, 'obs_type'] = 'jitter'
                
                dataframe.update(single_star_df_hires)
                break

    return dataframe

def plot(dataframe):
                       
    df = dataframe.sort_values(by='ra')
    fig = px.scatter(df, x='bjd', y='star_id', color='instrument', symbol='obs_type')
    fig['layout']['yaxis']['autorange'] = "reversed"
    fig.show()


if __name__ == "__main__":
    
   plot(obs_type(sql_df))
   
   
   
   