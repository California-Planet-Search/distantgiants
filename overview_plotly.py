"""
Module to make overview plot on Jump
"""

import numpy as np
import pandas as pd

import plotly.express as px
import plotly.graph_objects as go
from astropy.time import Time
pd.options.mode.chained_assignment = None


def obs_type(dataframe, hires_template_counts = 1e5, apf_template_counts = 9e7, radius = 5):
    """
    Takes a dataframe with observations of a set of stars and assigns each observation a 
    type based on iodine, counts, and instrument.
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
    
    ##########################
    # names_ra = dataframe[['star_id', 'ra']]
    #
    # display_window = 100
    # start_date = Time.now().jd-3*display_window/4
    # end_date = Time.now().jd+display_window/4
    #
    #
    # dataframe = dataframe[(start_date < dataframe['bjd'])&(dataframe['bjd']<end_date)]
    #
    # missing = names_ra[~names_ra['star_id'].isin(dataframe.star_id)]#.drop_duplicates(subset='star_id')
    #
    #
    # dataframe=dataframe.append(missing).fillna(1)
    #
    # # dataframe.drop_duplicates(subset='star_id').to_csv('whynot.csv')
    # # dataframe =
    ###########################
    
    for i in range(len(dataframe)):
        row = dataframe.iloc[i]
        if row['iodine_in'] == 't':
            dataframe['obs_type'][i] = 'RV'
        elif row['iodine_in'] == 'f':
            if row['instrument'] == 'HIRES':
                if row['counts'] >= hires_template_counts:
                    dataframe['obs_type'][i] = 'Template'
                else:
                    dataframe['obs_type'][i] = 'Recon'
            elif row['instrument'] == 'APF':
                if row['counts'] >= apf_template_counts:
                    dataframe['obs_type'][i] = 'Template'
                else:
                    dataframe['obs_type'][i] = 'Recon'

    
    # HIRES-only df for getting jitter
    for name in dataframe.drop_duplicates(subset='star_id')['star_id']:
        
        single_star_df_hires = dataframe.query('star_id == "{}" and instrument == "hires_j"'.format(name))
        single_star_df_hires['temp_index'] = [i for i in range(len(single_star_df_hires))]
        
        # For every observationm, check the one before and the after
        for k in range(len(single_star_df_hires))[1:-1]:

            # Look through each hires obs and see if it has an observation on either side within time radius
            diff_list = [abs(single_star_df_hires.query('temp_index == {}'.format(k))['bjd'].values\
                        -single_star_df_hires.query('temp_index == {}'.format(k-1+j))['bjd'].values)\
                         for j in range(3)]
            truth_list = [diff_list[j] < radius and\
                         single_star_df_hires.query('temp_index == {}'.format(k-1+j))['iodine_in'].values == 't' \
                         for j in range(3)]
            
            if truth_list.count(True) == 3:
                # jitter = True
                
                # These are the (consecutive) temporary indices for the jitter sequence
                jit_temp_ind = [k-1, k, k+1]
                
                # These are the (not necessarily consecutive) original indices of the jitter sequence
                jit_ind = single_star_df_hires[single_star_df_hires['temp_index'].isin(jit_temp_ind)].index
    
                for i in jit_ind:
                    single_star_df_hires.loc[i, 'obs_type'] = 'Jitter'
                
                dataframe.update(single_star_df_hires)
                # Once we find one jitter test, no need to keep looking at the same star
                break
                
    # dataframe.to_csv('test.csv')
    dataframe = dataframe.query('obs_type == "RV"')
    # dataframe = dataframe.sort_values(by='ra').reset_index(drop=True)
    return dataframe
    
def make_plot_df(dataframe, display_window=1200):
    
    now = Time.now().jd # Define a single "now" time
    
    # bjd_offset = 2450000
    start_date = now-display_window
    
    name_change_fn = lambda name: 'HD '+name if name[:1].isdigit() else name # The stars starting with numbers are HD
    dataframe['star_id'] = list(map(name_change_fn, dataframe.star_id))
    dataframe['star_id'] = dataframe.star_id.str.replace('T00', 'TOI-') # Change TOIs to have standard names
    
    
    df = dataframe[dataframe['bjd']>=start_date] # Only take observations past desired date
    # df = df.sort_values(by='ra', ascending=True)
    df = df[df['star_id'].isin(df.drop_duplicates(subset='star_id').star_id)]
    df = df.replace({'hires_j':'HIRES', 'apf':'APF'})
    df = obs_type(df) # Assign each observation a type
    
    ###
    # df = df[df['star_id']=='TOI-1669']
    ###
    df.to_csv('csv/overview_plot_df.csv')
    
    return df
    

def plot(df, display_window=826):
    
    now = Time.now().jd # Define a single "now" time


    # bjd_offset = 2450000
    start_date = now-display_window
    end_date = now+10
    

    y_length = np.arange(len(df.drop_duplicates(subset='star_id')['star_id']))[::-1]


    fig = go.Figure()
    
    ## Plot HIRES points
    df_hires = df.query('instrument=="HIRES"')
    fig.add_trace(go.Scatter(x=df_hires.bjd, y=df_hires.star_id,
                            mode='markers',
                            marker=dict(size=6,
                                        symbol='square'),
                            marker_color='rgb(228,26,28)',
                            name='HIRES RV'))

    ## Plot APF points
    df_apf = df.query('instrument=="APF"')
    fig.add_trace(go.Scatter(x=df_apf.bjd, y=df_apf.star_id,
                            mode='markers',
                            marker=dict(size=4,
                                        symbol='circle'),
                            marker_color='rgba(200,200,200,1)',
                            name='APF RV'))


    fig.update_traces(marker=dict(line=dict(width=0.75,
                                  color='DarkSlateGrey')),
                      selector=dict(mode='markers'))
    fig.update_yaxes(categoryorder='array', categoryarray= df.star_id.to_list())
    fig.update_layout(legend_traceorder="reversed")

    ## Change grid lines to black
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='rgba(100, 100, 100, 0.5)')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='rgba(100, 100, 100, 0.5)')


    # Make background transparent
    # fig.update_layout({'plot_bgcolor':'rgba(0,0,0,0)',
    #                    'paper_bgcolor':'rgba(0,0,0,0)'})


    # # Show target RA in the right margin
    # fig.add_annotation(dict(font=dict(color="black",size=12),
    #                         x=1.09,
    #                         y=0.955,
    #                         showarrow=False,
    #                         text="Right Ascension",
    #                         xref="paper",
    #                         yref="paper",
    #                         xanchor='center',
    #                         yanchor='bottom'
    #                        ))
    
    # x_loc = end_date+10
    # names_ra = df.drop_duplicates(subset='star_id')\
    #              .sort_values(by='ra', ascending=False)\
    #              .reset_index(drop=True)[['star_id', 'ra']]
    # names_ra_len = len(names_ra) # 47
    # for i in range(names_ra_len):
    #
    #     y_pos = 0.055 + (1-2*0.055)*i/46 # Hand-tweaked name position values
    #     ra_hr = names_ra['ra'][i]
    #     ra_str = str(int(ra_hr)).zfill(2) + 'h ' + str(int(np.round(ra_hr%1*60))).zfill(2) + 'm'
    #     fig.add_annotation(dict(font=dict(color="black",size=10),
    #                             x=1.02,
    #                             y=y_pos,
    #                             showarrow=False,
    #                             text=ra_str,
    #                             xref="paper",
    #                             yref="paper",
    #                             xanchor='left',
    #                             yanchor='middle'
    #                            ))
    
    # Add 1-month line segment for reference
    # Make y-coordinates blank strings so Plotly gives the line segment its own line without a label.
    # fig.add_trace(go.Scatter(x=[end_date-40, end_date-10],
    #                          y=['', ''], name='1 month', mode = 'lines+markers',
    #                       marker = {'color' : 'green'}))
    
    ## Human-readable dates on the x-axis
    ###############################################
    ## Earliest date and latest date
    ed_str = Time(start_date, format='jd').iso
    ld_str = Time(end_date, format='jd').iso


    # Pandas lets me create a list of the first day of the month in a set of months
    tick_timestamps = pd.date_range(ed_str, ld_str, freq='5M') - pd.offsets.MonthBegin(1)
    tick_strings = [tstamp.date().isoformat() for tstamp in tick_timestamps] # Convert pd object to string
    tick_values = [Time(string, format='iso').jd for string in tick_strings] # Convert strings to jds
    
    
    fig.update_xaxes(range=[start_date, end_date])
    fig.update_xaxes(
        ticktext=tick_strings,
        tickvals=tick_values,
    )
    ################################################
    
    fig['layout']['yaxis']['autorange'] = "reversed"
    fig.update_yaxes(type='category')
    height = len(df.drop_duplicates(subset='star_id'))*20
    # fig.update_layout(yaxis={"dtick":1},margin={"t":0,"b":0})#, height=height)

    fig.update_layout(margin={"t":0,"b":30,"r":150},
                      font=dict(size=10),
                      legend=dict(yanchor='bottom',
                                  y=-0.01,
                                  xanchor='left',
                                  x=1.0,
                                  font=dict(size=14)))

    # fig.show()
    # fig.write_image('abacus.png', scale=5, height=height, width=850)
    fig.write_image('abacus.png', scale=3, height=1500, width=1000)
    return


if __name__ == "__main__":
    display_window = 1095
    
    sql_df = pd.read_csv('csv/Distant_Giants_Observing_Requests.csv')
    make_plot_df(sql_df, display_window=display_window)
    
    plot_df = pd.read_csv('csv/overview_plot_df.csv')
    plot(plot_df, display_window=display_window)
   
   
   
   