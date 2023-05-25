"""
Judah Van Zandt
Python 3
This module is used to compile general information about the state of each target in the SC2A program.
Written: June 1, 2020
Last modified: June 24, 2021
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from git import repo

from astropy.time import Time
from astroplan import Observer, FixedTarget
import astropy.units as u
from astropy.coordinates import SkyCoord

import ephem
import observing as obs

import distantgiants

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)




def make_overview(plot = True, observability = False):
    """
    With a lot of intermediate steps, this function generates overview_df and an overview plot. They both contain important
    information about each target in SC2A.
    """
    
    distantgiants = pd.read_csv('csv/distantgiants.csv')
    
    gamma_df = pd.read_csv('csv/Distant_Giants_gamma.csv')

    sql_df = pd.read_csv('csv/Distant_Giants_Observing_Requests.csv')[['star_id', 'instrument', 'counts', 'iodine_in', 'bjd']]

    sql_df = pd.merge(distantgiants, sql_df, on = 'star_id', how = 'inner')
    
    ## Special case: 207897 has a long baseline, but it's sparse and not useful until June 1 2020. Cut it down so the code thinks it's not cooked
    ## Similar for 93963: not useful until Jul. 1 2020
    ## Similar for T001444: not useful until Aug. 1 2020
    iod_condition = (sql_df['iodine_in'] == 't')
    drop_inds_1 = sql_df[(sql_df['star_id'] == '207897') 
                       & (sql_df['bjd'] < 2459001.5) 
                       & iod_condition].index
                       
    drop_inds_2 = sql_df[(sql_df['star_id'] == '93963') 
                       & (sql_df['bjd'] < 2459031.5) 
                       & iod_condition].index
                       
    drop_inds_3 = sql_df[(sql_df['star_id'] == 'T001444') 
                       & (sql_df['bjd'] < 2459062.5) 
                       & iod_condition].index
    drop_inds = drop_inds_1.append([drop_inds_2, drop_inds_3])
    
    sql_df = sql_df.drop(drop_inds).sort_values(by='star_id')
    
    ## The Radvel Coordination Squad made a list of each target's preferred radvel setup file. Use them to get gdot/gddot
    rcs_df = pd.read_csv('csv/rcs_master.csv')[['CPS_NAME', 'Fit preferred (official file name without .py)']]\
                .rename(columns={'CPS_NAME':'star_id','Fit preferred (official file name without .py)':'fit_pref'})
    rcs_df['fit_pref'] = rcs_df['fit_pref'].fillna('default_cps') # Use default_cps if a target has no preferred fit


    # gam_pref is shorter than distantgiants, meaning that not every target has a preferred radvel fit. Not urgent, but address later
    gam_pref = gamma_df.merge(rcs_df, left_on = ['star_id', 'runname'], right_on = ['star_id', 'fit_pref'])\
                       .drop(['runname'], axis=1)\
                       .drop_duplicates(subset='star_id')

    sql_df = pd.merge(sql_df, gam_pref, on = 'star_id', how = 'left')


    # We want to find templates by excluding 10k recon spectra (on HIRES). Also take obs with counts < 0 b/c in May/June 2021 the exposure meter was broken, so all obs have exp_count = -1 (and I didn't request any recon spectra). At most 9 systems have such templates (I can find them later by removing the counts < 0 option).
    # 5.4E8 counts on APF is the same as 60k on HIRES, so a recon spectrum on APF would have ~9E7 counts
    template_apf = sql_df.query('instrument == "apf" and counts > 9e7 and iodine_in == "f"')\
                         .sort_values(by=['star_id', 'iodine_in', 'counts'])\
                         .drop_duplicates(subset=['star_id', 'instrument'], keep = 'last')\
                                                [['star_id', 'bjd', 'counts', 'iodine_in']]\
                         .rename(columns = {'bjd':'apf_template_bjd'}).replace(np.nan, 0)


    template_hires = sql_df.query('instrument == "hires_j" and (counts > 6e4 or counts < 0) and iodine_in == "f"')\
                           .sort_values(by=['star_id', 'iodine_in', 'counts'])\
                           .drop_duplicates(subset=['star_id', 'instrument'], keep = 'last')\
                                                  [['star_id', 'bjd', 'counts', 'iodine_in']]\
                           .rename(columns = {'bjd':'hires_template_bjd'}).replace(np.nan, 0)

    
    template_df = pd.merge(template_apf[['star_id', 'apf_template_bjd']],\
                           template_hires[['star_id', 'hires_template_bjd']], on = 'star_id', how = 'outer')
    
    distantgiants = pd.merge(distantgiants, template_df, on = 'star_id', how = 'left')
    
    # Assign 'no' if a target has neither a HIRES nor APF template.
    # Note that the 'no' string must match with the jitter if statement in target_request_generator.py (near line 130)
    template_lambda = lambda x, y: 'no' if np.isnan(x) and np.isnan(y) else 'YES'
    distantgiants['have_template'] = list(map(template_lambda, distantgiants['apf_template_bjd'],\
                                                               distantgiants['hires_template_bjd']))
    
    
    full_star_list = pd.merge(sql_df.drop_duplicates(subset = 'star_id'), \
                              distantgiants.drop_duplicates(subset = 'star_id'), \
                              on = 'star_id', how = 'outer')['star_id'].to_frame()
    
    # Stars that need recon
    sql_df['day'] = Time(sql_df.bjd, format='jd').iso
    sql_df_yesiod = sql_df.query("iodine_in == 't'").reset_index(drop = True)
    sql_df_noiod = sql_df.query("iodine_in == 'f'").reset_index(drop = True)

    recon_df = pd.merge(sql_df_noiod.query('counts >= 4500')[['star_id', 'instrument', 'bjd']].\
                        sort_values(['star_id', 'bjd']).drop_duplicates\
                        (subset = 'star_id', keep = 'first'), full_star_list['star_id'], \
                        on = 'star_id', how = 'outer').rename(columns = {'bjd':'bjd_recon'}).\
                        rename(columns = {'instrument':'instrument_recon'})

    recon_df['have_recon'] = list(map(lambda x: 0 if pd.isna(x) else 1, recon_df['bjd_recon']))
    
    
    # Stars that need jitter tests
    # Create a groupby so that size() gives the number of rows in each (star, day, instrument) group
    grouped_list = sql_df_yesiod.groupby(by=['star_id', 'day', 'instrument'])
    # has_jitter just asks if a single star got 3 observations in a single (UT) day.
    has_jitter = grouped_list.size()\
                .to_frame()\
                .reset_index()\
                .rename(columns = {0:'count'})\
                .query('count >= 3 and instrument == "hires_j"')\
                .sort_values(['star_id', 'day'])\
                .drop_duplicates(subset = 'star_id', keep = 'first')

    jitter_df = pd.merge(has_jitter, full_star_list, on = 'star_id', how = 'outer')
    
    # 1 if target has jitter, otherwise 0
    have_jitter_lambda = lambda x: 0 if pd.isna(x) else 1
    # NaN if object doesn't have jitter, otherwise bjd of jitter date
    bjd_jitter_lambda = lambda x: x if pd.isna(x) else Time(x, format = 'iso', out_subfmt = 'date').jd
    
    jitter_df['have_jitter'] = list(map(have_jitter_lambda, jitter_df['count']))
    jitter_df['bjd_jitter']  = list(map(bjd_jitter_lambda, jitter_df['day']))


    # dates_df gives every observation and its jd, sorted so that the most and least recent can be picked for each target
    # latest_obs finds the most recent observation with both HIRES and APF for each target
    dates_df = pd.merge(sql_df_yesiod[['star_id', 'instrument', 'bjd']]\
                 .sort_values(by = ['star_id', 'instrument', 'bjd'])\
                 .reset_index(drop = True), full_star_list, on = 'star_id', how = 'outer')
    
    # Find latest and earliest for each instrument
    latest_obs = dates_df.drop_duplicates(subset = ['star_id', 'instrument'], keep = 'last')
    earliest_obs = dates_df.drop_duplicates(subset = ['star_id', 'instrument'], keep = 'first')
    
    # Find latest and earliest between both instruments
    latest_obs_tot = dates_df.sort_values(by = ['star_id', 'bjd']).drop_duplicates(subset = ['star_id'], keep = 'last')
    earliest_obs_tot = dates_df.sort_values(by = ['star_id', 'bjd']).drop_duplicates(subset = ['star_id'], keep = 'first')

    # Adding columns for the most recent HIRES and APF iodine-in obs to overview_df. Displayed as days since last obs.
    last_obs_hires = pd.DataFrame({'star_id': latest_obs.query("instrument == 'hires_j'")['star_id'], \
                                   'last_obs_hires': -latest_obs.query("instrument == 'hires_j'")['bjd'].\
                                   subtract(Time(Time.now().iso, out_subfmt='date').jd)})

    last_obs_apf = pd.DataFrame({'star_id': latest_obs.query("instrument == 'apf'")['star_id'], \
                                 'last_obs_apf': -latest_obs.query("instrument == 'apf'")['bjd'].\
                                 subtract(Time(Time.now().iso, out_subfmt='date').jd)})
    
    #############################################################################################

    # Creating overview_df to show general information
    overview_df = recon_df[['star_id', 'have_recon']]
    overview_df = pd.merge(overview_df, jitter_df[['star_id', 'have_jitter']], how = 'outer', on = 'star_id')
    overview_df = pd.merge(overview_df, distantgiants[['star_id', 'have_template']], how = 'outer', on = 'star_id') # Stars that need templates 
    overview_df = pd.merge(overview_df, 
                           pd.DataFrame(sql_df_yesiod.groupby(by=['star_id', 'instrument']).size()).reset_index()\
                           .rename(columns = {0: 'tot_iodine_hires'})\
                           .query("instrument == 'hires_j'")\
                           .drop(columns = {'instrument'}), how = 'outer', on = 'star_id')

    overview_df = pd.merge(overview_df, 
                           pd.DataFrame(sql_df_yesiod.groupby(by=['star_id', 'instrument']).size()).reset_index()\
                           .rename(columns = {0: 'tot_iodine_apf'}).query("instrument == 'apf'")\
                           .drop(columns = {'instrument'}), how = 'outer', on = 'star_id')

    overview_df = pd.merge(overview_df, last_obs_hires, on = 'star_id', how = 'outer')

    overview_df = pd.merge(overview_df, last_obs_apf, on = 'star_id', how = 'outer')
    
    overview_df = pd.merge(overview_df, 
                           earliest_obs_tot.rename(columns={'bjd':'earliest_obs'}),
                           on='star_id', how='inner')
    
    
    ## Determine the baseline of each target
    baseline_hires = pd.DataFrame({'star_id': latest_obs.query("instrument == 'hires_j'")['star_id'].values,\
                               'baseline_hires': latest_obs.query("instrument == 'hires_j'")['bjd'].values
                                             - earliest_obs.query("instrument == 'hires_j'")['bjd'].values})

    baseline_apf = pd.DataFrame({'star_id': latest_obs.query("instrument == 'apf'")['star_id'].values,\
                               'baseline_apf': latest_obs.query("instrument == 'apf'")['bjd'].values\
                                           - earliest_obs.query("instrument == 'apf'")['bjd'].values})
                               
    baseline_tot = pd.DataFrame({'star_id': latest_obs_tot['star_id'].values,\
                               'baseline_tot': latest_obs_tot['bjd'].values\
                                           - earliest_obs_tot['bjd'].values})


    overview_df = pd.merge(overview_df, baseline_hires, on = 'star_id', how = 'outer')
    overview_df = pd.merge(overview_df, baseline_apf, on = 'star_id', how = 'outer')
    overview_df = pd.merge(overview_df, baseline_tot, on = 'star_id', how = 'outer')
    
    # Cleaning up overview_df
    overview_df.replace({'have_recon':{0:'no', 1:'YES'},
                        'have_jitter':{0:'no', 1:'YES'}}, 
                        inplace=True)

    overview_df.update(overview_df[['tot_iodine_hires', 'tot_iodine_apf']].fillna(0))
    overview_df.update(overview_df[['last_obs_hires', 'last_obs_apf']].fillna(np.inf))
    overview_df.update(overview_df[['baseline_hires', 'baseline_apf']].fillna('N/A'))
    
    overview_df['tot_iodine'] = overview_df['tot_iodine_hires']+overview_df['tot_iodine_apf']
    
    # Add a new column 'cooked?' to indicate any target with both a 3+ year baseline and 35+ observations. 
    cooked = []
    for i in range(len(overview_df)):
        if isinstance(overview_df['baseline_hires'][i], str):
            baseline_hires = 0
        elif isinstance(overview_df['baseline_hires'][i], float):
            baseline_hires = float(overview_df['baseline_hires'][i])
        else:
            print('Error: type is neither str nor float')
            print(type(overview_df['baseline_hires'][i]))
        
        if isinstance(overview_df['baseline_apf'][i], str):
            baseline_apf = 0
        elif isinstance(overview_df['baseline_apf'][i], float):
            baseline_apf = float(overview_df['baseline_apf'][i])
        else:
            print('Error: type is neither str nor float')
            print(type(overview_df['baseline_apf'][i]))
        
        if (overview_df['baseline_hires'][i] > 1200) and\
           (overview_df['tot_iodine'][i] >= 35):
           
           cooked.append('COOKED')
           #print(overview_df[i]['star_id'])
            
        else:
            cooked.append('cookin')
            # print(overview_df[i]['star_id'], baseline_tot[i]['baseline_tot'], overview_df[i]['tot_iodine'])
    overview_df['cooked?'] = cooked


    overview_cols = ['star_id', 'vmag', 'ra', 'dec', 'fit_pref', 
                     'dvdt', 'u_dvdt', 'curv', 'u_curv']
    overview_df = pd.merge(overview_df, 
                           sql_df.drop_duplicates(subset = 'star_id')[overview_cols], 
                           on = 'star_id')
        
    # Using strings instead of True/False bc TRUE vs. f is easier to compare in a csv file than TRUE and FALSE
    three_sig_lambda = lambda val, err: "TRUE" if abs(val) > 3*abs(err) else "f"
    overview_df['3sig_trend'] = list(map(three_sig_lambda, overview_df['dvdt'], overview_df['u_dvdt']))
    overview_df['3sig_curv'] = list(map(three_sig_lambda, overview_df['curv'], overview_df['u_curv']))
    
    overview_df['ra'] = overview_df['ra']*15
    overview_df = overview_df.rename(columns = {'ra':'ra_deg', 'dec':'dec_deg'}).sort_values(by='star_id')


    # Creating plot_df with all of the information to create an image with an overview for each target
    plot_df = pd.merge(overview_df, recon_df.drop(columns = 'have_recon'), on = 'star_id')
    plot_df = pd.merge(plot_df, template_df, on = 'star_id', how = 'left')

    plot_cols = ['star_id', 'instrument_recon', 'baseline_tot', 'tot_iodine', 'tot_iodine_hires', 'tot_iodine_apf', 
                'bjd_recon', 'bjd_jitter', 'apf_template_bjd', 'hires_template_bjd', 
                'have_template', 'cooked?', 'last_obs_hires', 'last_obs_apf', 'ra_deg', 
                'dec_deg', '3sig_trend', '3sig_curv']
    plot_df = pd.merge(plot_df, jitter_df.drop(columns = 'have_jitter'), on = 'star_id')[plot_cols]

    # Initializing variables for the plot

    y_length = np.arange(len(plot_df['star_id']))[::-1]

    y_increment = (y_length[1]-y_length[0])
    display_window = 100
    keck = Observer.at_site('W. M. Keck Observatory') # Site for astroplan to get sunrise/set
    observatory = obs.setupObservatory('keck') # Site for Ian's observing code

    start_date = Time.now().jd-2*display_window/3
    end_date = Time.now().jd+display_window/3

    sidereal_offset = 0.002731 # This is how much earlier a target rises each day (in days)

    # Dates to display on the x-axis. Roughly one tick per month by dividing the window by 30
    date_intervals_jd = np.linspace(start_date, end_date, int(display_window/30 + 1))
    
    date_intervals_iso = [str(Time(date_intervals_jd[i], format = 'jd').iso).split(' ')[0]
                                                 for i in range(len(date_intervals_jd))]
    

    recon_color = 'green'
    template_color = 'black'
    jitter_color = 'pink'
    color_hires = 'red'
    color_apf = 'blue'

    ######################
    # Making the plot
    ######################

    fig, ax = plt.subplots(figsize=(14, 8), dpi= 100, facecolor='w', edgecolor='k')

    z_order_list = z_order_list = ['observability', 'h_line_plot', 'rv_pts', 'recon_pts', 'template_pts', 'jitter_pts', 'recon_arrows', 'jitter_arrows']

    plt.hlines(y_length, start_date, end_date, colors = ['black', 'gray'], linewidths = 0.2, zorder = z_order_list.index('h_line_plot'))
    plt.grid(visible = True, which = 'major', axis = 'x', linewidth = .2)

    ax.set_xlim(start_date, end_date)

    plot_df = plot_df.sort_values(by = 'ra_deg', ignore_index = True)

    # # Plotting the dates of each star's recon and jitter
    # recon_pts = ax.scatter(plot_df['bjd_recon'], y_length, color = recon_color, s=6, zorder = z_order_list.index('recon_pts'))
    # jitter_pts = ax.scatter(plot_df['bjd_jitter'], y_length, color = jitter_color, s=6, zorder = z_order_list.index('jitter_pts'))
    # apf_template_pts = ax.scatter(plot_df['apf_template_bjd'], y_length, color = 'black', marker = '*', s=30, zorder = z_order_list.index('template_pts'))
    # hires_template_pts = ax.scatter(plot_df['hires_template_bjd'], y_length, color = 'black', marker = '*', s=30, zorder = z_order_list.index('template_pts'))


    ax.text(start_date-23, y_length[0]+1, 'J?', size = 12, horizontalalignment = 'right')
    ax.text(start_date-18, y_length[0]+1, 'T?', size = 12, horizontalalignment = 'right')
    ax.text(end_date+3, y_length[0]+1, 'Nobs', size = 12)
    ax.text(end_date+17, y_length[0]+1, 'Last Obs', size = 12)
    ax.text(end_date+37, y_length[0]+1, r'$\dot{\gamma}$', size = 14)
    ax.text(end_date+43, y_length[0]+1, r'$\ddot{\gamma}$', size = 14)

    for i in range(len(plot_df)):

        # Put star name on right and left vertical axes, with cooked targets faded
        if plot_df['cooked?'][i] == 'cookin':
            alpha = 1
        elif plot_df['cooked?'][i] == 'COOKED':
            alpha = 0.5
    
        if plot_df['have_template'][i] == 'YES':
            template_facecolors = 'k'
        else:
            template_facecolors = 'none'

        if not np.isnan(plot_df['bjd_jitter'][i]):
            jitter_facecolors = 'k'
        else:
            jitter_facecolors = 'none'

        ###############
        if plot_df['3sig_trend'][i] == "TRUE":
            trend_facecolors = 'k'
        else:
            trend_facecolors = 'none'

        if plot_df['3sig_curv'][i] == "TRUE":
            curv_facecolors = 'k'
            print(plot_df['star_id'][i], ' has a 3σ curvature term!')
        else:
            curv_facecolors = 'none'
        ################

        last_obs = np.min([plot_df.iloc[i]['last_obs_hires'], plot_df.iloc[i]['last_obs_apf']])

        gap_list = [20, 30, 60]
        color_list = ['gold', 'orange', 'red']

        last_obs_color = 'k'
        for j in range(len(gap_list)):
            if last_obs > gap_list[j]:
                last_obs_color = color_list[j]

        # Plotting the name of each star on the y-axis
        ax.text(start_date-1, y_length[i]+0.3*y_increment, plot_df['star_id'][i], size = 7, alpha = alpha, horizontalalignment = 'right')

        # Template and jitter
        ax.scatter(start_date-20, y_length[i], s=14, marker = 'o', color = 'k', facecolors = template_facecolors, clip_on=False)
        ax.scatter(start_date-25, y_length[i], s=14, marker = 'o', color = 'k', facecolors = jitter_facecolors, clip_on=False)

        # Number of iodine-in obs
        ax.text(end_date+3, y_length[i]-.3, s='{:.0f}'.format(plot_df.iloc[i]['tot_iodine']), size=12)
        # Days since last observation
        ax.text(end_date+17, y_length[i]-.3, s='{:.0f}'.format(last_obs), size=12, color=last_obs_color)

        # Trend or curvature detected?
        ax.scatter(end_date+38, y_length[i], s=14, marker = 'o', color = 'k', facecolors = trend_facecolors, clip_on=False)
        ax.scatter(end_date+44, y_length[i], s=14, marker = 'o', color = 'k', facecolors = curv_facecolors, clip_on=False)

        # # If a target got recon or jitter before the earliest displayed date, use an arrow
        # if plot_df['bjd_recon'][i] < start_date:
        #     ax.scatter(start_date+0.5, y_length[i] - 0.25*y_increment, marker = '<', color = recon_color, s=10, zorder = z_order_list.index('recon_arrows'))
        #
        # if plot_df['bjd_jitter'][i] < start_date:
        #     ax.scatter(start_date+0.5, y_length[i] + 0.25*y_increment, marker = '<', color = jitter_color, s=10, zorder = z_order_list.index('recon_arrows'))

        ### Observability ###
        # Plotting observability bars takes a while; pick how many stars you want bars for here.
        # Also I've tried plotting observability according to cpsutils as in target_request_generator. That way agrees better with Jump (no surprise there) but it takes ~twice as long to run. We're stuck with this for now I think.

        if observability == True:
    
            star_name = plot_df['star_id'][i]
            ra_deg = plot_df['ra_deg'][i]
            dec_deg = plot_df['dec_deg'][i]
            full_marker_list = []
            star_marker_list = []
            print(star_name)

            star_target = obs.setupTarget(ra_deg, dec_deg, star_name)

            # List of all days from start_date to end_date, incrementing by one sidereal day to get back to the same sidereal time. I believe this will avoid eventually missing a day, and will instead double count a day each time the sidereal offset adds up to a whole day.
            day_list = [start_date + i*(1-sidereal_offset) for i in range(int(end_date-start_date))]

            for day in day_list:
       
               # Find the next time the sun sets
               sunset_time = keck.sun_set_time(Time(day, format='jd'), which = 'next', horizon = -12*u.deg)
       
               # Find the sunrise that follows the next sunset
               sunrise_time = keck.sun_rise_time(Time(sunset_time, format='jd'), which = 'next', horizon = -12*u.deg)
       
               # Divide the night into ~minutes, erring on the side of finer-than-minute sampling
               night_division = int(np.ceil((sunrise_time - sunset_time).jd*24*60))

               time_list = np.linspace(sunset_time.jd, sunrise_time.jd, night_division)
       
               # vis is an output of Ian's code. It is an array of tuples of the form (date, bool), where bool expresses the observability at the jd date.
               vis = obs.isObservable(time_list, observatory, star_target)
       
       
        ########## The truth counter starts at 0 and increments for each minute of observability. After 30 consecutive minutes, the target is considered observable for that night.
               truth_counter = 0
               for truth_value in vis:
                   if truth_value:
                       truth_counter += 1
                       if truth_counter == 30:
                   
                           # On the first night, initialize old_observability
                           if day_list.index(day) == 0:
                               old_observability = True
                       
                           # After the first night, assign the value to current_observability instead
                           else:
                               current_observability = True
                           break
           
                   elif not truth_value:
                       truth_counter = 0
                       if day_list.index(day) == 0:
                           old_observability = False
                       else:
                           current_observability = False
                   
               if day_list.index(day) == 0:
                   star_marker_list.append((sunset_time.jd, old_observability))


               elif current_observability != old_observability:
                   star_marker_list.append((sunset_time.jd, current_observability))
                   old_observability = current_observability
    
            # star_marker_list gives the 'landmark' dates for a target, where it transitions from being observable to unobservable, or vice versa
            star_marker_list.append((end_date, False))
           #######
            # print(star_marker_list)
            # Break down star_marker_list into dates and boolean values for plotting
            date_list, truth_list = np.transpose(np.vstack(star_marker_list))[0], np.transpose(np.vstack(star_marker_list))[1]

            for j in range(len(truth_list)):
                if truth_list[j]:
                    plt.hlines(y_length[i], date_list[j], date_list[j+1], colors = ['gray'], 
                               linewidths = 4, alpha = 0.5, zorder = z_order_list.index('observability'))


    cadenced_hires_dates = []
    cadenced_hires_y = []

    cadenced_apf_dates = []
    cadenced_apf_y = []

    for i in range(len(sql_df_yesiod)):
       index_of_planet_in_plot_df = pd.Index(plot_df['star_id']).get_loc(sql_df_yesiod['star_id'][i])


       if sql_df_yesiod['instrument'][i] == 'hires_j':
           cadenced_hires_dates.append(sql_df_yesiod['bjd'][i])
           cadenced_hires_y.append(y_length[index_of_planet_in_plot_df])

       elif sql_df_yesiod['instrument'][i] == 'apf':
           cadenced_apf_dates.append(sql_df_yesiod['bjd'][i])
           cadenced_apf_y.append(y_length[index_of_planet_in_plot_df])


    cadenced_hires_rvs = ax.scatter(cadenced_hires_dates, cadenced_hires_y - 0.15*(y_length[1]-y_length[0]), marker = 'v', 
                             color = color_hires, s = 10, zorder=z_order_list.index('rv_pts'))

    cadenced_apf_rvs = ax.scatter(cadenced_apf_dates, cadenced_apf_y - 0.15*(y_length[1]-y_length[0]), marker = 'v', 
                             color = color_apf, s = 10, zorder=z_order_list.index('rv_pts'))


    ############
    observing_schedule_df = pd.read_csv('../jump-config/allocations/hires_j/hires_schedule_2022A.csv')[['Date', 'start', 'stop']]

     # The schedule has multiple rows for some nights because the nights were paid for by 2 programs. This chunk combines the duplicate lines and their night fractions. It can NOT account for non-contiguous observing periods. For example, if CPS has the first 1/4 of the night, then we hand off for the second 1/4, then we get it back for the second half, this chunk will think we have the whole night.

    for i in observing_schedule_df.drop_duplicates(subset='Date')['Date']:
        if len(observing_schedule_df.query('Date == "{}"'.format(i))) > 1:
            index_list = observing_schedule_df.query('Date == "{}"'.format(i)).index
            date = observing_schedule_df['Date'][index_list[0]]
            start = observing_schedule_df['start'][index_list[0]]
            stop = observing_schedule_df['stop'][index_list[-1]]

            observing_schedule_df = observing_schedule_df.drop(index_list).reset_index(drop=True)

            observing_schedule_df = observing_schedule_df.append(pd.DataFrame({'Date':[date], 'start':[start], 'stop':[stop]})).reset_index(drop=True)

    observing_schedule_df = observing_schedule_df.sort_values(by='Date').reset_index(drop=True)

    # The dates in the schedule are given for Hawaii time at midnight that morning. If we start observing Jan 1 at 6 pm Hawaii time, then the JD is Jan 2 at 5 am. 6 pm is early, but we don't need to be too precise because we are going to find the next sunset time anyway.
    observing_dates = Time(observing_schedule_df['Date'].values.tolist(), format='iso').jd + 1 + 5/24


    ## Uses the fact that dates are in chronological order, so the min index corresponds to the earliest date
    # index_of_next_date = min([np.where(observing_dates == i) for i in observing_dates if i > Time.now().jd])[0][0]
    #
    #
    line_today, = ax.plot((Time.now().jd, Time.now().jd), (y_length[0], y_length[-1]), c = 'gray', linestyle = 'dashed')
    line_15, = ax.plot((Time.now().jd-15, Time.now().jd-15), (y_length[0], y_length[-1]), c = 'gray', linestyle = 'dotted')
    line_25, = ax.plot((Time.now().jd-25, Time.now().jd-25), (y_length[0], y_length[-1]), c = 'black', linestyle = 'dashdot')
    #
    # next_obs = [observing_dates[index_of_next_date+n] for n in range(1)]
    # # next_obs = [Time('2021-02-21', format = 'iso').jd]
    # for i in next_obs:
    #     line_next_obs,= ax.plot((i, i), (y_length[0], y_length[-1]), linewidth = 1, c = 'red')

    plt.xticks(date_intervals_jd, date_intervals_iso, fontsize = 8)


    plt.yticks([], [], size = 14)
    # ax.legend([recon_pts, jitter_pts, hires_template_pts, apf_template_pts,
    #            cadenced_hires_rvs, cadenced_apf_rvs, line_today, line_25, line_15],
    #            ['Recon', 'Jitter', 'HIRES Template', 'APF Template', 'HIRES',
    #             'APF', 'Today', '25', '15'], loc = (0.10,1.02), prop = {'size':10}, ncol = 4);
    ax.legend([ 
               cadenced_hires_rvs, cadenced_apf_rvs, line_today, line_25, line_15], 
               ['HIRES', 'APF', 'Today', '25', '15'], loc = (0.10,1.02), prop = {'size':10}, ncol = 4);
    fig.tight_layout()
    plt.savefig('csv/overview_plot.jpg', dpi = 500, pil_kwargs={'quality':95})
    
    if plot: 
        plt.show()
        
    
    return overview_df

def update_overview(overview_df):
    
    ov_cols = ['star_id', 'vmag', 'ra_deg', 'dec_deg', 'fit_pref', 'dvdt', 'u_dvdt', 'curv', 'u_curv', 'baseline_tot',
               'have_recon', 'have_jitter', 'have_template',
               'tot_iodine', 'tot_iodine_hires', 'tot_iodine_apf',
               'cooked?', 'last_obs_hires', 'last_obs_apf',
               'baseline_hires', 'baseline_apf', 'earliest_obs']
    
    overview_df[ov_cols].to_csv('csv/overview_df.csv', index = False)
    
    print('Updated overview_df.csv')
    
def update_trends(overview_df):
    
    out_file = open('../jump-config/programs/tks_distantgiants_trends.txt', 'w+')

    trend_count = 0 
    for i in range(len(overview_df)):
        
        if overview_df['3sig_trend'][i] == "TRUE" or overview_df['3sig_curv'][i] == "TRUE":
            
            star_name = overview_df['star_id'][i]
            out_file.write(star_name+'\n')
            trend_count += 1


    out_file.close()

    
    jump_config_path = r'../jump-config'
    commit_message = 'Updated tks_distantgiants_trends.txt'
    
    my_repo = repo.Repo(jump_config_path) # Path to jump_config repo
    my_repo.git.add(update=True) # Adds updates to existing files to index, rather than index.add, which I think works for                              adding *new* files
    my_repo.index.commit(commit_message) # Add commit message to the index (staging area between working dir and repo)
    origin = my_repo.remote(name='origin') # Specify where to push changes
    origin.pull()
    origin.push()
    print('Updated tks_distantgiants_trends.txt ({} targets)'.format(trend_count))   
    
    
    

if __name__ == "__main__":
    
    ov = make_overview(plot = False, observability = False)
    
    update_overview(ov)
    update_trends(ov)
    
    
    
    
    
    
    
    
    
    
    