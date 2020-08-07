"""
Judah Van Zandt
Python 3

This module is used to compile general information about the state of each target in the SC2A program.

Written: June 1, 2020
Last modified: July 7, 2020
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.time import Time
from astroplan import Observer, FixedTarget
import astropy.units as u
from astropy.coordinates import SkyCoord

import ephem
import observing as obs

import distantgiants




def make_overview(plot = True, observability = False):
    """
    With a lot of intermediate steps, this function generates overview_df and an overview plot. They both contain important
    information about each target in SC2A.
    """
    # Need jump_df for template information
    jump_df = pd.read_csv('csv/candidates.csv')
    
    jump_df = jump_df.drop_duplicates(subset='Name')
    jump_df['rp'].replace(np.nan, -100, inplace = True)
    jump_df = jump_df[['Name', 'rp', 'have_template_hires_j', 'have_template_apf', 'mass', 'radius', 'svalue', 'ruwe', 'per']]
    # jump_df = jump_df.drop(index = jump_df.query("Name == 'T001244'").index.values) # Remove T001244, which has vmag = 11.9

    # Creating a new column to tell whether a target has a template on either HIRES or APF
    jump_df['have_template_hires_j'].replace(np.nan, 0, inplace = True)
    jump_df['have_template_apf'].replace(np.nan, 0, inplace = True)
    jump_df['have_template'] = list(map(lambda x, y: 'NO' if x+y < 1 else 'YES', jump_df['have_template_hires_j'], jump_df['have_template_apf']))
    
    
    distantgiants_spec = pd.read_csv('csv/distantgiants_spec.csv')
    distantgiants_spec = pd.merge(distantgiants_spec, jump_df, left_on = 'star_id', right_on = 'Name')
    
    sql_df = pd.read_csv('csv/Distant_Giants_Observing_Requests.csv')

    sql_df = pd.merge(distantgiants_spec[['star_id', 'vmag', 'ra', 'dec']], sql_df, on = 'star_id', how = 'inner')


    full_star_list = pd.merge(sql_df.drop_duplicates(subset = 'star_id'), \
                              distantgiants_spec.drop_duplicates(subset = 'star_id'), \
                              on = 'star_id', how = 'outer')['star_id'].to_frame()
    
    # Stars that need recon
    sql_df['day'] = Time(sql_df['bjd'], format='jd', scale='utc', out_subfmt='date').iso
    sql_df_yesiod = sql_df.query("iodine_in == 't'").reset_index(drop = True)
    sql_df_noiod = sql_df.query("iodine_in == 'f'").reset_index(drop = True)

    recon_df = pd.merge(sql_df_noiod.query('counts >= 4500')[['star_id', 'instrument', 'bjd']].\
                        sort_values(['star_id', 'bjd']).drop_duplicates\
                        (subset = 'star_id', keep = 'first'), full_star_list['star_id'], \
                        on = 'star_id', how = 'outer').rename(columns = {'bjd':'bjd_recon'}).\
                        rename(columns = {'instrument':'instrument_recon'})

    recon_df['have_recon'] = list(map(lambda x: 0 if pd.isna(x) else 1, recon_df['bjd_recon']))
         
    # Stars that need jitter tests

    grouped_list = sql_df_yesiod.groupby(by=['star_id', 'day', 'instrument'])
    has_jitter = grouped_list.size().to_frame().reset_index().rename(columns = {0:'count'}).query('count >= 3').sort_values(['star_id', 'day']).drop_duplicates(subset = 'star_id', keep = 'first')
    jitter_df = pd.merge(has_jitter, full_star_list, on = 'star_id', how = 'outer')
    jitter_df['have_jitter'] = list(map(lambda x: 0 if pd.isna(x) else 1, jitter_df['count']))
    jitter_df['bjd_jitter'] = list(map(lambda x: x if pd.isna(x) else Time(x, format = 'iso', out_subfmt = 'date').jd, jitter_df['day']))



    # dates_df gives every observation and its jd, sorted so that the most and least recent can be picked for each target
    # latest_obs finds the most recent observation with both HIRES and APF for each target
    dates_df = pd.merge(sql_df_yesiod[['star_id', 'instrument', 'bjd']].sort_values(by = ['star_id', 'instrument', 'bjd']).reset_index(drop = True),\
                        full_star_list, on = 'star_id', how = 'outer')

    latest_obs = dates_df.drop_duplicates(subset = ['star_id', 'instrument'], keep = 'last')
    earliest_obs = dates_df.drop_duplicates(subset = ['star_id', 'instrument'], keep = 'first')

    # Adding columns for the most recent HIRES and APF iodine-in obs to overview_df. Displayed as days since last obs.
    last_obs_hires = pd.DataFrame({'star_id': latest_obs.query("instrument == 'hires_j'")['star_id'], \
                                   'last_obs_hires': -latest_obs.query("instrument == 'hires_j'")['bjd'].\
                                   subtract(Time(Time.now().iso, out_subfmt='date').jd)})

    last_obs_apf = pd.DataFrame({'star_id': latest_obs.query("instrument == 'apf'")['star_id'], \
                                 'last_obs_apf': -latest_obs.query("instrument == 'apf'")['bjd'].\
                                 subtract(Time(Time.now().iso, out_subfmt='date').jd)})



    # Creating overview_df to show general information
    overview_df = recon_df[['star_id', 'have_recon']]
    overview_df = pd.merge(overview_df, jitter_df[['star_id', 'have_jitter']], how = 'outer', on = 'star_id')
    overview_df = pd.merge(overview_df, distantgiants_spec[['star_id', 'have_template']], how = 'outer', on = 'star_id') # Stars that need templates 
    overview_df = pd.merge(overview_df, 
                           pd.DataFrame(sql_df_yesiod.groupby(by=['star_id', 'instrument']).size()).reset_index().\
                           rename(columns = {0: 'tot_iodine_hires'}).query("instrument == 'hires_j'").drop(columns = {'instrument'}),
                           how = 'outer', on = 'star_id')

    overview_df = pd.merge(overview_df, 
                           pd.DataFrame(sql_df_yesiod.groupby(by=['star_id', 'instrument']).size()).reset_index().\
                           rename(columns = {0: 'tot_iodine_apf'}).query("instrument == 'apf'").drop(columns = {'instrument'}),
                           how = 'outer', on = 'star_id')

    overview_df = pd.merge(overview_df, last_obs_hires, on = 'star_id', how = 'outer')

    overview_df = pd.merge(overview_df, last_obs_apf, on = 'star_id', how = 'outer')

    # Creating HIRES and APF baselines and adding them to overview_df
    baseline_hires = pd.DataFrame({'star_id': latest_obs.query("instrument == 'hires_j'")['star_id'].reset_index(drop = True),\
                               'baseline_hires': latest_obs.query("instrument == 'hires_j'")['bjd'].reset_index(drop = True) - \
                               earliest_obs.query("instrument == 'hires_j'")['bjd'].reset_index(drop = True)})

    baseline_apf = pd.DataFrame({'star_id': latest_obs.query("instrument == 'apf'")['star_id'].reset_index(drop = True),\
                               'baseline_apf': latest_obs.query("instrument == 'apf'")['bjd'].reset_index(drop = True) - \
                               earliest_obs.query("instrument == 'apf'")['bjd'].reset_index(drop = True)})

    overview_df = pd.merge(overview_df, baseline_hires, on = 'star_id', how = 'outer')

    overview_df = pd.merge(overview_df, baseline_apf, on = 'star_id', how = 'outer')
    
    # Cleaning up overview_df
    overview_df.replace({'have_recon':{0:'NO', 1:'YES'},
                        'have_jitter':{0:'NO', 1:'YES'}}, 
                        inplace=True)

    overview_df.update(overview_df[['tot_iodine_hires', 'tot_iodine_apf']].fillna(0))
    overview_df.update(overview_df[['last_obs_hires', 'last_obs_apf']].fillna('NEVER'))
    overview_df.update(overview_df[['baseline_hires', 'baseline_apf']].fillna('N/A'))
    
    # Add a new column 'cooked?' to indicate any target with both a 3+ year baseline and 40+ observations. 
    cooked = []
    for i in range(len(overview_df)):
        if type(overview_df['baseline_hires'][i]) == str:
            baseline_hires = 0
        elif type(overview_df['baseline_hires'][i]) == float:
            baseline_hires = float(overview_df['baseline_hires'][i])
        else:
            print(type(overview_df['baseline_hires'][i]))
        
        if type(overview_df['baseline_apf'][i]) == str:
            baseline_apf = 0
        elif type(overview_df['baseline_apf'][i]) == float:
            baseline_apf = float(overview_df['baseline_apf'][i])
        else:
            print(type(overview_df['baseline_apf'][i]))
        
        if (baseline_hires > 1095 or baseline_apf >= 1095) and\
        (overview_df['tot_iodine_hires'][i] + overview_df['tot_iodine_apf'][i] >= 40):
            cooked.append('COOKED')
        else:
            cooked.append('cookin')
    overview_df['cooked?'] = cooked
    
    # I've commented this line out to keep from cutting on membership in distantgiants
    overview_df = pd.merge(overview_df, sql_df.drop_duplicates(subset = 'star_id')[['star_id', 'vmag', 'ra', 'dec']], on = 'star_id')
    overview_df = overview_df.rename(columns = {'ra':'ra_deg', 'dec':'dec_deg'})

    if plot == True:
        # Creating plot_df with all of the information to create an image with an overview for each target
        plot_df = pd.merge(overview_df, recon_df.drop(columns = 'have_recon'), on = 'star_id')
        plot_df = pd.merge(plot_df, jitter_df.drop(columns = 'have_jitter'), on = 'star_id')[['star_id', 'instrument_recon', 'bjd_recon', 'bjd_jitter', 'have_template','cooked?', 'ra_deg', 'dec_deg']]

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
        date_intervals_iso = Time(date_intervals_jd, format = 'jd', out_subfmt = 'date').iso

        recon_color = 'green'
        template_color = 'black'
        jitter_color = 'pink'
        rv_color_hires = 'red'
        rv_color_apf = 'blue'

        ######################
        # Making the plot
        ######################

        fig, ax = plt.subplots(figsize=(10, 8), dpi= 100, facecolor='w', edgecolor='k')

        z_order_list = z_order_list = ['observability', 'h_line_plot', 'rv_pts', 'recon_pts', 'jitter_pts', 'recon_arrows', 'jitter_arrows']

        plt.hlines(y_length, start_date, end_date, colors = ['black', 'gray'], linewidths = 0.2, zorder = z_order_list.index('h_line_plot'))
        plt.grid(b = True, which = 'major', axis = 'x', linewidth = .2)

        ax.set_xlim(start_date, end_date)

        plot_df = plot_df.sort_values(by = 'ra_deg', ignore_index = True)

        # Plotting the dates of each star's recon and jitter
        recon_pts = ax.scatter(plot_df['bjd_recon'], y_length, color = recon_color, s=6, zorder = z_order_list.index('recon_pts'))
        jitter_pts = ax.scatter(plot_df['bjd_jitter'], y_length, color = jitter_color, s=6, zorder = z_order_list.index('jitter_pts'))


        ax.text(end_date+10, y_length[0]+1, 'Template?', size = 8)
        ax.text(start_date-15, y_length[0]+1, 'Jitter?', size = 8)
        cooked = [] # List determining which stars already have significant obs and shouldn't be observed
        
        for i in range(len(plot_df)):
            
            if plot_df['have_template'][i] == 'YES':
                template_facecolors = 'k'
            elif plot_df['have_template'][i] == 'NO':
                template_facecolors = 'none'

            if str(plot_df['bjd_jitter'][i]) != 'nan':
                jitter_facecolors = 'b'
            else: 
                jitter_facecolors = 'none'


            # Open or closed circle to show template status
            ax.scatter(end_date+14.5, y_length[i], s=10, marker = 'o', color = 'k', facecolors = template_facecolors, clip_on=False)

            # Open or closed circle to show jitter status
            ax.scatter(start_date-12.5, y_length[i], s=10, marker = 'o', color = 'k', facecolors = jitter_facecolors, clip_on=False)

            # Put star name on right and left vertical axes, with cooked targets faded
            if plot_df['cooked?'][i] == 'cookin':
                alpha = 1
            elif plot_df['cooked?'][i] == 'COOKED':
                alpha = 0.5

            # Plotting the name of each star on the y-axis
            ax.text(start_date-1, y_length[i]+0.3*y_increment, plot_df['star_id'][i], size = 7, alpha = alpha, horizontalalignment = 'right')
            ax.text(end_date+1, y_length[i]+0.3*y_increment, plot_df['star_id'][i], size = 7, alpha = alpha)

            # If a target got recon or jitter before the earliest displayed date, use an arrow
            if plot_df['bjd_recon'][i] < start_date:
                ax.scatter(start_date+0.5, y_length[i] - 0.25*y_increment, marker = '<', color = recon_color, s=10, zorder = z_order_list.index('recon_arrows'))
    
            if plot_df['bjd_jitter'][i] < start_date:
                ax.scatter(start_date+0.5, y_length[i] + 0.25*y_increment, marker = '<', color = jitter_color, s=10, zorder = z_order_list.index('recon_arrows'))
            
            ### Observability ###
            # Plotting observability bars takes a while; pick how many stars you want bars for here.
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
                   
                   # vis is an output of Ian's code. It is an array of tuple of the form (date, bool), where bool expresses the observability at the jd date.
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
                print(star_marker_list)
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
                                 color = 'red', s = 10, zorder=z_order_list.index('rv_pts'))

        cadenced_apf_rvs = ax.scatter(cadenced_apf_dates, cadenced_apf_y - 0.15*(y_length[1]-y_length[0]), marker = 'v', 
                                 color = 'blue', s = 10, zorder=z_order_list.index('rv_pts'))
        
        
        ############
        observing_schedule_df = pd.read_csv('../jump-config/allocations/hires_j/hires_schedule_2020B.csv')[['Date', 'start', 'stop']]
        
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
    
    
        # Uses the fact that dates are in chronological order, so the min index corresponds to the earliest date
        index_of_next_date = min([np.where(observing_dates == i) for i in observing_dates if i > Time.now().jd])[0][0]
        
        
        next_obs = [observing_dates[index_of_next_date+n] for n in range(10)]
        # next_obs = [Time('2020-08-04', format = 'iso').jd]
        ############


        line_today, = ax.plot((Time.now().jd, Time.now().jd), (y_length[0], y_length[-1]), c = 'gray', linestyle = 'dashed')
        line_15, = ax.plot((Time.now().jd-15, Time.now().jd-15), (y_length[0], y_length[-1]), c = 'gray', linestyle = 'dotted')
        line_25, = ax.plot((Time.now().jd-25, Time.now().jd-25), (y_length[0], y_length[-1]), c = 'black', linestyle = 'dashdot')
        for i in next_obs:
            line_next_obs, = ax.plot((i, i), (y_length[0], y_length[-1]), linewidth = 1, c = 'red')

        plt.xticks(date_intervals_jd, date_intervals_iso, fontsize = 8)


        plt.yticks([], [], size = 14)
        ax.legend([recon_pts, jitter_pts, cadenced_hires_rvs, cadenced_apf_rvs, line_today, line_25, line_15, line_next_obs], ['Recon', 'Jitter', 'HIRES', 'APF', 'Today', '25', '15', 'Future Nights'], loc = (0.24,1.02), prop = {'size':10}, ncol = 4);
        fig.tight_layout()
        plt.savefig('csv/overview_plot.jpg', dpi = 500, quality=95)
        plt.show()
        
    
    return overview_df

def update_overview(overview_df):
    
    overview_df.to_csv('csv/overview_df.csv', index = False)
    
    print('Updated overview_df.csv')

if __name__ == "__main__":
    
   update_overview(make_overview(plot = True, observability = False))






