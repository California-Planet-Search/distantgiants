"""
Judah Van Zandt
Python 3

This module generates a list of observing requests based on the current status of each star in the SC2A program.

Written: June 1, 2020
Last modified: August 27, 2020
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tks_target_list_gen.sim.hires.exposure as exp
import astropy.units as u
from astropy.time import Time

import warnings
from astropy.utils.exceptions import AstropyWarning
# warnings.simplefilter('ignore', category=AstropyWarning)

from cpsutils.obsplanning import Star
from cpsutils.obsplanning import times

from overview import make_overview
from tks_distantgiants import make_distantgiants

from astroplan import download_IERS_A


def init_overview(iers = False):
    
    if iers == True:
        download_IERS_A()
    
    overview_df = make_overview(plot=False) 
    
    return overview_df


def obs_request_list_gen(overview_df):
    """
    This function uses the current status of each target in overview_df (jitter, template, recon status as well as 
    time since last RV obs) to create a list of stars that need each type of observation. These are  fed to the generator
    function to be turned into actual script lines.
    """
    overview_df = overview_df.sort_values(by = 'ra_deg').reset_index(drop = True)
    
    observing_schedule_df = pd.read_csv('../jump-config/allocations/hires_j/hires_schedule_2020B.csv')[['Date', 'start', 'stop']].sort_values(by=['Date', 'start'])
    
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
    
    
    # Uses the fact that the dates are in chronological order, so the min index corresponds to the earliest date
    
    index_of_next_date = min([np.where(observing_dates == i) for i in observing_dates if i > Time.now().jd])[0][0]
    
    next_observing_date = observing_dates[index_of_next_date]
    # next_observing_date = Time('2020-10-30', format='iso').jd
    
    
    time_gap = next_observing_date - Time.now().jd
    
    
    start = observing_schedule_df['start'][index_of_next_date]
    stop = observing_schedule_df['stop'][index_of_next_date]
    # start, stop = 0, 1
    
    request_list = [[], [], [], []]
    obs_list = ['have_recon', 'have_jitter', 'have_template']
    obs_type = ['recon', 'jitter', 'template', 'rv']
    obs_prio = [2, 3, 3]
    
    for i in range(len(overview_df)):
        
        if overview_df['cooked?'][i] == 'COOKED':
            continue
        
        star_name = overview_df['star_id'][i]
        vmag = overview_df['vmag'][i]
        RA_deg = overview_df['ra_deg'][i]
        Dec_deg = overview_df['dec_deg'][i]
        
        if Dec_deg >= 0:
                sign = '+'
        else:
            sign = '-'
        
        RA_new = str(int(RA_deg/15)).zfill(2) +'h'+ str(int((RA_deg/15%1)*60)).zfill(2) +'m'+ str(np.round((((RA_deg/15%1)*60)%1)*60, 1)).zfill(4)+'s'
        Dec_new = sign + str(int(abs(Dec_deg))).zfill(2) +'d'+ str(int((abs(Dec_deg) - int(abs(Dec_deg)))*60)).zfill(2) +'m'+ str(int(np.round((((abs(Dec_deg) - int(abs(Dec_deg)))*60)%1)*60, 0))).zfill(2)+'s'
       
        star_object = Star.star(star_name, RA = RA_new, Dec = Dec_new)
        
        observer_times = times.ObserverTimes(utc_date = Time(next_observing_date, format='jd').iso.split(' ')[0], night_kind = [start, stop])
        
        visibility = star_object.visibility(observer_times, verbose = False)
        
        visible_time = max([i[2] for i in star_object.visibility(observer_times, verbose = False)])
        if star_name == 'T001691':
            print(visible_time)
        if visible_time < 0.5*u.h:
            continue

        for j in range(len(obs_type)):
            # Creating requests for recon, jitter, and templates
            if j <= 2:
                if overview_df[obs_list[j]][i] == 'NO':
                    # Jitter test stars must pass the extra criterion of being visible for ~3 hours. Skipping 191939 at Howard's suggestion.
                    if (obs_list[j] == 'have_jitter' and visible_time < 3*u.h) or star_name == '191939' :
                        continue
                    # Templates will only be requested if the star already has 3 or more iodine-in RVs
                    if (obs_list[j]) == 'have_template_hires_j' and overview_df['tot_iodine_hires'][i] + overview_df['tot_iodine_apf'][i] < 3:
                        continue
                        
                    request_list[j].append((overview_df['star_id'][i], obs_type[j], obs_prio[j], vmag, RA_deg, Dec_deg))
           
            # Creating requests for cadence RVs  
            elif obs_type[j] == 'rv':
                if overview_df['last_obs_hires'][i] == 'NEVER':
                    hires_never = True
                else:
                    hires_never = False
                    hires_days = float(overview_df['last_obs_hires'][i]) + time_gap # Computes days that will have passed since last obs on the actual day of observation
    
                if overview_df['last_obs_apf'][i] == 'NEVER':
                    apf_never = True
                else:
                    apf_never = False
                    apf_days = float(overview_df['last_obs_apf'][i])
                if (hires_never or hires_days > 25)\
                and (apf_never or apf_days > 25):
                    prio = 1
                elif (hires_never or hires_days > 20)\
                and (apf_never or apf_days > 20):
                    prio = 0
                elif (hires_never or hires_days > 15)\
                and (apf_never or apf_days > 15):
                    prio = 0
                # elif (hires_never or hires_days > 0)\
 #                and (apf_never or apf_days > 0):
 #                    prio = 4
                else:
                    prio = 0
              
                
                request_list[j].append((overview_df['star_id'][i], 'rv', prio, vmag, RA_deg, Dec_deg))
        
    return request_list


# print(obs_request_list_gen(init_overview(iers=False)))
    
    
def generator(star_requests):
    """
    The list of star_requests is expected as a list of lists: [[recon_requests], [jitter_requests], ...],
    where for example, [recon_requests] looks like [(star_name_1, 'recon'), (star_name_2, 'recon'), ...].
    star_name is the CPS ID of the star as it appears in the all_TOIs spreadsheet, and the obs_type is
    an element of the list ['recon', 'jitter', 'template', 'rv'].
    """
    
    date_list = Time.now().iso.split(' ')[0]
    
    out_file = open('observing_requests/'+date_list+'.txt', 'w+')
    
    list_of_line_lists = []
    
    total_time = 0
    
    for obs_type_list in star_requests:
        line_list = []
        for j in obs_type_list:
            star_name = j[0]
            obs_type = j[1]
            prio = j[2]
            v_mag = j[3]
            RA_deg = j[4]
            Dec_deg = j[5]

            if prio == 0:
                continue
            


            RA_new = str(int(RA_deg/15)).zfill(2) +' '+ str(int((RA_deg/15%1)*60)).zfill(2) +' '+ str(np.round((((RA_deg/15%1)*60)%1)*60, 1)).zfill(4)

            if Dec_deg >= 0:
                sign = '+'
            else:
                sign = '-'

            Dec_new = sign + str(int(abs(Dec_deg))).zfill(2) +' '+ str(int((abs(Dec_deg) - int(abs(Dec_deg)))*60)).zfill(2) +' '+ str(int(np.round((((abs(Dec_deg) - int(abs(Dec_deg)))*60)%1)*60, 0))).zfill(2)
            
            
            if obs_type == 'recon':
                iodine = 'out'
                iod_status = False
                
                if v_mag <= 10:
                    decker = 'B1'
                elif v_mag > 10:
                    decker = 'B3'
                
                counts = 10
                n_shots = '1x'
                initials = 'DG'
                string = 'Recon for TKS Distant Giants'

            elif obs_type == 'jitter':
                iodine = 'in'
                iod_status = True
                
                if v_mag <= 10:
                    decker = 'B5'
                elif v_mag > 10:
                    decker = 'C2'
                
                counts = 60
                n_shots = '1x'
                initials = 'DG'
                string = '** Jitter test'
                # v_mag = 0
            

            elif obs_type == 'template':
                iodine = 'out'
                iod_status = False
                
                if v_mag <= 10:
                    decker = 'B1'
                elif v_mag > 10:
                    decker = 'B3'
                # v_mag = 0
                
                counts = 125
                n_shots = '1x'
                initials = 'DG'
                string = '** template please add B-stars'
                

            elif obs_type == 'rv':
                iodine = 'in'
                iod_status = True
                
                if v_mag <= 10:
                    decker = 'B5'
                elif v_mag > 10:
                    decker = 'C2'
                # v_mag = 0
                    
                counts = 60
                n_shots = '1x'
                initials = 'DG'
                string = 'Cadence RV for TKS Distant Giants'
                
            t_exp = int(np.round(exp.exposure_time(v_mag, counts, iod = iod_status)))
            if obs_type == 'jitter':
                total_time += 3*t_exp
            else:
                total_time += t_exp

            t_max_list = np.array([500, 600, 900, 1200, 1500, 1800, 2700, 3600, 5400, 6000])
            t_max = str(2*t_exp + min([i for i in t_max_list - 2*t_exp if i>0]))

            line = star_name.rjust(15)+' '+RA_new+' '+Dec_new+' '+'2000'+' '+('vmag='+str(np.round(v_mag, 1))).rjust(9)+' '+\
                    str(t_exp).rjust(4)+'/'+t_max.rjust(4)+str(counts).rjust(4)+'k'+' '+decker+' '+n_shots+' '+\
                    iodine.rjust(3)+' '+'p'+str(prio)+' '+initials+' '+string+'\n'
            line_list.append(line)
        list_of_line_lists.append(line_list) 
    
    time_string = 'Total requested time = {} seconds = {:.0f} hr {:.0f} min {:.1f} sec'.format(total_time, np.floor(total_time/3600), np.floor((total_time/3600%1)*60), total_time/3600%1*60%1*60)+'\n \n'
    out_file.write(time_string)
        
    for line_list in list_of_line_lists:
        for line in line_list:
            out_file.write(line)
    print('Did you remember to update jump_df and sql_df?')
    print('New observation request list generated')
    

if __name__ == '__main__':
    
    request_list = obs_request_list_gen(init_overview(iers=False))
    
    generator(request_list)
    
    
    
    
    
    
    