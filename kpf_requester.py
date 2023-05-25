import numpy as np
import pandas as pd

import tks_target_list_gen.sim.hires.exposure as exp



def df_maker():
    
    overview_df = pd.read_csv('csv/overview_df.csv')
    
    ## Name
    star_name = overview_df.star_id
    
    ## RA columns
    RA_deg = np.array(overview_df.ra_deg)
    
    if RA_deg.max() < 25:
        print('WARNING: RA is probably provided in Hr instead of deg')
    
    ra_hr = (RA_deg/15).astype(int)
    ra_min = ((RA_deg/15%1)*60).astype(int)
    ra_sec = np.round((((RA_deg/15%1)*60)%1)*60, 1)
    
    
    ## Dec columns
    DEC_deg = np.array(overview_df.dec_deg)
    
    dec_deg = DEC_deg.astype(int)
    dec_min = (( abs(DEC_deg) - (abs(DEC_deg)).astype(int) )*60).astype(int)
    dec_sec = (np.round((((abs(DEC_deg) - (abs(DEC_deg)).astype(int) )*60)%1)*60, 1))
    
    ## Epoch
    epoch = [2000 for i in range(len(overview_df))]
    
    ## vmag= dummy column
    vmag_str = ['vmag=' for i in range(len(overview_df))]

    ## Vmag
    vmag = overview_df.vmag
    
    ## Counts (in thousands)
    counts = np.array([60 for i in range(len(overview_df))]) # 60k for all DG targets
    
    ## exp_time
    t_exp = (np.round(exp.exposure_time(vmag, counts, iod = True))).astype(int)
    
    ## max_exp_time
    t_max_list = np.array([500, 600, 900, 1200, 1500, 1800, 2700, 3600, 5400, 6000])
    t_max = [(2*t_exp[j] + min([i for i in t_max_list - 2*t_exp[j] if i>0])) for j in range(len(t_exp))]

    ## Decker
    decker = ['B5' if vmag[i]<= 10 else 'C2' for i in range(len(overview_df))]
    
    ## Desired obs
    cadence = 40
    
    # Choose the maximum of the following 2 options:
    #######################
    # Baseline-limited
    N_baseline = np.ceil((3*365 - overview_df.baseline_hires)/cadence) # Specified cadence until 3-year baseline reached
    # Observation-limited
    N_iod = 30 - overview_df.tot_iodine_hires
    #######################
    
    N_needed = np.amax([N_baseline, N_iod], axis=0)
    N_needed = np.maximum(N_needed, 0) # Cooked targets give negative obs_needed. Take 0 instead

    
    # Number of possible obs per semester depends on cadence
    N_sem = int(180/cadence)
    print('With the given cadence, we get {} obs per semester'.format(N_sem))
    N_obs = np.minimum(N_needed, N_sem) # Only ask for at most 7 obs per semester
    
    ## Number of observations acquired
    N_acq = [1 for i in range(len(overview_df))]
    
    
    ## Iodine
    iodine = ['in' for i in range(len(overview_df))]
    
    ## Priority
    prio = ['p1' for i in range(len(overview_df))]
    
    ## Program code
    prog_code = ['DG' for i in range(len(overview_df))]
    
    ## Telescope proposal code
    tel_code = ['U062' for i in range(len(overview_df))]
    
    
    # print(star_name, ra_hr, ra_min, ra_sec)
    # sfd
    
    kpf_df = pd.DataFrame({'N_obs':N_obs, 
                          'Starname':star_name, 
                          'RAH':ra_hr,
                          'RAM':ra_min,
                          'RAS':ra_sec,
                          'DECD':dec_deg,
                          'DECM':dec_min,
                          'DECS':dec_sec,
                          'epoch':epoch,
                          'vmag_str':vmag_str,
                          'vmag':vmag,
                          'T_exp':t_exp,
                          'T_max':t_max,
                          'counts':counts,
                          'decker':decker,
                          'N_acq':N_acq,
                          'iodine':iodine,
                          'prio':prio,
                          'prog_code':prog_code,
                          'tel_code':tel_code})
    # print(kpf_df.head(5))
    
    kpf_df.to_csv('kpf_df.csv')


if __name__ == "__main__":
    df_maker()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    