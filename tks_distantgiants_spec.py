"""
Judah Van Zandt
Python 3

Module to generate the tks_distgiants target list, which comprises all stars from tks_distgiants_pool that pass
spectroscopic vetting. These are then written out to the jump-config/programs folder on Github.

Written: June 1, 2020
Last modified: July 8, 2020
"""


import numpy as np
import pandas as pd
import sys
from astropy.time import Time
import matplotlib.pyplot as plt
from git import repo

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def make_distantgiants_spec(distantgiants_photo):
    """
    Creates spec_cuts, the dataframe of all SC2A stars on Jump that pass spectroscopic cuts. Then writes spec_cuts out to
    the jump-config Github repo.
    """
    
    # manual_cuts_1 is the manually-compiled set of targets that pass MES>12 and NOT qlp-only. manual_cuts_2 is the result of enforcing another cut on close companions using Gaia DR2 data
    manual_cuts_2 = pd.read_csv('csv/manual_cuts_2.csv')
    print('{} targets from distantgiants_photo pass manual cuts'.format(len(manual_cuts_2)))
    
    ###################
    ## Reactivate manual cuts after All-Hands Meeting
    distantgiants_photo = pd.merge(distantgiants_photo, manual_cuts_2['cps'], on = 'cps')
    ###################
    # Spec
    spec_values = pd.read_csv('csv/TKS_-_Distant_Giants_Spectroscopic_Properties.csv')
    
    # Take the logrhk and vsini/teff values that correspond to the highest-SNR observation
    logrhk_sval_df = spec_values.dropna(subset=['logrhk', 'sval'])\
                                .sort_values(by=['star_id', 'counts'])\
                                .drop_duplicates(subset='star_id', keep='last')\
                                [['star_id', 'observation_id', 'logrhk', 'sval']]
                                
    vsini_teff_df = spec_values.dropna(subset=['teff', 'vsini'])\
                               .sort_values(by=['star_id', 'counts'])\
                               .drop_duplicates(subset='star_id', keep = 'last')[['star_id', 'vsini', 'teff']]
                               
    
    distantgiants_spec = pd.merge(distantgiants_photo.rename(columns = {'Vmag':'vmag'}),\
                                  logrhk_sval_df, left_on = 'cps', right_on = 'star_id', how='left')\
                                      .drop(columns='star_id').rename(columns={'cps':'star_id'})
    distantgiants_spec = pd.merge(distantgiants_spec, vsini_teff_df, on = 'star_id', how='left')

    ## Removing cuts on logr'hk to avoid losing targets.
    # distantgiants_spec['vsini'].replace(np.nan, -100, inplace = True)
    # distantgiants_spec['logrhk'].replace(np.nan, -100, inplace = True)
    # distantgiants_spec['teff'].replace(np.nan, -100, inplace = True)
    # Remove cuts on spec parameters. We set the 47-target sample below, and some targets fluctuate above -4.7 logrhk.
    #distantgiants_spec = distantgiants_spec.query('teff < 6250 and vsini <= 5.00 and logrhk < -4.7').reset_index(drop = True)
    
    # # samp47 is a list of the 44 stars assigned to SC2A by Ashley's code, one target (T002088) assigned to another program, and 2 cooked targets not assigned to SC2A (219134, 75732).
    # Because samp47 has the final say on spectroscopic cuts, this module is basically obsolete. I keep it for historical purposes.
    # Some targets, (T001235, HIP70705, and T001772) were not selected by Ashley's code despite passing all cuts.
    
    ## These lines are unrelated to spec. They just update vmag, ra, and dec with values from Jump.
    ## The original values are from tois_perfect and are inaccurate in some cases.
    dg_updates = pd.read_csv('csv/Distant_Giants_Update_Params.csv')
    distantgiants_spec = distantgiants_spec.merge(dg_updates, on='star_id')
    distantgiants_spec['vmag'] = distantgiants_spec['vmag_new'].fillna(distantgiants_spec['vmag'])
    distantgiants_spec['ra'] = distantgiants_spec['ra_new'].fillna(distantgiants_spec['ra'])
    distantgiants_spec['dec'] = distantgiants_spec['dec_new'].fillna(distantgiants_spec['dec'])
    distantgiants_spec = distantgiants_spec.drop(columns=['ra_new', 'dec_new', 'vmag_new'])

    samp47 = pd.read_csv('csv/samp47.csv')
    distantgiants_spec = distantgiants_spec.merge(samp47, on = 'star_id', how='inner')
    print('47-target sample enforced in tks_distantgiants_spec.py')
    # print(distantgiants_spec.sort_values(by='star_id', ignore_index=True).star_id)
    # print(samp47.sort_values(by='star_id', ignore_index=True).star_id)
    # dfd
    
    distantgiants_spec.sort_values(by='toi').to_csv('csv/distantgiants_spec.csv', index = False)

  
    return distantgiants_spec
    
    
def update_distantgiants_spec(distantgiants_spec):
    out_file = open('../jump-config/programs/tks_distantgiants_spec.txt', 'w+')

    for star_name in distantgiants_spec.sort_values(by = 'star_id').star_id.values:
        out_file.write(star_name+'\n')


    out_file.close()

    
    jump_config_path = r'../jump-config'
    commit_message = 'Updated distantgiants_spec.txt'
    
    my_repo = repo.Repo(jump_config_path) # Path to jump_config repo
    my_repo.git.add(update=True) # Adds updates to existing files to index, rather than index.add, which I think works for adding *new* files
    my_repo.index.commit(commit_message) # Add commit message to the index (staging area between working dir and repo)
    origin = my_repo.remote(name='origin') # Specify where to push changes
    origin.pull()
    origin.push()
    print('Updated distantgiants_spec.txt ({} targets)'.format(len(distantgiants_spec)))


if __name__ == "__main__":
    
    distantgiants_photo = pd.read_csv('csv/distantgiants_photo.csv')

    update_distantgiants_spec(make_distantgiants_spec(distantgiants_photo))
    
    # distantgiants_photo = pd.read_csv('/Users/judahvz/research/code/all_hands_meetings/january_2021/dg_photo_hypothetical.csv')
    # print(len(distantgiants_photo))
    # make_distantgiants_spec(distantgiants_photo)
    
    
