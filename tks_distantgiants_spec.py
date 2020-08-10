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

def make_distantgiants_spec():
    """
    Creates spec_cuts, the dataframe of all SC2A stars on Jump that pass spectroscopic cuts. Then writes spec_cuts out to
    the jump-config Github repo.
    """
    
    # manual_cuts_2 is essentially distantgiants_photo with the added qlp, MES, and close companion cuts
    manual_cuts_2 = pd.read_csv('/Users/judahvz/research/code/GitHub/distantgiants/csv/manual_cuts_2.csv')
    print('{} targets from distantgiants_photo pass manual cuts'.format(len(manual_cuts_2)))
    

    # Spec
    spec_values = pd.read_csv('/Users/judahvz/research/code/GitHub/distantgiants/csv/TKS_-_Distant_Giants_Spectroscopic_Properties.csv')
    
    # Take the logrhk and vsini/teff values that correspond to the highest-SNR observation
    logrhk_df = spec_values.dropna(subset=['logrhk']).sort_values(by=['star_id', 'counts']).drop_duplicates(subset='star_id', keep='last')[['star_id', 'logrhk']]
    vsini_teff_df = spec_values.dropna(subset=['teff', 'vsini']).sort_values(by=['star_id', 'counts']).drop_duplicates(subset='star_id', keep = 'last')[['star_id', 'vsini', 'teff']]
    
    # print('{} targets have APF/HIRES observations, {} have logrhk values, and {} have vsini/teff values'.format(len(spec_values.drop_duplicates(subset='star_id')), len(logrhk_df), len(vsini_teff_df)))
    

    distantgiants_spec = pd.merge(manual_cuts_2.rename(columns = {'Vmag':'vmag'}), logrhk_df, left_on = 'cps', right_on = 'star_id', how='left').drop(columns='star_id').rename(columns={'cps':'star_id'})
    distantgiants_spec = pd.merge(distantgiants_spec, vsini_teff_df, on = 'star_id', how='left')
    
    # I have removed these for now to have a robust, 'for sure' list. If we want to add in targets that don't have spec properties later, I can use it
    # distantgiants_spec['vsini'].replace(np.nan, -100, inplace = True)
 #    distantgiants_spec['logrhk'].replace(np.nan, -100, inplace = True)
 #    distantgiants_spec['teff'].replace(np.nan, -100, inplace = True)
    
    distantgiants_spec = distantgiants_spec.query('teff < 6250 and vsini <= 5.00 and logrhk < -4.7').reset_index(drop = True)
    
    
    distantgiants_spec.sort_values(by='toi').to_csv('/Users/judahvz/research/code/GitHub/distantgiants/csv/distantgiants_spec.csv', index = False)
   
  
    return distantgiants_spec
    
    
def update_distantgiants_spec(distantgiants_spec):
    out_file = open('../jump-config/programs/tks_distantgiants_spec.txt', 'w+')

    for star_name in distantgiants_spec.sort_values(by = 'star_id').star_id.values:
        out_file.write(star_name+'\n')


    out_file.close()

    
    jump_config_path = r'../jump-config'
    commit_message = 'Updated distantgiants_spec.txt'
    
    my_repo = repo.Repo(jump_config_path) # Path to jump_config repo
    my_repo.git.add(update=True) # Adds updates to existing files to index, rather than index.add, which I think works for                              adding *new* files
    my_repo.index.commit(commit_message) # Add commit message to the index (staging area between working dir and repo)
    origin = my_repo.remote(name='origin') # Specify where to push changes
    origin.pull()
    origin.push()
    print('Updated distantgiants_spec.txt ({} targets)'.format(len(distantgiants_spec)))


if __name__ == "__main__":
    
    update_distantgiants_spec(make_distantgiants_spec())
    
    
    
    
