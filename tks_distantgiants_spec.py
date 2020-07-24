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
    the jump-config Github repo. SC2A stars are defined by the output of tks_distgiants_pool.py.
    Then pushes the changes made to tks_distgiants.txt to the jump-config Github.
    """
    
    jump_df = pd.read_csv('csv/candidates.csv')
    distantgiants_photo = pd.read_csv('csv/distantgiants_photo.csv')
    
    jump_df = jump_df.drop_duplicates(subset='Name')
    jump_df['vsini'].replace(np.nan, -100, inplace = True)
    jump_df['logrhk'].replace(np.nan, -100, inplace = True)
    jump_df['rp'].replace(np.nan, -100, inplace = True)
    jump_df = jump_df.drop(index = jump_df.query("Name == 'T001244'").index.values) # Remove T001244, which has vmag = 11.9
    
    
    distantgiants_photo = pd.merge(jump_df.drop(columns = ['vmag']), distantgiants_photo[['cps', 'Vmag']], left_on = 'Name', right_on = 'cps', how = 'right').rename(columns = {'Vmag':'vmag'})
    
    # Creating a new column to tell whether a target has a template on either HIRES or APF
    distantgiants_photo['have_template_hires_j'].replace(np.nan, 0, inplace = True)
    distantgiants_photo['have_template_apf'].replace(np.nan, 0, inplace = True)
    distantgiants_photo['have_template'] = list(map(lambda x, y: 'NO' if x+y < 1 else 'YES', distantgiants_photo['have_template_hires_j'], distantgiants_photo['have_template_apf']))

    distantgiants_spec = distantgiants_photo.query('vsini <= 5.00 and logrhk < -4.7').reset_index(drop = True)
    distantgiants_spec.rename(columns = {'Name':'star_id'}, inplace = 'True')
    # spec_cuts.at[pd.Index(spec_cuts['star_id']).get_loc('55CNC'), 'star_id'] = '75732'
    
    no_no = ['1300', 'T001538', '97658', '110067', '192279', 'T000561', 'T001386', 'T001391', 
             'T001504', 'T001537', 'T001609', 'T001648', 'T001655', 'T001699', 'T001711']
    distantgiants_spec = distantgiants_spec.drop(distantgiants_spec.index[distantgiants_spec['star_id'].isin(no_no)])
    
    distantgiants_spec.to_csv('csv/distantgiants_spec.csv', index = False)
    
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
    
    
    
    
