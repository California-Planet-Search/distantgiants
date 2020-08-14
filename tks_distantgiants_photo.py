"""
Judah Van Zandt
Python 3

Module to generate the tks_distgiants_pool target list, which comprises all stars from selected_TOIs that pass
pre-spectroscopic vetting. These are then written out to the jump-config/programs folder on Github.

Written: June 1, 2020
Last modified: July 8, 2020
"""

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import gspread
from git import repo

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


def make_distantgiants_photo():
    """
    Accesses the selected_TOIs spreadsheet so that SC2A can be regularly updated.
    
    Performs photometric cuts on selected_TOIs and writes out the stars that pass to tks_distgiants_pool.txt, which is in the
    jump-config repo.
    
    Pushes the changes made to tks_distgiants_photo.txt to the jump-config Github.
    """
    tois_perfect_location = '../tks_target_list_gen/prioritization/info/TOIs_perfect.csv'
    column_changes = {'cps_name':'cps', 'm_s':'smass', 'r_s':'sradius'}
    tois_perfect = pd.read_csv(tois_perfect_location).rename(columns = column_changes)
    tois_perfect = tois_perfect[['tic','toi','cps','ra','dec','vmag','smass','sradius','evol','ruwe','rp', 'period']]

    tks_drop = pd.read_csv('csv/tks_drop.csv')
    
    # Turn these columns from strings to floats
    for column in ['dec', 'smass', 'ruwe', 'vmag', 'rp']:
        tois_perfect['{}'.format(column)] = tois_perfect['{}'.format(column)].astype(float)
        
   

    ## Photometric Cuts
    tois_perfect = tois_perfect.drop_duplicates(subset = 'cps')
    tois_perfect = tois_perfect.query('dec > 0 and evol=="MS" and ruwe < 1.3 and vmag <= 12 and rp < 10 and smass > 0.5 and smass < 1.5')
    tois_perfect = tois_perfect.rename(columns = {'vmag' : 'Vmag', 'sradius' : 'Rs', 'rp' : 'Rp'}).reset_index(drop=True)
    tois_perfect.at[pd.Index(tois_perfect['cps']).get_loc('T001290'), 'cps'] = 'K00246'
    # T001443 has a close companion listed on Jump; Erik says no reason to drop
    # tois_perfect = tois_perfect.drop(tois_perfect[tois_perfect['cps'] == 'T001443'].index)
    # tois_perfect.at[pd.Index(tois_perfect['cps']).get_loc('T001823'), 'cps'] = 'TIC142381532'
    
    distantgiants_photo = tois_perfect[['cps','toi','tic','ra','dec','Vmag','Rs', 'smass', 'Rp', 'ruwe', 'evol', 'period']]
    distantgiants_photo = distantgiants_photo[~distantgiants_photo['cps'].isin(tks_drop['Name'])].reset_index(drop = True) # Drop any stars that are in tks_drop
    
    
    distantgiants_photo.to_csv('csv/distantgiants_photo.csv', index = False)
    
    return distantgiants_photo

def update_distantgiants_photo(distantgiants_photo):
    out_file = open('../jump-config/programs/tks_distantgiants_photo.txt', 'w+')

    for star_name in distantgiants_photo.sort_values(by = 'cps').cps.values:
        out_file.write(star_name+'\n')


    out_file.close()

    
    jump_config_path = r'../jump-config'
    commit_message = 'Updated distantgiants_photo.txt'
    
    my_repo = repo.Repo(jump_config_path) # Path to jump_config repo
    my_repo.git.add(update=True) # Adds updates to existing files to index, rather than index.add, which I think works for                              adding *new* files
    my_repo.index.commit(commit_message) # Add commit message to the index (staging area between working dir and repo)
    origin = my_repo.remote(name='origin') # Specify where to push changes
    origin.pull()
    origin.push()
    print('Updated distantgiants_photo.txt ({} targets)'.format(len(distantgiants_photo)))


if __name__ == "__main__":
    
  update_distantgiants_photo(make_distantgiants_photo())



