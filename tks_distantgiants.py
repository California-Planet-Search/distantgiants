"""
Judah Van Zandt
Python 3

Module to update tks_distantgiants.txt on GitHub and generate the distantgiants.csv file, both of which contain all targets
that have been chosen by Ashley Chontos's prioritization code and can therefore be requested for observations.

Written: June 19, 2020
Last modified: July 8, 2020
"""

import numpy as np
import pandas as pd
from git import repo
import tks_distantgiants_spec as spec

def make_distantgiants():
    
    distantgiants_spec = pd.read_csv('csv/distantgiants_spec.csv')

    # observing_priorities = pd.read_csv('../tks_target_list_gen/data/observing_priorities.csv')
    observing_priorities = pd.read_csv('csv/observing_priorities.csv')
    
    overlap = [observing_priorities['toi'][i] for i in range(len(observing_priorities)) if 'SC2A' in observing_priorities['programs'][i]]
    print('There are {} targets assigned to SC2A.'.format(len(overlap)))
    # print('number of SC2A in Ashley list:', len(overlap))
    # The following are columns that I want displayed in the distantgiants paper
    relevant_columns = ['tic', 'star_id', 'toi', 'ra', 'dec', 'vmag', 'teff', 'mass', 'radius', 'vsini', 'logrhk', 'sval', 'ruwe', 'rp', 'period'] # Add in per and svalue
    
    ###################
    ## I used to compose DG of only targets assigned by Ashley's code. I'm taking this out because it only removes cooked targets, which I actually want to include in final analysis.
    # distantgiants = pd.merge(observing_priorities['tic'], distantgiants_spec, how = 'inner', on = 'tic').rename(columns={'Rp':'rp', 'Rs':'radius', 'smass':'mass'})[relevant_columns].sort_values(by='toi')
    
    distantgiants = distantgiants_spec
    ###################
    
    distantgiants.to_csv('csv/distantgiants.csv', index = False)
    
    return distantgiants
    
    
def update_distantgiants(distantgiants):
    
    out_file = open('../jump-config/programs/tks_distantgiants.txt', 'w+')

    for star_name in distantgiants.sort_values(by = 'star_id').star_id.values:
        out_file.write(star_name+'\n')


    out_file.close()

    
    jump_config_path = r'../jump-config'
    commit_message = 'Updated distantgiants.txt'
    
    my_repo = repo.Repo(jump_config_path) # Path to jump_config repo
    my_repo.git.add(update=True) # Adds updates to existing files to index, rather than index.add, which I think works for                              adding *new* files
    my_repo.index.commit(commit_message) # Add commit message to the index (staging area between working dir and repo)
    origin = my_repo.remote(name='origin') # Specify where to push changes
    origin.pull()
    origin.push()
    print('Updated distantgiants.txt ({} targets)'.format(len(distantgiants)))
    

if __name__ == '__main__':
    
    update_distantgiants(make_distantgiants())