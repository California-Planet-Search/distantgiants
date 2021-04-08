"""
Judah Van Zandt
Python 3

This module updates the tks_distantgiants_pool.txt and tks_distantgiants.txt files on GitHub, as well as creating the 
overview.csv file and updating the list of observing requests for the SC2A program.

Written: June 1, 2020
Last modified: July 8, 2020
"""

import tks_distantgiants_photo as tks_photo
import tks_distantgiants_spec as tks_spec
import tks_distantgiants as tks

import overview as ov
from target_request_generator import generator, obs_request_list_gen, init_overview
from astropy.time import Time


tks_photo.update_distantgiants_photo(tks_photo.make_distantgiants_photo())
tks_spec.update_distantgiants_spec(tks_spec.make_distantgiants_spec(tks_photo.make_distantgiants_photo()))
tks.update_distantgiants(tks.make_distantgiants())

ov.update_overview(ov.make_overview(plot = True, observability = False))
generator(obs_request_list_gen(init_overview(iers=False)))


