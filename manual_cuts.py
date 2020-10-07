import numpy as np
import pandas as pd

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

def update_manual_cuts_2(verbose = False):
    # Set the search radius - I've noticed that if you set it too small, it misses some targets, weirdly?
    # So I usually set it large at the beginning then cut the list down
    radius = 300     # units of arcsec
    width = u.Quantity(radius, u.arcsec)
    height = u.Quantity(radius, u.arcsec)


    distantgiants_photo = pd.read_csv('/Users/judahvz/research/code/GitHub/distantgiants/csv/distantgiants_photo.csv')

    # These are all photo targets, but with MES and qlp information included. This df is made MANUALLY from distantgiants_photo, so it will not update automatically if new TOIs are added to tois_perfect.csv
    manual_cuts_1 = pd.read_csv('/Users/judahvz/research/code/GitHub/distantgiants/csv/manual_cuts_1.csv')
    manual_cuts_1 = manual_cuts_1.query("`qlp_only?`== False and MES >= 12.0")

    # Merge to get ra and dec
    manual_cuts_1 = pd.merge(manual_cuts_1[['cps', 'MES']], distantgiants_photo, on = 'cps')
    data = manual_cuts_1

    # List of TIC IDs and RA / Dec coordinates

    all_companions = []
    for star_index in range(len(data)):
    
        # ra and dec values according to Jump
        ra = float(data["ra"][star_index])
        dec = float(data["dec"][star_index])
        vmag = float(data["Vmag"][star_index])
        tic = str(data["tic"][star_index])
        star_id = str(data["cps"][star_index])

    
        # Search Gaia database for companions near the target
        coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
        r = Gaia.query_object_async(coordinate=coord, width=width, height=height, verbose=False)[['ra', 'dec', 'phot_g_mean_mag']]
        r = pd.DataFrame(np.transpose([r['ra'], r['dec'], r['phot_g_mean_mag']]), columns = r.columns).sort_values(by='phot_g_mean_mag')
    
        # There is often significant disagreement between Jump and Gaia, so we look for the primary within a 100" radius
        # and 1 magnitude difference
        mag_condition = abs(vmag-r['phot_g_mean_mag']) < 1
        dist_condition = np.sqrt((ra - r['ra'])**2 + (dec-r['dec'])**2)*3600 < 100
    
        primary_index_list = r[dist_condition & mag_condition].index
    
        # If no star in the field matches my values for the primary, make a note of that with the print statement and then just assign the first target as the primary. Shaky solution but it prevents erroring out.
        if len(primary_index_list) == 0:
            if verbose:
                print('limits too strict!!!!')
            primary_index_list = [0]
        
        increment = 0.1
        sep_cutoff = 10 # arcsec
        while len(primary_index_list) > 1:
            if verbose:
                print('bumpin it down a lil')

            cutoff -= increment
            dist_condition = np.sqrt((ra - r['ra'])**2 + (dec-r['dec'])**2)*3600 < cutoff
        
        
            primary_index_list = r[ mag_condition & dist_condition ].index
        
        primary_star_index = primary_index_list[0]
        # ra and dec values according to Gaia
        primary_ra = r['ra'][primary_star_index]
        primary_dec = r['dec'][primary_star_index]
    
        # This number is the distance cutoff to look for companions in arcsec (currently 5 arcsec)
        cutoff = 5   # arcsec
        line_list = []

        for companion_index in range(len(r)):
        
            comp_ra = r['ra'][companion_index]
            comp_dec = r['dec'][companion_index]
            comp_mag = r['phot_g_mean_mag'][companion_index]
        
            distance = np.sqrt((primary_ra - comp_ra)**2 + (primary_dec - comp_dec)**2)*3600
            distance = distance
    #         print('distance is ',distance)
        
            if distance < cutoff:
                if verbose:
                    print('made the 5 cut')
                line = pd.DataFrame([[star_id, comp_ra, comp_dec, distance, comp_mag]], columns=['star_id', 'ra', 'dec', 'distance', 'gmag'])
                line_list.append(line)
            
        if len(line_list) > 1:
            companion_df = pd.concat([i for i in line_list], ignore_index=True)
            all_companions.append(companion_df.sort_values(by='gmag'))

    companion_list = [i['star_id'][0] for i in all_companions]
    print('{} targets were cut on close companions'.format(len(companion_list)))

    # manual_cuts_2 has undergone cuts on qlp/MES (1) and on Gaia separation (2) 
    manual_cuts_2 = manual_cuts_1.drop(manual_cuts_1[manual_cuts_1['cps'].isin(companion_list)].index).reset_index(drop=True)

    manual_cuts_2.to_csv('/Users/judahvz/research/code/GitHub/distantgiants/csv/manual_cuts_2.csv', index=False)
    print('Updated manual_cuts_2 ({} targets)'.format(len(manual_cuts_2)))
    
if __name__ == '__main__':
    update_manual_cuts_2(verbose = True)
