import astropy.io.ascii as ascii
import astropy.table
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS
from astropy import units as u
from astropy.table import Table
from astropy.io.ascii import FixedWidth
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from astropy.coordinates import FK5

def get_aspecs_radec(frame='fk5'):
    coords = []
    freqs = []

    with open(os.path.join("data", "ASPECS_lines.txt")) as data_file:
        # Read in the locations and
        for line in data_file:
            no_circle = line.split("(")[1]
            first_split = no_circle.split(",")
            ra = float(first_split[0])
            dec = float(first_split[1])
            coords.append(SkyCoord(ra*u.deg, dec*u.deg, frame=frame))
            frequency = first_split[2].split("{")[1]
            frequency = float(frequency.split("}")[0])
            freqs.append(frequency)

    coords = SkyCoord(coords)
    return coords, freqs

def build_aspecs_catalog(initial_catalog=None, dec_key='dec', ra_key='ra', frame='fk5'):
    hdu_list = fits.open(initial_catalog)
    initial_catalog = hdu_list[1].data
    ra_dec = SkyCoord(initial_catalog[ra_key] * u.deg, initial_catalog[dec_key] * u.deg, frame=frame)

    coords = []
    freqs = []

    with open(os.path.join("data", "ASPECS_lines.txt")) as data_file:
        # Read in the locations and
        for line in data_file:
            no_circle = line.split("(")[1]
            first_split = no_circle.split(",")
            ra = float(first_split[0])
            dec = float(first_split[1])
            coords.append(SkyCoord(ra*u.deg, dec*u.deg, frame=frame))
            frequency = first_split[2].split("{")[1]
            frequency = float(frequency.split("}")[0])
            freqs.append(frequency)

    coords = SkyCoord(coords)
    idx, d2d, d3d = coords.match_to_catalog_sky(ra_dec)
    print("\n----------------- Number of Matches: " + str(len(idx)) + "/" + str(len(initial_catalog[idx])))
    print("Distances: ")
    num_in_close = 0
    catalog_ids = []
    for index, id in enumerate(idx):
        if coords[index].separation(ra_dec[id]).arcsecond < 1.0:
            num_in_close += 1
            print("\nMatch: " + str(index))
            print("Distance: " + str(coords[index].separation(ra_dec[id]).arcsecond))
            print("Location (RA): " + str(coords[index].ra.hms))
            print("Location (Dec): " + str(coords[index].dec.hms))
            print("Location (Deg): " + str(coords[index]))
            try:
                print("Catalog ID: " + str(initial_catalog[id]['id']))
                catalog_ids.append((initial_catalog[id]['id'],index))
                print("In Skelton et al. Catalog: " + str(initial_catalog[id]['flag_3dh']))
            except:
                try:
                    catalog_ids.append((initial_catalog[id]['unique_id'],index))
                    print("Catalog ID: " + str(initial_catalog[id]['unique_id']))
                    print("Skelton et al. ID: " + str(initial_catalog[id]['skelton_id']))
                except:
                    continue
    print("Number Close to Catalog: ", num_in_close)
    print("Catalog IDs")
    print(catalog_ids)

    # Get the IDs of the matched values
    #catalog_ids = initial_catalog[idx]['id']

def build_aspecs_catalog_ascii(initial_catalog=None, dec_key='dec', ra_key='ra', frame='fk5'):
    hdu_list = fits.open(initial_catalog)
    initial_catalog = hdu_list[1].data
    ra_dec = SkyCoord(initial_catalog[ra_key] * u.deg, initial_catalog[dec_key] * u.deg, frame=frame)

    aspecs_lines = Table.read("data/line_search_P3_wa_crop.out", format="ascii")

def compare_catalog_locations(roberto_catalog, initial_catalog, ra_key='ra', dec_key='dec', frame='fk5'):
    hdu_list = fits.open(initial_catalog)
    initial_catalog = hdu_list[1].data
    emline = fits.open("data/MW_44fields_emline_table_v1.0.fits")[1].data
    diff_ids = []
    for id in initial_catalog['unique_id']:
        if id in emline['UNIQUE_ID']:
            diff_ids.append(id)
    rows_to_use = []
    for index, row in enumerate(initial_catalog):
        if row['unique_id'] in diff_ids:
            rows_to_use.append(index)
    initial_catalog = initial_catalog[rows_to_use]

    ra_dec = SkyCoord(initial_catalog[ra_key] * u.deg, initial_catalog[dec_key] * u.deg, frame=frame)
    hdu_list = fits.open(roberto_catalog)
    rob_cat = hdu_list[1].data
    roberto_ra_dec = SkyCoord(rob_cat['ra'] * u.deg, rob_cat['dc'] * u.deg, frame=frame)

    # Now compare the two catalogs to find matches
    idx, d2d, d3d = roberto_ra_dec.match_to_catalog_sky(ra_dec)
    # Matches less than 0.5 arc seconds
    less_than_5 = []
    # Less than 0.25 arc seconds
    less_than_25 = []
    # Less than 1 arc second
    less_than_1 = []
    for index, id in enumerate(idx):
        if roberto_ra_dec[index].separation(ra_dec[id]).arcsecond < 1:
            less_than_1.append(id)
        if roberto_ra_dec[index].separation(ra_dec[id]).arcsecond < 0.5:
            less_than_5.append(id)
        if roberto_ra_dec[index].separation(ra_dec[id]).arcsecond < 0.25:
            less_than_25.append(id)
            if False:
                print("\nMatch: " + str(index))
                print("Location: " + str(roberto_ra_dec[index]))
                print("Distance: " + str(roberto_ra_dec[index].separation(ra_dec[id]).arcsecond))
                try:
                    print("MUSE Catalog Sep: " + str(initial_catalog[id]['skelton_sep']))
                    print("Difference Between MUSE and Roberto: " + str(initial_catalog[id]['skelton_sep'] - roberto_ra_dec[index].separation(ra_dec[id]).arcsecond))
                    print("Catalog ID: " + str(rob_cat[index]['id']))
                    print("In Skelton et al. Catalog: " + str(rob_cat[index]['flag_3dh']))
                    print("Skelton et al. ID: " + str(initial_catalog[id]['skelton_id']))
                    print("MUSE Catalog ID: " + str(initial_catalog[id]['unique_id']))
                except:
                    continue

    print("Less than 1.0 arcseconds: " + str(len(np.unique(less_than_1))))
    print("Less than 0.5 arcseconds: " + str(len(np.unique(less_than_5))))
    print("Less than 0.25 arcseconds: " + str(len(np.unique(less_than_25))))

    roberto_catalog = Table.read("roberto_catalog_muse_skelton.fits", format='fits')
    # Now add those matches within a given constraint to the FITS file to get add another Z fit
    roberto_catalog['muse_wide_z'] = np.zeros(len(rob_cat['ra']))
    roberto_catalog['muse_wide_z_err'] = np.zeros(len(rob_cat['ra']))
    roberto_catalog['muse_id'] = np.zeros(len(rob_cat['ra']))
    roberto_catalog['muse_quality'] = np.zeros(len(rob_cat['ra']))
    roberto_catalog['skelton_id'] = np.zeros(len(rob_cat['ra']))
    print(idx)
    total_matches = 0
    not_in_error_z = []
    for index, galaxy in enumerate(idx):
        # Use
        # Change this if want better/worse coverage
        if roberto_ra_dec[index].separation(ra_dec[galaxy]).arcsecond < 0.5:
            total_matches += 1
            roberto_catalog[index]['muse_wide_z'] = initial_catalog[galaxy]['z']
            roberto_catalog[index]['muse_wide_z_err'] = initial_catalog[galaxy]['z_err']
            roberto_catalog[index]['muse_id'] = initial_catalog[galaxy]['unique_id']
            roberto_catalog[index]['muse_quality'] = initial_catalog[galaxy]['confidence']
            roberto_catalog[index]['skelton_id'] = initial_catalog[galaxy]['skelton_id']
            low_z_range = np.round(roberto_catalog[index]['muse_wide_z'] - roberto_catalog[index]['muse_wide_z_err'], decimals=3)-0.15
            high_z_range = np.round(roberto_catalog[index]['muse_wide_z'] + roberto_catalog[index]['muse_wide_z_err'], decimals=3)+0.15

            if roberto_catalog[index]['z'] < low_z_range or roberto_catalog[index]['z'] > high_z_range:
                # Offset not consistent
                not_in_error_z.append(index)

    print(not_in_error_z)
    muse_error_ids = []
    for index in not_in_error_z:
        muse_error_ids.append(roberto_catalog[index]['muse_id'])
        print("\nMUSE ID: " + str(roberto_catalog[index]['muse_id']))
        print("Roberto Best Z: " + str(roberto_catalog[index]["z"]))
        print("MUSE Best Z: " + str(np.round(roberto_catalog[index]["muse_wide_z"], decimals=3)) +"+-" + str(np.round(roberto_catalog[index]["muse_wide_z_err"], decimals=3)))
        print("Difference: " + str(np.round(roberto_catalog[index]['muse_wide_z'] - roberto_catalog[index]['z'], decimals=3)))
    print("Number of inconsistent matches: " + str(len(not_in_error_z)) + "/" + str(total_matches))
    roberto_magphys = Table.read("/home/jacob/Research/magphys_in_latest.fits", format='fits')
    print(roberto_magphys.columns)
    print(len(roberto_magphys))

    # Check the same Z redshifts vs both default and muse wide Z
    same_as_muse = 0
    same_as_overall = 0
    not_same = 0
    for index, galaxy in enumerate(roberto_magphys):
        # Get the ID
        id_galaxy = galaxy['id']

        id_mask = (roberto_catalog['id'] == id_galaxy)
        # Check if it is same Z
        gal = roberto_catalog[id_mask]
        # Same ID
        if np.isclose(galaxy['z'],gal['muse_wide_z']):
            #print("\nSame As MUSE Z ID: " + str(id_galaxy))
            same_as_muse += 1
        if np.isclose(galaxy['z'], gal['z']):
            #print("\nSame Overall Z ID: " + str(id_galaxy))
            same_as_overall += 1
        if not np.isclose(galaxy['z'], gal['z']) and not np.isclose(galaxy['z'], gal['muse_wide_z']):
            not_same += 1

    print("Same Overall: " + str(same_as_overall))
    print("Same as MUSE: " + str(same_as_muse))
    print("Not Same: " + str(not_same))
    print(muse_error_ids)

    #roberto_catalog.write("roberto_catalog_muse_skelton_matched.fits", format='fits')
    return roberto_catalog

def compare_open_catalog_locations(roberto_catalog, initial_catalog, ra_key='ra', dec_key='dec', frame='fk5'):
    ra_dec = SkyCoord(initial_catalog[ra_key] * u.deg, initial_catalog[dec_key] * u.deg, frame=frame)
    roberto_ra_dec = SkyCoord(roberto_catalog['ra'] * u.deg, roberto_catalog['dc'] * u.deg, frame=frame)

    # Now compare the two catalogs to find matches
    idx, d2d, d3d = roberto_ra_dec.match_to_catalog_sky(ra_dec)
    # Matches less than 0.5 arc seconds
    less_than_5 = []
    # Less than 0.25 arc seconds
    less_than_25 = []
    # Less than 1 arc second
    less_than_1 = []
    for index, id in enumerate(idx):
        if roberto_ra_dec[index].separation(ra_dec[id]).arcsecond < 1:
            less_than_1.append([index, id])
        if roberto_ra_dec[index].separation(ra_dec[id]).arcsecond < 0.5:
            less_than_5.append([index, id])
        if roberto_ra_dec[index].separation(ra_dec[id]).arcsecond < 0.25:
            less_than_25.append([index, id])
            if False:
                print("\nMatch: " + str(index))
                print("Location: " + str(roberto_ra_dec[index]))
                print("Distance: " + str(roberto_ra_dec[index].separation(ra_dec[id]).arcsecond))
                try:
                    print("MUSE Catalog Sep: " + str(initial_catalog[id]['skelton_sep']))
                    print("Difference Between MUSE and Roberto: " + str(initial_catalog[id]['skelton_sep'] - roberto_ra_dec[index].separation(ra_dec[id]).arcsecond))
                    print("Catalog ID: " + str(roberto_catalog[index]['id']))
                    print("In Skelton et al. Catalog: " + str(roberto_catalog[index]['flag_3dh']))
                    print("Skelton et al. ID: " + str(initial_catalog[id]['skelton_id']))
                    print("MUSE Catalog ID: " + str(initial_catalog[id]['unique_id']))
                except:
                    continue

    print("Less than 1.0 arcseconds: " + str(len(less_than_1)))
    print("Less than 0.5 arcseconds: " + str(len(less_than_5)))
    print("Less than 0.25 arcseconds: " + str(len(less_than_25)))


if __name__ == "__main__":

    ids = [13077, 9228, 14566, 14890,
           17087, 17425, 18199, 18470,
           19407, 21424, 22074, 22619,
           23331, 24915, 28752, 28802,
           51344, 51553, 56747, 57167,
           57179, 61484]

    rest_co = [377,182,102,111,235,249,178,311,246,280,143,209,219,184,212,234,212,241,420,246,243,365,246]
    ones = {"1-0":3, "2-1":15, "3-2": 3, "4-3": 1}

    one_sets = [3,15,3,1]

    ones_expanded = [0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,3]

    plt.hist(ones_expanded, bins=4)

    locs, labels = plt.xticks()
    plt.xticks(locs, ['', '1-0', '2-1', '3-2', '4-3'])
    plt.show()

    zs = [3.11, 0.925, 0.085,
          0.163, 1.535,1.685,
          0.872, 2.32, 1.665,
          1.675, 0.355, 0.96,
          1.094, 0.738, 1.019,
          1.191, 1.269, 1.537,
          3.496, 1.599, 1.597,
          2.835, 1.599]

    #build_aspecs_catalog(os.path.join("data", "jacob_aspecs_catalog_fixed_magphys_jcb3.fits"), dec_key='dc')
    #build_aspecs_catalog(os.path.join("data", "MW_44fields_main_table_v1.0.fits"))
    #build_aspecs_catalog(os.path.join("roberto_catalog_muse.fits"), dec_key='dc')



    #compare_catalog_locations(os.path.join("data", "jacob_aspecs_catalog_fixed_magphys_jcb3.fits"), os.path.join("data", "MW_44fields_main_table_v1.0.fits"))

    roberto_catalog = Table.read("roberto_catalog_muse_skelton_matched.fits", format='fits')
    magphys = Table.read("/home/jacob/Development/Wide_ASPECS/mapghys_in_nov2018_all.fits", format='fits')

    print(len(magphys))
    same_as_muse = 0
    different_z = 0
    same_z = 0
    different_from_muse = 0
    big_diff = 0
    total = 0

    for x, id in enumerate(ids):
        mask = (magphys['id'] == id)
        print(magphys[id])
        gal = magphys[mask]
        print("Zs: {} Cat: {}".format(zs[x], gal['z']))
    for galaxy in magphys:
        gal_id = galaxy['id']
        total += 1
        mask = roberto_catalog['id'] == gal_id
        masked_big = roberto_catalog[mask]
        if np.isclose(np.round(masked_big['z'],3), np.round(galaxy['z'],3)):
            same_z += 1
        if np.isclose(np.round(masked_big['muse_wide_z'],3), np.round(galaxy['z'],3)):
            same_as_muse += 1
        if not np.isclose(np.round(masked_big['z'],3), np.round(galaxy['z'],3)):
            if 1.01*np.round(masked_big['z'],3) < galaxy['z'] or 0.99*np.round(masked_big['z'],3) > galaxy['z']:
                different_z += 1
                #print("Different ID: {} \nMy Z: {}\n Magphys Z: {}\n".format(masked_big['id'][0], np.round(masked_big['muse_wide_z'],3), galaxy['z']))
        if not np.isclose(np.round(masked_big['muse_wide_z'],3), np.round(galaxy['z'],3)):
            if masked_big['muse_wide_z'] > 0:
                different_from_muse += 1
                if 1.15*np.round(masked_big['muse_wide_z'],3) < galaxy['z'] or 0.85*np.round(masked_big['muse_wide_z'],3) > galaxy['z']:
                    print("ID: {} \nMUSE Z: {}\n Catalog Z: {}\n".format(masked_big['id'][0], np.round(masked_big['muse_wide_z'],3), galaxy['z']))
                    big_diff += 1

    print("Smae As Muse: " + str(same_as_muse))
    print("Same Overall " + str(same_z))
    print("different Overall " + str(different_z))
    print("different From Muse " + str(different_from_muse))
    print("Large Diff From Muse " + str(big_diff))
    print("Total: " + str(total))



    catalog_goodss = fits.open("/home/jacob/Research/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS")
    skelton_goodss = catalog_goodss[1].data

    muse_hdu = fits.open(os.path.join("data", "MW_44fields_main_table_v1.0.fits"))
    muse_catalog = muse_hdu[1].data

    # Add Skelton IDs to the Roberto Catalog with MUSE
    muse_roberto_catalog = fits.open("roberto_catalog_muse.fits")[1].data

    #ra_dec = SkyCoord(skelton_goodss['ra'] * u.deg, skelton_goodss['dec'] * u.deg, frame='fk5')
    roberto_ra_dec = SkyCoord(muse_roberto_catalog['ra'] * u.deg, muse_roberto_catalog['dc'] * u.deg, frame='fk5')
    #idx, d2d, d3d = roberto_ra_dec.match_to_catalog_sky(ra_dec)

    # Matches less than 0.5 arc seconds
    less_than_5 = []
    # Less than 0.25 arc seconds
    less_than_25 = []
    # Less than 1 arc second
    less_than_1 = []
    '''
    roberto_catalog = Table.read("roberto_catalog_muse.fits", format='fits')

    num_less_than_15 = 0
    roberto_catalog['skelton_id'] = np.zeros(len(roberto_ra_dec))

    for index, id in enumerate(idx):
        if roberto_ra_dec[index].separation(ra_dec[id]).arcsecond < 0.1:
            num_less_than_15 += 1
            roberto_catalog[index]['skelton_id'] = skelton_goodss[id]['id']


    print("Less than 0.1 arcseconds: " + str(num_less_than_15))
    roberto_catalog.write("roberto_catalog_muse_skelton.fits", format='fits')
    '''

    # Now have skelton IDs, can match directly with MUSE matchings
    # Only MUSE ones without skelton IDs are the separations used
    roberto_catalog = Table.read("roberto_catalog_muse_skelton.fits", format='fits')
    roberto_catalog['muse_wide_z'] = np.zeros(len(roberto_catalog['ra']))
    roberto_catalog['muse_wide_z_err'] = np.zeros(len(roberto_catalog['ra']))
    roberto_catalog['muse_id'] = np.zeros(len(roberto_catalog['ra']))
    roberto_catalog['muse_quality'] = np.zeros(len(roberto_catalog['ra']))

    # Now have that, match by skelton_id, then if id not zero, match by Sky
    for index, row in enumerate(roberto_catalog):
        skelton_id = row['skelton_id']
        if skelton_id != 0:
            muse_mask = (np.isclose(muse_catalog['skelton_id'], skelton_id))
            if len(muse_catalog[muse_mask]) > 1:
                print(muse_catalog[muse_mask])
            if len(muse_catalog[muse_mask]) == 1:
                roberto_catalog[index]['muse_wide_z'] = muse_catalog[muse_mask]['z']
                roberto_catalog[index]['muse_wide_z_err'] = muse_catalog[muse_mask]['z_err']
                roberto_catalog[index]['muse_id'] = muse_catalog[muse_mask]['unique_id']
                roberto_catalog[index]['muse_quality'] = muse_catalog[muse_mask]['confidence']
                roberto_catalog[index]['z'] = muse_catalog[muse_mask]['z']
            if len(muse_catalog[muse_mask]) == 2:
                roberto_catalog[index]['muse_wide_z'] = muse_catalog[muse_mask][0]['z']
                roberto_catalog[index]['muse_wide_z_err'] = muse_catalog[muse_mask][0]['z_err']
                roberto_catalog[index]['muse_id'] = muse_catalog[muse_mask][0]['unique_id']
                roberto_catalog[index]['muse_quality'] = muse_catalog[muse_mask][0]['confidence']
                roberto_catalog[index]['z'] = muse_catalog[muse_mask][0]['z']

    # Now only those with no skelton id match:
    muse_ra_dec = SkyCoord(muse_catalog['ra'] * u.deg, muse_catalog['dec'] * u.deg, frame='fk5')

    idx, d2d, d3d = roberto_ra_dec.match_to_catalog_sky(muse_ra_dec)
    num_less_than_15 = 0
    less_than_15 = []
    for index, id in enumerate(idx):
        if np.isclose(roberto_catalog[index]['muse_id'], 0):
            if np.isclose(muse_catalog[id]['skelton_id'], 0):
                # No match and not already in Skelton
                if roberto_ra_dec[index].separation(muse_ra_dec[id]).arcsecond < 1:
                    less_than_1.append([index, id])
                    if roberto_ra_dec[index].separation(muse_ra_dec[id]).arcsecond < 0.5:
                        # Check if already matched to Skelton object
                        if muse_catalog[id]['unique_id'] not in roberto_catalog['muse_id']:
                            # Now do the by hand adding
                            less_than_5.append([roberto_catalog[index]['id'], muse_catalog[id]['unique_id']])
                            roberto_catalog[index]['muse_wide_z'] = muse_catalog[id]['z']
                            roberto_catalog[index]['muse_wide_z_err'] = muse_catalog[id]['z_err']
                            roberto_catalog[index]['muse_id'] = muse_catalog[id]['unique_id']
                            roberto_catalog[index]['muse_quality'] = muse_catalog[id]['confidence']
                        else:
                            print("ID Already Exists with Roberto ID ")
                        if roberto_ra_dec[index].separation(muse_ra_dec[id]).arcsecond < 0.25:
                            less_than_25.append([index, id])
                            if roberto_ra_dec[index].separation(muse_ra_dec[id]).arcsecond < 0.1:
                                less_than_15.append([index, id])



    print("Less than 1.0 arcseconds: " + str(len(np.unique(less_than_1))) + "/" + str(len(less_than_1)))
    print("Less than 0.5 arcseconds: " + str(len(np.unique(less_than_5))) + "/" + str(len(less_than_5)))
    print("Less than 0.25 arcseconds: " + str(len(np.unique(less_than_25))) + "/" + str(len(less_than_25)))
    print("Less than 0.1 arcseconds: " + str(len(np.unique(less_than_15))) + "/" + str(len(less_than_15)))
    print(less_than_15)
    print("----------------  Now Less Than 0.5 arcseconds")
    print(less_than_5)

    roberto_catalog.write("roberto_catalog_muse_skelton_matched_manFix.fits", format='fits')






"""
full_goodss = fits.open("data/jacob_aspecs_catalog_fixed_magphys_jcb3.fits")
full_goodss = full_goodss[1].data
roberto_catalog = fits.open("data/magphys_in.fits")
spec_ids = full_goodss['id']
diff_ids = []
roberto_catalog = roberto_catalog[1].data
for id in roberto_catalog['id']:
    if id in spec_ids:
        diff_ids.append(id)
# now create smaller catalog with full catalog info

rows_to_use = []
for index, row in enumerate(full_goodss):
    if row['id'] in diff_ids:
        rows_to_use.append(index)

smaller_catalog = full_goodss[rows_to_use]

initial_catalog = fits.open(os.path.join("data", "MW_44fields_main_table_v1.0.fits"))[1].data

compare_open_catalog_locations(smaller_catalog, initial_catalog)
"""

"""

Vastly different Roberto and MUSE redshifts

[146041313.0, 146069355.0, 145049108.0, 146052337.0, 145006019.0, 145006019.0, 146082368.0, 145022065.0, 145034089.0, 
145019060.0, 145005017.0, 146018244.0, 146044319.0, 145038096.0, 146061346.0, 145050109.0, 144049136.0, 144024074.0, 
144025075.0, 143031111.0, 143042127.0, 144050137.0, 143040124.0, 143009021.0, 144057147.0, 143037118.0, 143002006.0, 
107020110.0, 142044159.0, 141028124.0, 141052165.0, 141011083.0, 141016089.0, 142002064.0, 141039152.0, 141002074.0, 
141003075.0, 141003075.0, 141040153.0, 141006078.0, 141004076.0, 131021114.0, 132020035.0, 141017090.0, 132015026.0, 
132009020.0, 131051155.0, 132048102.0, 131011094.0, 132003006.0, 130033059.0, 131016105.0, 132045098.0, 131048152.0, 
131046150.0, 137070146.0, 130015022.0, 130001001.0, 131041145.0, 131014103.0, 129017133.0, 131053157.0, 137053126.0, 
130006006.0, 128064270.0, 130038064.0, 129019135.0, 131013100.0, 136012137.0, 130003003.0, 130009010.0, 137083161.0, 
128045247.0, 136007129.0, 137012028.0, 128044246.0, 135037227.0, 135038228.0, 135009176.0, 136044198.0, 129002081.0, 
135015184.0, 128032224.0, 137043105.0, 128007194.0, 136045199.0, 128055260.0, 128040241.0, 135005170.0, 134006016.0, 
128025217.0, 128038236.0, 136037187.0, 134027046.0, 140002014.0, 128019210.0, 128019210.0, 134009019.0, 140024077.0, 
134036056.0, 122017082.0, 135035225.0, 140018054.0, 139013229.0, 134037057.0, 135002046.0, 134011023.0, 134002002.0, 
134002002.0, 140005024.0, 134025044.0, 140012041.0, 134030050.0, 139046300.0, 139019244.0, 135026216.0, 139051305.0, 
133020056.0, 139055311.0, 134038060.0, 135041231.0, 139032271.0, 133033082.0, 134018030.0, 135046236.0, 133029074.0, 
134020036.0, 134003004.0, 134003003.0, 133028073.0, 139054310.0, 139054310.0, 133017053.0, 139044297.0, 133032080.0, 
125007023.0, 125006022.0, 125059136.0, 125001001.0, 125001001.0, 125052125.0, 125026077.0, 125037106.0, 125049122.0, 
125009025.0, 125009025.0, 125040110.0, 126015035.0, 126047134.0, 126061172.0, 126042110.0, 126049137.0] 

"""