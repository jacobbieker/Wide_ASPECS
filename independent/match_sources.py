"""

This is focused on matching sources in the catalog to those detected in the cubes

"""
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS, Distance
from astropy.coordinates import match_coordinates_sky, search_around_sky
from stats import comoving_volume, get_kms, has_spec_z, get_co_z, convert_deltaZ_to_kms

transitions = {"1-0": [0.0030, 0.3694, 115.271, 0.2801, 89],
               "2-1": [1.0059, 1.7387, 230.538, 1.4277, 1920],
               "3-2": [2.0088, 3.1080, 345.796, 2.6129, 3363],
               "4-3": [3.0115, 4.4771, 461.041, 3.8030, 4149], }


def convert_observed_line_to_restframe():
    return NotImplementedError


def calculate_delta_z():
    return NotImplementedError


def estimate_redshift():
    return NotImplementedError


def match_lines_to_catalog(lines, catalog, max_redshift=0.3, snr_limit=6., max_sep=1.0, method='closest'):
    aspecs_table = Table(names=(
    'RA (J2000)', 'DEC (J2000)', 'Roberto ID', 'Roberto RA', 'Roberto DEC', 'Observed CO (GHz)', 'Restframe CO (GHz)',
    'Transition', 'Z (Matched)', 'Z (CO)',
    'Spec Z', 'Delta Z', 'Delta V (Km/s)', 'Km/s', 'Separation (Arcsecond)', 'S/N', 'Flux Density at Peak (Jy/beam)',
    'Integrated Flux (Jy km/s)', 'Width (Channels)', 'Cosmic Volume (Mpc^3)', 'Log(M*)', 'Error Log(M*)', 'Log(SFR)',
    'Error Log(SFR)'),
                         dtype=(
                         'f8', 'f8', 'int32', 'f8', 'f8', 'f4', 'f4', 'U6', 'f4', 'f4', 'bool', 'f4', 'f8', 'f8', 'f4',
                         'f4', 'f4', 'f4', 'int8', 'f4', 'f4', 'f4', 'f4', 'f4'))

    """
    Steps to do so:
    
    Find separations between line coordinates and catalog coordinates
    
    For those that are within the arcsecond limit, see if the galactic redshift is within the range that ASPECS can find
    
        If so, then get the difference in delta_z to see if that is within the range allowed
            If so, then get the properties and put together a whole entry on it
        If not, see if line matches to a different CO line within that range
            If so, save it out
    If not within range, see which line it could go to and use that one
    
    
    """

    # first step is to do is get the SkyCoords

    catalog_ra = 'ra_2'
    catalog_dec = 'dc'

    # Only choose ones above SN limit
    lines = lines[lines['rsnrrbin'] >= snr_limit]

    line_skycoords = make_skycoords(lines, ra='rra', dec='rdc')
    catalog_skycoords = make_skycoords(catalog, ra=catalog_ra, dec=catalog_dec)

    catalog_ids = []

    # Second step is to calculate the catalog matches
    if method == 'all_closest':
        # This is for getting all the matches, and only keeping the one with the closest redshift
        # Do it where it goes through all matches within a given radius
        idxc, idxcatalog, d2d, d3d = search_around_sky(line_skycoords, catalog_skycoords, max_sep * u.arcsecond)

        # Many to many is way too large to work, so doing it one by one
        print("Matching done")
        print(len(idxc))

        # Get the set of chosen lines, all not chosen ones are sent to the other thing
        chosen_lines = set(idxc)
        full_set = set([i for i in range(len(lines))])
        non_matched_set_indexes = full_set - chosen_lines

        closest_redshift_and_catalog_id = [(-9999, -9999, -99999) for i in range(len(full_set))]

        for index, separation in enumerate(d2d):
            matched_line = lines[idxc[index]]
            if separation.arcsecond < max_sep:
                # Could be a match!
                # Get the catalog match
                matched_galaxy = catalog[idxcatalog[index]]  # index is the index in line_skycoord matched
                # idx[index] is then the index into catalog that is matched to this one
                for key, values in transitions.items():
                    if (values[0] - max_redshift) < matched_galaxy['z_1'] < (values[1] + max_redshift):
                        # Now within range of this transition
                        delta_z, matched_key = get_delta_z(matched_galaxy['z_1'], matched_line['rfreq'] * u.GHz)
                        if np.abs(delta_z) <= max_redshift:  # Checks that delta z within the range
                            # Now check with offset if the z is within the range
                            if matched_galaxy['z_1'] + delta_z < (0.4) or (1.1) <= matched_galaxy['z_1'] + delta_z <= (
                                    1.8) or (2.2) < matched_galaxy['z_1'] + delta_z < (4.4):
                                if np.abs(closest_redshift_and_catalog_id[idxc[index]][0]) > np.abs(delta_z):
                                    closest_redshift_and_catalog_id[idxc[index]] = (
                                    delta_z, separation.arcsecond, catalog[idxcatalog[index]]['id'])
                                elif np.isclose(np.abs(closest_redshift_and_catalog_id[idxc[index]][0]),
                                                np.abs(delta_z)):
                                    # Closest so based off separation
                                    if closest_redshift_and_catalog_id[idxc[index]][1] > separation.arcsecond:
                                        closest_redshift_and_catalog_id[idxc[index]] = (
                                        delta_z, separation.arcsecond, catalog[idxcatalog[index]]['id'])
                                # so with offset, the galaxy is now within the range, is above SNR, and have a transition
                                # Now get the KMS, if there is a Spec Z, Comoving volume, etc. and add to the table
                                volume = comoving_volume(values[0], values[1], 52.5)
                                spec_z = has_spec_z(matched_galaxy)
                                kms = get_kms(matched_line['width'], matched_line['rfreq'])
                                rest_frame_ghz = convert_to_rest_frame_ghz(matched_galaxy['z_1'],
                                                                           matched_line['rfreq'] * u.GHz)
                                co_z = get_co_z(matched_line['rfreq'], matched_key)
                                delta_v = convert_deltaZ_to_kms(delta_z)
                                new_row = (np.round(matched_line['rra'], 6),
                                           np.round(matched_line['rdc'], 6),
                                           np.int(matched_galaxy['id']),
                                           np.round(matched_galaxy[catalog_ra], 6),
                                           np.round(matched_galaxy[catalog_dec], 6),
                                           matched_line['rfreq'],
                                           rest_frame_ghz,
                                           matched_key,
                                           matched_galaxy['z_1'],
                                           co_z,
                                           spec_z,
                                           delta_z,
                                           delta_v,
                                           kms,
                                           np.round(separation.arcsecond, 4),
                                           matched_line['rsnrrbin'],
                                           matched_line['rpeak'],
                                           matched_line['rflux'],
                                           matched_line['width'],
                                           np.round(volume, 3),
                                           matched_galaxy['Mstar_50_1'],
                                           matched_galaxy['Mstar_84_1'] - matched_galaxy['Mstar_50_1'],
                                           matched_galaxy['SFR_50_1'],
                                           matched_galaxy['SFR_84_1'] - matched_galaxy['SFR_50_1'])
                                aspecs_table.add_row(new_row)
                            else:  # The matched_galaxy + delta_z not within range so move on to trying to do the CO line matching
                                table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                                aspecs_table.add_row(table_input)
                        else:  # The matched_galaxy not within range so move on to trying to do the CO line matching
                            table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                            aspecs_table.add_row(table_input)
                    else:  # The separation is larger than the max allowed, so move on to trying to do the CO line matching
                        table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                        aspecs_table.add_row(table_input)
            else:  # The separation is larger than the max allowed, so move on to trying to do the CO line matching
                table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                aspecs_table.add_row(table_input)
        # Now have to do it for the non-matched ones
        for index in non_matched_set_indexes:
            matched_line = lines[index]
            table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
            aspecs_table.add_row(table_input)

        # Now need to only get the catalog ids that are relevant, so not -99999
        catalog_ids = [i[2] for i in closest_redshift_and_catalog_id if i[2] > -999]

    if method == 'all':
        # Do it where it goes through all matches within a given radius
        idxc, idxcatalog, d2d, d3d = search_around_sky(line_skycoords, catalog_skycoords, max_sep * u.arcsecond)

        # Many to many is way too large to work, so doing it one by one
        print("Matching done")
        print(len(idxc))

        # Get the set of chosen lines, all not chosen ones are sent to the other thing
        chosen_lines = set(idxc)
        full_set = set([i for i in range(len(lines))])
        non_matched_set_indexes = full_set - chosen_lines

        for index, separation in enumerate(d2d):
            matched_line = lines[idxc[index]]
            if separation.arcsecond < max_sep:
                # Could be a match!
                # Get the catalog match
                matched_galaxy = catalog[idxcatalog[index]]  # index is the index in line_skycoord matched
                # idx[index] is then the index into catalog that is matched to this one
                for key, values in transitions.items():
                    if (values[0] - max_redshift) < matched_galaxy['z_1'] < (values[1] + max_redshift):
                        # Now within range of this transition
                        delta_z, matched_key = get_delta_z(matched_galaxy['z_1'], matched_line['rfreq'] * u.GHz)
                        if np.abs(delta_z) <= max_redshift:  # Checks that delta z within the range
                            # Now check with offset if the z is within the range
                            if matched_galaxy['z_1'] + delta_z < (0.4) or (1.1) <= matched_galaxy['z_1'] + delta_z <= (
                                    1.8) or (2.2) < matched_galaxy['z_1'] + delta_z < (4.4):
                                catalog_ids.append(catalog[idxcatalog[index]]['id'])
                                # so with offset, the galaxy is now within the range, is above SNR, and have a transition
                                # Now get the KMS, if there is a Spec Z, Comoving volume, etc. and add to the table
                                volume = comoving_volume(values[0], values[1], 52.5)
                                spec_z = has_spec_z(matched_galaxy)
                                kms = get_kms(matched_line['width'], matched_line['rfreq'])
                                rest_frame_ghz = convert_to_rest_frame_ghz(matched_galaxy['z_1'],
                                                                           matched_line['rfreq'] * u.GHz)
                                co_z = get_co_z(matched_line['rfreq'], matched_key)
                                delta_v = convert_deltaZ_to_kms(delta_z)
                                new_row = (np.round(matched_line['rra'], 6),
                                           np.round(matched_line['rdc'], 6),
                                           np.int(matched_galaxy['id']),
                                           np.round(matched_galaxy[catalog_ra], 6),
                                           np.round(matched_galaxy[catalog_dec], 6),
                                           matched_line['rfreq'],
                                           rest_frame_ghz,
                                           matched_key,
                                           matched_galaxy['z_1'],
                                           co_z,
                                           spec_z,
                                           delta_z,
                                           delta_v,
                                           kms,
                                           np.round(separation.arcsecond, 4),
                                           matched_line['rsnrrbin'],
                                           matched_line['rpeak'],
                                           matched_line['rflux'],
                                           matched_line['width'],
                                           np.round(volume, 3),
                                           matched_galaxy['Mstar_50_1'],
                                           matched_galaxy['Mstar_84_1'] - matched_galaxy['Mstar_50_1'],
                                           matched_galaxy['SFR_50_1'],
                                           matched_galaxy['SFR_84_1'] - matched_galaxy['SFR_50_1'])
                                aspecs_table.add_row(new_row)
                            else:  # The matched_galaxy + delta_z not within range so move on to trying to do the CO line matching
                                table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                                aspecs_table.add_row(table_input)
                        else:  # The matched_galaxy not within range so move on to trying to do the CO line matching
                            table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                            aspecs_table.add_row(table_input)
                    else:  # The separation is larger than the max allowed, so move on to trying to do the CO line matching
                        table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                        aspecs_table.add_row(table_input)
            else:  # The separation is larger than the max allowed, so move on to trying to do the CO line matching
                table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                aspecs_table.add_row(table_input)
        # Now have to do it for the non-matched ones
        for index in non_matched_set_indexes:
            matched_line = lines[index]
            table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
            aspecs_table.add_row(table_input)

    if method == 'closest':
        idx, d2d, d3d = match_coordinates_sky(line_skycoords, catalog_skycoords)
        # So now, idx is the index into catalog_skycoords to get the matched coordinate for line_skycoords, shape=line_skycoord
        # d2d is on sky separation between line_skycoords and its closest match
        # So the catalog_skycoord[idx] is the match for the line_skycoord

        for index, ident in enumerate(idx):
            matched_line = lines[index]
            matched_to_galaxy = False
            # Check if above the SNR limit
            if d2d[index].arcsecond < max_sep:
                # Could be a match!
                # Get the catalog match
                matched_galaxy = catalog[ident]  # index is the index in line_skycoord matched
                # idx[index] is then the index into catalog that is matched to this one
                for key, values in transitions.items():
                    if (values[0] - max_redshift) < matched_galaxy['z_1'] < (values[1] + max_redshift):
                        # Now within range of this transition
                        rest_frame_ghz = convert_to_rest_frame_ghz(matched_galaxy['z_1'],
                                                                   matched_line['rfreq'])
                        delta_z, matched_key = get_delta_z(matched_galaxy['z_1'], rest_frame_ghz)
                        if np.abs(delta_z) <= max_redshift:  # Checks that delta z within the range
                            # Now check with offset if the z is within the range
                            if matched_galaxy['z_1'] + delta_z < (0.4) or (1.1) <= matched_galaxy['z_1'] + delta_z <= (
                                    1.8) or (2.2) < matched_galaxy['z_1'] + delta_z < (4.4):
                                matched_to_galaxy = True
                                catalog_ids.append(matched_galaxy['id'])
                                # so with offset, the galaxy is now within the range, is above SNR, and have a transition
                                # Now get the KMS, if there is a Spec Z, Comoving volume, etc. and add to the table
                                volume = comoving_volume(values[0], values[1], 52.5)
                                spec_z = has_spec_z(matched_galaxy)
                                kms = get_kms(matched_line['width'], matched_line['rfreq'])
                                co_z = get_co_z(matched_line['rfreq'], matched_key)
                                delta_v = convert_deltaZ_to_kms(delta_z)
                                new_row = (np.round(matched_line['rra'], 6),
                                           np.round(matched_line['rdc'], 6),
                                           np.int(matched_galaxy['id']),
                                           np.round(matched_galaxy[catalog_ra], 6),
                                           np.round(matched_galaxy[catalog_dec], 6),
                                           matched_line['rfreq'],
                                           rest_frame_ghz,
                                           matched_key,
                                           matched_galaxy['z_1'],
                                           co_z,
                                           spec_z,
                                           delta_z,
                                           delta_v,
                                           kms,
                                           np.round(d2d[index].arcsecond, 4),
                                           matched_line['rsnrrbin'],
                                           matched_line['rpeak'],
                                           matched_line['rflux'],
                                           matched_line['width'],
                                           np.round(volume, 3),
                                           matched_galaxy['Mstar_50_1'],
                                           matched_galaxy['Mstar_84_1'] - matched_galaxy['Mstar_50_1'],
                                           matched_galaxy['SFR_50_1'],
                                           matched_galaxy['SFR_84_1'] - matched_galaxy['SFR_50_1'])
                                aspecs_table.add_row(new_row)
            if not matched_to_galaxy:
                table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                if table_input is not None:
                    aspecs_table.add_row(table_input)

    # now have the catalog matches:
    aspecs_catalog = catalog[catalog_ids]

    return aspecs_table, aspecs_catalog


def match_to_co_line(single_line, max_redshift=0.3):
    """
    Match a single line to a CO line transition
    :param single_line:
    :return: Data to then be added to the Table
    """
    estimated_transition, estimated_z = get_estimated_z(single_line['rfreq'])
    if estimated_z < 0.4 or 1.1 <= estimated_z <= 1.8 or 2.2 < estimated_z < 4.4:
        # TODO Get the MS, rest_ghs, delta_z etc.
        rest_frame_ghz = convert_to_rest_frame_ghz(estimated_z, single_line['rfreq'])
        delta_z, matched_key = get_delta_z(estimated_z, rest_frame_ghz)
        if np.abs(delta_z) <= max_redshift:
            spec_z = False
            volume = comoving_volume(transitions[estimated_transition][0], transitions[estimated_transition][1], 52.5)
            kms = get_kms(single_line['width'], single_line['rfreq'])
            co_z = get_co_z(single_line['rfreq'], matched_key)
            delta_v = convert_deltaZ_to_kms(delta_z)
            new_row = (np.round(single_line['rra'], 6),
                       np.round(single_line['rdc'], 6),
                       -999,
                       -999,
                       -999,
                       single_line['rfreq'],
                       rest_frame_ghz,
                       matched_key,
                       estimated_z,
                       co_z,
                       spec_z,
                       delta_z,
                       delta_v,
                       kms,
                       -999,
                       single_line['rsnrrbin'],
                       single_line['rpeak'],
                       single_line['rflux'],
                       single_line['width'],
                       np.round(volume, 3),
                       -999,
                       -999,
                       -999,
                       -999)

            return new_row
        else:
            return None
    else:
        return None


def make_skycoords(source, ra='ra', dec='dec', distance=None):
    """
    Makes and returns a SkyCoord array from given source
    :param source: Source with information
    :param ra: Key for RA
    :param dec: Key for Dec
    :return: SkyCoord list
    """
    if distance is None:
        skycoords = SkyCoord(source[ra] * u.deg, source[dec] * u.deg, frame='fk5')
    else:
        distances = Distance(z=source[distance])
        skycoords = SkyCoord(source[ra] * u.deg, source[dec] * u.deg, distance=distances, frame='fk5')

    return skycoords


def get_observed_ghz(z, transition):
    """
    Get the observed GHz for given redshift based on transition, from Wide ASPECS paper
    :param z: Z to calculate for
    :param transition:
    :return:
    """
    emitted_ghz = transitions[transition][2] * u.GHz

    observed_ghz = emitted_ghz / (z + 1)

    return observed_ghz


def get_delta_z(z, rest_ghz):
    """
    Take a measured GHz value, and calculates the restframe GHz value based on the given z of the matched galaxy
    :param z:
    :param ghz:
    :return:
    """

    # First step is to convert to nm rom rest frame GHz
    set_zs = []
    for key, values in transitions.items():
        if values[0] - 0.3 <= z <= values[1] + 0.3:
            sghz = values[2] * u.GHz  # Gets the GHz of the CO line
            rest_ghz /= (z + 1)
            set_z = np.round((sghz - rest_ghz) / rest_ghz, 3)  # (Freq_emitted - Freq_obs)/ Freq_obs = z
            set_z = z - set_z
            rest_ghz *= (z + 1)
            set_zs.append((key, set_z))
    set_z = np.min([np.abs(i[1]) for i in set_zs])
    for element in set_zs:
        if np.isclose(np.abs(element[1]), set_z):
            return element[1], element[0]


def convert_to_rest_frame_ghz(z, ghz):
    """
    Take a measured GHz value, and calculates the restframe GHz value based on the given z of the matched galaxy
    :param z:
    :param ghz:
    :return:
    """

    observed_ghz = ghz * u.GHz

    emitted = observed_ghz * (z + 1)

    return emitted


def get_estimated_z(ghz):
    """
    Estimate the CO line based on Wide-ASPECS one, (3-2), z > 2 or higher J, calculate possible Z's and find which Z is closest
    to the <z> value from Wide ASPECS
    :param ghz:
    :return: transition, estimated_z
    """
    differences = []
    for key, values in transitions.items():
        # Convert ghz to rest_ghz of the Z value, otherwise always closest to lowest one
        sghz = convert_to_rest_frame_ghz(values[3], ghz)
        delta_z, matched_key = get_delta_z(values[3], sghz)
        differences.append((matched_key, delta_z))

    min_diff = np.min([np.abs(i[1]) for i in differences])

    for index, element in enumerate(differences):
        if np.isclose(np.abs(element[1]), min_diff):
            return element[0], transitions[element[0]][3]
