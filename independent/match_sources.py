"""

This is focused on matching sources in the catalog to those detected in the cubes

"""
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle, SkyOffsetFrame, ICRS, Distance
>>> from astropy.coordinates import match_coordinates_sky



transitions = {"1-0": [0.0030, 0.3694, 115.271, 0.2801, 89],
               "2-1": [1.0059, 1.7387, 230.538, 1.4277, 1920],
               "3-2": [2.0088, 3.1080, 345.796, 2.6129, 3363],
               "4-3": [3.0115, 4.4771, 461.041, 3.8030, 4149],}

def convert_observed_line_to_restframe():
    return NotImplementedError

def calculate_delta_z():
    return NotImplementedError

def estimate_redshift():
    return NotImplementedError

def match_lines_to_catalog(lines, catalog, max_redshift=0.3, snr_limit=6., max_sep=1.0):
    aspecs_table = Table(names=('RA (J2000)', 'DEC (J2000)', 'Roberto ID', 'Roberto RA', 'Roberto DEC',  'Observed CO (GHz)', 'Restframe CO (GHz)', 'Transition', 'Z (Matched)', 'Z (CO)',
                                'Spec Z', 'Delta Z', 'Km/s', 'Separation (Arcsecond)', 'S/N', 'Flux Density at Peak (Jy/beam)',
                                'Integrated Flux (Jy km/s)', 'Width (Channels)', 'Cosmic Volume (Mpc^3)', 'Log(M*)', 'Error Log(M*)', 'Log(SFR)', 'Error Log(SFR)'),
                         dtype=('f8', 'f8', 'int32', 'f8', 'f8', 'f4', 'f4', 'U6', 'f4', 'f4', 'bool', 'f4', 'f8', 'f4', 'f4', 'f4', 'f4', 'int8', 'f4', 'f4', 'f4', 'f4', 'f4'))

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

    line_skycoords = make_skycoords(lines, ra='rra', dec='rdc')
    catalog_skycoords = make_skycoords(catalog, ra='ra', dec='dc')

    # Second step is to calculate the catalog matches

    idx, d2d, d3d = match_coordinates_sky(line_skycoords, catalog_skycoords)
    # TODO This only matches to the closest match though... farther matches in the same area might not be shown
    # So now, idx is the index into catalog_skycoords to get the matched coordinate for line_skycoords, shape=line_skycoord
    # d2d is on sky separation between line_skycoords and its closest match
    # So the catalog_skycoord[idx] is the match for the line_skycoord

    for index, separation in enumerate(d2d):
        matched_line = lines[index]
        # Check if above the SNR limit
        if matched_line['rsnrrbin'] >= snr_limit:
            if separation.arcsecond < max_sep:
                # Could be a match!
                # Get the catalog match
                matched_galaxy = catalog_skycoords[idx[index]] # index is the index in line_skycoord matched
                # idx[index] is then the index into catalog that is matched to this one
                for key, values in transitions.items():
                    if (values[0] - max_redshift) < matched_galaxy['z_1'] < (values[1] + max_redshift):
                        # Now within range of this transition
                        delta_z, matched_key = get_delta_z(matched_galaxy['z_1'], matched_line['rfreq'] * u.GHz)
                        if np.abs(delta_z) < max_redshift: # Checks that delta z within the range
                            # Now check with offset if the z is within the range
                            if matched_galaxy['z_1'] + delta_z < (0.4) or (1.1) <= matched_galaxy['z_1'] + delta_z <= (1.8) or (2.2) < matched_galaxy['z_1'] + delta_z < (4.4):
                                # so with offset, the galaxy is now within the range, is above SNR, and have a transition
                                # Now get the KMS, if there is a Spec Z, Comoving volume, etc. and add to the table
                                # TODO Fill that OUT
                                return NotImplementedError
                            else: # The matched_galaxy + delta_z not within range so move on to trying to do the CO line matching
                                table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                                aspecs_table.add_row(table_input)
                    else: # The matched_galaxy not within range so move on to trying to do the CO line matching
                        table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                        aspecs_table.add_row(table_input)
            else: # The separation is larger than the max allowed, so move on to trying to do the CO line matching
                table_input = match_to_co_line(matched_line, max_redshift=max_redshift)
                aspecs_table.add_row(table_input)

    return aspecs_table

def match_to_co_line(single_line, max_redshift=0.3):
    """
    Match a single line to a CO line transition
    :param single_line:
    :return: Data to then be added to the Table
    """
    estimated_transition, estimated_z = get_estimated_z(single_line['rfreq'])
    if estimated_z < 0.4 or 1.1 <= estimated_z <= 1.8 or 2.2 < estimated_z < 4.4:
        #TODO Get the MS, rest_ghs, delta_z etc.
        rest_frame_ghz = convert_to_rest_frame_ghz(estimated_z, single_line['rfreq'])
        delta_z, matched_key = get_delta_z(estimated_z, rest_frame_ghz)
        if np.abs(delta_z) <= max_redshift:
            has_spectroscopic = False
            comoving_volume =
    return NotImplementedError


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

def get_delta_z(z, observed_ghz):
    """
    Take a measured GHz value, and calculates the restframe GHz value based on the given z of the matched galaxy
    :param z: Z of the matched galaxy
    :param observed_ghz: The frequency of the observed line in GHz
    :return: Separation, Transition
    """

    set_zs = []
    for key, values in transitions.items():
        if values[0] <= z <= values[1]:
            emitted_ghz = values[2] * u.GHz # Gets the GHz of the CO line
            set_z = np.round((emitted_ghz - observed_ghz)/ observed_ghz, 3) # (Freq_emitted - Freq_obs)/ Freq_obs = z
            set_z = z - set_z
            set_zs.append((key, set_z))
    set_z = np.min([np.abs(i[1]) for i in set_zs])
    for element in set_zs:
        if np.isclose(np.abs(element[1]),set_z):
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