"""

This is focused on data subsets, so creating cuts in the data, selecting z-selected sections, etc.

"""

from astropy.table import Table, hstack, join


def perform_cuts(catalog):
    """
    Perform quality cuts like expected
    :param catalog:
    :return:
    """
    quality_cuts = (catalog["Q0_1"] < 2.0) & (catalog["Q2_1"] < 1.0) & (
            (catalog["Mstar_50_1"] - catalog["Mstar_16_1"]) < 0.5) & \
                   ((catalog["Mstar_84_1"] - catalog["Mstar_50_1"]) < 0.5) & \
                   ((catalog["SFR_50_1"] - catalog["SFR_16_1"]) < 0.5) & \
                   ((catalog["SFR_84_1"] - catalog["SFR_16_1"]) < 0.5)

    return catalog[quality_cuts]


def select_spectroscopic_sources(catalog):
    return NotImplementedError


def load_table(ascii_table, header=0, start=1):
    ascii_table_data = Table.read(ascii_table, format="ascii", header_start=header, data_start=start)
    return ascii_table_data

def load_catalog(catalog):
    return NotImplementedError


def save_catalog(catalog):
    return NotImplementedError


def save_ascii(catalog):
    return NotImplementedError


def combine_catalogs(catalog_one, catalog_two):

    combined = join(catalog_one, catalog_two, keys='id')
    return combined

def create_points_and_error(column_base_name, initial_catalog):
    centerpoints = initial_catalog[str(column_base_name + "_50_1")]
    lower_error = initial_catalog[str(column_base_name + "_16_1")]
    upper_error = initial_catalog[str(column_base_name + "_84_1")]
    centerpoints = np.nan_to_num(centerpoints)
    zero_mask = centerpoints != 0.0
    centerpoints = centerpoints[zero_mask]
    lower_error = centerpoints - lower_error[zero_mask]
    upper_error = upper_error[zero_mask] - centerpoints
    error_bars = [lower_error, upper_error]

    return centerpoints, error_bars


def create_points_and_error_by_z(column_base_name, initial_catalog, low_z, high_z):
    z_mask = (initial_catalog["z_1"] >= low_z) & (initial_catalog["z_1"] <= high_z)
    centerpoints = initial_catalog[z_mask][str(column_base_name + "_50_1")]
    lower_error = initial_catalog[z_mask][str(column_base_name + "_16_1")]
    upper_error = initial_catalog[z_mask][str(column_base_name + "_84_1")]
    z_values = initial_catalog[z_mask]["z_1"]
    centerpoints = np.nan_to_num(centerpoints)
    lower_error = np.nan_to_num(lower_error)
    upper_error = np.nan_to_num(upper_error)
    zero_mask = centerpoints != 0.0
    centerpoints = centerpoints[zero_mask]
    lower_error = centerpoints - lower_error[zero_mask]
    upper_error = upper_error[zero_mask] - centerpoints
    z_values = z_values[zero_mask]
    error_bars = [lower_error, upper_error]
    return centerpoints, error_bars, z_values