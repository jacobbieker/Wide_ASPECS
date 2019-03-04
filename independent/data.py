"""

This is focused on data subsets, so creating cuts in the data, selecting z-selected sections, etc.

"""

from astropy.table import Table


def perform_cuts():
    return NotImplementedError


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
    return NotImplementedError
