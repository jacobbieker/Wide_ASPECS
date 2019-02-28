"""

This is focused on matching sources in the catalog to those detected in the cubes

"""
import numpy as np

def convert_observed_line_to_restframe():
    return NotImplementedError

def calculate_delta_z():
    return NotImplementedError

def estimate_redshift():
    return NotImplementedError

def match_lines_to_catalog(lines, catalog, max_redshift=0.3):
    return NotImplementedError

