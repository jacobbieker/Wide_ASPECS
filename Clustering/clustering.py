import numpy as np
import matplotlib.pyplot as plt


def angular_distance():
    raise NotImplementedError


def generate_random_catalog(number_of_points, mask):
    """
    Generates random number of points in position of Wide ASPECS  and within the mask



    :param number_of_points:
    :param mask:
    :return:
    """
    raise NotImplementedError


def angular_correlation_function(data_catalog, random_catalog):
    """
    Calculates the arrays for the data, random, and data_random for w(theta)

    :param data_catalog:
    :param random_catalog:
    :return:
    """
    distance_bins = np.logspace(0,np.log10(50), 10)



    raise NotImplementedError



def xi_r(data_array, data_random_array, random_array):
    """

    :param data_array:
    :param data_random_array:
    :param random_array:
    :return:
    """
    data_array /= np.sum(data_array)
    data_random_array /= np.sum(data_random_array)
    random_array /= np.sum(random_array)
    return data_array / random_array - 2 * data_random_array / random_array + 1


def xi_r_error(omega_theta, data_array):
    """
    Data array is not normalized
    :param omega_theta:
    :param data_array:
    :return:
    """

    return (1 + omega_theta) / np.sqrt(data_array)
