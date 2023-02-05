import nmrglue as ng
import numpy as np


# Helper function: get noise from region without peaks ##FALG->DOCSTRING##
def get_noise(dic, data, width=1):
    """ 
    Calculate noise from edge of spectrum


    PARAMETERS
    ----------
    dic (dictionary): dictionary of NMR parameters in NMRPipe format
    data (numpy.ndarray): array of NMR data in NMRPipe format
    width (float): width of region for noise determination in ppm

    RETURNS
    -------
    noise (float): noise of spectrum as standard deviation
    """
    # get unit converter object
    uc = get_uc(dic, data)

    # get spectrum limits
    l_ppm, r_ppm = get_limits(uc)

    noise = np.std(data[uc.i(l_ppm, 'ppm'):uc.i(l_ppm-width, 'ppm')].real)

    return noise


# Helper function: create unit converter object
def get_uc(dic, data):
    """
    Creates an unit conversion object from NMR data


    PARAMETERS
    ----------
    dic (dictionary): dictionary of NMR parameters in NMRPipe format
    data (numpy.ndarray): array of NMR data in NMRPipe format

    RETRUNS
    -------
    uc (unit conversion object): unit conversion object (see nmrglue)

    """
    uc = ng.pipe.make_uc(dic, data)

    return uc

# Helper function: get spectrum limits


def get_limits(uc):
    """
    Get spectrum limits via unit conversion object


    PARAMTERS
    ---------
    uc (unit conversion object): unit conversion object (see nmrglue)

    RETURNS
    -------
    l_ppm (float): left limit of spectrum in ppm
    r_ppm (float): right limit of spectrum in ppm

    """
    # Get limits in ppm
    l_ppm, r_ppm = uc.ppm_limits()

    return l_ppm, r_ppm
