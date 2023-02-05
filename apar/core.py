import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from .helpers import get_noise, get_uc, get_limits


# Load NMR data ##FLAG->execption##
def load_data(path, format):
    """
    Load raw data from directory

    Note: see docs of nmrglue for more options


    Parameters
    ----------
    path (str): path to directory with raw data
    format (str): format of raw data, resp., instrument used
        supported formats: bruker

    Return
    ------
    dic (dictionary): dictionary of NMR parameters
    data (numpy.ndarray): array of NMR data

    """
    if format == "bruker":
        # read data
        dic, data = ng.bruker.read(path)
    else:
        print("format not supported")
        # throw execpt => format not supported

    return dic, data


# Convert data to NMRPipe format
def convert_data(dic, data, format="bruker"):
    """
    Convert data to NMRPipe format


    PARAMETERS
    ----------
    dic (dictionary): dictionary of NMR parameters
    data (numpy.ndarray): array of NMR data
    format (str): format of raw data, resp., instrument used
        supported formats: bruker

    RETURNS
    -------
    converted_dic (dictionary): dictionary of NMR parameters
    converted_data (numpy.ndarray): array of NMR data in NMRPipe format

    """
    if format == "bruker":
        # remove digital filter
        data = ng.bruker.remove_digital_filter(dic, data)

        # convert to NMRPipe
        C = ng.convert.converter()
        C.from_bruker(dic, data)
        converted_dic, converted_data = C.to_pipe()
    else:
        print("format not supported")
        # throw execpt => format not supported

    return converted_dic, converted_data


# Pre-processing raw data ##FLAG->update processing##
def process_data(dic, data):
    """
    Pre-processing of NMR raw data

    Note: optimized for 1H NMR recorded on a BRUKER AVANCE III HD spectrometer 
    equipped with 5 mm BBO Prodigy CryoProbe and perfect echo W5 WATERGATE 
    solvent suppression pulse sequence


    PARAMETERS
    ----------
    dic (dictionary): dictionary of NMR parameters in NMRPipe format
    data (numpy.ndarray): array of NMR data in NMRPipe format

    RETURNS
    -------
    processed_dic (dictionary): dictionary of NMR parameters
    processed_data (numpy.ndarray): array of processed NMR data in NMRPipe format

    """
    # zero filling
    processed_data = ng.proc_base.zf(data)
    # bdata = ng.proc_base.zf_size(bdata, 32768)          # zero fill to 32768 points
    # dic, data = ng.pipe_proc.apod(dic, data, qName='EM')        # Apodization

    # ????
    processed_dic, processed_data = ng.pipe_proc.em(
        dic, processed_data, lb=0.5)

    # Fourier transformation
    processed_data = ng.proc_base.fft(processed_data)

    # automatic linear phase correction
    processed_data = ng.proc_autophase.autops(processed_data, fn="acme")

    # reverse data
    processed_data = ng.proc_base.rev(processed_data)

    # smooth data
    processed_data = ng.proc_base.smo(processed_data, 2)

    # adata = ng.proc_base.fft(adata)                       # Fourier transform
    # automatic linear phase correction
    # adata = ng.proc_autophase.autops(adata, fn="acme")
    # adata = ng.proc_base.rev(adata)                       # reverse the data
    # adata = ng.proc_base.smo(adata, 2)

    return processed_dic, processed_data


# Get parameters of internal standard used
def get_standards():
    """
    Get a dictionary of internal standards with chemical shift positions


    Returns
    -------
    standards (dictionary): dictionary of internal standards
        abbreviation, chemical shift in ppm, number of hydrogen atoms

    """
    standards = {"DMSO": (2.60, 6), "DMSO2": (2.96, 6)}

    return standards


# Identification of internal standard: find peak position
def find_standard(dic, data, limit=0.05, threshold=1e5, standard="DMSO"):
    """
    Find peak position of internal standards within limits


    Parameters
    ----------
    dic (dictionary): dictionary of NMR parameters in NMRPipe format
    data (numpy.ndarray): array of processed NMR data
    limit (float): limits for peak identification region in ppm
    threshold (float): intensity threshold for automated peak picking
    standard (str): name of internal standard 

    Returns
    -------
    real_shift (numpy.float64): peak position of standard in ppm
    expected_shift (numpy.float64): expected peak posistion of standard in ppm

    """
    # get unit converter object
    uc = get_uc(dic, data)

    # list of internal standards and their expected chemical shift
    standards = get_standards()
    expected_shift = standards[standard][0]

    # define upper and lower limit
    upper_limit = expected_shift + limit
    lower_limit = expected_shift - limit

    # slice array to region of standard
    sliced_data = data.copy()
    sliced_data[0:uc.i(upper_limit, 'ppm')], sliced_data[uc.i(
        lower_limit, 'ppm'):] = 0, 0

    # detect all peaks above threshold
    peaks = ng.peakpick.pick(
        sliced_data, pthres=threshold, algorithm="connected", table=True)

    # peak position of standard (chemical shift) in ppm assuming highest intensity
    real_shift = uc.ppm(peaks[0][0])

    return real_shift, expected_shift


# Chemical shift referencing
def referencing(dic, data, real_shift, expected_shift):
    """ 
    Chemical shift referencing based on resonance of internal standard


    PARAMETERS
    ----------
    dic (dictionary): dictionary of NMR parameters in NMRPipe format
    data (numpy.ndarray): array of NMR data in NMRPipe format
    real_shift (numpy.float64): peak position of standard in ppm
    expected_shift (numpy.float64): expected peak posistion of standard in ppm

    RETURNS
    -------
    referenced_data (numpy.ndarray): array of NMR data referenced on resonance 
    of internal standard

    """
    # get unit converter object
    uc = get_uc(dic, data)

    # get spectrum limits
    l_ppm, r_ppm = get_limits(uc)

    # calculate difference of real and expected chemical shifts and convert to pts
    shift = uc.i((l_ppm - (real_shift - expected_shift)), 'ppm')

    # perform frequency shift
    referenced_data = ng.process.proc_base.roll(data, pts=shift)

    return referenced_data


# Peak identification
def peak_identification(dic, data, SNR=4):
    """ 
    Peak identification with given SNR as intensity threshold


    PARAMETERS
    ----------
    dic (dictionary): dictionary of NMR parameters in NMRPipe format
    data (numpy.ndarray): array of NMR data in NMRPipe format
    SNR (float): signal-to-noise ratio used as intensity threshold

    RETURNS
    -------
    peaks (numpy.recarray): array of identified peaks
        consists of peak position (X_AXIS), cluster number (cID), 
        estimated peak scales (linewidths, X_LW), and estimated
        peak amplitudes (VOL)

    """
    # set threshold based on SNR
    threshold = SNR * get_noise(dic, data)

    # detect all peaks within threshold
    peaks = ng.peakpick.pick(data, pthres=threshold, algorithm="connected")

    return peaks


# Peak integration
def peak_integration(dic, data, peaks, water_shift=4.7, water_width=0.2, integration_limit=0.05):
    """
    Peak integration of identified peaks with fixed integration limits


    PARAMETERS
    ----------
    dic (dictionary): dictionary of NMR parameters in NMRPipe format
    data (numpy.ndarray): array of NMR data in NMRPipe format
    peaks (numpy.recarray): array of identified peaks
    water_shift (float): expected peak position of water in ppm
    water_width (float): expected width of water peak in ppm
    integration_limit (float): fixed integration limits in ppm

    RETURNS
    -------
    integration_results (numpy.ndarray): array of integration results
        peak position, lower and upper limit of integration, and area


    """
    # get unit converter object
    uc = get_uc(dic, data)

    # define water region; upper and lower limits
    water_limits = (water_shift+water_width/2, water_shift-water_width/2)

    # integrate all identified peaks
    integration_results = np.zeros(4)
    for n, peak in enumerate(peaks):

        # get peak position and convert to ppm
        pps = uc.ppm(peak[0])

        # tuple of defined integration limits
        limits = (pps+integration_limit, pps-integration_limit)

        # integration with nmrglue function
        area = ng.analysis.integration.integrate(data.real, uc, limits)

        # append integration results to int_results array if not in water region
        if water_limits[0] > pps or pps > water_limits[1]:
            integration_result = [pps, limits[0], limits[1], area[0]]
            integration_results = np.vstack(
                [integration_results, integration_result])

    return integration_results


# Get product information
def get_products():
    """
    Get a dictionary of eCO2RR product information


    Returns
    -------
    products (pandas.DataFrame): array of products
        product name, chemical shift in ppm, multiplicity, number of 
        hydrogen atoms contributing to peak, and number of transferred 
        electrons

    """

    products = pd.DataFrame(
        [["methanol", 3.23, "s", 3, 8], ["formate", 8.33, "s", 1, 2]])

    return products


# Peak assignment to products
def product_assignment(integration_results, products, threshold=0.02):
    """
    Assign product names to identified peaks


    PARAMETERS
    ----------
    integration_results (numpy.ndarray): array of integration results
    products (pandas.DataFrame): array of products
    threshold (float): deviation from expected chemical shift of product

    RETURNS
    -------
    identified_products (pandas.DataFrame): array of products
        product name, expected chemical shift in ppm, multiplicity, number of 
        hydrogen atoms contributing to peak, number of transferred electrons,
        observed chemical shift, and peak area     

    """
    # construct new data frame
    identified_products = products.copy()

    # assign identified peaks to products
    for n, peak in enumerate(integration_results[:, 0]):
        for m, product in enumerate(products.loc[:, 1]):
            if product-threshold < peak < product+threshold:
                # add peak position to data frame
                identified_products.loc[m, 6] = peak

                # add area (integration result) to data frame
                identified_products.loc[m, 7] = integration_results[n, 3]

    return identified_products


# Product quantification
def product_quantification(products, integration_results, concentration, standard="DMSO"):
    """
    Calculation of product concentrations based on internal standard


    PARAMETERS
    ----------
    products (pandas.DataFrame): array of identified products with 
        integration results
    concentration (float): concentration of internal standard
    standard (str): name of internal standard 

    RETURNS
    -------
    products (pandas.DataFrame): array of products
        product name, expected chemical shift in ppm, multiplicity, number of 
        hydrogen atoms contributing to peak, number of transferred electrons,
        observed chemical shift, peak area, and concentration based on standard   

    """

    # get number of hydrogen atoms and chemical shift of standard
    standards = get_standards()
    shift_standard = standards[standard][0]
    atoms_standard = standards[standard][1]

    # get area of internal standard
    for n, peak in enumerate(integration_results[:, 0]):
        if shift_standard == round(peak, 2):
            area_standard = integration_results[n, 3]

    # calculate concentrations and add to data frame
    for n, product in products.iterrows():
        product_concentration = (
            concentration/area_standard)*(atoms_standard/product[3])*product[7]
        products.loc[n, 8] = product_concentration

    return products


# Calculate Faradaic efficiencies
def product_faradaic_efficiency(products, charge, volume):
    """
    Calculate Faradaic efficienies of identified products from electrocatalytic data

    Note: concentration of products in µmol/L


    PARAMETERS
    ----------
    products (pandas.DataFrame): array of identified products
    charge (float): total charge passed during experiment in A*s
    volume (float): volume of electrolyte used in mL

    RETURNS
    -------
    products (pandas.DataFrame): array of products
        product name, expected chemical shift in ppm, multiplicity, number of 
        hydrogen atoms contributing to peak, number of transferred electrons,
        observed chemical shift, peak area, concentration based on standard, 
        and FEs of products 

    """

    # define constants
    F = 96485.33212  # (C/mol; A*s/mol)

    for n, product in products.iterrows():
        # calculation of product mols from concentration and volume
        mol = product[8]*(volume*10e-3)

        # convert µmol to mol and calculate charge with Faradaic constant
        product_charge = (mol*10e-6)*product[4]*F

        # calculate Faradaic efficiency from total charge passed
        FE = abs(product_charge/charge)

        # add FE to product data frame
        products.loc[n, 9] = FE

    return products
