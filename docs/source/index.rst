.. Automated Product Analysis Routine (APAR) documentation master file, created by
   sphinx-quickstart on Wed May  3 15:44:58 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2

   index

Welcome to APAR's documentation!
=====================================================================

**APAR** is an Automated Product Analysis Routine written in Python 
and designed to automate the determination of concentrations 
and Faradaic efficiencies of electrochemical CO\ :sub:`2` \ reduction 
reaction (eCO\ :sub:`2`\ RR) products from nuclear magnetic resonance 
(NMR) raw data.


Functionalities
---------------

The main functions are related to six consecutive steps described briefly:
(**1**) a series of pre-processing steps are performed (e.g., apodization, zero-filling, phase correction, etc.) to increase the sing-to-noise ratio (SNR),
(**2**) chemical shift referencing based on detection of userdefined internal standard,
(**3**) peak identification with adjustable threshold based on SNR and peak integration,
(**4**) product assignment using chemical shift positions and coupling constants reported in Preikschas, P. et al. *under review* (2023),
(**5**) calculation of product concentrations based on an internal standard and
(**6**) Faradaic efficiencies of (eCO\ :sub:`2`\ RR) products.


Installation
------------

APAR requires the following open-source Python packages:
(**1**) `Nmrglue <https://www.nmrglue.com/>`_ for reading NMR raw data and processing,
(**2**) `NumPy <https://numpy.org/>`_ and (**3**) `SciPy <https://scipy.org/>`_ for handling multidimensional arrays and mathematical operations,
(**4**) `pandas <https://pandas.pydata.org/>`_ for data export, and 
(**5**) `matplotlib <https://matplotlib.org/>`_ for data visualization (optional).

The current version and required packages can be installed directly from GitHub.

.. code-block:: python
   :linenos:

   git clone https://github.com/philpreikschas/apar
   pip install -r requirements.txt


Usage
-----

Raw NMR data can be analyzed with APAR's built-in functions following the above described six consecutive steps after importing APAR's core module.

.. code-block:: python
   :linenos:

   import apar.core as apar


1. Raw data import and processing

.. code-block:: python
   :linenos:

   # Load NMR data
   raw_dic, raw_data = apar.load_data("example-data/example-spectrum/1", "bruker")

   # Convert data to NMRPipe format
   converted_dic, converted_data = apar.convert_data(raw_dic, raw_data, "bruker")

   # Pre-processing raw data
   processed_dic, processed_data = apar.process_data(converted_dic, converted_data)


2. Identification of internal standard and chemical shift referencing

.. code-block:: python
   :linenos:

   # Identification of internal standard: find peak position
   real_shift, expected_shift = apar.find_standard(processed_dic, processed_data)


   # Chemical shift referencing
   referenced_data = apar.referencing(
         processed_dic, processed_data, real_shift, expected_shift)

3. Peak identification based on SNR and integration

.. code-block:: python
   :linenos:

   # Get all peaks from shift-corrected data
   peaks = apar.peak_identification(processed_dic, referenced_data, SNR=8)

   # Peak integration
   integration_results = apar.peak_integration(processed_dic, referenced_data, peaks)

4. Product assignment based on tabular data

.. code-block:: python
   :linenos:

   # Get product information
   products = apar.get_products()

   # Product assignment
   integration_results = apar.peak_integration(
         processed_dic, referenced_data, peaks)

5. Product quantification

.. code-block:: python
   :linenos:

   # Product quantification
   identified_products = apar.product_quantification(identified_products,
         integration_results, concentration=50, standard="DMSO")

6. Calculation of Faradaic efficienies

.. code-block:: python
   :linenos:

   # Calculate Faradaic efficiencies
   identified_products = apar.product_faradaic_efficiency(
         identified_products, charge=-434.272, volume=40)

An example analysis scheme and NMR data are provided on GitHub.

Main functions
--------------

load_data()
~~~~~~~~~~~~~~~~~~~~~~~
Load raw data from directory. 

.. note:: See docs of nmrglue for more options.

.. code-block:: python

   load_data(path, format)


**Parameters**

* path (str): path to directory with raw data
* format (str): format of raw data, resp., instrument used supported formats: bruker

**Returns**

* dic (dictionary): dictionary of NMR parameters
* data (numpy.ndarray): array of NMR data

convert_data()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Convert data to NMRPipe format.

.. code-block:: python

   convert_data(dic, data, format="bruker")

**Parameters**

* dic (dictionary): dictionary of NMR parameters
* data (numpy.ndarray): array of NMR data
* format (str): format of raw data, resp., instrument used; supported formats: bruker


**Returns**

* converted_dic (dictionary): dictionary of NMR parameters
* converted_data (numpy.ndarray): array of NMR data in NMRPipe format


process_data()
~~~~~~~~~~~~~~~~~~~~~~~
Pre-processing of NMR raw data.

.. note:: optimized for 1H NMR recorded on a BRUKER AVANCE III HD spectrometer equipped with 5 mm BBO Prodigy CryoProbe and perfect echo W5 WATERGATE solvent suppression pulse sequence.

.. code-block:: python

   process_data(dic, data)

**Parameters**

* dic (dictionary): dictionary of NMR parameters in NMRPipe format
* data (numpy.ndarray): array of NMR data in NMRPipe format

**Returns**

* processed_dic (dictionary): dictionary of NMR parameters
* processed_data (numpy.ndarray): array of processed NMR data in NMRPipe format


get_standards()
~~~~~~~~~~~~~~~
Get a dictionary of internal standards with chemical shift positions.

.. code-block:: python

   get_standards()


**Return**

* standards (dictionary): dictionary of internal standards (abbreviation, chemical shift in ppm, number of hydrogen atoms)

find_standard()
~~~~~~~~~~~~~~~
Find peak position of internal standards within limits.

.. code-block:: python

   find_standard(dic, data, limit=0.05, threshold=1e5, standard="DMSO")

**Parameters**

* dic (dictionary): dictionary of NMR parameters in NMRPipe format
* data (numpy.ndarray): array of processed NMR data
* limit (float): limits for peak identification region in ppm
* threshold (float): intensity threshold for automated peak picking
* standard (str): name of internal standard 


**Returns**

* real_shift (numpy.float64): peak position of standard in ppm
* expected_shift (numpy.float64): expected peak posistion of standard in ppm

referencing()
~~~~~~~~~~~~~
Chemical shift referencing based on resonance of internal standard.

.. code-block:: python

   referencing(dic, data, real_shift, expected_shift)

**Parameters**

* dic (dictionary): dictionary of NMR parameters in NMRPipe format
* data (numpy.ndarray): array of NMR data in NMRPipe format
* real_shift (numpy.float64): peak position of standard in ppm
* expected_shift (numpy.float64): expected peak posistion of standard in ppm


**Return**

* referenced_data (numpy.ndarray): array of NMR data referenced on resonance of internal standard


peak_identification()
~~~~~~~~~~~~~~~~~~~~~
Peak identification with given SNR as intensity threshold.

.. code-block:: python

   peak_identification(dic, data, SNR=4)

**Parameters**

* dic (dictionary): dictionary of NMR parameters in NMRPipe format
* data (numpy.ndarray): array of NMR data in NMRPipe format
* SNR (float): signal-to-noise ratio used as intensity threshold

**Return**

* peaks (numpy.recarray): array of identified peaks
      consists of peak position (X_AXIS), cluster number (cID), estimated peak scales (linewidths, X_LW), and estimated peak amplitudes (VOL)

peak_integration()
~~~~~~~~~~~~~~~~~~
Peak integration of identified peaks with fixed integration limits.

.. code-block:: python

   peak_integration(dic, data, peaks, water_shift=4.7, water_width=0.2, integration_limit=0.05)

**Parameters**

* dic (dictionary): dictionary of NMR parameters in NMRPipe format
* data (numpy.ndarray): array of NMR data in NMRPipe format
* peaks (numpy.recarray): array of identified peaks
* water_shift (float): expected peak position of water in ppm
* water_width (float): expected width of water peak in ppm
* integration_limit (float): fixed integration limits in ppm


**Return**

* integration_results (numpy.ndarray): array of integration results
        peak position, lower and upper limit of integration, and area

get_products()        
~~~~~~~~~~~~~~
Get a dictionary of eCO2RR product information

.. code-block:: python

   get_products()  

**Return**

* products (pandas.DataFrame): array of products
      product name, chemical shift in ppm, multiplicity, number of 
      hydrogen atoms contributing to peak, and number of transferred 
      electrons

product_assignment()
~~~~~~~~~~~~~~~~~~~~
Assign product names to identified peaks.

.. code-block:: python

   product_assignment(integration_results, products, threshold=0.02)

**Parameters**

* integration_results (numpy.ndarray): array of integration results
* products (pandas.DataFrame): array of products
* threshold (float): deviation from expected chemical shift of product


**Return**

* identified_products (pandas.DataFrame): array of products
      product name, expected chemical shift in ppm, multiplicity, number of 
      hydrogen atoms contributing to peak, number of transferred electrons,
      observed chemical shift, and peak area     

product_quantification()
~~~~~~~~~~~~~~~~~~~~~~~~
Calculation of product concentrations based on internal standard.

.. code-block:: python

   product_quantification(products, integration_results, concentration, standard="DMSO")

**Parameters**

* products (pandas.DataFrame): array of identified products with 
      integration results
* concentration (float): concentration of internal standard
* standard (str): name of internal standard 


**Return**

* products (pandas.DataFrame): array of products
      product name, expected chemical shift in ppm, multiplicity, number of 
      hydrogen atoms contributing to peak, number of transferred electrons,
      observed chemical shift, peak area, and concentration based on standard   

product_faradaic_efficiency()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Calculate Faradaic efficienies of identified products from electrocatalytic data.

.. note:: Concentrations of products in µmol/L.

.. code-block:: python

   product_faradaic_efficiency(products, charge, volume)

**Parameters**

* products (pandas.DataFrame): array of identified products
* charge (float): total charge passed during experiment in A*s
* volume (float): volume of electrolyte used in mL


**Return**

* products (pandas.DataFrame): array of products
      product name, expected chemical shift in ppm, multiplicity, number of 
      hydrogen atoms contributing to peak, number of transferred electrons,
      observed chemical shift, peak area, concentration based on standard, 
      and FEs of products