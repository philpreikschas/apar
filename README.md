![GitHub](https://img.shields.io/github/license/philpreikschas/apar)
[![documentation](https://img.shields.io/badge/docs-apar.readthedocs.io-lightgrey)](https://apar.readthedocs.io)
[![DOI article](https://img.shields.io/badge/DOI-10.1038/s42004--023--00948--9-red)](https://dx.doi.org/10.1038/s42004-023-00948-9)
[![DOI code](https://img.shields.io/badge/DOI_Code-10.5281/zenodo.8070371-blue)](https://zenodo.org/badge/latestdoi/587696258)

# Automated Product Analysis Routine (APAR)
APAR is an Automated Product Analysis Routine written in Python and designed to automate the determination of concentrations and Faradaic efficiencies of electrochemical CO2 reduction reaction (eCO2RR) products from nuclear magnetic resonance (NMR) raw data.

# Functionalities
The main functions are related to six consecutive steps described briefly: (1) a series of pre-processing steps are performed (e.g., apodization, zero-filling, phase correction, etc.) to increase the sing-to-noise ratio (SNR), (2) chemical shift referencing based on detection of user-defined internal standard, (3) peak identification with adjustable threshold based on SNR and peak integration, (4) product assignment using chemical shift positions and coupling constants reported in Preikschas, P. et al. under review (2023), (5) calculation of product concentrations based on an internal standard and (6) Faradaic efficiencies of (eCO2RR) products.

# Citing APAR
**NMR-based quantification of liquid products in CO2 electroreduction on phosphate-derived nickel catalysts**.

P. Preikschas, A.J. Martín, B.S. Yeo, and J. Pérez-Ramírez, [_Commun. Chem._ **2023**, in press](https://dx.doi.org/10.1038/s42004-023-00948-9).
