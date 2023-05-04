# Automated Product Analysis Routine (APAR)
APAR is an Automated Product Analysis Routine written in Python and designed to automate the determination of concentrations and Faradaic efficiencies of electrochemical CO2 reduction reaction (eCO2RR) products from nuclear magnetic resonance (NMR) raw data.

# Functionalities
The main functions are related to six consecutive steps described briefly: (1) a series of pre-processing steps are performed (e.g., apodization, zero-filling, phase correction, etc.) to increase the sing-to-noise ratio (SNR), (2) chemical shift referencing based on detection of userdefined internal standard, (3) peak identification with adjustable threshold based on SNR and peak integration, (4) product assignment using chemical shift positions and coupling constants reported in Preikschas, P. et al. under review (2023), (5) calculation of product concentrations based on an internal standard and (6) Faradaic efficiencies of (eCO2RR) products.
