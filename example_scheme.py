import apar.core as apar

# suppress warnings
import warnings
warnings.filterwarnings("ignore")

# Load NMR data
raw_dic, raw_data = apar.load_data("example-data/example-spectrum/1", "bruker")

# Convert data to NMRPipe format
converted_dic, converted_data = apar.convert_data(raw_dic, raw_data, "bruker")


# Pre-processing raw data ##FLAG->update processing##
processed_dic, processed_data = apar.process_data(
    converted_dic, converted_data)


# Identification of internal standard: find peak position
real_shift, expected_shift = apar.find_standard(processed_dic, processed_data)


# Chemical shift referencing
referenced_data = apar.referencing(
    processed_dic, processed_data, real_shift, expected_shift)


# Get all peaks from shift-corrected data
peaks = apar.peak_identification(processed_dic, referenced_data, SNR=8)


# Peak integration
integration_results = apar.peak_integration(
    processed_dic, referenced_data, peaks)


# Get product information
products = apar.get_products()

# Product assignment
identified_products = apar.product_assignment(integration_results, products)

# Product quantification
identified_products = apar.product_quantification(identified_products,
                                                  integration_results, concentration=50, standard="DMSO")

# Calculate Faradaic efficiencies
identified_products = apar.product_faradaic_efficiency(
    identified_products, charge=-434.272, volume=40)

print("results:")
print(identified_products)
