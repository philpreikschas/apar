"""
v0.1

things to do:
- list of internal standards in separate function
- calculate relative areas based on compound list
- solvent peak picking: add case of multiple hits
- calculate peak threshold based on SNR
- remove debugs
"""

from pickle import TRUE
import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt

# PP: add pandas for csv export
import pandas as pd

# PP: suppress warnings
import warnings
warnings.filterwarnings("ignore")

""" 
USER INPUT 
"""
# read data
bdic, bdata = ng.bruker.read('real-nmr-data/PP22043C/1')

# internal standard for shift referencing and quantification
IS = 'DMSO'


# remove the digital filter
bdata = ng.bruker.remove_digital_filter(bdic, bdata)

# process the spectrum
# bdata = ng.proc_base.zf_size(bdata, 32768)          # zero fill to 32768 points
""" bdata = ng.proc_base.zf(bdata)                      # zero fill by padding with zeros. """

# bdata = ng.proc_base.em(bdata, lb=0.01, inv=False, rev=False)             # Apodization
# bdata = ng.proc_base.gmb(bdata)             # Apodization
# bdata = ng.proc_base.em(bdata, lb=0.00021)             # Apodization
# bdic, bdata = ng.pipe_proc.em(bdic, bdata)    # Apodization

""" bdata = ng.proc_base.fft(bdata)                     # Fourier transform """
# bdata = ng.proc_base.ps(bdata, p0=-50.0)          # phase correction
""" bdata = ng.proc_autophase.autops(bdata, fn="acme")  # automatic linear phase correction
bdata = ng.proc_base.di(bdata)                      # discard the imaginaries
bdata = ng.proc_base.rev(bdata)                     # reverse the data """



C = ng.convert.converter()
C.from_bruker(bdic, bdata)
dic, data = C.to_pipe()
ng.pipe.write("test.fid", dic, data, overwrite=TRUE)

data = ng.proc_base.zf(data)                        # zero fill by padding with zeros.
# bdata = ng.proc_base.zf_size(bdata, 32768)          # zero fill to 32768 points
# dic, data = ng.pipe_proc.apod(dic, data, qName='EM')        # Apodization

adic, adata = ng.pipe_proc.em(dic, data, lb=0.5)

data = ng.proc_base.fft(data)                       # Fourier transform
data = ng.proc_autophase.autops(data, fn="acme")    # automatic linear phase correction
data = ng.proc_base.rev(data)                       # reverse the data


adata = ng.proc_base.fft(adata)                       # Fourier transform
adata = ng.proc_autophase.autops(adata, fn="acme")    # automatic linear phase correction
adata = ng.proc_base.rev(adata)                       # reverse the data
adata = ng.proc_base.smo(adata, 2)                       
uc_adata = ng.pipe.make_uc(adic, adata)               

uc = ng.pipe.make_uc(dic, data)


# get limits in ppm
l_ppm, r_ppm = uc.ppm_limits()

""" 
DEBUG: determine limits
"""
# print("DEBUG: limits")
# print(l_ppm)
# print(r_ppm)
# print(type(l_ppm))

# PP: internal standard identification
def findIS(data, limit=0.05, IS=IS):
    """ 
    find peak position of internal standard within limits
    """

    # list of internal standards and their expected chemical shift
    ISs = {'DMSO': 2.60}

    # intensity threshold for automated peak picking
    threshold = 1e5

    # slice array to region of IS
    sdata = data.copy()
    sdata[0:uc.i(ISs[IS]+limit, 'ppm')], sdata[uc.i(ISs[IS]-limit, 'ppm'):] = 0, 0
    
    # detect all peaks with a threshold
    s_pps = ng.peakpick.pick(sdata, pthres=threshold, algorithm="connected")

    return uc.ppm(s_pps[0]["X_AXIS"])   

s_pps = findIS(adata)

""" 
DEBUG: print IS chemical shift
"""
# print("DEBUG: IS chemical shift")
# print(s_pps)
# print(type(s_pps))

# PP: chemical shift referencing
def referencing(data, s_pps, IS='DMSO'):
    """ 
    chemical shift referencing based on resonance of internal standard
    """

    # list of internal standards and their expected chemical shift
    ISs = {'DMSO': 2.60}

    # calculate difference of chemical shifts and convert to pts
    shift = uc.i((l_ppm - (s_pps - ISs[IS])), 'ppm')

    # perform frequency shift
    rdata = ng.process.proc_base.roll(data, pts=shift)

    return rdata, shift

rdata, shift = referencing(adata, s_pps, IS='DMSO')
# uc2 = ng.pipe.make_uc(rdic, rdata)

# PP: helper functions
def noise(data, limits=[l_ppm, l_ppm-1]):
    """ 
    calculate noise from area wihtout products
    """

    std = np.std(data[uc.i(limits[0], 'ppm'):uc.i(limits[1], 'ppm')].real)

    return std

# PP: peak identification and integration
def identify_peaks(data, SNR=4):
    """ 
    peak identification with given SNR as threshold
    """

    # set threshold based on SNR
    threshold = SNR * noise(data)

    # detect all peaks with a threshold
    peaks = ng.peakpick.pick(rdata, pthres=threshold, algorithm="connected")

    return peaks

# PP: add to parameter lists
SNR_int = 8

# PP: get all peaks from shift-corrected data
peaks = identify_peaks(rdata, SNR=SNR_int)



# integrate all identified peaks
int_results = np.zeros(4)
for n, peak in enumerate(peaks):
    
    # define water region & threshold
    pps_water = 4.70
    threshold_water = 0.1

    # get peak position
    pps = uc.ppm(peak["X_AXIS"])

    # define integration limits
    limits = (pps+0.05, pps-0.05)

    # integration with nmrglue function
    area = ng.analysis.integration.integrate(rdata.real, uc, limits)

    # append integration results to int_results array if not in water region
    if pps_water-threshold_water > pps or pps > pps_water+threshold_water:
        int_result = [pps, limits[0], limits[1], area[0]]
        int_results = np.vstack([int_results, int_result])

""" 
DEBUG: integration results
"""
print("DEBUG: integration results")
print(int_results)

# plot and indicate all peaks
fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
""" ax.plot(uc.ppm_scale(), data.real) """

# ax.plot(uc.ppm_scale(), rdata.real)
ax.plot(uc_adata.ppm_scale(), rdata.real, linewidth=2.0, color="#00A7E9")

# add markers for peak positions
for n, peak in enumerate(peaks):

    height = rdata[int(peak["X_AXIS"])]
    ppm = uc.ppm(peak["X_AXIS"])
    if ppm < 4.6 or ppm > 4.8:
        ax.scatter(ppm, height, marker="o", color="#FA5649", s=100, alpha=0.5)
        ax.text(ppm, height + 5e5, round(ppm, 2), ha="center", va="center")

ax.hlines(SNR_int * noise(rdata), *uc.ppm_limits(), linestyle="--", color="#3B3838")
ax.axvspan(4.6, 4.8, color='#F37E27', alpha=0.5)
ymin, ymax = ax.get_ylim()
ax.text(4.70, ymax, "water region", ha="center", va="bottom")
ax.text(12-0.05, (SNR_int * noise(rdata))+(SNR_int * noise(rdata) * 0.1), "integration threshold", ha="left", va="bottom")

# ax.set_ylim(top=1.2e7)
# ax.set_xlim(200, 0)

plt.xlim(12, 0)

plt.show()

# PP: export raw data to csv
df1 = pd.DataFrame({'Chemical Shift': uc.ppm_scale(), 'Intensity': data.real, 'Intensity Corrected': rdata.real}, columns=['Chemical Shift', 'Intensity', 'Intensity Corrected'], )
df2 = pd.DataFrame({'Peaks': peaks}, columns=['Peaks'], )
df3 = pd.DataFrame(int_results[1:,0:], columns=['Peak Position', 'Lower Limit', 'Upper Limit', 'Area'])
with pd.ExcelWriter('export/export_nmr-integration.xlsx') as writer:  
    df1.to_excel(writer, sheet_name='spectrum', index=False, header=True)
    """ df2.to_excel(writer, sheet_name='peaks', index=False, header=True) """
    df3.to_excel(writer, sheet_name='integration', index=False, header=True)



""" 
DEBUG: random output
"""
print("DEBUG: random output")

std = np.std(rdata[uc.i(l_ppm, 'ppm'):uc.i(l_ppm-1, 'ppm')].real)
mean = np.average(rdata[uc.i(l_ppm, 'ppm'):uc.i(l_ppm-1, 'ppm')].real)

print(std)
print(mean)
print(std/mean)
print(l_ppm)
print(data[uc("2.6 ppm")])
print(data[0].real)
udic = ng.pipe.guess_udic(dic,data)
print(udic[0]['time'])
print(udic)