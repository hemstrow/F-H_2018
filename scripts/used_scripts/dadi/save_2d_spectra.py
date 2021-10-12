from __future__ import division
import sys

if len(sys.argv) != 7:
  print("Usage: save_2d_spectra.py\tparameters\tmodel\tinfile\tpops\tprojection\tpolarized")
  print("This script generates and saves site frequency spectra for a given model for both observed data and estimated demographic parameters. The file adjusted_2d_portik_models.py containing models must be located in the directory from which this script is run.\n")
  print("Arguments:\n")
  print("\tparameters: demographic parameters, in the order expected by the model, as [parm1, parm2, ...]\n")
  print("\tmodel: name of the selected model, matching that in adjusted_2p_portik_models.py\n")
  print("\tinfile: path to the dadi formated input file for the observed data\n")
  print("\tpops: names of the populations to use, as [A, B]\n")
  print("\tprojection: projection to use, as [10, 10]\n")
  print("\tpolarized: boolean (True or False), specifies if the spectra should be polarized.\n")
  print("Provided", len(sys.argv), "arguments, terminating.\n")
  sys.exit()


import time
import numpy
import scipy
import matplotlib
import dadi
import os
import array
import math
import adjusted_2p_portik_models
import matplotlib.pyplot as pyplot

parms = sys.argv[1]
selected_mod = sys.argv[2]
datafile = sys.argv[3]
pops = sys.argv[4]
projection = sys.argv[5]
polarized = sys.argv[6]

# import the parameters
parms_mapped = list(map(float, parms.strip('[]').split(",")))

# import the projection
projection = list(map(int, projection.strip('[]').split(",")))

# import pops
pops = pops.strip('[]').split(',')

# specify number of grid points
pts_l = [200,220,240]

# Make the extrapolating version of our demographic model function of choice.
func_ex = dadi.Numerics.make_extrap_log_func(eval('adjusted_2p_portik_models' + '.' + selected_mod))

# grab the number of samples in the original spectra
## import the data
dd = dadi.Misc.make_data_dict(datafile)
if polarized == "True":
  data = dadi.Spectrum.from_data_dict(dd, pops, projection, polarized=True)
else:
  data = dadi.Spectrum.from_data_dict(dd, pops, projection, polarized=False)
## get number of samples
ns = data.sample_sizes


# generate the extrapolated sfs given our parameters
model = func_ex(numpy.float64(parms_mapped), ns, pts_l)

# plot
# pyplot.figure()
# dadi.Plotting.plot_2d_comp_Poisson(model, data)

# generate residuals between observed and estimated spectra
masked_model, masked_data = dadi.Numerics.intersect_masks(model, data)
min_toplot = min(masked_model.min(), masked_data.min())
resid = dadi.Inference.linear_Poisson_residual(model, data, mask=min_toplot)

# save both spectra and residual
model.to_file(fname = "modeled_spectra.txt")
data.to_file(fname = "obs_spectra.txt")
resid.to_file(fname = "residual.txt")
