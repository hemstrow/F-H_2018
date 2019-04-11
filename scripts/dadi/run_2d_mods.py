import sys

if len(sys.argv) != 14:
  print("Usage: run_3d_mods.py\tpath_to_dir_with_data\tinfile\toutfile\tmodel\tmaxiters\tpops\tfold_spectra\tfold_parms\tinitial_parms\tupper_bounds\tlower_bounds\tprojection\toptimizer\n")
  print("This script runs one of several 3 pop demographic models in dadi given the correct SNP format input data.\n")
  print("Arguments:")
  print("\tpath_to_dir_with_data: The path to the directory containing input data.")
  print("\tinfile: Input file name.")
  print("\toutfile: Name for output file, contained in the directory from which this script is run.")
  print("\tmodel: Name of the model to run. For options, see adjusted_2p_portik_models.py")
  print("\tmaxiters: Maximum number of iterations overwhich to optimize.")
  print("\tpops: list of populations, in the correct order for the model.")
  print("\tpolarized: Is the spectra polarized. True or False.")
  print("\tfold_parms: How strongly should the initial parameters be perturbed?")
  print("\tinitial_parms: List of initial parameter values.")
  print("\tupper_bounds: List of upper bounds for parameters.")
  print("\tlower_bounds: List of lower bounds for parameters.")
  print("\tprojection: How should the populations be projected?")
  print("\toptimizer: Which optimizer should be used? Options:")
  print("\t\tlog: optimize_log.")
  print("\t\tfmin: optimize_log_fmin")
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
from __future__ import division

os.chdir(sys.argv[1])



# import the data
dd = dadi.Misc.make_data_dict(sys.argv[2])

# Extract the spectrum from data for the populations listed
# projected down to requested number of samples per population (150, 15, 15 for NA, GUA, HAW is good).
pops = sys.argv[6]
pops = pops.strip('[]').split(',')
projection = sys.argv[12]
projection = map(int, projection.strip('[]').split(","))
if sys.argv[7] == "True":
  data = dadi.Spectrum.from_data_dict(dd, pops, projection, polarized=True)
else:
  data = dadi.Spectrum.from_data_dict(dd, pops, projection, polarized=False)
ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.
pts_l = [40,50,60]

# Now let's optimize parameters for this model.

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.
upper_bound = sys.argv[10]
upper_bound = map(float, upper_bound.strip('[]').split(','))
lower_bound = sys.argv[11]
lower_bound = map(float, lower_bound.strip('[]').split(','))

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = sys.argv[9]
p0 = map(float, p0.strip('[]').split(','))

# Make the extrapolating version of our demographic model function of choice.
func_ex = dadi.Numerics.make_extrap_log_func(eval('adjusted_2p_portik_models' + '.' + sys.argv[4]))

# Perturb our parameters before optimization. This does so by taking each
# parameter a up to a factor of two up or down.
pfold = map(float, sys.argv[8])
p0 = dadi.Misc.perturb_params(p0, fold = pfold, upper_bound=upper_bound,
                              lower_bound=lower_bound)

print("\nReady to begin run. Parameters:\n")

print("Model:\n" + sys.argv[4])
print("Maxiters:\n" + sys.argv[5])
print("Pops:")
print(pops)
print("Projection:")
print(projection)
print("Number of segregating sites in projected sfs:")
print(data.S())
print("Polarized spectra:")
print(sys.argv[7])
print("Level of parm folding:\n" + sys.argv[8])
print("Optimizer:\n" + sys.argv[13])
print("p0:")
print(p0)
print("Upper Bounds:")
print(upper_bound)
print("Lower Bounds:")
print(lower_bound)



# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters.
# The maxiter argument restricts how long the optimizer will run. For real 
# runs, you will want to set this value higher (at least 10), to encourage
# better convergence. You will also want to run optimization several times
# using multiple sets of intial parameters, to be confident you've actually
# found the true maximum likelihood parameters.

print('\nBeginning optimization ************************************************')
if sys.argv[13] == "log":
  print("Optimize log.")
  popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=map(int, sys.argv[5]))
elif sys.argv[13] == "fmin":
  print("optimize_log_fmin.")
  popt = dadi.Inference.optimize_log_fmin(p0, data, func_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound, verbose=1, maxiter=map(int, sys.argv[5]))
else:
  print("Unaccepted optimizer\n")
  sys.exit()

print('\nFinshed optimization **************************************************')

# plot and calculate theta0
print("Finalizing model and calculating stats...")
model = func_ex(popt, ns, pts_l)

# theta 0
t0 = dadi.Inference.optimal_sfs_scaling(model, data)

# log likelihood
ll = dadi.Inference.ll_multinom(model, data)

# AIC
AIC = (-2*(float(ll))) + (2*len(p0))

#print results
f = open(sys.argv[3],'a')

f.write("model:\t" + sys.argv[4] + "\tpops:\t" + ' '.join(map(str, pops)) + "\ttheta:\t" + str(t0) + "\tll:\t" + str(numpy.around(ll, 4)) + "\tAIC:\t" + str(numpy.around(AIC, 4)) + "\toptimal_parameters:\t" + ' '.join(map(str, numpy.around(popt, 9))) + "\n")

print("Finished.")
