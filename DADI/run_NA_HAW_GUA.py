import numpy
import dadi
import os
import array
import pylab
import matplotlib
os.chdir("../DADI")

# define the model

def three_seq_estab((nu2B, nu3B, nu2F, nu3F, ts, tp, m21, m31, m32, m23), (n1,n2,n3), pts):
    """
    Model with growth, split, bottleneck in pop2, exp recovery, migration
    
    nu2B: HAW pop size after founding.
    nu3B: GUA pop size after founding
    nu2F: final HAW pop size.
    nu3F: final GUA pop size.
    ts: Time from start (HAW split) to GUA split.
    tp: Time from GUA founding to present.
    m21: Migration into HAW from NA.
    m31: Migration into GUA from NA.
    m32: Migration into GUA from HAW.
    m23: Migration in HAW from GUA.
    
    n1,n2,n3: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    # Define the grid we'll use
    xx = yy = zz = dadi.Numerics.default_grid(pts)
    
    # print("Initializing pop.\n")
    # phi for the equilibrium ancestral population (NA)
    phi = dadi.PhiManip.phi_1D(xx)
    
    # print("Diverging HAW.\n")
    # NA to HAW divergence.
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # print("Growing HAW with migration.\n")
    # We need to define a function to describe the non-constant population 2
    # size. lambda is a convenient way to do so. Lambda just makes a function with parameter t, always just does an expression. In this case, given a time (t), what is the pop size?
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/(ts+tp))
  
    # Now we move forward in time until the GUA split (ts) with migration between the two.
    phi = dadi.Integration.two_pops(phi, xx, ts, nu2=nu2_func, m21=m21)
    
    # print("Splitting GUA.")
    # Now split off GUA from HAW
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    # print("Growning HAW and GUA with migration.\n")
    # Grow the two split populations until present
    # Need an equation for nu3 and to redefine nu2, since we are starting at a new time period.
    nu3_func = lambda t: nu3B*(nu3F/nu3B)**(t/tp)
    nu2_func2 = lambda t: nu2B*(nu2F/nu2B)**((ts+t)/(ts+tp))
    
    # Move forward in time till present.
    phi = dadi.Integration.three_pops(phi, xx, tp, nu2=nu2_func2, nu3=nu3_func, m21=m21, m31=m31, m32=m32, m23=m23)
    
    # print("Finishing.\n")
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,yy,zz))
    return(sfs)

# import the data

dd = dadi.Misc.make_data_dict('dadi_snps_10kgap.txt')

# Extract the spectrum for NAM, GUA, and HAW from that dictionary,
# projected down to 20 samples per population.
data= dadi.Spectrum.from_data_dict(dd, ['NAM','GUA','HAW'], [20,20,20], polarized=True)
ns = data.sample_sizes

# These are the grid point settings will use for extrapolation.
pts_l = [40,50,60]

# Now let's optimize parameters for this model.

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.
# Parameters are: (nu2B, nu3B, nu2F, nu3F, ts, tp, m21, m31, m32, m23)
upper_bound = [10, 10, 10, 10, 20, 20, 10, 10, 10, 10]
lower_bound = [1e-3, 1e-3, 1e-3, 1e-3, 0, 0, 0, 0, 0, 0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [0.1,0.1,2,1,2,1,0.2,0.1,0.2,0.1]
# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(three_seq_estab)

# Perturb our parameters before optimization. This does so by taking each
# parameter a up to a factor of two up or down.
p0 = dadi.Misc.perturb_params(p0, fold = 1, upper_bound=upper_bound,
                              lower_bound=lower_bound)
# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters.
# The maxiter argument restricts how long the optimizer will run. For real 
# runs, you will want to set this value higher (at least 10), to encourage
# better convergence. You will also want to run optimization several times
# using multiple sets of intial parameters, to be confident you've actually
# found the true maximum likelihood parameters.


print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=20, maxiter=50)
print('Finshed optimization **************************************************')

# plot and calculate theta0
model = func_ex(popt, ns, pts_l)

# theta 0
t0 = dadi.Inference.optimal_sfs_scaling(model, data)

#print results
print("Unfolded optimal parms:\n")
print(popt)
print("\ntheta:\n")
print(t0)
print("\n")

###############################################
# folded
dataf= dadi.Spectrum.from_data_dict(dd, ['NAM','GUA','HAW'], polarized=False)
nsf = data.sample_sizes

print('Beginning optimization folded************************************************')
popt_f = dadi.Inference.optimize_log(p0, dataf, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=20, maxiter=50)
print('Beginning optimization folded************************************************')

modelf = func_ex(popt_f, nsf, pts_l)

# theta 0
t0f = dadi.Inference.optimal_sfs_scaling(modelf, data)

print("Folded optimal parms:\n")
print(popt_f)
print("\ntheta:\n")
print(t0f)
print("\n")


