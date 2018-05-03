import sys

if len(sys.argv) != 14:
  print("Usage: run_3d_mods.py\tpath_to_dir_with_data\tinfile\toutfile\tmodel\tmaxiters\tpops\tfold_spectra\tfold_parms\tinitial_parms\tupper_bounds\tlower_bounds\tprojection\toptimizer\n")
  print("This script runs one of several 3 pop demographic models in dadi given the correct SNP format input data.\n")
  print("Arguments:")
  print("\tpath_to_dir_with_data: The path to the directory containing input data.")
  print("\tinfile: Input file name.")
  print("\toutfile: Name for output file, contained in the directory with the input data this script is run.")
  print("\tmodel: Name of the model to run. Options:")
  print("\t\tcgrowth\n\t\tlgrowth_both\n\t\tlgrowth_p3\n\t\tlgrowth_p2\n\t\tlgb_p1tb_2f\n\t\t2p_cgrowth\n\t\t2p_lgrowth_both\n\t\t2p_lgrowth_1\n\t\t2p_lgrowth_2")
  print("\t\tNote: all of the 2p_* models are for two populations only!")
  print("\tmaxiters: Maximum number of iterations overwhich to optimize.")
  print("\tpops: list of populations, in the correct order for the model.")
  print("\tfold_spectra: Should the spectra be folded? True or False.")
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
import dadi
import os
import array
import math
import matplotlib

os.chdir(sys.argv[1])

# Models for three pops
# define the models. Some of these should be run both with the order of GUA and HAW flipped!
def cgrowth((nu2B, nu3B, nu2F, nu3F, ts, tp, m21, m31, m32, m23), (n1,n2,n3), pts):
    """
    Models p2 establishment followed by p3 establishment. Migration between p3 and p2 as well as to both of those from p1.
    Assumes constant, exponential growth in both p2 and p3. Assumes equal growth rates in both areas.
    
    nu2B: p2 pop size after founding.
    nu3B: p3 pop size after founding
    nu2F: final p2 pop size.
    nu3F: final p3 pop size.
    ts: Time from start (p2 split) to p3 split.
    tp: Time from p3 founding to present.
    m21: Migration into p2 from p1.
    m31: Migration into p3 from p1.
    m32: Migration into p3 from p2.
    m23: Migration in p2 from p3.
    
    n1,n2,n3: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    start = time.time()
    # Define the grid we'll use
    xx = yy = zz = dadi.Numerics.default_grid(pts)
    
    # print("Initializing pop.\n")
    # phi for the equilibrium ancestral population (p1)
    phi = dadi.PhiManip.phi_1D(xx)
    
    # print("Diverging p2.\n")
    # p1 to p2 divergence.
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # print("Growing p2 with migration.\n")
    # We need to define a function to describe the non-constant population 2
    # size. lambda is a convenient way to do so. Lambda just makes a function with parameter t, always just does an expression. In this case, given a time (t), what is the pop size?
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/(ts+tp))
  
    # Now we move forward in time until the p3 split (ts) with migration between the two.
    phi = dadi.Integration.two_pops(phi, xx, ts, nu2=nu2_func, m21=m21)
    
    # print("Splitting p3.")
    # Now split off p3 from p2
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    # print("Growing p2 and p3 with migration.\n")
    # Grow the two split populations until present
    # Need an equation for nu3 and to redefine nu2, since we are starting at a new time period.
    nu3_func = lambda t: nu3B*(nu3F/nu3B)**(t/tp)
    nu2_func2 = lambda t: nu2B*(nu2F/nu2B)**((ts+t)/(ts+tp))
    
    # Move forward in time till present.
    phi = dadi.Integration.three_pops(phi, xx, tp, nu2=nu2_func2, nu3=nu3_func, m21=m21, m31=m31, m32=m32, m23=m23)
    
    # print("Finishing=======================================================\n")
    end = time.time()
    print("Iter time: " + str(end - start))
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,yy,zz))
    return(sfs)

def lgrowth_both((nu2B, nu3B, K2, K3, ts, tp, m21, m31, m32, m23, r2, r3), (n1,n2,n3), pts):
    """
    Models pop2 establishment followed by pop3 establishment. Migration between p1 and p2 as well as to both of those from p1.
    Assumes logistic growth to different K values in both p2 and p3.
    
    nu2B: p2 pop size after founding.
    nu3B: p3 pop size after founding
    k2: p2 carrying capacity.
    k3: p3 carrying capacity.
    ts: Time from start (p2 split) to p3 split.
    tp: Time from p3 founding to present.
    m21: Migration into p2 from p1.
    m31: Migration into p3 from p1.
    m32: Migration into p3 from p2.
    m23: Migration in p2 from p3.
    r2: intrinsic growth rate of p2
    r3: intrinsic growth rate of p3
    
    n1,n2,n3: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    start = time.time()
    
    # Define the grid we'll use
    xx = yy = zz = dadi.Numerics.default_grid(pts)
    
    # print("Initializing pop.\n")
    # phi for the equilibrium ancestral population (p1)
    phi = dadi.PhiManip.phi_1D(xx)
    
    # p1 to p2 divergence.
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # print("Growing p2 with migration.\n")
    # We need to define a function to describe the non-constant population 2
    # size. lambda is a convenient way to do so. Lambda just makes a function with parameter t, always just does an expression. In this case, given a time (t), what is the pop size?
    nu2_func = lambda t: (K2*nu2B*math.exp(r2*t))/(K2 + nu2B*(math.exp(r2*t) - 1))
  
    # Now we move forward in time until the p3 split (ts) with migration between the two.
    phi = dadi.Integration.two_pops(phi, xx, ts, nu2=nu2_func, m21=m21)
    
    # Now split off p3 from p2
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    # print("Growing p2 and p3 with migration.\n")
    # Grow the two split populations until present
    # Need an equation for nu3 and to redefine nu2 where the time elapsed is t plus the time from founding until the split.
    nu3_func = lambda t: (K3*nu3B*math.exp(r3*t))/(K3 + nu3B*(math.exp(r3*t) - 1))
    nu2_func2 = lambda t: (K2*nu2B*math.exp(r2*(t+ts)))/(K2 + nu2B*(math.exp(r2*(t+ts)) - 1))
    
    # Move forward in time till present.
    phi = dadi.Integration.three_pops(phi, xx, tp, nu2=nu2_func2, nu3=nu3_func, m21=m21, m31=m31, m32=m32, m23=m23)
    
    end = time.time()
    
    # print("Finishing=======================================================\n")
    print("Iter time: " + str(end - start))
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,yy,zz))
    return(sfs)

def lgrowth_p3((nu2B, nu3B, nu2F, K3, ts, tp, m21, m31, m32, m23, r3), (n1,n2,n3), pts):
    """
    Models p2 establishment followed by p3 establishment. Migration between p3 and p2 as well as to both of those from p1.
    Assumes logistic growth to different K values in both p2 and p3.
    
    nu2B: p2 pop size after founding.
    nu3B: p3 pop size after founding
    k2: p2 carrying capacity.
    k3: p3 carrying capacity.
    ts: Time from start (p2 split) to p3 split.
    tp: Time from p3 founding to present.
    m21: Migration into p2 from p1.
    m31: Migration into p3 from p1.
    m32: Migration into p3 from p2.
    m23: Migration in p2 from p3.
    r2: intrinsic growth rate of p2
    r3: intrinsic growth rate of p3
    
    n1,n2,n3: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    start = time.time()
    # Define the grid we'll use
    xx = yy = zz = dadi.Numerics.default_grid(pts)
    
    # print("Initializing pop.\n")
    # phi for the equilibrium ancestral population (p1)
    phi = dadi.PhiManip.phi_1D(xx)
    
    # print("Diverging p2.\n")
    # p1 to p2 divergence.
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # print("Growing p2 with migration.\n")
    # We need to define a function to describe the non-constant population 2
    # size. lambda is a convenient way to do so. Lambda just makes a function with parameter t, always just does an expression. In this case, given a time (t), what is the pop size?
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/(ts+tp))
  
    # Now we move forward in time until the p3 split (ts) with migration between the two.
    phi = dadi.Integration.two_pops(phi, xx, ts, nu2=nu2_func, m21=m21)
    
    # print("Splitting p3.")
    # Now split off p3 from p2
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    # print("Growning p2 and p3 with migration.\n")
    # Grow the two split populations until present
    # Need an equation for nu3 and to redefine nu2 where the time elapsed is t plus the time from founding until the split.
    nu3_func = lambda t: (K3*nu3B*math.exp(r3*t))/(K3 + nu3B*(math.exp(r3*t) - 1))
    nu2_func2 = lambda t: nu2B*(nu2F/nu2B)**((ts+t)/(ts+tp))
    
    # Move forward in time till present.
    phi = dadi.Integration.three_pops(phi, xx, tp, nu2=nu2_func2, nu3=nu3_func, m21=m21, m31=m31, m32=m32, m23=m23)
    
    # print("Finishing.\n")
    end = time.time()
    print("Iter time: " + str(end - start))
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,yy,zz))
    return(sfs)

def lgrowth_p2((nu2B, nu3B, nu3F, K2, ts, tp, m21, m31, m32, m23, r2), (n1,n2,n3), pts):
    """
    Models p2 establishment followed by p3 establishment. Migration between p3 and p2 as well as to both of those from p1.
    Assumes logistic growth to different K values in both p2 and p3.
    
    nu2B: p2 pop size after founding.
    nu3B: p3 pop size after founding
    k2: p2 carrying capacity.
    k3: p3 carrying capacity.
    ts: Time from start (p2 split) to p3 split.
    tp: Time from p3 founding to present.
    m21: Migration into p2 from p1.
    m31: Migration into p3 from p1.
    m32: Migration into p3 from p2.
    m23: Migration in p2 from p3.
    r2: intrinsic growth rate of p2
    r3: intrinsic growth rate of p3
    
    n1,n2,n3: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    start = time.time()
    # Define the grid we'll use
    xx = yy = zz = dadi.Numerics.default_grid(pts)
    
    # print("Initializing pop.\n")
    # phi for the equilibrium ancestral population (p1)
    phi = dadi.PhiManip.phi_1D(xx)
    
    # print("Diverging p2.\n")
    # p1 to p2 divergence.
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # print("Growing p2 with migration.\n")
    # We need to define a function to describe the non-constant population 2
    # size. lambda is a convenient way to do so. Lambda just makes a function with parameter t, always just does an expression. In this case, given a time (t), what is the pop size?
    nu2_func = lambda t: (K2*nu2B*math.exp(r2*t))/(K2 + nu2B*(math.exp(r2*t) - 1))
  
    # Now we move forward in time until the p3 split (ts) with migration between the two.
    phi = dadi.Integration.two_pops(phi, xx, ts, nu2=nu2_func, m21=m21)
    
    # print("Splitting p3.")
    # Now split off p3 from p2
    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    # print("Growning p2 and p3 with migration.\n")
    # Grow the two split populations until present
    # Need an equation for nu3 and to redefine nu2 where the time elapsed is t plus the time from founding until the split.
    nu3_func = lambda t: nu3B*(nu3F/nu3B)**(t/tp)
    nu2_func2 = lambda t: (K2*nu2B*math.exp(r2*(t+ts)))/(K2 + nu2B*(math.exp(r2*(t+ts)) - 1))
    
    # Move forward in time till present.
    phi = dadi.Integration.three_pops(phi, xx, tp, nu2=nu2_func2, nu3=nu3_func, m21=m21, m31=m31, m32=m32, m23=m23)
    
    # print("Finishing.\n")
    end = time.time()
    print("Iter time: " + str(end - start))
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,yy,zz))
    return(sfs)

def lgb_p1tb_2f((nu2B, nu3B, K2, K3, ts, tp, m21, m31, m32, m23, r2, r3), (n1,n2,n3), pts):
    """
    Models p2 and p3 establishment seperately from NA. p2 founded first. Migration between p3 and p2 as well as to both of those from p1.
    Assumes logistic growth to different K values in both p2 and p3.
    Reverse order of p2 and p3 to test the opposite introduction order.
    
    nu2B: p2 pop size after founding.
    nu3B: p3 pop size after founding
    k2: p2 carrying capacity.
    k3: p3 carrying capacity.
    ts1: Time from start (p2 split) to p3 split.
    ts2: Time from p3 founding to present.
    m21: Migration into p2 from p1.
    m31: Migration into p3 from p1.
    m32: Migration into p3 from p2.
    m23: Migration in p2 from p3.
    r2: intrinsic growth rate of p2
    r3: intrinsic growth rate of p3
    
    n1,n2,n3: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    start = time.time()
    # Define the grid we'll use
    xx = yy = zz = dadi.Numerics.default_grid(pts)
    
    # print("Initializing pop.\n")
    # phi for the equilibrium ancestral population (p1)
    phi = dadi.PhiManip.phi_1D(xx)
    
    # print("Diverging p2.\n")
    # p1 to p2 divergence.
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # print("Growing p2 with migration.\n")
    # We need to define a function to describe the non-constant population 2
    # size. lambda is a convenient way to do so. Lambda just makes a function with parameter t, always just does an expression. In this case, given a time (t), what is the pop size?
    nu2_func = lambda t: (K2*nu2B*math.exp(r2*t))/(K2 + nu2B*(math.exp(r2*t) - 1))
  
    # Now we move forward in time until the p3 split (ts) with migration between the two.
    phi = dadi.Integration.two_pops(phi, xx, ts, nu2=nu2_func, m21=m21)
    
    # print("Splitting p3.")
    # Now split off p3 from p2
    phi = dadi.PhiManip.phi_2D_to_3D_split_1(xx, phi)
    
    # print("Growning p2 and p3 with migration.\n")
    # Grow the two split populations until present
    # Need an equation for nu3 and to redefine nu2 where the time elapsed is t plus the time from founding until the split.
    nu3_func = lambda t: (K3*nu3B*math.exp(r3*t))/(K3 + nu3B*(math.exp(r3*t) - 1))
    nu2_func2 = lambda t: (K2*nu2B*math.exp(r2*(t+ts)))/(K2 + nu2B*(math.exp(r2*(t+ts)) - 1))
    
    # Move forward in time till present.
    phi = dadi.Integration.three_pops(phi, xx, tp, nu2=nu2_func2, nu3=nu3_func, m21=m21, m31=m31, m32=m32, m23=m23)
    
    # print("Finishing.\n")
    end = time.time()
    print("Iter time: " + str(end - start))
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2,n3), (xx,yy,zz))
    return(sfs)



# Models for two pops
def p2_cgrowth((nu1B, nu2B, nu1f, nu2f, ts, tp, m12, m21), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then exponential growth in each.
    
    nu1B: p1 pop size after founding bottleneck
    nu2B: p2 pop size after founding.
    nu1f: Final pop size for p1.
    nu2f: Final p2 pop size.
    ts: Time from start to p2 split.
    tp: Time from p2 founding to present.
    m21: Migration into p2 from p1.
    m12: Migration into p1 from p2.
    
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    start = time.time()
    xx = yy = dadi.Numerics.default_grid(pts)
    t_tot = tp + ts #total time.
    
    # phi for the equilibrium ancestral population (NA)
    phi = dadi.PhiManip.phi_1D(xx)
    
    # print("First growth...")
    # Pop 1 bottleneck followed by growth over time.
    nu1_func_pre = lambda t: nu1B*(nu1f/nu1B)**(t/t_tot)
    phi = dadi.Integration.one_pop(phi, xx, ts, nu=nu1_func_pre)
    
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # print("Second growth...")
    # growth functions:
    nu1_func = lambda t: nu1B*(nu1f/nu1B)**((t+ts)/t_tot)
    nu2_func = lambda t: nu2B*(nu2f/nu2B)**(t/tp)
    
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    end = time.time()
    print("Iter time: " + str(end - start))
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)

def p2_lgrowth_both((nu1B, nu2B, K1, K2, ts, tp, m12, m21, r1, r2), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then logistic growth in each.
    
    nu1B: p1 pop size after founding bottleneck
    nu2B: p2 pop size after founding.
    K1: Carrying capacity in pop 1.
    K2: Carrying capacity in pop 2.
    ts: Time from start to p2 split.
    tp: Time from p2 founding to present.
    m21: Migration into p2 from p1.
    m12: Migration into p1 from p2.
    r1: growth rate, pop 1
    r2: growth rate, pop 2

    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    start = time.time()
    
    # Define grid
    start = time.time()
    xx = yy = dadi.Numerics.default_grid(pts)
    
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Pop 1 bottleneck followed by growth over time.
    nu1_func_pre = lambda t: (K1*nu1B*math.exp(r1*t))/(K1 + nu1B*(math.exp(r1*t) - 1))
    phi = dadi.Integration.one_pop(phi, xx, ts, nu=nu1_func_pre)
    
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    # growth functions:
    nu1_func = lambda t: (K1*nu1B*math.exp(r1*(t+ts)))/(K1 + nu1B*(math.exp(r1*(t+ts)) - 1))
    nu2_func = lambda t: (K2*nu2B*math.exp(r2*t))/(K2 + nu2B*(math.exp(r2*t) - 1))
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)

def p2_lgrowth_1((nu1B, nu2B, K1, nu2f, ts, tp, m12, m21, r1), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then logistic growth in only pop 1.
    
    nu1B: p1 pop size after founding bottleneck
    nu2B: p2 pop size after founding.
    K1: Carrying capacity in pop 1.
    nu2f: final pop 2 pop size
    ts: Time from start to p2 split.
    tp: Time from p2 founding to present.
    m21: Migration into p2 from p1.
    m12: Migration into p1 from p2.
    r1: growth rate, pop 1
    
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    start = time.time()
    
    # Define grid
    start = time.time()
    xx = yy = dadi.Numerics.default_grid(pts)
    
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    # Pop 1 bottleneck and growth over time.
    nu1_func_pre = lambda t: (K1*nu1B*math.exp(r1*t))/(K1 + nu1B*(math.exp(r1*t) - 1))
    phi = dadi.Integration.one_pop(phi, xx, ts, nu=nu1_func_pre)
    
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # growth functions:
    nu1_func = lambda t: (K1*nu1B*math.exp(r1*(t+ts)))/(K1 + nu1B*(math.exp(r1*(t+ts)) - 1))
    nu2_func = lambda t: nu2B*(nu2f/nu2B)**(t/tp)
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)

def p2_lgrowth_2((nu1B, nu2B, nu1f, K2, ts, tp, m12, m21, r2), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then logistic growth in only pop 1.
    
    nu1B: p1 pop size after founding bottleneck
    nu1ag: p1 pop size after initial growth period
    nu2B: p2 pop size after founding.
    nu1f: final pop 1 pop size
    K2: Carrying capacity in pop 1.
    ts: Time from start to p2 split.
    tp: Time from p2 founding to present.
    m21: Migration into p2 from p1.
    m12: Migration into p1 from p2.
    r1: growth rate, pop 1
    
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    start = time.time()
    
    # Define grid
    start = time.time()
    xx = yy = dadi.Numerics.default_grid(pts)
    t_tot = tp + ts #total time.
    
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Pop 1 bottleneck and growth over time.
    nu1_func_pre = lambda t: nu1B*(nu1f/nu1B)**(t/t_tot)
    phi = dadi.Integration.one_pop(phi, xx, ts, nu=nu1_func_pre)
    
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # growth functions:
    nu1_func = lambda t: nu1B*(nu1f/nu1B)**((t+ts)/t_tot)
    nu2_func = lambda t: (K2*nu2B*math.exp(r2*t))/(K2 + nu2B*(math.exp(r2*t) - 1))
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)




# import the data
dd = dadi.Misc.make_data_dict(sys.argv[2])

# Extract the spectrum from data for the populations listed
# projected down to requested number of samples per population (150, 15, 15 for NA, GUA, HAW is good).
pops = sys.argv[6]
pops = pops.strip('[]').split(',')
projection = sys.argv[12]
projection = map(int, projection.strip('[]').split(","))
data = dadi.Spectrum.from_data_dict(dd, pops, projection, polarized=sys.argv[7])
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
func_ex = dadi.Numerics.make_extrap_log_func(eval(sys.argv[4]))

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
f.write("model:\t" + sys.argv[4] + "\tpops:\t" + ' '.join(map(str, pops)) + "\ttheta:\t" + str(t0) + "\tll:\t" + str(numpy.around(ll, 4)) + "\tAIC:\t" + str(numpy.around(AIC, 4)) + "\toptimal_parameters:\t" + ' '.join(map(str, popt)) + "\n")

print("Finished.")
