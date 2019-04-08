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

def vic_no_mig(params, ns, pts):
    """
    Split into two populations, no migration. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented 
    by nuA*(1-s).
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ti: Time in the past of split (in units of 2*Na generations) 
    """
    start = time.time()
    nuA, nu1, nu2, Ti, s = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu1 = nuA*(1-s)
    nu2 = nuA*s

    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx)
    return(fs)
    
def vic_anc_asym_mig(params, ns, pts):
    """
    Split with asymmetric migration followed by isolation. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented 
    by nuA*(1-s).
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the ancient migration (in units of 2*Na generations).
    T2: The scaled time between the ancient migration and present.
    """
    start = time.time()
    nuA, nu1, nu2, m12, m21, T1, T2, s = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu1 = nuA*(1-s)
    nu2 = nuA*s
   
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=m12, m21=m21)
    
    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def vic_sec_contact_asym_mig(params, ns, pts):
    """
    Split with no gene flow, followed by period of asymmetrical gene flow. Populations are 
    fractions of ancient population, where population 2 is represented by nuA*(s), and 
    population 1 is represented by nuA*(1-s).
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    T1: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    T2: The scaled time between the secondary contact and present.
    """
    start = time.time()
    nuA, nu1, nu2, m12, m21, T1, T2, s = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu1 = nuA*(1-s)
    nu2 = nuA*s

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=m12, m21=m21)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)
 
def founder_nomig(params, ns, pts):
    """
    Split into two populations, with no migration. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant. 
		nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    Ti: Time in the past of split (in units of 2*Na generations) 
    """
    start = time.time()
    nuA, nu1, nu2, Ti, s = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)
 
def founder_sym(params, ns, pts):
    """
    Split into two populations, with one migration rate. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant. 
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    Ti: Time in the past of split (in units of 2*Na generations) 
    m: Migration (2*Na*m12)
    """
    start = time.time()
    nuA, nu1, nu2, m, Ti, s = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_func, m12=m, m21=m)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_asym(params, ns, pts):
    """
    Split into two populations, with two migration rates. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant. 
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    Ti: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    start = time.time()
    nuA, nu1, nu2, m12, m21, Ti, s = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_func, m12=m12, m21=m21)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def vic_no_mig_admix_early(params, ns, pts):
    """
    Split into two populations, no migration but a discrete admixture event from pop 1 into
    pop 2 occurs. Populations are fractions of ancient population, where population 2 is 
    represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ti: Time in the past of split (in units of 2*Na generations) 
    f: Fraction of updated population 2 to be derived from population 1. 
    """
    start = time.time()
    nuA, nu1, nu2, Ti, s, f = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    nu1 = nuA*(1-s)
    nu2 = nuA*s

    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def vic_no_mig_admix_late(params, ns, pts):
    """
    Split into two populations, no migration but a discrete admixture event from pop 1 into
    pop 2 occurs. Populations are fractions of ancient population, where population 2 is 
    represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Ti: Time in the past of split (in units of 2*Na generations) 
    f: Fraction of updated population 2 to be derived from population 1. 
    """
    start = time.time()
    nuA, nu1, nu2, Ti, s, f = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu1 = nuA*(1-s)
    nu2 = nuA*s

    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def vic_two_epoch_admix(params, ns, pts):
    """
    Split with no gene flow, followed by no migration but a discrete admixture 
    event from pop 1 into pop 2 occurs. Populations are fractions of ancient population, where population 2 is 
    represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T1: The scaled time between the split and admixture event (in units of 2*Na generations).
    T2: The scaled time between the admixture event and present.
    f: Fraction of updated population 2 to be derived from population 1. 
    """
    start = time.time()
    nuA, nu1, nu2, T1, T2, s, f = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu1 = nuA*(1-s)
    nu2 = nuA*s

    phi = Integration.two_pops(phi, xx, T1, nu1, nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)

    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_early(params, ns, pts):
    """
    Split into two populations, with no migration. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant.
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    Ti: Time in the past of split (in units of 2*Na generations) 
    f: Fraction of updated population 2 to be derived from population 1.
    """
    start = time.time()
    nuA, nu1, nu2, Ti, s, f = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_late(params, ns, pts):
    """
    Split into two populations, with no migration. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant.
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    Ti: Time in the past of split (in units of 2*Na generations) 
    f: Fraction of updated population 2 to be derived from population 1.
    """
    start = time.time()
    nuA, nu1, nu2, Ti, s, f = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_two_epoch(params, ns, pts):
    """
    Split into two populations, with no migration. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant.
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T1: Time in the past of split (in units of 2*Na generations)
    T2: The scaled time between the admixture event and present.
    f: Fraction of updated population 2 to be derived from population 1.
    """
    start = time.time()
    nuA, nu1, nu2, T1, T2, s, f = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)


##########################same founder models, but with logistic growth

def founder_nomig_logistic(params, ns, pts):
    """
    Split into two populations, with no migration. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant. 
		nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    K2: Carrying capacity for pop 2
    r2: Logistic growth rate for pop 2
    nu1: Final size of pop 1.
    Ti: Time in the past of split (in units of 2*Na generations) 
    """
    start = time.time()
    nuA, nu1, nu2, Ti, s, K2, r2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: (K2*(nuA*s)*math.exp(r2*t))/(K2 + (nuA*s)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)
 
def founder_sym_logistic(params, ns, pts):
    """
    Split into two populations, with one migration rate. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant. 
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    K2: Carrying capacity for pop 2
    r2: Logistic growth rate for pop 2
    nu1: Final size of pop 1.
    Ti: Time in the past of split (in units of 2*Na generations) 
    m: Migration (2*Na*m12)
    """
    start = time.time()
    nuA, nu1, nu2, m, Ti, s, K2, r2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: (K2*(nuA*s)*math.exp(r2*t))/(K2 + (nuA*s)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_func, m12=m, m21=m)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_asym_logistic(params, ns, pts):
    """
    Split into two populations, with two migration rates. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant. 
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    K2: Carrying capacity for pop 2
    r2: Logistic growth rate for pop 2
    nu1: Final size of pop 1.
    Ti: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    start = time.time()
    nuA, nu1, nu2, m12, m21, Ti, s, K2, r2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: (K2*(nuA*s)*math.exp(r2*t))/(K2 + (nuA*s)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_func, m12=m12, m21=m21)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_early_logistic(params, ns, pts):
    """
    Split into two populations, with no migration. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant.
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    K2: Carrying capacity for pop 2
    r2: Logistic growth rate for pop 2
    nu1: Final size of pop 1.
    Ti: Time in the past of split (in units of 2*Na generations) 
    f: Fraction of updated population 2 to be derived from population 1.
    """
    start = time.time()
    nuA, nu1, nu2, Ti, s, f, K2, r2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: (K2*(nuA*s)*math.exp(r2*t))/(K2 + (nuA*s)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_late_logistic(params, ns, pts):
    """
    Split into two populations, with no migration. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant.
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    K2: Carrying capacity for pop 2
    r2: Logistic growth rate for pop 2
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    Ti: Time in the past of split (in units of 2*Na generations) 
    f: Fraction of updated population 2 to be derived from population 1.
    """
    start = time.time()
    nuA, nu1, nu2, Ti, s, f, K2, r2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: (K2*(nuA*s)*math.exp(r2*t))/(K2 + (nuA*s)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_two_epoch_logistic(params, ns, pts):
    """
    Split into two populations, with no migration. Populations are fractions of ancient
    population, where population 2 is represented by nuA*(s), and population 1 is represented by nuA*(1-s).
    Population two undergoes an exponential growth event, while population one is constant.
    nuA: Ancient population size
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T1: Time in the past of split (in units of 2*Na generations)
    T2: The scaled time between the admixture event and present.
    f: Fraction of updated population 2 to be derived from population 1.
    """
    start = time.time()
    nuA, nu1, nu2, T1, T2, s, f, K2, r2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T1)
    nu2_func = lambda t: (K2*(nuA*s)*math.exp(r2*t))/(K2 + (nuA*s)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    nu2_func = lambda t: (K2*(nuA*s)*math.exp(r2*(t+T1)))/(K2 + (nuA*s)*(math.exp(r2*(t+T1)) - 1))
    
    phi = Integration.two_pops(phi, xx, T2, nu1, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)



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
