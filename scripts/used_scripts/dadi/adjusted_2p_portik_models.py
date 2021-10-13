from __future__ import division
import numpy
import time
import math
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

######################################growth subfunctions#######################################
# nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
# 
# nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
# 
# nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
# 
# nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))

#######################################vicariance models, no growth#############################

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
    nuA, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu1 = nuA*(1-s)
    nu2 = nuA*s
    
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
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
    nuA, m12, m21, T1, T2, s = params
    
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
    nuA, m12, m21, T1, T2, s = params
    
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
    nuA, Ti, s, f = params
    
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
    nuA, Ti, s, f = params
    
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
    nuA, T1, T2, s, f = params
    
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

##########################founder models, exponential growth###################################

#########growth pop 2################
def founder_nomig_growth_pop_2(params, ns, pts):
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
    nuA, nu2, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_exp_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_sym_growth_pop_2(params, ns, pts):
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
    nuA, nu2, m, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_exp_func, m12=m, m21=m)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_asym_growth_pop_2(params, ns, pts):
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
    nuA, nu2, m12, m21, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_exp_func, m12=m12, m21=m21)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_early_growth_pop_2(params, ns, pts):
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
    nuA, nu2, Ti, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_exp_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_late_growth_pop_2(params, ns, pts):
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
    nuA, nu2, Ti, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_exp_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_two_epoch_growth_pop_2(params, ns, pts):
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
    nuA, nu2, T1, T2, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/(T1+T2))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2_exp_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    
    nu2_0 = nu2_exp_func(T1)
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**((t+T1)/(T1+T2))
    phi = Integration.two_pops(phi, xx, T2, nu1, nu2_exp_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)


#########growth pop 1################
def founder_nomig_growth_pop_1(params, ns, pts):
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
    nuA, nu1, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_exp_func, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_sym_growth_pop_1(params, ns, pts):
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
    nuA, nu1, m, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_exp_func, nu2, m12=m, m21=m)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_asym_growth_pop_1(params, ns, pts):
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
    nuA, nu1, m12, m21, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_exp_func, nu2, m12=m12, m21=m21)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_early_growth_pop_1(params, ns, pts):
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
    nuA, nu1, Ti, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_exp_func, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_late_growth_pop_1(params, ns, pts):
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
    nuA, nu1, Ti, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_exp_func, nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_two_epoch_growth_pop_1(params, ns, pts):
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
    nuA, nu1, T1, T2, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/(T1 + T2))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, T1, nu1_exp_func, nu1, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    nu1_0 = nu1_exp_func(T1)
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**((t+T1)/(T1 + T2))
    phi = Integration.two_pops(phi, xx, T2, nu1_exp_func, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

##########growth both pops############
def founder_nomig_growth_both(params, ns, pts):
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
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_exp_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_exp_func, nu2_exp_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_sym_growth_both(params, ns, pts):
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
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_exp_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_exp_func, nu2_exp_func, m12=m, m21=m)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_asym_growth_both(params, ns, pts):
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
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_exp_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_exp_func, nu2_exp_func, m12=m12, m21=m21)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_early_growth_both(params, ns, pts):
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
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_exp_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_exp_func, nu2_exp_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_late_growth_both(params, ns, pts):
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
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/Ti)
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/Ti)
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_exp_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_exp_func, nu2_exp_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_two_epoch_growth_both(params, ns, pts):
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
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/(T1 + T2))
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**(t/(T1 + T2))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_exp_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, T1, nu1_exp_func, nu2_exp_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    
    nu1_0 = nu1_exp_func(T1)
    nu2_0 = nu1_exp_func(T1)
    nu1_exp_func = lambda t: nu1_0 * (nu1/nu1_0)**((t+T1)/(T1 + T2))
    nu2_exp_func= lambda t: nu2_0 * (nu2/nu2_0)**((t+T1)/(T1 + T2))
    phi = Integration.two_pops(phi, xx, T2, nu1_exp_func, nu2_exp_func, m12=0, m21=0)

    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)



##########################same founder models, but with logistic growth
#########growth pop 2################
def founder_nomig_logistic_pop_2(params, ns, pts):
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
    nuA, K2, r2, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_logistic_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_sym_logistic_pop_2(params, ns, pts):
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
    nuA, K2, r2, m, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_logistic_func, m12=m, m21=m)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_asym_logistic_pop_2(params, ns, pts):
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
    nuA,  K2, r2, m12, m21, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_logistic_func, m12=m12, m21=m21)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_early_logistic_pop_2(params, ns, pts):
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
    nuA, K2, r2, Ti, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_logistic_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_late_logistic_pop_2(params, ns, pts):
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
    nuA, K2, r2, Ti, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1, nu2_logistic_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_two_epoch_logistic_pop_2(params, ns, pts):
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
    nuA, K2, r2, T1, T2, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1 = nuA*(1-s)
    nu2_0 = nuA*s
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, T1, nu1, nu2_logistic_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    
    nu2_0 = nu2_logistic_func(T1)
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    phi = Integration.two_pops(phi, xx, T2, nu1, nu2_logistic_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)


#########growth pop 1################
def founder_nomig_logistic_pop_1(params, ns, pts):
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
    nuA, K1, r1, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_logistic_func, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_sym_logistic_pop_1(params, ns, pts):
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
    nuA, K1, r1, m, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_logistic_func, nu2, m12=m, m21=m)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_asym_logistic_pop_1(params, ns, pts):
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
    nuA, K1, r1, m12, m21, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_logistic_func, nu2, m12=m12, m21=m21)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_early_logistic_pop_1(params, ns, pts):
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
    nuA, K1, r1, Ti, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_logistic_func, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_late_logistic_pop_1(params, ns, pts):
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
    nuA, K1, r1, Ti, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_logistic_func, nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_two_epoch_logistic_pop_1(params, ns, pts):
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
    nuA, K1, r1, T1, T2, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, T1, nu1_logistic_func, nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    nu1_0 = nu1_logistic_func(T1)
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    phi = Integration.two_pops(phi, xx, T2, nu1_logistic_func, nu2, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

##########growth both pops############
def founder_nomig_logistic_both(params, ns, pts):
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
    nuA, K1, K2, r1, r2, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_logistic_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_logistic_func, nu2_logistic_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_sym_logistic_both(params, ns, pts):
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
    nuA, K1, K2, r1, r2, m, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_logistic_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_logistic_func, nu2_logistic_func, m12=m, m21=m)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_asym_logistic_both(params, ns, pts):
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
    nuA, K1, K2, r1, r2, m12, m21, Ti, s = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_logistic_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_logistic_func, nu2_logistic_func, m12=m12, m21=m21)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_early_logistic_both(params, ns, pts):
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
    nuA, K1, K2, r1, r2, Ti, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_logistic_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_logistic_func, nu2_logistic_func, m12=0, m21=0)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_late_logistic_both(params, ns, pts):
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
    nuA, K1, K2, r1, r2, Ti, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_logistic_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, Ti, nu1_logistic_func, nu2_logistic_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_nomig_admix_two_epoch_logistic_both(params, ns, pts):
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
    nuA, K1, K2, r1, r2, T1, T2, s, f = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx, nu=nuA)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    nu1_0 = nuA*(1-s)
    nu2_0 = nuA*s
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    #note, the nu2_0 can be eliminated and the function can appear as:
    #nu2_logistic_func = lambda t: (nuA*(1-s)) * (nu2/(nuA*(1-s)))**(t/Ti)
    phi = Integration.two_pops(phi, xx, T1, nu1_logistic_func, nu2_logistic_func, m12=0, m21=0)
    phi = PhiManip.phi_2D_admix_1_into_2(phi, f, xx,xx)
    
    
    nu1_0 = nu1_logistic_func(T1)
    nu2_0 = nu1_logistic_func(T1)
    nu1_logistic_func = lambda t: (K1*(nu1_0)*math.exp(r1*t))/(K1 + (nu1_0)*(math.exp(r1*t) - 1))
    nu2_logistic_func = lambda t: (K2*(nu2_0)*math.exp(r2*t))/(K2 + (nu2_0)*(math.exp(r2*t) - 1))
    phi = Integration.two_pops(phi, xx, T2, nu1_logistic_func, nu2_logistic_func, m12=0, m21=0)

    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

#######################founder models, but with historic growth prior to split######################
def founder_asym_hist_igrowth_p2(params, ns, pts):
    """

    
    Population undergoes historic growth some time prior to split, p2  grows after split, p1 stays
    constant
    
    nuA: Ancient population size
    nuG: Size after historic growth
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu2F: Final size of pop 2.
    Tg: Time between the historic growth event and split.
    Ts: Time between split and second growth egent (in units of 2*Na generations)
    Tg2: Time to between second growth event and present
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    start = time.time()
    nuA, nuG, nu2F, m12, m21, Tg, Ts, Tg2, s = params
    
    # initialize
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx, nu=nuA)
    
    # first growth, runs for Tg
    T = Tg
    phi = Integration.one_pop(phi, xx, T, nuG)
    
    # split
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu1 = nuG*(1-s)
    nu2 = nuG*s
    
    # time passes, runs for Ts
    T = Ts
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    
    # pop 2 grows, run the rest of the time
    T = Tg2
    phi = Integration.two_pops(phi, xx, T, nu1, nu2F, m12=m12, m21=m21)

    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)

def founder_asym_hist_3epoch_exp_growth_p1(params, ns, pts):
    """
   
    Population undergoes historic growth some time prior to split, then later some exponential
    growth, then splits, p2 grows after split, p1 grows/shrinks again.

    nuA: Ancient population size
    nuG: Size after historic growth
    s: Fraction of nuA that goes to pop2. (Pop 1 has size nuA*(1-s).)
    nu2F: Final size of pop 2.
    nu1F: Final size of pop 1.
    nuG2: Size of pop 1 at the split.
    Tg: Time between historic growth and second growth phase
    Tg2: Time between second growth phase start and split
    Ts: Time between split and instant growth of p2
    Tg3: Time between instant growth of p2 and present
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    """
    start = time.time()
    nuA, nuG, nu1F, nuG2, nu2F, m12, m21, Tg, Tg2, Ts, Tg3, s = params
    
    # initialize
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx, nu=nuA)
    
    # first growth
    T = Tg
    phi = Integration.one_pop(phi, xx, T, nuG)
    
    # second growth
    T = Tg2
    nu1_exp_func = lambda t: nuG * (nuG2/nuG)**(t/T)
    phi = Integration.one_pop(phi, xx, T, nu1_exp_func)
    
    # split
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    nu1 = nuG2*(1-s)
    nu2 = nuG2*s
    
    # time passes, p1 starts final growth phase, p2 is constant.
    T = Ts
    T_exp_func = Ts + Tg3 # since we are only doing the start of the growth, we set the full growth period as T for the growth function
    nu1_exp_func = lambda t: nu1 * (nu1F/nu1)**(t/T_exp_func)
    phi = Integration.two_pops(phi, xx, T, nu1_exp_func, nu2, m12=m12, m21=m21)
    
    # pop 2 grows, run the rest of the time
    T = Ts + Tg3
    init_T = Ts
    phi = Integration.two_pops(phi, xx, T, nu1_exp_func, nu2F, m12=m12, m21=m21, initial_t=init_T)

    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return(fs)