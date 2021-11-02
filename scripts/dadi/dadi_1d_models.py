from __future__ import division
import numpy
import time
import math
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum


#############Neutral/instantanious change##############
def neutral(notused, ns, pts):
    """
    Standard neutral model. No params. pop size is theta
    """
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

def two_epoch(params, ns, pts):
    """
    Instantaneous size change some time ago.

    nu: Size after change.
    T: Time in the past at which size change happened (in units of 2*Na 
       generations) 
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,T = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    
    phi = Integration.one_pop(phi, xx, T, nu)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

#############exponential growth########################

def bottleneck_exp_growth(params, ns, pts):
    """
    Models a bottleneck at Tb, followed by recovery to present 
    nuB: bottleneck size
    nuF: final size
    Tb: Time since start of bottleneck
    Tg: Time since start of growth
    Note, ancestral size is theta
    """
    
    nuB,nuF,Tb,Tg = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    
    # bottleneck
    Tb = Tb - Tg
    phi = Integration.one_pop(phi, xx, Tb, nuB)
    
    # recovery
    nu_func = lambda t: nuB * (nuF/nuB)**(t/Tg)
    phi = Integration.one_pop(phi, xx, Tg, nu_func)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

def instant_bottleneck_exp_growth(params, ns, pts):
    """
    Models a bottleneck at Tb, followed by recovery to present 
    nuB: bottleneck size
    nuF: final size
    Tg: Time since bottleneck and growth start
    Note, ancestral size is theta
    """
    
    nuB,nuF,Tg = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    
    # recovery
    nu_func = lambda t: nuB * (nuF/nuB)**(t/Tg)
    phi = Integration.one_pop(phi, xx, Tg, nu_func)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

def one_phase_exp_growth(params, ns, pts):
    """
    Models exponential growth/decline starting at T until present
    nuA: ancestral pop size
    nuF: current pop size
    T: time since start of growth
    """

    nuA,nuF,T = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    nu_func = lambda t: nuA * (nuF/nuA)**(t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs
    
def two_phase_exp_growth(params, ns, pts):
    """
    Models exponential growth/decline starting at Tg, then a different growth/decline starting at Td
    Often, growth then decline or vice versa. Be careful with bounds on decline!
    nuA: ancestral pop size
    nuG: max size after phase 1
    nuF: current pop size
    Tg: time since start of phase 1
    Td: time since start of phase 2
    """

    nuA,nuG,nuF,Tg,Td = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    
    # first phase, start at nuA, run till nuG over Tg
    nuS = nuA
    nuE = nuG
    T = Tg - Td
    nu_func = lambda t: nuS * (nuE/nuS)**(t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)

    # second phase, start at nuG
    nuS = nuG
    nuE = nuF
    T = Td
    nu_func = lambda t: nuS * (nuE/nuS)**(t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)

    
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

#############logistic growth##########################

def bottleneck_logistic_growth(params, ns, pts):
    """
    Models a bottleneck at Tb, followed by recovery to present 
    nuB: bottleneck size
    K: carrying capacity
    Tb: Time since start of bottleneck
    Tg: Time since start of growth
    r: growth rate
    Note, ancestral size is theta
    """
    
    nuB,K,Tb,Tg,r = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    
    # bottleneck
    Tb = Tb - Tg
    phi = Integration.one_pop(phi, xx, Tb, nuB)
    
    # recovery
    nu_logistic_func = lambda t: (K*(nuB)*math.exp(r*t))/(K + (nuB)*(math.exp(r*t) - 1))
    phi = Integration.one_pop(phi, xx, Tg, nu_logistic_func)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

def instant_bottleneck_logistic_growth(params, ns, pts):
    """
    Models a bottleneck at Tb, followed by recovery to present 
    nuB: bottleneck size
    K: carrying capacity
    Tg: Time since bottleneck and growth start
    r: growth rate
    Note, ancestral size is theta
    """
    
    nuB,K,Tg,r = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    
    # recovery
    nu_logistic_func = lambda t: (K*(nuB)*math.exp(r*t))/(K + (nuB)*(math.exp(r*t) - 1))
    phi = Integration.one_pop(phi, xx, Tg, nu_logistic_func)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

def one_phase_logistic_growth(params, ns, pts):
    """
    Models exponential growth/decline starting at T until present
    nuA: ancestral pop size
    K: carrying capacity
    T: time since start of growth
    r: growth rate
    """

    nuA,K,T,r = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    nu_logistic_func = lambda t: (K*(nuA)*math.exp(r*t))/(K + (nuA)*(math.exp(r*t) - 1))
    phi = Integration.one_pop(phi, xx, T, nu_func)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs
    
def two_phase_logistic_growth_vary_all(params, ns, pts):
    """
    Models exponential growth/decline starting at Tg, then a different growth/decline starting at Td
    r can differ between periods
    Often, growth then decline or vice versa. Be careful with bounds on decline!
    nuA: ancestral pop size
    K1: Carrying capacity for phase 1
    K2: Carrying capacity for phase 2
    Tg: time since start of phase 1
    Td: time since start of phase 2
    r1: growth rate, phase 1
    r2: growth rate, phase 2
    """

    nuA,K1,K2,Tg,Td,r1,r2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    
    # first phase, start at nuA, run till nuG over Tg
    nuS = nuA
    T = Tg - Td
    nu_logistic_func = lambda t: (K1*(nuS)*math.exp(r1*t))/(K1 + (nuS)*(math.exp(r1*t) - 1))
    phi = Integration.one_pop(phi, xx, T, nu_logistic_func)

    # second phase, start at ending point of first period
    nuS = nu_logistic_func(T)
    T = Td
    nu_logistic_func = lambda t: (K2*(nuS)*math.exp(r2*t))/(K2 + (nuS)*(math.exp(r2*t) - 1))
    phi = Integration.one_pop(phi, xx, T, nu_logistic_func)

    
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

def two_phase_logistic_growth_fixed_r(params, ns, pts):
    """
    Models logistic growth/decline starting at Tg, then a different growth/decline starting at Td.
    r is the same in both
    Often, growth then decline or vice versa. Be careful with bounds on decline!
    nuA: ancestral pop size
    K1: Carrying capacity for phase 1
    K2: Carrying capacity for phase 2
    Tg: time since start of phase 1
    Td: time since start of phase 2
    r: growth rate
    """

    nuA,K1,K2,Tg,Td,r = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    
    # first phase, start at nuA, run till nuG over Tg
    nuS = nuA
    T = Tg - Td
    nu_logistic_func = lambda t: (K1*(nuS)*math.exp(r*t))/(K1 + (nuS)*(math.exp(r*t) - 1))
    phi = Integration.one_pop(phi, xx, T, nu_logistic_func)

    # second phase, start at the ending pop size of the first period
    nuS = nu_logistic_func(T)
    T = Td
    nu_logistic_func = lambda t: (K2*(nuS)*math.exp(r*t))/(K2 + (nuS)*(math.exp(r*t) - 1))
    phi = Integration.one_pop(phi, xx, T, nu_logistic_func)

    
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs
