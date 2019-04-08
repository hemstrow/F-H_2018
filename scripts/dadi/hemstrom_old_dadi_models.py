# Models for three pops
# define the models. Some of these should be run both with the order of GUA and HAW flipped!
def p3_cgrowth((nu2B, nu3B, nu2F, nu3F, ts, tp, m21, m31, m32, m23), (n1,n2,n3), pts):
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

def p3_lgrowth_both((nu2B, nu3B, K2, K3, ts, tp, m21, m31, m32, m23, r2, r3), (n1,n2,n3), pts):
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

def p3_lgrowth_p3((nu2B, nu3B, nu2F, K3, ts, tp, m21, m31, m32, m23, r3), (n1,n2,n3), pts):
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

def p3_lgrowth_p2((nu2B, nu3B, nu3F, K2, ts, tp, m21, m31, m32, m23, r2), (n1,n2,n3), pts):
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

def p3_lgb_p1tb_2f((nu2B, nu3B, K2, K3, ts, tp, m21, m31, m32, m23, r2, r3), (n1,n2,n3), pts):
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


# Models for two pops. Initial bottlneck in one, then growth going forward. A second bottleneck at the split for pop two, then growth.
def ib_cgrowth((nu1B, nu2B, nu1f, nu2f, ts, tp, m12, m21), (n1,n2), pts):
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

def ib_lgrowth_both((nu1B, nu2B, K1, K2, ts, tp, m12, m21, r1, r2), (n1,n2), pts):
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

def ib_lgrowth_1((nu1B, nu2B, K1, nu2f, ts, tp, m12, m21, r1), (n1,n2), pts):
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

def ib_lgrowth_2((nu1B, nu2B, nu1f, K2, ts, tp, m12, m21, r2), (n1,n2), pts):
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

# Two pops, no initial bottleneck or growth. Just an initial pop and then a split w bottleneck in pop 2.
def nig_cgrowth((nu1B, nu2B, nu1f, nu2f, tp, m12, m21), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then exponential growth in each.
    
    nu1B: p1 pop size
    nu2B: p2 pop size after founding.
    nu1f: Final pop size for p1.
    nu2f: Final p2 pop size.
    tp: Time from p2 founding to present.
    m21: Migration into p2 from p1.
    m12: Migration into p1 from p2.
    
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    
    start = time.time()
    xx = yy = dadi.Numerics.default_grid(pts)
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    # growth functions:
    nu1_func = lambda t: nu1B*(nu1f/nu1B)**(t/tp)
    nu2_func = lambda t: nu2B*(nu2f/nu2B)**(t/tp)
    
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    end = time.time()
    print("Iter time: " + str(end - start))
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)

def nig_lgrowth_both((nu1B, nu2B, K1, K2, tp, m12, m21, r1, r2), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then logistic growth in each.
    
    
    nu1B: p1 pop size
    nu2B: p2 pop size after founding.
    K1: Carrying capacity in pop 1.
    K2: Carrying capacity in pop 2.
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
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    # growth functions:
    nu1_func = lambda t: (K1*nu1B*math.exp(r1*t))/(K1 + nu1B*(math.exp(r1*t) - 1))
    nu2_func = lambda t: (K2*nu2B*math.exp(r2*t))/(K2 + nu2B*(math.exp(r2*t) - 1))
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)

def nig_lgrowth_1((nu1B, nu2B, K1, nu2f, tp, m12, m21, r1), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then logistic growth in only pop 1.
    
    
    nu1B: p1 pop size
    nu2B: p2 pop size after founding.
    K1: Carrying capacity in pop 1.
    nu2f: final pop 2 pop size
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
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # growth functions:
    nu1_func = lambda t: (K1*nu1B*math.exp(r1*t))/(K1 + nu1B*(math.exp(r1*t) - 1))
    nu2_func = lambda t: nu2B*(nu2f/nu2B)**(t/tp)
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)

def nig_lgrowth_2((nu1B, nu2B, nu1f, K2, tp, m12, m21, r2), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then logistic growth in only pop 1.
    
    
    nu1B: p1 pop size
    nu2B: p2 pop size after founding.
    nu1f: final pop 1 pop size
    K2: Carrying capacity in pop 1.
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
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # growth functions:
    nu1_func = lambda t: nu1B*(nu1f/nu1B)**(t/tp)
    nu2_func = lambda t: (K2*nu2B*math.exp(r2*t))/(K2 + nu2B*(math.exp(r2*t) - 1))
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)

# Two pops, intial population growth followed by divergence with a bottleneck.
def onegrow_cgrowth((nu1ag, nu2B, nu1f, nu2f, ts, tp, m12, m21), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then exponential growth in each.
    
    nu1ag: p1 pop size after initial growth period
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
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Pop 1 growth over time.
    phi = dadi.Integration.one_pop(phi, xx, ts, nu=nu1ag)
    
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    # growth functions:
    nu1_func = lambda t: nu1ag*(nu1f/nu1ag)**(t/tp)
    nu2_func = lambda t: nu2B*(nu2f/nu2B)**(t/tp)
    
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    end = time.time()
    print("Iter time: " + str(end - start))
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)

def onegrow_lgrowth_both((nu1ag, nu2B, K1, K2, ts, tp, m12, m21, r1, r2), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then logistic growth in each.
    
    
    nu1ag: p1 pop size after initial growth period
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
    # Pop 1 growth over time.
    phi = dadi.Integration.one_pop(phi, xx, ts, nu=nu1ag)
    
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    
    # growth functions:
    nu1_func = lambda t: (K1*nu1ag*math.exp(r1*t))/(K1 + nu1ag*(math.exp(r1*t) - 1))
    nu2_func = lambda t: (K2*nu2B*math.exp(r2*t))/(K2 + nu2B*(math.exp(r2*t) - 1))
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)

def onegrow_lgrowth_1((nu1ag, nu2B, K1, nu2f, ts, tp, m12, m21, r1), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then logistic growth in only pop 1.
    
    
    nu1ag: p1 pop size after initial growth period
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
    # Pop 1 growth over time.
    phi = dadi.Integration.one_pop(phi, xx, ts, nu=nu1ag)
    
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # growth functions:
    nu1_func = lambda t: (K1*nu1ag*math.exp(r1*t))/(K1 + nu1ag*(math.exp(r1*t) - 1))
    nu2_func = lambda t: nu2B*(nu2f/nu2B)**(t/tp)
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)

def onegrow_lgrowth_2((nu1ag, nu2B, nu1f, K2, ts, tp, m12, m21, r2), (n1,n2), pts):
    """
    Models growth in pop 1 until a split, then logistic growth in only pop 1.
    
    
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
    
    
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Pop 1 growth over time.
    phi = dadi.Integration.one_pop(phi, xx, ts, nu=nu1ag)
    
    
    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # growth functions:
    nu1_func = lambda t: nu1ag*(nu1f/nu1ag)**(t/tp)
    nu2_func = lambda t: (K2*nu2B*math.exp(r2*t))/(K2 + nu2B*(math.exp(r2*t) - 1))
    phi = dadi.Integration.two_pops(phi, xx, tp, nu1=nu1_func, nu2=nu2_func, 
                                    m12=m12, m21=m21)
    
    # print("Finishing=======================================================\n")
    
    end = time.time()
    print("Iter time: " + str(end - start))
    
    # Finally, calculate the spectrum. n1, n2, and n3 are the sample sizes to take from the populations. xx, yy, and zz are the grid sizes.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return(sfs)
