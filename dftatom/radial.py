"""
This module contains helper functions for all features, that are needed for
atomic DFT calculations. The low-level things are written in fortran in
rdirac.f, high-level things are in python.
"""

from math import pi, sqrt, log, atan

from numpy import array, zeros

def solve_radial_eigenproblem(n,l,u,r,Emin=-3000,Emax=-1000,Z=82,
        eps=1e-5,rel=0):
    """Solves the radial Schrodinger or Dirac (depending on "rel") equation.

    Finds the wavefunction with defined "n" and "l". The potential is "u".
    rel ... 0 nonrelat
            1 semirelat
            2 relat - spin up
            3 relat - spin down
           -1 original

    for rel=2,3 the grid "r" can be arbitrary
    for rel=-1,0, it must be hyperbolic
    rel=1 is not yet implemented

    Emin,Emax gives the window, where to look at the solution by halving the
      interval. if there is no solution in the interval, the interval is
      automatically adjusted, so you don't have to worry about it. The
      calculation is a little faster, if you set Emin and Emax near the correct
      energy, but it's not necessary.
    eps ... the solver halves the difference between Emin and Emax until
            |Emax-Emin|<eps
    r ... the grid on which to solve the equation, must be hyperbolic (the
        solver counts with it)

    returns (E,y), where E is the energy and y is the numpy.array, the
        wavefunction
    """
    assert rel in [-1,0,1,2,3]
    assert n>0
    assert 0<=l and l<n
    import rdirac
    return rdirac.solve_radial_eigenproblem(n,l,Emin,Emax-Emin,eps,u,r,Z,rel)

def create_log_grid(Z=82.,rmin=0.00625,rmax=20.,inc=1.0247):
    r=[rmin/Z]
    while r[-1]<rmax:
        r.append(r[-1]*inc)
    return array(r)

def create_log_grid2(rmin=1,rmax=1,N=1):
    r=[rmin]
    while r[-1]<rmax:
        n = len(r)
        r.append(rmin*(rmax/rmin)**(float(n)/N))
    return array(r)

def create_hyperbolic_grid(Z=82.,rmin=0.00625,rmax=20.,k0=500):
    """r=ap*j/(jm-j) for j=1,2,3... and k0 is going to be the len(r)"""
    R1=rmin/Z
    R2=rmax
    jm=int((R1-R2)/(R1-R2/k0)+0.5)
    ap=R1*(jm-1)
    return array([ap*j/(jm-j) for j in range(1,k0+1)])

def vxc_from_n(n,relat=2):
    if relat in [1,2]:
        r = 1
    else:
        assert relat in [-1,0]
        r = 0
    import rdirac
    return rdirac.getvxc(n,r)

    #equivalent code in python:
    if n == 0:
        return 0
    def Y(y):
        return y**2+b*y+c
    y0 = -0.10498
    b = 3.72744
    c = 12.93532
    A = 0.0621814
    Q = sqrt(4*c-b**2)
    rs = (3/(4*pi*n))**(1./3)
    y = sqrt(rs)
    e = A/2 * (log(y**2/Y(y))+2*b/Q * atan(Q/(2*y+b)) - \
        b*y0/Y(y0)*(log((y-y0)**2/Y(y))+2*(b+2*y0)/Q*atan(Q/(2*y+b))))
    Vc = e - A/6 * (c*(y-y0)-b*y0*y)/((y-y0)*Y(y))
    ex = -3/(4*pi) * (3*pi**2*n)**(1./3)

    if r:
        c = 137.036
        beta = 1./c * (3*pi**2*n)**(1./3)
        mu = sqrt(1+beta**2)
        A = (beta*mu-log(beta+mu))/beta**2
        R = 1-3./2*A**2
        dAdbeta = 2/mu - 2*A/beta
        dRdbeta = -3*A*dAdbeta
        Vx = 4*ex/3*R+ex*dRdbeta*beta/3
    else:
        Vx = 4*ex/3
    return Vx+Vc

def getVxc(density,R,relat=True):
    """Input: density on R, density = n, "R" is the grid
    Output: return Vxc potential

    """
    return array([vxc_from_n(n/(4*pi),relat) for n,r in zip(density,R)])

def integratepoisson(density,r):
    """density is "n", it integrates the poisson eq.,
    see the docstring in the rdirac.integrate_radial_poisson for more info."""
    import rdirac
    return rdirac.integrate_radial_poisson(density,r,0,0)

def getVh(density,r,Z):
    """density is "n", it integrates the poisson eq.

    we have these conditions on the solution U:
      * U[-1] == Z/r[-1]
      * derivative of U at r[0] is zero.

    This determines the solution uniquely.

    """
    U=integratepoisson(density,r)
    return U+Z/r[-1]-U[-1]

def KS_construct_density(eigenvalues,R,Z):
    """Takes KS energies and wavefunctions and constructs
    n = sum y^2
    
    returns n on the grid
    """
    Es,data = eigenvalues
    density = zeros(len(R), "d")
    for E,index in Es[:Z]:
        y = data[index][0]
        density = density + y**2
    return density

def integrate(func,r):
    """Returns an integral of a function 'func' on r.
    """
    import rdirac
    return rdirac.integrate(r,func)

def KS_solve(u,r,Z,nl=None, relat = 2, eps=1e-5):
    """Solves Dirac equation on the radial grid "r", potential "u" and
    a nucleus charge "Z". 

    nl is a list of (n,l), for which you want to calculate the energies and
    wavefunctions. If None, it is automatically determined
    from Z, so that the lowest Z eigenvalues are present (plus some more).

    relat .....
        2 ... Dirac equation, distinguish cases j=l+1/2 and j=l-1/2
        0 ... Schrodinger equation, spin independent (thus 2x less
            calculations)
        1 ... scalar relativistic, this means, that we neglect spin, thus we
            have the same amount of calculation as in relat=0, but we want to
            add correcting terms to the Schrodinger equation in order to get as
            precise result as possible, while neglecting spins.
    
    Returns (Es, data), where:
      Es ...... is a list of (E,index), E ... energy for data[index]
      data .... is a list of (y,n,l,p), y ... wavefunction for "n" and "l"
         and "p", where p=+1 (j=l+1/2) or p=-1 (j=l-1/2)
    """
    assert relat in [-1,0,1,2]
    if nl == None:
        nl=[]
        s=0
        n=1
        while s<Z:
            for l in range(n):
                nl.append((n,l))
                s+=2*(2*l+1)
            n+=1
        for l in range(n):
            nl.append((n,l))

    Es=[]
    data=[]
    for n,l in nl:
        if relat==2:
            #spin up
            fn=(2*l+1+1)
            E,y = solve_radial_eigenproblem(n,l,u,r,Z=Z,rel=2,eps=eps)
            data.append([y,n,l,1])
            Es.extend([(E,len(data)-1)]*fn)

            #spin down
            fn=(2*l-1+1)
            if l!=0:
                E,y = solve_radial_eigenproblem(n,l,u,r,Z=Z,rel=3,eps=eps)
                p = -1
            else:
                #use the result from the spin up
                p = 1
            data.append([y,n,l,p])
            Es.extend([(E,len(data)-1)]*fn)
        else:
            E,y = solve_radial_eigenproblem(n,l,u,r,Z=Z,rel=relat,eps=eps)
            fn=2*(2*l+1)
            data.append([y,n,l,1])
            Es.extend([(E,len(data)-1)]*fn)
    Es.sort(key=lambda x: x[0])
    return Es,data
