import numpy as np

# TODO
# 1. Get bspline_derivs working
# 2. Optimize derivs to work with imin, imax

def find_left(t, x):
    if x == t[0]:
        left = t.searchsorted(x+1.0e-9)-1
    else:
        left = t.searchsorted(x)-1
    return left

def find_nsplines(nknots, order):
    return nknots-order-1

def multiple_bsplines(t, order, x, sparse=False):
    nknots = t.size
    nsplines_min = nknots-order-1
    left = find_left(t, x)
    b = np.zeros((order+1,nknots-1), dtype='float64')
    b[0,left] = 1.0
    for n in range(1,order+1):
        imin = max(0,left-n)
        imax = min(left,nknots-n-2)
        for i in range(imin, imax+1):
            bi = b[n-1,i]
            bi1 = b[n-1,i+1]
            ti = t[i]
            tin = t[i+n]
            tin1 = t[i+n+1]
            ti1 = t[i+1]
            if bi == 0.0:
                result = 0.0
            else:
                result = bi*(x-ti)/(tin-ti)
            if bi1 != 0.0:
                result += bi1*(tin1-x)/(tin1-ti1)
            b[n][i] = result
    if sparse:
        return b[:,imin:imax+1], imin, imax
    else:
        return b[:,0:nsplines_min]
    
def single_bspline(t, order, x, sparse=False):
    nknots = t.size
    nsplines = nknots-order-1
    left = find_left(t, x)
    b = np.zeros(nknots-1, dtype='float64')
    b[left] = 1.0
    for n in range(1,order+1):
        imin = max(0,left-n)
        imax = min(left,nknots-n-2)
        for i in range(imin, imax+1):
            bi = b[i]
            bi1 = b[i+1]
            ti = t[i]
            tin = t[i+n]
            tin1 = t[i+n+1]
            ti1 = t[i+1]
            if bi == 0.0:
                result = 0.0
            else:
                result = bi*(x-ti)/(tin-ti)
            if bi1 != 0.0:
                result += bi1*(tin1-x)/(tin1-ti1)
            b[i] = result
    if sparse:
        return b[imin:imax+1], imin, imax
    else:
        return b[0:nsplines]

def deriv_iter(t, n, b):
    db = np.zeros_like(b)
    nsplines = find_nsplines(t.size, n)
    for i in range(nsplines-1):
        bi = b[i]
        bi1 = b[i+1]
        tin = t[i+n]
        ti = t[i]
        ti1 = t[i+1]
        tin1 = t[i+n+1]
        if bi == 0.0:
            result = 0.0
        else:
            result = n*bi/(tin-ti)
        if bi1 != 0.0:
            result -= n*bi1/(tin1-ti1)
        db[i] = result
        # db[i] = n*b[i]/(t[i+n]-t[i]) - n*b[i+1]/(t[i+n+1]-t[i+1])
    # The last spline only has the first term
    i = nsplines-1
    bi = b[i]
    tin = t[i+n]
    ti = t[i]
    if bi != 0.0:
        db[i] = n*bi/(tin-ti)
    return db

def bspline_deriv(t, order, deriv, x):
    init_order = order-deriv
    b = single_bspline(t, init_order, x)
    for d in range(1, deriv+1):
        b = deriv_iter(t, init_order+d, b)
    return b[:-deriv]

def bspline_derivs(t, order, deriv, x):
    nknots = t.size
    nsplines = find_nsplines(nknots, order)
    b = multiple_bsplines(t, order, x, sparse=False)
    for n in range(1,order):
        for new_order in range(n+1, order+1):
            b[n][:] = deriv_iter(t, new_order, b[n][:])

    # b[1][:] = deriv_iter(t, 2, b[1][:])
    # b[1][:] = deriv_iter(t, 3, b[1][:])
    # b[2][:] = deriv_iter(t, 3, b[2][:])
    return b

def plot_bsplines(t, order, npoints=100):
    import pylab
    xmin = t[0]
    xmax = t[-1]
    xvalues = np.linspace(xmin, xmax, npoints)
    yvalues = np.zeros((npoints, t.size-order-1), dtype='float64')
    for i, x in enumerate(xvalues):
        b = single_bspline(t, order, x)
        yvalues[i,:] = b
    pylab.plot(xvalues, yvalues)
    return xvalues, yvalues

def plot_dbsplines(t, order, deriv, npoints=100):
    import pylab
    xmin = t[0]
    xmax = t[-1]
    xvalues = np.linspace(xmin, xmax, npoints)
    yvalues = np.zeros((npoints, t.size-order-1), dtype='float64')
    for i, x in enumerate(xvalues):
        b = bspline_deriv(t, order, deriv, x)
        yvalues[i,:] = b
    pylab.plot(xvalues, yvalues)
    return xvalues, yvalues

