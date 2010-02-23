import numpy as np

def find_nsplines(n_knots, order):
    return n_knots-order-1

class BSplineGenerator(object):
    
    def __init__(self, knotset):
        self.knotset = knotset
        self.knots = self.knotset.knots_as_array()
        self._n_knots = len(self.knots)
    
    def find_left(self, x):
        if x == self.knots[0]:
            left = self.knots.searchsorted(x+1.0e-9)-1
        else:
            left = self.knots.searchsorted(x)-1
        return left
    
    def n_splines(self, order):
        return find_nsplines(self.n_knots, order)
    
    @property
    def n_knots(self):
        return self._n_knots
    
    def single_bspline(self, x, order, sparse=False):
        _n_splines = find_nsplines(self.n_knots, order)
        left = self.find_left(x)
        b = np.zeros(self.n_knots-1, dtype='float64')
        b[left] = 1.0
        for n in range(1, order+1):
            imin = max(0, left-n)
            imax = min(left, self.n_knots-n-2)
            for i in range(imin, imax+1):
                bi = b[i]
                bi1 = b[i+1]
                ti = self.knots[i]
                tin = self.knots[i+n]
                tin1 = self.knots[i+n+1]
                ti1 = self.knots[i+1]
                if bi == 0.0:
                    result = 0.0
                else:
                    result = bi*(x-ti)/(tin-ti)
                if bi1 != 0.0:
                    result += bi1*(tin1-x)/(tin1-ti1)
                b[i] = result
        if sparse:
            return b[imin:imax+1], _n_splines, imin, imax
        else:
            return b[0:_n_splines]

    def multiple_bsplines(self, x, order, sparse=False):
        _n_splines = find_nsplines(self.n_knots, order)
        left = self.find_left(x)
        b = np.zeros((order+1,self.n_knots-1), dtype='float64')
        b[0,left] = 1.0
        for n in range(1, order+1):
            imin = max(0, left-n)
            imax = min(left, self.n_knots-n-2)
            for i in range(imin, imax+1):
                bi = b[n-1,i]
                bi1 = b[n-1,i+1]
                ti = self.knots[i]
                tin = self.knots[i+n]
                tin1 = self.knots[i+n+1]
                ti1 = self.knots[i+1]
                if bi == 0.0:
                    result = 0.0
                else:
                    result = bi*(x-ti)/(tin-ti)
                if bi1 != 0.0:
                    result += bi1*(tin1-x)/(tin1-ti1)
                b[n][i] = result
        if sparse:
            return b[:,imin:imax+1], _n_splines, imin, imax
        else:
            return b[:,0:_n_splines]
    
    def plot_bsplines(self, order, npoints=100):
        import pylab
        xmin = self.knots[0]
        xmax = self.knots[-1]
        _n_splines = find_nsplines(self.n_knots, order)
        xvalues = np.linspace(xmin, xmax, npoints)
        yvalues = np.zeros((npoints, _n_splines), dtype='float64')
        for i, x in enumerate(xvalues):
            b = self.single_bspline(x, order)
            yvalues[i,:] = b
        return pylab.plot(xvalues, yvalues)

    def deriv_iter(self, n, b):
        print "deriv_inter.n:", n
        db = np.zeros_like(b)
        _nsplines = find_nsplines(self.n_knots, n)
        # for i in range(_nsplines-1):
        for i in range(_nsplines):
            bi = b[i]
            bi1 = b[i+1]
            tin = self.knots[i+n]
            ti = self.knots[i]
            ti1 = self.knots[i+1]
            tin1 = self.knots[i+n+1]
            if bi == 0.0:
                result = 0.0
            else:
                result = n*bi/(tin-ti)
            if bi1 != 0.0:
                result -= n*bi1/(tin1-ti1)
            db[i] = result
            # db[i] = n*b[i]/(t[i+n]-t[i]) - n*b[i+1]/(t[i+n+1]-t[i+1])
        # Sometimes, the last spline only has the first term
        # i = _nsplines-1
        # bi = b[i]
        # tin = self.knots[i+n]
        # ti = self.knots[i]
        # if bi != 0.0:
        #     db[i] = n*bi/(tin-ti)
        return db
    
    def bspline_deriv(self, x, order, deriv):
        init_order = order-deriv
        print "init_order:", init_order
        b = self.single_bspline(x, init_order)
        for d in range(1, deriv+1):
            b = self.deriv_iter(init_order+d, b)
        return b[:-deriv]
    
    def plot_dbsplines(self, order, deriv, npoints=100):
        import pylab
        xmin = self.knots[0]
        xmax = self.knots[-1]
        _n_splines = find_nsplines(self.n_knots, order)
        xvalues = np.linspace(xmin, xmax, npoints)
        yvalues = np.zeros((npoints, _n_splines), dtype='float64')
        for i, x in enumerate(xvalues):
            b = self.bspline_deriv(x, order, deriv)
            yvalues[i,:] = b
        pylab.plot(xvalues, yvalues)
        return xvalues, yvalues