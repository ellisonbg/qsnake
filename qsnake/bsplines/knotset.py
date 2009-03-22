import numpy as np


class KnotSetError(Exception):
    pass


class KnotSet(object):
    """
    A set of knot points that define a set of B-splines.
    
    >>> ks = KnotSet()
    """   

    def __init__(self):        
        """
        >>> ks = KnotSet()
        >>> ks.knots()
        []
        """
        self._knots = {}

    def __len__(self):
        """
        >>> ks = KnotSet()
        >>> len(ks)
        0
        >>> ks.add_knot(0.0)
        >>> len(ks)
        1
        """
        return sum(self._knots.itervalues())

    def knots_as_dict(self):
        """
        Return the knots as a dict of positions and repititions.
        
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.knots_as_dict()
        {0.0: 1, 1.0: 1}
        """
        return self._knots

    def iter_unique_knots(self):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.repeat_knot(0.0)
        >>> list(ks.iter_unique_knots())
        [0.0, 1.0]
        """
        keys = self._knots.keys()
        keys.sort() # do we need this?
        for x in keys:
            yield x

    def iter_knots(self):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.repeat_knot(0.0)
        >>> list(ks.iter_knots())
        [0.0, 0.0, 1.0]
        
        """
        for x in self.iter_unique_knots():
            for i in range(self._knots[x]):
                yield x

    def knots(self):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.repeat_knot(0.0)
        >>> ks.knots()
        [0.0, 0.0, 1.0]
        """
        return list(self.iter_knots())

    def unique_knots(self):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.repeat_knot(0.0)
        >>> ks.unique_knots()
        [0.0, 1.0]
        """
        return list(self.iter_unique_knots())

    def unique_knots_as_array(self):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.repeat_knot(0.0)
        >>> ks.unique_knots_as_array()
        array([ 0.,  1.])
        """
        return np.array(self.unique_knots(), dtype='float64')

    def knots_as_array(self):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.repeat_knot(0.0)
        >>> ks.knots_as_array()
        array([ 0.,  0.,  1.])
        """
        return np.array(self.knots(), dtype='float64')

    @property
    def xmin(self):
        return min(self._knots.iterkeys())

    @property
    def xmax(self):
        return max(self._knots.iterkeys())

    @property
    def n_knots(self):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.repeat_knot(0.0)
        >>> ks.n_knots
        3
        """
        return len(self)

    @property
    def n_unique_knots(self):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.repeat_knot(0.0)
        >>> ks.n_unique_knots
        2
        """
        return len(self._knots)

    def add_knot(self, x):
        if x not in self._knots:
            self._knots[float(x)] = 1

    def delete_knot(self, x):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(1.0)
        >>> ks.add_knot(0.0)
        >>> ks.repeat_knot(0.0)
        >>> ks.knots()
        [0.0, 0.0, 1.0]
        >>> ks.delete_knot(0.0)
        >>> ks.knots()
        [1.0]
        >>> ks.delete_knot(1.0)
        >>> ks.knots()
        []
        """
        x = float(x)
        if x in self._knots:
            del self._knots[x]

    def repeat_knot(self, x, n=1):
        x = float(x)
        if x in self._knots:
            for i in range(n):
                self._knots[x] += 1
        else:
            raise KnotSetError('knot must be present to repeat: %r' % x)
    
    def repeat_xmin(self, n=1):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.repeat_xmin(3)
        >>> ks.knots()
        [0.0, 0.0, 0.0, 0.0, 1.0]
        """
        xmin = self.xmin
        self.repeat_knot(xmin, n)
    
    def repeat_xmax(self, n=1):
        """
        >>> ks = KnotSet()
        >>> ks.add_knot(0.0)
        >>> ks.add_knot(1.0)
        >>> ks.repeat_xmin(3)
        >>> ks.knots()
        [0.0, 0.0, 0.0, 0.0, 1.0]
        >>> ks.repeat_xmax(4)
        >>> ks.knots()
        [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        """
        xmax = self.xmax
        self.repeat_knot(xmax, n)
    
    def _build(self):
        """Build the actual knot points.
        
        Subclasses should define this method.  The method should:
        
        1. Create a python list (self.knots) with the self.npoints 
           unique knot points
        2. Call self._add_repeated_knots()
        
        """
        raise NotImplementedError("implement in a subclass")
    
    def plot(self):
        """Plot the knot points using matplotlib."""
        import pylab
        xvals = self.unique_knots_as_array()
        yvals = self.n_unique_knots*[0]
        return pylab.plot(xvals, yvals, 'bo')
    
    # def bracketing_knots(self, i):
    #     """Returns the bracketing knot points for interval i.
    #     
    #     The first interval is labeled by 0.
    #     """
    #     return (self.knots[self.order - 1 + i], self.knots[self.order + i])



class UniformKnotSet(KnotSet):
    """A KnotSet with equal spaced knot points."""
    
    def __init__(self, xmin, xmax, n_knots):
        """
        >>> ks = UniformKnotSet(0.0, 10.0, 5)
        >>> ks.knots()
        [0.0, 2.5, 5.0, 7.5, 10.0]
        >>> len(ks)
        5
        """
        super(UniformKnotSet, self).__init__()
        self._build(xmin, xmax, n_knots)
        
    def _build(self, xmin, xmax, n_knots):
        """Create the knot points on a linear mesh.
        """
        
        dl = (xmax - xmin)/(n_knots-1)
        
        for i in range(n_knots):
            self._knots[xmin + i*dl] = 1


class PowerKnotSet(KnotSet):
    """A KnotSet with knot points spaced by powers."""
    
    def __init__(self, xmin, xmax, n_knots, power):
        """
        >>> ks = PowerKnotSet(0.0, 10.0, 6, 2.0)
        >>> ks.knots()
        [0.0, 0.8571428571428571, 2.2857142857142856, 4.2857142857142856, 6.8571428571428568, 9.9999999999999982]
        """
        super(PowerKnotSet, self).__init__()
        self.power = power
        self._build(xmin, xmax, n_knots)
        
    def _build(self, xmin, xmax, n_knots):
        knots = [pow(x, self.power) for x in range(1, n_knots + 1)]
        scale = (xmax - xmin)/(knots[-1] - knots[0])
        knots = [x*scale for x in knots]
        shift = xmin - knots[0]
        knots = [x + shift for x in knots]
        for x in knots:
            self.add_knot(x)

__all__ = [KnotSetError, KnotSet, UniformKnotSet, PowerKnotSet]

