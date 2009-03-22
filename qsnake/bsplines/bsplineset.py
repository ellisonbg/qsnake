class BSplineSet(object):
    """An object representing a B-spline basis set.
    
    *Attributes*:
    
        self.knotset : {KnotSet}
            A KnotSet subclass instance used to define the splines.
        self.n_splines : {int}
            The total number of splines.
        self.order : {int}
            The order of splines.
        self.mesh_points : {ndarray}
            An array of the mesh points, the splines are defined on.
        self.mesh_weights : {ndarray}
            An array of Guass-Leg. integration weights corresponding
            to the mesh_points.
        self.b : {ndarray}
            An array of shape (self.n_mesh_points, self.n_splines), where
            b[p,i] contains the value of the ith spline at
            self.mesh_points[p].
        self.db : {ndarray}
            An array of shape (self.n_mesh_points, self.n_splines), where
            db[p,i] contains the value of the first derivative of the
            ith spline at self.mesh_points[p].
        self.ddb : {ndarray}
            An array of shape (self.n_mesh_points, self.n_splines), where
            ddb[p,i] contains the value of the second derivative of the
            ith spline at self.mesh_points[p].
        """
    
    def __init__(self, knotset, bcmin=0, bcmax=0):
        """Create a set of B-splines from a KnotSet and BCs.
        
        *Parameters*:
            
            knotset : {KnotSet}
                The KnotSet subclass instance used to define the B-splines.
            bcmin : {0,1}
                The boundary condition to use at the minimum of the interval.
                0 -> function vanishes
                1 -> derivative of the function vanishes
            bcmax : {0,1}
                The boundary conditionto use at the max of the interval.
        """
        self.knotset = knotset
        self.bcmin = bcmin
        self.bcmax = bcmax
        self.n_splines = self.knotset.n_splines
        self.order = self.knotset.order
        self.n_legendre = self.order + 1
        self._build()
    
    # Methods that build the splines 
    
    def _build(self):
        """Build the B-splines at the mesh points."""
        self._build_intmesh()
        self._build_splines_on_intmesh()
        self._build_bracketing_index_list()
        self._build_overlapping_list()
        
    def _build_intmesh(self):
        """Build the mesh_points and mesh_weights."""
        
        unique_points = self.knotset.knots[self.order - 1: self.n_splines + 1]
        im = GLIntMesh(unique_points, self.n_legendre)
        self.mesh_points = im.mesh_points
        self.mesh_weights = im.mesh_weights
        self.n_mesh_points = im.n_mesh_points
        
    def _build_splines_on_intmesh(self):
        """Build the splines at the integration mesh points."""
        
        self.b = np.zeros((self.n_mesh_points, self.n_splines),'float64')
        self.db = np.zeros((self.n_mesh_points, self.n_splines),'float64')
        self.ddb = np.zeros((self.n_mesh_points, self.n_splines),'float64')        

        bsf = BSplineFunction(self.knotset)
        
        # Build the splines at each integration mesh point
        for index, p in enumerate(self.mesh_points):
            (tempb, tempdb, tempddb) = bsf.all(p)
            self.b[index,:] = tempb
            self.db[index,:] = tempdb
            self.ddb[index,:] = tempddb
       
        # Slice the arrays to impose boundary conditions
        # Now the derivative vanishing
        if self.bcmin == 1:
            self.b[:,1] = (self.b[:,0] + self.b[:,1])/2.0
            self.db[:,1] = (self.db[:,0] + self.db[:,1])/2.0
            self.ddb[:,1] = (self.ddb[:,0] + self.ddb[:,1])/2.0
        if self.bcmax == 1:
            self.b[:,-2] = (self.b[:,-1] + self.b[:,-2])/2.0
            self.db[:,-2] = (self.db[:,-1] + self.db[:,-2])/2.0
            self.ddb[:,-2] = (self.ddb[:,-1] + self.ddb[:,-2])/2.0

        # First the wavefunction vanishing case
        if self.bcmin == 0 or self.bcmin == 1:
            self.b = self.b[:,1:]
            self.db = self.db[:,1:]
            self.ddb = self.ddb[:,1:]
            self.n_splines = self.n_splines - 1
        if self.bcmax == 0 or self.bcmax == 1:
            self.b = self.b[:,:-1]
            self.db = self.db[:,:-1]
            self.ddb = self.ddb[:,:-1]
            self.n_splines = self.n_splines - 1
       
    def _build_bracketing_index_list(self):
        self.bracketing_index_list = []
        for i in range(self.n_splines):
            for j in range(self.b.shape[0]):
                if not self.b[j][i]==0.0:
                    guess_min_index = j
                    break
            for j in range(guess_min_index, self.b.shape[0]):
                if not self.b[self.b.shape[0]-1][i]==0.0:
                    guess_max_index = self.b.shape[0] - 1
                    break
                elif self.b[j][i]==0.0:
                    guess_max_index = j - 1
                    break
            self.bracketing_index_list.append((guess_min_index,
                guess_max_index))
        return self.bracketing_index_list
  
    def _build_overlapping_list(self):
        self.op_min_index = np.zeros((self.n_splines,self.n_splines))
        self.op_max_index = np.zeros((self.n_splines,self.n_splines))
        for p in range(self.n_splines):
            for q in range(self.n_splines):
                min_and_max = self.overlap_indices(p,q)
                if min_and_max:
                    self.op_min_index[p,q] = min_and_max[0]
                    self.op_max_index[p,q] = min_and_max[1]
                else:
                    self.op_min_index[p,q] = 1
                    self.op_max_index[p,q] = -1
                
    # Code for computing if two splines overlap
                                    
    def overlapping_spline_list(self, i):
        """Returns a list of the splines j (with j < i) that overlap 
        the ith spline.
        
        See full_overlapping_spline_list for the complete list.
        """
        return range(max(i-self.order+1,0),i+1)
        
    def overlapping_spline_indices(self, i):
        return (max(i-self.order+1,0),i+1)
        
    def full_overlapping_spline_list(self, i):
        """Returns a list of all the splines that overlap the ith spline."""
        if i >= 0 or i < self.n_splines:
            return range(max(i-self.order+1,0),
                min(self.n_splines,i + self.order))
        else:
            return []
        
    # I need to precompute all thee indices and store them in arrays!!!
    
    def interval_indices(self,i):
        """Returns the indices of the knots that bracket the ith spline."""
        return self.bracketing_index_list[i]
                        
    def overlap_indices(self, i, j):
        """Returns the indices of the integration points where 
        the i and jth splines overlap."""
        (mina,maxa) = self.interval_indices(i)
        (minb,maxb) = self.interval_indices(j)
        if minb > maxa or mina > maxb:
            return None
        else:
            return (max(mina,minb), min(maxa, maxb))
            
    def overlapping_point_list(self, i, j):
        """Returns a list of the integration points where two splines
        overlap."""
        min_and_max = self.overlap_indices(i,j)
        if min_and_max:
            return range(min_and_max[0], min_and_max[1] + 1)
        else:
            return []
        
    def overlapping_point_indices(self, i, j):
        min_and_max = self.overlap_indices(i,j)
        if min_and_max:
            return min_and_max[0], min_and_max[1] + 1
        else:
            return -1,-1
        
    def plot_b(self):
        """Plot the B-splines as a function of x."""
        import pylab
        return pylab.plot(self.mesh_points, self.b)
        
    def plot_db(self):
        """Plot the 1st deriv. of the B-splines as a function of x."""
        import pylab
        return pylab.plot(self.mesh_points, self.db)

    def plot_ddb(self):
        """Plot the 2nd deriv. of the B-splines as a function of x."""
        import pylab
        return pylab.plot(self.mesh_points, self.ddb)
    
    def plot_pattern(self):
        """Plot the sparsity pattern of self.b."""
        import pylab
        pylab.spy(self.b, aspect='auto')
        
    def plot_contour(self):
        """Plot the B-splines as a filled contour plot."""
        import pylab
        pylab.contourf(self.b)
        
class RangeError(Exception):
    pass

class BSplineFunction(object):

    def __init__(self, knotset):
        self.knotset = knotset
        self.n_splines = self.knotset.n_splines
        self.order = self.knotset.order
        self.min = self.knotset.min
        self.max = self.knotset.max
        
    def b(self, x):
        if x < self.min or x > self.max:
            raise RangeError("Value %f is outside range [%f, %f]" % (x, self.min, self.max))
        # deboor.interv has problems in x is too close to self.max
        if abs(self.max - x) < 1.0e-12:
            x = self.max - 1.0e-12
        bvals = np.zeros(self.n_splines,'float64')
        nderiv = 1
        (left,error) = deboor.interv(self.knotset.knots,x)
        nonzero_atx = deboor.bsplvd(self.knotset.knots,
            self.order,x,left,nderiv)
        for i in range(self.order):
            bvals[left - self.order + i] = nonzero_atx[i]
        return bvals
                     
    def db(self, x):
        if x < self.min or x > self.max:
            raise RangeError("Value %f is outside range [%f, %f]" % (x, self.min, self.max))

        # deboor.interv has problems in x is too close to self.max
        if abs(self.max - x) < 1.0e-12:
            x = self.max - 1.0e-12   
        dbvals = np.zeros(self.n_splines,'float64')
        nderiv = 2
        (left,error) = deboor.interv(self.knotset.knots,x)
        nonzero_atx = deboor.bsplvd(self.knotset.knots,
            self.order,x,left,nderiv)
        for i in range(self.order):
            dbvals[left - self.order + i] = nonzero_atx[i][1]
        return dbvals
    
    def ddb(self, x):
        if x < self.min or x > self.max:
            raise RangeError("Value %f is outside range [%f, %f]" % (x, self.min, self.max))

        # deboor.interv has problems in x is too close to self.max
        if abs(self.max - x) < 1.0e-12:
            x = self.max - 1.0e-12   
        ddbvals = np.zeros(self.n_splines,'float64')
        nderiv = 3
        (left,error) = deboor.interv(self.knotset.knots,x)
        nonzero_atx = deboor.bsplvd(self.knotset.knots,
            self.order,x,left,nderiv)
        for i in range(self.order):
            ddbvals[left - self.order + i] = nonzero_atx[i][2]
        return ddbvals
    
    def all(self, x):
        if x < self.min or x > self.max:
            raise RangeError("Value %f is outside range [%f, %f]" % (x, self.min, self.max))

        # deboor.interv has problems in x is too close to self.max
        if abs(self.max - x) < 1.0e-12:
            x = self.max - 1.0e-12    
        bvals = np.zeros(self.n_splines,'float64')
        dbvals = np.zeros(self.n_splines,'float64')
        ddbvals = np.zeros(self.n_splines,'float64')
        nderiv = 3
        (left,error) = deboor.interv(self.knotset.knots,x)
        nonzero_atx = deboor.bsplvd(self.knotset.knots,
            self.order,x,left,nderiv)
        for i in range(self.order):
            bvals[left - self.order + i] = nonzero_atx[i][0]
            dbvals[left - self.order + i] = nonzero_atx[i][1]
            ddbvals[left - self.order + i] = nonzero_atx[i][2]
        return (bvals, dbvals, ddbvals)
    
    def plot_b(self, step_size=0.1):
        import pylab
        plot_points = np.arange(self.min, self.max, step_size)
        bvals = [self.b(p) for p in plot_points]
        pylab.plot(plot_points, bvals)
        
    def plot_db(self, step_size=0.1):
        import pylab
        plot_points = np.arange(self.min, self.max, step_size)
        bvals = [self.db(p) for p in plot_points]
        pylab.plot(plot_points, bvals)
        
    def plot_ddb(self, step_size=0.1):
        import pylab
        plot_points = np.arange(self.min, self.max, step_size)
        bvals = [self.ddb(p) for p in plot_points]
        pylab.plot(plot_points, bvals)