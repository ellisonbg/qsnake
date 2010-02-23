import numpy as np
import scipy as sp
from scipy import special

class IntMesh(object):
    
    def __init__(self, grid):
        self.grid = np.array(grid, dtype='float64')
        self._mesh_points = None
        self._mesh_weights = None
    
    def _build(self):
        """
        This must set self._mesh_points and self._mesh_weights
        """
        raise NotImplementedError("implement in a subclass")
    
    def mesh_points(self):
        return self._mesh_points
    
    def mesh_weights(self):
        return self._mesh_weights
    
    def integrate(self, f):
        # vf = np.vectorize(f)
        return np.sum(self._mesh_weights*f(self._mesh_points))
    
    @property
    def n_mesh_points(self):
        return len(self._mesh_points)
    
    @property
    def n_intervals(self):
        return len(self.grid)-1
        
        

class GLIntMesh(IntMesh):
    """
    Builds an integration mesh an weights for Gauss-Legendre quadrature.
    """
    
    def __init__(self, grid, order):
        """Create the integration mesh.

        This class creates the integration points and weights needed
        to perform Gauss-Legendre quadrature over a predefined grid.
        
        The integral of a function f can be computed using::
        
            self.mesh_weights*f(self.mesh_points)
        
        *Parameters*:
        
            grid : {array_like}
                A sequence of unique points defining the mesh.
            n_legendre: {int}
                The order the GL quadrature.
        
        Here is a basic example::
        
            >>> im = GLIntMesh(range(10), 3)
            >>> im.mesh_points()
            array([ 0.11270167,  0.5       ,  0.88729833,  1.11270167,  1.5       ,
                    1.88729833,  2.11270167,  2.5       ,  2.88729833,  3.11270167,
                    3.5       ,  3.88729833,  4.11270167,  4.5       ,  4.88729833,
                    5.11270167,  5.5       ,  5.88729833,  6.11270167,  6.5       ,
                    6.88729833,  7.11270167,  7.5       ,  7.88729833,  8.11270167,
                    8.5       ,  8.88729833])
            >>> im.mesh_weights()
            array([ 0.27777778,  0.44444444,  0.27777778,  0.27777778,  0.44444444,
                    0.27777778,  0.27777778,  0.44444444,  0.27777778,  0.27777778,
                    0.44444444,  0.27777778,  0.27777778,  0.44444444,  0.27777778,
                    0.27777778,  0.44444444,  0.27777778,  0.27777778,  0.44444444,
                    0.27777778,  0.27777778,  0.44444444,  0.27777778,  0.27777778,
                    0.44444444,  0.27777778])
            >>> im.n_intervals
            9
            >>> im.order
            3
        
        Now, let's actually try to integrate a basic function::
        
            >>> import numpy as np
            >>> grid = np.linspace(0.0, 2.0*np.pi, 16)
            >>> im = GLIntMesh(grid, 3)
            >>> np.sum(im.mesh_weights()*np.cos(im.mesh_points()))
            -5.2735593669694936e-16
        
        Or, we can just use the builtin integrate method::
        
            >>> im.integrate(np.cos)
            -5.2735593669694936e-16
            >>> grid = np.linspace(0.0, 1.0, 16)
            >>> im = GLIntMesh(grid, 3)
            >>> im.integrate(lambda x: x)
            0.50000000000000022
        """
        IntMesh.__init__(self, grid)
        self.order = order
        self._build()

    def _build(self):
        _n_mesh_points = self.n_intervals*self.order
        self._mesh_points = np.zeros(_n_mesh_points, 'float64')
        self._mesh_weights = np.zeros(_n_mesh_points, 'float64')
        
        [p_points, p_weights] = special.orthogonal.p_roots(self.order)
        p_points = np.real(p_points)
        
        for i in range(self.n_intervals):
            (xmin, xmax) = (self.grid[i], self.grid[i + 1])
            scale = (xmax - xmin)/2.0
            scaled_p_points = scale*p_points
            shifted_p_points = scaled_p_points + xmin + scale
            lspp = len(shifted_p_points)
            self._mesh_points[i * lspp : (i + 1)*lspp] = shifted_p_points
            self._mesh_weights[i * lspp : (i + 1)*lspp] = p_weights*scale

__all__ = [IntMesh, GLIntMesh]