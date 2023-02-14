import numpy as np

from plate import Plate


class Load:
    """Base class for loads on a rectangular plate.
    
    Parameters
    ----------
    p0 : float, optional
        Load magnitude. Defaults to 1.

    plate : Plate
        Plate on which the load is being applied.
    """

    def __init__(self, **kwargs):

        # Store parameters
        self.p0 = kwargs.get("p0", 1.0)
        self.plate = kwargs.get("plate", Plate())


    def p(self, x, y):
        return 0.0


    def p_from_series(self, x, y, m_terms=10, n_terms=10):
        p = 0.0
        for mi in self.m[:min(m_terms, len(self.m))]:
            for ni in self.n[:min(n_terms,len(self.n))]:
                p += self.P_mn(mi, ni)*np.sin(mi*np.pi*x/self.plate.a)*np.sin(ni*np.pi*y/self.plate.b)

        return p


class UniformLoad(Load):
    """Class for a uniform load.
    
    Parameters
    ----------
    p0 : float, optional
        Load magnitude. Defaults to 1.

    plate : Plate
        Plate on which the load is being applied.
    """


    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Initialize Fourier coeff calcs
        self.m = np.array(range(1, 99, 2))
        self.n = np.array(range(1, 99, 2))


    def p(self, x, y):
        return self.p0

    
    def P_mn(self, m, n):
        return 16.0*self.p0/(np.pi**2*m*n)


class SinusoidalLoad(Load):
    """Class for a sinusoidal load.
    
    Parameters
    ----------
    p0 : float, optional
        Load magnitude. Defaults to 1.

    plate : Plate
        Plate on which the load is being applied.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Initialize Fourier coefs calcs
        self.m = np.array([1], dtype=int)
        self.n = np.array([1], dtype=int)


    def p(self, x, y):
        return self.p0*np.sin(np.pi*x/self.plate.a)*np.sin(np.pi*y/self.plate.b)


    def P_mn(self, m, n):
        if m == 1 and n == 1:
            return self.p0
        else:
            return 0.0


class HydrostaticLoad(Load):
    """Class for a hydrostatic load.
    
    Parameters
    ----------
    p0 : float, optional
        Load magnitude. Defaults to 1.

    plate : Plate
        Plate on which the load is being applied.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Initialize Fourier coefficient calcs
        self.m = np.array(range(1, 50))
        self.n = np.array(range(1, 99, 2))


    def p(self, x, y):
        return self.p0*x/self.plate.b


    def P_mn(self, m, n):
        return 8.0*self.p0/(np.pi**2*m*n)*-1**(m+1)


class PatchLoad(Load):
    """Class for a patch load.
    
    Parameters
    ----------
    p0 : float, optional
        Load magnitude. Defaults to 1.

    c : float
        Width of patch in x-direction.

    d : float
        Width of patch in y-direction.

    x : float
        Center of patch in x-direction.

    y : float
        Center of patch in y-direction.

    plate : Plate
        Plate on which the load is being applied.
    """


    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Store patch information
        self.c = kwargs.get("c")
        self.d = kwargs.get("d")
        self.x0 = kwargs.get("x")
        self.y0 = kwargs.get("y")

        # Initialize Fourier coef calcs
        self.m = np.array(range(1,50))
        self.n = np.array(range(1,50))


    def p(self, x, y):
        if x < self.x0 + self.c and x > self.x0 - self.c and y < self.y0 + self.d and y > self.y0 - self.d:
            return self.p0
        else:
            return 0.0


    def P_mn(self, m, n):
        x = 4.0*self.p0/(np.pi**2*m*n*self.c*self.d)
        x *= np.sin(m*np.pi*self.x0/self.plate.a)
        x *= np.sin(n*np.pi*self.y0/self.plate.b)
        x *= np.sin(m*np.pi*self.c/self.plate.a)
        x *= np.sin(n*np.pi*self.d/self.plate.b)
        return x