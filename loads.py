import numpy

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


    def p(x, y):
        return 0.0


class UniformLoad(Load):
    """Class for a uniform load.
    
    Parameters
    ----------
    p0 : float, optional
        Load magnitude. Defaults to 1.

    plate : Plate
        Plate on which the load is being applied.
    """

    def p(x, y):
        return self.p0


class SinusoidalLoad(Load):
    """Class for a sinusoidal load.
    
    Parameters
    ----------
    p0 : float, optional
        Load magnitude. Defaults to 1.

    plate : Plate
        Plate on which the load is being applied.
    """

    def p(x, y):
        return self.p0*np.sin(np.pi*x/self.plate.a)*np.sin(np.pi*y/self.plate.b)


class HydrostaticLoad(Load):
    """Class for a hydrostatic load.
    
    Parameters
    ----------
    p0 : float, optional
        Load magnitude. Defaults to 1.

    plate : Plate
        Plate on which the load is being applied.
    """

    def p(x, y):
        return self.p0*x/self.plate.a