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

    x1 : float
        Center of patch in x-direction.

    x2 : float
        Center of patch in y-direction.

    plate : Plate
        Plate on which the load is being applied.
    """


    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Store patch information
        self.c = kwargs.get("c")
        self.d = kwargs.get("d")
        self.x1 = kwargs.get("x1")
        self.y1 = kwargs.get("y1")


    def p(x, y):
        if x < self.x1 + self.c and x > self.x1 - self.c and y < self.y1 + self.d and y > self.y1 - self.d:
            return self.p0
        else:
            return 0.0