import numpy as np


class Plate:
    """Defines a plate for analysis using the Navier or Levy solution methods.
    
    Parameters
    ----------
    a : float, optional
        Plate dimension in the x-direction. Defaults to 1.

    b : float, optional
        Plate dimension in the y-direction. Defaults to 1.

    h : float, optional
        Plate thickness. Defaults to 0.01.

    E : float, optional
        Young's modulus. Defaults to 30 Msi.

    v : float, optional
        Poisson's ratio. Defaults to 0.33.
    """

    def __init__(self, **kwargs):

        # Store parameters
        self.a = kwargs.get("a", 1)
        self.b = kwargs.get("b", 1)
        self.h = kwargs.get("h", 0.01)
        self.E = kwargs.get("E", 30.0e6)
        self.v = kwargs.get("v", 0.33)