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
        self.D = self.E*self.h**3/(12.0*(1.0-self.v**2))


class OrthotropicPlate:
    """Defines an orthotropic plate for analysis using the Navier or Levy solution methods.
    
    Parameters
    ----------
    D : ndarray
        3x3 D matrix.
    """

    def __init__(self):
        pass


class ZeroNinetyZeroPlate(OrthotropicPlate):
    """A laminated plate having a 0/90/0 lamination. Material properties are assumed the same between layers.
    
    Parameters
    ----------
    a : float, optional
        Plate dimension in the x-direction. Defaults to 1.

    b : float, optional
        Plate dimension in the y-direction. Defaults to 1.

    t : float, optional
        Lamina thickness. Defaults to 0.01.

    E1 : float
        Young's modulus in the 1-direction.

    E2 : float
        Young's modulus in the 2-direction.

    v12 : float
        Poisson's ratio.

    G12 : float
        Shear modulus.
    """

    def __init__(self, **kwargs):

        # Store parameters
        self.a = kwargs.get("a", 1)
        self.b = kwargs.get("b", 1)
        self.t = kwargs.get("t", 0.01)
        self.h = self.t*3.0
        self.E1 = kwargs.get("E1")
        self.E2 = kwargs.get("E2")
        self.v12 = kwargs.get("v12")
        self.G12 = kwargs.get("G12")

        # Calculate derived properties
        self.v21 = self.E2*self.v12/self.E1
        self.calc_D()

    
    def calc_D(self):
        # Calculates the elements of the D matrix

        # Get Q matrices
        denom = 1.0 - self.v12*self.v21
        Q11_TB = self.E1/denom
        Q22_TB = self.E2/denom
        Q12_TB = self.v21*self.E1/denom
        Q66_TB = self.G12
        Q11_mid = Q22_TB
        Q22_mid = Q11_TB
        Q12_mid = Q12_TB
        Q66_mid = Q66_TB

        # Get D
        f_TB = 2.0*((0.5*self.h)**3 - (0.5*self.t)**3)/3.0
        f_mid = 2.0*(0.5*self.t)**3/3.0
        self.D11 = f_TB*Q11_TB + f_mid*Q11_mid
        self.D22 = f_TB*Q22_TB + f_mid*Q22_mid
        self.D12 = f_TB*Q12_TB + f_mid*Q12_mid
        self.D66 = f_TB*Q66_TB + f_mid*Q66_mid