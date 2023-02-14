import matplotlib.pyplot as plt

from plate import Plate
from loads import *


class Solver:
    """Generic plate self. class.
    
    Parameters
    ----------
    plate : Plate
        Plate to solve the deflection for.
    
    load : Load
        Load applied.

    m_max : int
        Maximum m index for sine series.

    n_max : int
        Maximum n index for sine series.
    """


    def __init__(self, plate, load, m_max, n_max=1):
        
        # Initialize
        self.plate = plate
        self.load = load
        self.m_max = m_max
        self.n_max = n_max


    def get_deflection_field(self, Nx, Ny):
        # Returns the deflections across the plate

        # Get deflections at various points
        Nx = 50
        Ny = 50
        X = np.linspace(0.0, self.plate.a, Nx)
        Y = np.linspace(0.0, self.plate.b, Ny)
        w = np.zeros((Nx, Ny))
        for i, xi in enumerate(X):
            for j, yj in enumerate(Y):
                w[i,j] = self.w(xi, yj)

        return X, Y, w


    def plot_field(self, X, Y, value, label=''):

        # Plot deflection
        fig, ax = plt.subplots()
        contour_plot = ax.contourf(X, Y, value, 100)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.invert_yaxis()
        cbar = fig.colorbar(contour_plot)
        cbar.ax.set_title(label)
        plt.show()


    def plot_deflection_field(self, Nx=50, Ny=50):

        X, Y, w = self.get_deflection_field(Nx, Ny)
        self.plot_field(X, Y, w, 'Deflection')


    def compare_to_analytic_sinusoidal_solution(self):

        # Get deflections at various points
        Nx = 50
        Ny = 50
        X = np.linspace(0.0, self.plate.a, Nx)
        Y = np.linspace(0.0, self.plate.b, Ny)
        w = np.zeros((Nx, Ny))
        w_anl = np.zeros_like(w)
        for i, xi in enumerate(X):
            for j, yj in enumerate(Y):
                w[i,j] = self.w(xi, yj)
                w_anl[i,j] = self.load.p0*np.sin(np.pi*xi/self.plate.a)*np.sin(np.pi*yj/self.plate.b)
        w_anl /= np.pi**4*self.plate.D*(self.plate.a**-2 + self.plate.b**-2)**2

        # Plot deflection
        fig, ax = plt.subplots()
        contour_plot = ax.contourf(X, Y, w-w_anl, 100)
        cbar = fig.colorbar(contour_plot)
        plt.show()


class NavierSolver(Solver):
    """Solves for the deflection using the Navier solution method.
    
    Parameters
    ----------
    plate : Plate
        Plate to solve the deflection for.
    
    load : Load
        Load applied.

    m_max : int
        Maximum m index for sine series.

    n_max : int
        Maximum n index for sine series.
    """


    def __init__(self, plate, load, m_max, n_max):
        super().__init__(plate, load, m_max, n_max)


    def W_mn(self, m, n):
        # Returns the desired Fourier deflection coefficient
        P_mn = self.load.P_mn(m, n)
        return P_mn/(((m/self.plate.a)**2 + (n/self.plate.b)**2)**2 * (np.pi**4*self.plate.D))


    def w(self, x, y):
        # Returns the plate deflection at the given location

        w = 0.0
        for mi in self.load.m[:self.m_max]:
            for ni in self.load.n[:self.n_max]:
                Sx = np.sin(mi*np.pi*x/self.plate.a)
                Sy = np.sin(ni*np.pi*y/self.plate.b)
                w += self.W_mn(mi, ni)*Sx*Sy

        return w


class LevySolver(Solver):
    """Solves for the deflection using the Levy solution method.
    
    Parameters
    ----------
    plate : Plate
        Plate to solve the deflection for.
    
    load : Load
        Load applied.

    BC : str
        Boundary conditions to be applied.

    m_max : int
        Maximum m index for sine series.
    """


    def __init__(self, plate, load, BC, m_max):
        super().__init__(plate, load, m_max)

        self.BC = BC