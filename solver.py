import matplotlib.pyplot as plt

from plate import Plate
from loads import *


class Solver:
    """Generic plate solver class.
    
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
        contour_plot = ax.contourf(X, Y, value.T, 100)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_aspect('equal')
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

    
    def sigma_x(self, x, y, z):
        # Returns the x bending stress at the given location

        P = self.d2w_dx2(x, y) + self.plate.v*self.d2w_dy2(x, y)
        return -self.plate.E*z/(1.0-self.plate.v**2)*P

    
    def sigma_y(self, x, y, z):
        # Returns the y bending stress at the given location

        P = self.d2w_dy2(x, y) + self.plate.v*self.d2w_dx2(x, y)
        return -self.plate.E*z/(1.0-self.plate.v**2)*P


    def tau_xy(self, x, y, z):
        # Returns the plane shear stress at the given location

        return -self.plate.E*z/(1.0+self.plate.v)*self.d2w_dxdy(x, y)


    def tau_xz(self, x, y, z):
        # Returns the transverse shear stress (from equilibrium) at the given location

        P = self.d3w_dx3(x, y) + self.d3w_dxdy2(x,y)
        Z = 0.25**self.plate.h**2 - z**2
        return self.plate.E*Z*P/(2.0*(1.0-self.plate.v**2))


    def tau_yz(self, x, y, z):
        # Returns the transverse shear stress (from equilibrium) at the given location

        P = self.d3w_dx2dy(x, y) + self.d3w_dy3(x,y)
        Z = 0.25**self.plate.h**2 - z**2
        return self.plate.E*Z*P/(2.0*(1.0-self.plate.v**2))


    def corner_reactions(self):
        # Returns the corner reaction forces
        pass


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
        Pmn = self.load.Pmn(m, n)
        return Pmn/(((m/self.plate.a)**2 + (n/self.plate.b)**2)**2 * (np.pi**4*self.plate.D))


    def w(self, x, y):
        # Returns the plate deflection at the given location

        w = 0.0
        for mi in self.load.m[:self.m_max]:
            for ni in self.load.n[:self.n_max]:
                Sx = np.sin(mi*np.pi*x/self.plate.a)
                Sy = np.sin(ni*np.pi*y/self.plate.b)
                w += self.W_mn(mi, ni)*Sx*Sy

        return w


    def d2w_dx2(self, x, y):
        # Returns the second derivative of w wrt x

        d2w_dx2 = 0.0
        for mi in self.load.m[:self.m_max]:
            for ni in self.load.n[:self.n_max]:
                Sx = np.sin(mi*np.pi*x/self.plate.a)
                Sy = np.sin(ni*np.pi*y/self.plate.b)
                d2w_dx2 += -(mi*np.pi/self.plate.a)**2*self.W_mn(mi, ni)*Sx*Sy

        return d2w_dx2


    def d3w_dx3(self, x, y):
        # Returns the third derivative of w wrt x

        d3w_dx3 = 0.0
        for mi in self.load.m[:self.m_max]:
            for ni in self.load.n[:self.n_max]:
                Cx = np.cos(mi*np.pi*x/self.plate.a)
                Sy = np.sin(ni*np.pi*y/self.plate.b)
                d3w_dx3 += -(mi*np.pi/self.plate.a)**3*self.W_mn(mi, ni)*Cx*Sy

        return d3w_dx3


    def d2w_dy2(self, x, y):
        # Returns the second derivative of w wrt y

        d2w_dy2 = 0.0
        for mi in self.load.m[:self.m_max]:
            for ni in self.load.n[:self.n_max]:
                Sx = np.sin(mi*np.pi*x/self.plate.a)
                Sy = np.sin(ni*np.pi*y/self.plate.b)
                d2w_dy2 += -(ni*np.pi/self.plate.b)**2*self.W_mn(mi, ni)*Sx*Sy

        return d2w_dy2


    def d3w_dy3(self, x, y):
        # Returns the third derivative of w wrt y

        d3w_dy3 = 0.0
        for mi in self.load.m[:self.m_max]:
            for ni in self.load.n[:self.n_max]:
                Sx = np.sin(mi*np.pi*x/self.plate.a)
                Cy = np.cos(ni*np.pi*y/self.plate.b)
                d3w_dy3 += -(ni*np.pi/self.plate.b)**3*self.W_mn(mi, ni)*Sx*Cy

        return d3w_dy3


    def d2w_dxdy(self, x, y):
        # Returns the second mixed derivative of W

        d2w_dxdy = 0.0
        for mi in self.load.m[:self.m_max]:
            for ni in self.load.n[:self.n_max]:
                Cx = np.cos(mi*np.pi*x/self.plate.a)
                Cy = np.cos(ni*np.pi*y/self.plate.b)
                d2w_dxdy += (mi*np.pi/self.plate.a)*(ni*np.pi/self.plate.b)*self.W_mn(mi, ni)*Cx*Cy

        return d2w_dxdy


    def d3w_dx2dy(self, x, y):
        # Returns the third derivative of w wrt x squared and y

        d3w_dx2dy = 0.0
        for mi in self.load.m[:self.m_max]:
            for ni in self.load.n[:self.n_max]:
                Sx = np.sin(mi*np.pi*x/self.plate.a)
                Cy = np.cos(ni*np.pi*y/self.plate.b)
                d3w_dx2dy += -(mi*np.pi/self.plate.a)**2*(ni*np.pi/self.plate.b)*self.W_mn(mi, ni)*Sx*Cy

        return d3w_dx2dy


    def d3w_dxdy2(self, x, y):
        # Returns the third derivative of w wrt x and y squared

        d3w_dxdy2 = 0.0
        for mi in self.load.m[:self.m_max]:
            for ni in self.load.n[:self.n_max]:
                Cx = np.cos(mi*np.pi*x/self.plate.a)
                Sy = np.sin(ni*np.pi*y/self.plate.b)
                d3w_dxdy2 += -(mi*np.pi/self.plate.a)*(ni*np.pi/self.plate.b)**2*self.W_mn(mi, ni)*Cx*Sy

        return d3w_dxdy2


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
        self.symmetric = self.BC[1] == self.BC[3]


    def w(self, x, y):
        # Returns the deflection at the given location

        if self.symmetric:
            y -= 0.5*self.plate.b

        w = 0.0
        for m in self.load.m[:self.m_max]:

            Am, Bm, Cm, Dm, km = self.coefs(m)
            beta_m = m*np.pi/self.plate.a
            Shy = np.sinh(beta_m*y)
            Chy = np.cosh(beta_m*y)
            yShy = y*np.sinh(beta_m*y)
            yChy = y*np.cosh(beta_m*y)
            Sx = np.sin(beta_m*x)
            w += (Am*Shy + Bm*Chy + Cm*yShy + Dm*yChy + km)*Sx

        return w


    def km(self, m):
        # Returns the km coefficient

        return self.load.Pm(m)*(self.plate.a/(m*np.pi))**4/self.plate.D
    

    def coefs(self, m):
        # Returns Am, Bm, Cm, Dm, and km

        # Get k and beta
        km = self.km(m)
        beta_m = m*np.pi/self.plate.a

        # Get alpha based on symmetry
        if self.symmetric:
            alpha_m = 0.5*beta_m*self.plate.b
        else:
            alpha_m = beta_m*self.plate.b

        # Get intermediate calcs
        Tha = np.tanh(alpha_m)
        Cha = np.cosh(alpha_m)
        
        # Calculate coefs based on boundary conditions

        # All simply-supported
        if self.BC == "SSSS":

            # Get coefficients
            Am = 0.0
            Bm = -km*(2.0+alpha_m*Tha)/(2.0*Cha)
            Cm = km*beta_m/(2.0*Cha)
            Dm = 0.0

        # Clamped on y-faces
        elif self.BC == "SCSC":

            # Get coefficients
            denom = 1.0/Cha*(Tha + alpha_m*(1.0-Tha**2))
            Am = 0.0
            Bm = -km*(alpha_m+Tha)*denom
            Cm = km*beta_m*Tha*denom
            Dm = 0.0

        # Clamped at y=0, free at y=b
        elif self.BC == "SCSF":

            # Get matrix elements
            a11 = -beta_m*((1.0+self.plate.v)*Tha + alpha_m*(1.0-self.plate.v))
            a12 = alpha_m*(1.0-self.plate.v)*Tha + 2.0
            a21 = -beta_m*(2.0 + alpha_m*(self.plate.v-1.0)*Tha)
            a22 = (1.0+self.plate.v)*Tha + alpha_m*(self.plate.v-1.0)
            b1 = km*beta_m*(self.plate.v/Cha + 1.0 - self.plate.v)
            b2 = km*beta_m*(self.plate.v-1.0)*Tha

            # Get solution
            denom = a12*a21 - a11*a22
            Am = (b2*a12 - b1*a22)/denom
            Cm = (b1*a21 - b2*a11)/denom
            Bm = -km
            Dm = -beta_m*Am

        else:
            IOError("{0} is not an allowable BC set for the Levy solution! Quitting...")

        return Am, Bm, Cm, Dm, km

    
    def d2w_dx2(self, x, y):
        # Returns the second derivative of w wrt x at the given location

        if self.symmetric:
            y -= 0.5*self.plate.b

        d2w_dx2 = 0.0
        for m in self.load.m[:self.m_max]:

            Am, Bm, Cm, Dm, km = self.coefs(m)
            beta_m = m*np.pi/self.plate.a
            Shy = np.sinh(beta_m*y)
            Chy = np.cosh(beta_m*y)
            yShy = y*np.sinh(beta_m*y)
            yChy = y*np.cosh(beta_m*y)
            Sx = np.sin(beta_m*x)

            d2w_dx2 -= (Am*Shy + Bm*Chy + Cm*yShy + Dm*yChy + km)*beta_m**2*Sx

        return d2w_dx2

    
    def d3w_dx3(self, x, y):
        # Returns the third derivative of w wrt x at the given location

        if self.symmetric:
            y -= 0.5*self.plate.b

        d3w_dx3 = 0.0
        for m in self.load.m[:self.m_max]:

            Am, Bm, Cm, Dm, km = self.coefs(m)
            beta_m = m*np.pi/self.plate.a
            Shy = np.sinh(beta_m*y)
            Chy = np.cosh(beta_m*y)
            yShy = y*np.sinh(beta_m*y)
            yChy = y*np.cosh(beta_m*y)
            Cx = np.cos(beta_m*x)

            d3w_dx3 -= (Am*Shy + Bm*Chy + Cm*yShy + Dm*yChy + km)*beta_m**3*Cx

        return d3w_dx3

    
    def d2w_dxdy(self, x, y):
        # Returns the second derivative of w wrt x and y at the given location

        if self.symmetric:
            y -= 0.5*self.plate.b

        d2w_dxdy = 0.0
        for m in self.load.m[:self.m_max]:

            Am, Bm, Cm, Dm, km = self.coefs(m)
            beta_m = m*np.pi/self.plate.a
            Shy = np.sinh(beta_m*y)
            Chy = np.cosh(beta_m*y)
            yShy = y*np.sinh(beta_m*y)
            yChy = y*np.cosh(beta_m*y)
            Cx = np.cos(beta_m*x)

            d2w_dxdy += (Am*beta_m*Chy + Bm*beta_m*Shy + Cm*(Shy + beta_m*yChy) + Dm*(Chy + beta_m*yShy))*beta_m*Cx

        return d2w_dxdy


    def d2w_dy2(self, x, y):
        # Returns the second derivative of w wrt y at the given location

        if self.symmetric:
            y -= 0.5*self.plate.b

        d2w_dy2 = 0.0
        for m in self.load.m[:self.m_max]:

            Am, Bm, Cm, Dm, km = self.coefs(m)
            beta_m = m*np.pi/self.plate.a
            Shy = np.sinh(beta_m*y)
            Chy = np.cosh(beta_m*y)
            yShy = y*np.sinh(beta_m*y)
            yChy = y*np.cosh(beta_m*y)
            Sx = np.sin(beta_m*x)

            d2w_dy2 += (Am*beta_m**2*Shy + Bm*beta_m**2*Chy + Cm*(2.0*beta_m*Chy + beta_m**2*yShy) + Dm*(2.0*beta_m*Shy + beta_m**2*yChy))*Sx

        return d2w_dy2


    def d3w_dy3(self, x, y):
        # Returns the third derivative of w wrt y at the given location

        if self.symmetric:
            y -= 0.5*self.plate.b

        d3w_dy3 = 0.0
        for m in self.load.m[:self.m_max]:

            Am, Bm, Cm, Dm, km = self.coefs(m)
            beta_m = m*np.pi/self.plate.a
            Shy = np.sinh(beta_m*y)
            Chy = np.cosh(beta_m*y)
            yShy = y*np.sinh(beta_m*y)
            yChy = y*np.cosh(beta_m*y)
            Sx = np.sin(beta_m*x)

            d3w_dy3 += (Am*beta_m**3*Chy + Bm*beta_m**3*Shy + Cm*(3.0*beta_m**2*Shy + beta_m**3*yChy) + Dm*(3.0*beta_m**2*Chy + beta_m**3*yShy))*Sx

        return d3w_dy3


    def d3w_dx2dy(self, x, y):
        # Returns the third derivative of w wrt x squared and y at the given location

        if self.symmetric:
            y -= 0.5*self.plate.b

        d3w_dx2dy = 0.0
        for m in self.load.m[:self.m_max]:

            Am, Bm, Cm, Dm, km = self.coefs(m)
            beta_m = m*np.pi/self.plate.a
            Shy = np.sinh(beta_m*y)
            Chy = np.cosh(beta_m*y)
            yShy = y*np.sinh(beta_m*y)
            yChy = y*np.cosh(beta_m*y)
            Sx = np.sin(beta_m*x)

            d3w_dx2dy -= (Am*beta_m*Chy + Bm*beta_m*Shy + Cm*(Shy + beta_m*yChy) + Dm*(Chy + beta_m*yShy))*beta_m**2*Sx

        return d3w_dx2dy


    def d3w_dxdy2(self, x, y):
        # Returns the third derivative of w wrt x and y squared at the given location

        if self.symmetric:
            y -= 0.5*self.plate.b

        d3w_dxdy2 = 0.0
        for m in self.load.m[:self.m_max]:

            Am, Bm, Cm, Dm, km = self.coefs(m)
            beta_m = m*np.pi/self.plate.a
            Shy = np.sinh(beta_m*y)
            Chy = np.cosh(beta_m*y)
            yShy = y*np.sinh(beta_m*y)
            yChy = y*np.cosh(beta_m*y)
            Cx = np.cos(beta_m*x)

            d3w_dxdy2 += (Am*beta_m**2*Shy + Bm*beta_m**2*Chy + Cm*(2.0*beta_m*Chy+beta_m**2*yShy) + Dm*(2.0*beta_m*Shy + beta_m**2*yChy))*beta_m*Cx

        return d3w_dxdy2