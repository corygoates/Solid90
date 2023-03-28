import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d


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
        Young's modulus. Defaults to 30e6.

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
    BC : str, optional
        Boundary conditions. May be 'SSSS', 'SSSF', 'SCSC', or 'SSSC'. Defaults to 'SSSS'.
    """

    def __init__(self, **kwargs):

        # Store
        self.BC = kwargs.get('BC', 'SSSS')


    def get_alpha_m(self, m):
        # Returns alpha_m

        return m*np.pi/self.a


    def lambda_1_squared(self, w, m):
        # Returns the first root
        alpha_m = self.get_alpha_m(m)
        return alpha_m**2*self.D_hat/self.D22 + np.sqrt(self.D_hat**2/self.D22**2*alpha_m**4 + self.J0*w**2/self.D22 - self.D11/self.D22*alpha_m**4)


    def lambda_2_squared(self, w, m):
        # Returns the second root
        alpha_m = self.get_alpha_m(m)
        return -alpha_m**2*self.D_hat/self.D22 + np.sqrt(self.D_hat**2/self.D22**2*alpha_m**4 + self.J0*w**2/self.D22 - self.D11/self.D22*alpha_m**4)


    def high_freq_characteristic(self, w, m):
        # Returns the result of the high-frequency characteristic equation

        # Get prereqs
        alpha_m = self.get_alpha_m(m)
        l1_2 = self.lambda_1_squared(w, m)
        l2_2 = self.lambda_2_squared(w, m)
        l1 = np.sqrt(l1_2)
        l2 = np.sqrt(l2_2)
        Ch = np.cosh(l1*self.b)
        Sh = np.sinh(l1*self.b)
        C = np.cos(l2*self.b)
        S = np.sin(l2*self.b)

        if self.BC == 'SSSS':

            return -S*Sh

        elif self.BC == 'SSSF':

            F11 = l1_2 - self.D12/self.D22*alpha_m**2
            F12 = l2_2 + self.D12/self.D22*alpha_m**2
            F21 = l1_2 - self.D_bar/self.D22*alpha_m**2
            F22 = l2_2 + self.D_bar/self.D22*alpha_m**2

            return l2*F11*F22*Sh*C - l1*F21*F12*Ch*S

        elif self.BC == 'SCSC':

            return 2.0*(1.0 - Ch*C) + (l1/l2 - l2/l1)*Sh*S

        elif self.BC == 'SCSS':

            return l1*Ch*S - l2*Sh*C

    
    def low_freq_characteristic(self, w, m):
        # Returns the result of the low-frequency characteristic equation

        # Get roots
        l1 = np.sqrt(self.lambda_1_squared(w, m))
        l2 = np.sqrt(self.lambda_2_squared(w, m))

        if self.BC == 'SSSS':
            raise RuntimeError("No low-frequencies for SSSS!")

        elif self.BC == 'SSSF':
            raise RuntimeError("No low-frequencies for SSSF!")

        elif self.BC == 'SCSC':
            raise RuntimeError("No low-frequencies for SCSC!")


    def get_freq_transition_point(self, m):
        # Returns the transition between low and high frequencies
        alpha_m = m*np.pi/self.a
        return np.sqrt(self.D11*alpha_m**4/self.J0).item()


    def get_natural_freqs(self, m, N_freq):
        # Returns the first N natural frequencies for the given m

        # Initialize storage
        ws = np.zeros(N_freq)

        # Get low natural frequencies
        N_low = 0

        # Get high natural frequencies
        ws[N_low:] = self.get_high_natural_freqs(m, N_freq-N_low)

        return ws


    def get_high_natural_freqs(self, m, N_freq):
        # Returns the first N natural frequencies above the transition point for the given m

        # Get transition point
        w_trans = self.get_freq_transition_point(m)

        # Initialize storage
        ws = np.zeros(N_freq)

        # Find first natural frequency
        mult = 1.0
        while True:
            mult *= 1.1
            try:
                ws[0] = self.find_high_natural_freq_w_secant(w_trans*mult, m, w_min=w_trans)
            except RuntimeError:
                continue
            else:
                break

        # Find the rest
        for i in range(1, N_freq):

            # Find where we change sign
            f1 = self.high_freq_characteristic(ws[i-1]+1.0, m)
            w2 = ws[i-1]+2.0
            f2 = self.high_freq_characteristic(w2, m)
            while np.sign(f1) == np.sign(f2):
                f1 = f2
                w2 += 1.0
                f2 = self.high_freq_characteristic(w2, m)

            # Refine with secant method
            ws[i] = self.find_high_natural_freq_w_secant(w2, m, w_min=1.1*ws[i-1])

        return ws


    def find_high_natural_freq_w_secant(self, w0, m, w_min=0.0, i_max=1000):
        # Uses the secant method to find a high natural frequency close to w0 for the given m

        # Get transition frequency
        w_trans = self.get_freq_transition_point(m)

        # Initialize secant method
        f0 = self.high_freq_characteristic(w0, m)
        w1 = w0*1.1
        f1 = self.high_freq_characteristic(w1, m)

        # Loop
        i = 0
        while abs(w1-w0) > 1.0e-12 and i < i_max:

            i += 1

            # Get new root location
            w2 = w1 - f1*(w1-w0)/(f1-f0)
            if w2 < w_min:
                w2 = w_min + w_trans + np.random.random(1).item()

            # Update for next iteration
            f0 = f1
            w0 = w1
            w1 = w2
            f1 = self.high_freq_characteristic(w1, m)

        if i > i_max:
            raise RuntimeError("Secant method did not converge.")

        return w1


    def get_mode_coefficients(self, m, n):
        # Returns the coefficients for the given mode shape

        # Find the natural frequency for this mode
        ws = self.get_natural_freqs(m, n)
        w = ws[n-1]

        # Get roots
        alpha_m = self.get_alpha_m(m)
        l1_2 = self.lambda_1_squared(w, m)
        l2_2 = self.lambda_2_squared(w, m)
        l1 = np.sqrt(l1_2)
        l2 = np.sqrt(l2_2)

        if self.BC == 'SSSS':
            A = 0.0
            B = 0.0
            C = 0.0
            D = 1.0

        elif self.BC == 'SSSF':
            F11 = l1_2 - self.D12/self.D22*alpha_m**2
            F12 = l2_2 + self.D12/self.D22*alpha_m**2
            A = 0.0
            B = F12*np.sin(l2*self.b)
            C = 0.0
            D = F11*np.sinh(l1*self.b)

        elif self.BC == 'SCSC':
            K11 = np.cosh(l1*self.b) - np.cos(l2*self.b)
            K12 = np.sinh(l1*self.b) - l1/l2*np.sin(l2*self.b)
            A = K12
            B = -K11
            C = -K12
            D = K11*l1/l2

        elif self.BC == 'SCSS':
            A = 0.0
            B = np.sin(l2*self.b)
            C = 0.0
            D = -np.sinh(l1*self.b)

        return A, B, C, D, w, l1, l2


    def plot_mode_shape(self, m, n):
        """Plots the given mode shape for the plate."""

        # Get range of x and y
        x = np.linspace(0.0, self.a, 20*m)
        y = np.linspace(0.0, self.b, 20*n)

        # Get coefficients
        A, B, C, D, w, l1, l2 = self.get_mode_coefficients(m, n)

        # Determine displacements
        X, Y = np.meshgrid(x, y)
        W = (A*np.cosh(l1*Y) + B*np.sinh(l1*Y) + C*np.cos(l2*Y) + D*np.sin(l2*Y))*np.sin(m*np.pi*X/self.a)
        W /= np.max(np.max(np.abs(W))).item()

        # Plot
        ax = plt.figure().add_subplot(projection='3d')
        ax.plot_surface(X, Y, W, cmap='coolwarm')
        ax.set(xlabel='x', ylabel='y', zlabel='w', zlim=(2.0, -2.0), xlim=(self.a, 0.0))
        plt.show()


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

    BC : str, optional
        Boundary conditions. May be 'SSSS', 'SSSF', 'SCSC', or 'SSSC'. Defaults to 'SSSS'.
    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        # Store parameters
        self.a = kwargs.get("a", 1)
        self.b = kwargs.get("b", 1)
        self.tl = kwargs.get("tl", 0.01)
        self.h = self.tl*3.0
        self.E1 = kwargs.get("E1")
        self.E2 = kwargs.get("E2")
        self.v12 = kwargs.get("v12")
        self.G12 = kwargs.get("G12")
        self.rho = kwargs.get("rho")

        # Calculate derived properties
        self.v21 = self.E2*self.v12/self.E1
        self.h = self.tl*3.0
        self.J0 = self.rho*self.h
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
        f_TB = 2.0*((0.5*self.h)**3 - (0.5*self.tl)**3)/3.0
        f_mid = 2.0*(0.5*self.tl)**3/3.0
        self.D11 = f_TB*Q11_TB + f_mid*Q11_mid
        self.D22 = f_TB*Q22_TB + f_mid*Q22_mid
        self.D12 = f_TB*Q12_TB + f_mid*Q12_mid
        self.D66 = f_TB*Q66_TB + f_mid*Q66_mid

        # Calculate D_hat and D_bar
        self.D_hat = self.D12 + 2.0*self.D66
        self.D_bar = self.D12 + 4.0*self.D66