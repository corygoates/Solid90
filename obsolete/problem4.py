import numpy as np
import matplotlib.pyplot as plt


class ZeroNinetyZeroPlate:
    """Defines a 0-90-0 laminate plate that is SSSS. Used for vibrations."""


    def __init__(self, **kwargs):

        # Store
        self.E1 = kwargs.get("E1")
        self.E2 = kwargs.get("E2")
        self.G12 = kwargs.get("G12")
        self.v12 = kwargs.get("v12")
        self.a = kwargs.get("a")
        self.b = kwargs.get("b")
        self.tl = kwargs.get("tl")
        self.rho = kwargs.get("rho")

        # Calc derived properties
        self.v21 = self.E2*self.v12/self.E1
        self.h = self.tl*3.0
        self.J0 = self.rho*self.h

        self.calc_D()


    def calc_D(self):
        # Calculates the D matrix for the laminate plate

        d = 1.0-self.v12*self.v21

        # Top-bottom Q matrix
        Q_TB = np.zeros((3,3))
        Q_TB[0,0] = self.E1/d
        Q_TB[0,1] = Q_TB[0,0]*self.v21
        Q_TB[1,0] = Q_TB[0,1]
        Q_TB[1,1] = self.E2/d
        Q_TB[2,2] = self.G12

        # Mid Q matrix
        Q_mid = np.zeros((3,3))
        Q_mid[0,0] = Q_TB[1,1]
        Q_mid[1,1] = Q_TB[0,0]
        Q_mid[0,1] = Q_TB[0,1]
        Q_mid[1,0] = Q_TB[1,0]
        Q_mid[2,2] = Q_TB[2,2]

        # Combine
        f_tb = ((-0.5*self.tl)**3 - (-1.5*self.tl)**3)*2.0/3.0
        f_mid = ((0.5*self.tl)**3 - (-0.5*self.tl)**3)/3.0
        self.D = Q_TB*f_tb + Q_mid*f_mid
        self.D_hat = self.D[0,1] + 2.0*self.D[2,2]


    def get_alpha_m(self, m):
        # Returns alpha_m
        return m*np.pi/self.a


    def lambda_1_squared(self, w, m):
        # Returns the first root
        alpha_m = self.get_alpha_m(m)
        return alpha_m**2*self.D_hat/self.D[1,1] + np.sqrt(self.D_hat**2/self.D[1,1]**2*alpha_m**4 + self.J0*w**2/self.D[1,1] - self.D[0,0]/self.D[1,1]*alpha_m**4)


    def lambda_2_squared(self, w, m):
        # Returns the second root
        alpha_m = self.get_alpha_m(m)
        return -alpha_m**2*self.D_hat/self.D[1,1] + np.sqrt(self.D_hat**2/self.D[1,1]**2*alpha_m**4 + self.J0*w**2/self.D[1,1] - self.D[0,0]/self.D[1,1]*alpha_m**4)


    def high_freq_characteristic(self, w, m):
        # Returns the result of the high-frequency characteristic equation

        # Get roots
        l1 = np.sqrt(self.lambda_1_squared(w, m))
        l2 = np.sqrt(self.lambda_2_squared(w, m))

        # Get sines
        S = np.sin(l2*self.b)
        Sh = np.sinh(l1*self.b)

        return -S*Sh


    def get_freq_transition_point(self, m):
        # Returns the transition between low and high frequencies
        alpha_m = m*np.pi/self.a
        return np.sqrt(self.D[0,0]*alpha_m**4/self.J0).item()


    def get_high_natural_freqs(self, m, N_freq):
        # Returns the first N natural frequencies above the transition point for the given m

        # Get transition point
        w_trans = self.get_freq_transition_point(m)

        # Initialize storage
        ws = np.zeros(N_freq)

        # Find first natural frequency
        ws[0] = self.find_high_natural_freq_w_secant(w_trans*1.1, m, w_min=w_trans*1.05)

        # Find the rest
        for i in range(1, N_freq):

            # Find where we're approaching another root
            f0 = 0.0
            f1 = self.high_freq_characteristic(ws[i-1]+1.0, m)
            w2 = ws[i-1]+2.0
            f2 = self.high_freq_characteristic(w2, m)
            while np.sign(f1-f0) == np.sign(f2-f1):
                f0 = f1
                f1 = f2
                w2 += 1.0
                f2 = self.high_freq_characteristic(w2, m)

            # Find where we change sign
            while np.sign(f1) == np.sign(f2):
                f1 = f2
                w2 += 1.0
                f2 = self.high_freq_characteristic(w2, m)

            # Refine with secant method
            ws[i] = self.find_high_natural_freq_w_secant(w2, m, w_min=1.1*ws[i-1])

        return ws


    def find_high_natural_freq_w_secant(self, w0, m, w_min=0.0):
        # Uses the secant method to find a high natural frequency close to w0 for the given m

        # Get transition frequency
        w_trans = self.get_freq_transition_point(m)

        # Initialize secant method
        f0 = self.high_freq_characteristic(w0, m)
        w1 = w0*1.1
        f1 = self.high_freq_characteristic(w1, m)

        # Loop
        while abs(w1-w0) > 1.0e-6:

            # Get new root location
            w2 = w1 - f1*(w1-w0)/(f1-f0)
            if w2 < w_min:
                w2 = w_min + w_trans + np.random.random(1).item()

            # Update for next iteration
            f0 = f1
            w0 = w1
            w1 = w2
            f1 = self.high_freq_characteristic(w1, m)
        return w1


if __name__=="__main__":

    # initialize plate
    plate = ZeroNinetyZeroPlate(a=2.0, b=1.0, E1=155.0e9, E2=12.1e9, v12=0.248, G12=4.4e9, tl=1.5e-3, rho=7860.0)

    w1s = plate.get_high_natural_freqs(1, 10)
    print("w1s: ", w1s)
    w2s = plate.get_high_natural_freqs(2, 10)
    print("w2s: ", w2s)