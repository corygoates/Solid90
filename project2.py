import numpy as np
import matplotlib.pyplot as plt
from plate import ZeroNinetyZeroPlate


if __name__=="__main__":

    # Loop through boundary conditions
    w1s = []
    w2s = []
    for BC in ['SSSS', 'SSSF', 'SCSC', 'SCSS']:

        # initialize plate
        plate = ZeroNinetyZeroPlate(a=1.0, b=2.0, E1=155.0e9, E2=12.1e9, v12=0.248, G12=4.4e9, tl=1.5e-3, rho=2000.0, BC=BC)

        if False:#BC == 'SCSC':
            m = 1
            w_trans = plate.get_freq_transition_point(m)
            W = np.logspace(np.log10(w_trans), 3, 1000)
            plt.figure()
            plt.plot(W, abs(plate.high_freq_characteristic(W, m)))
            plt.xscale('log')
            plt.yscale('log')
            plt.show()

        # Get natural frequencies
        w1s.append(plate.get_natural_freqs(1, 5))
        w2s.append(plate.get_natural_freqs(2, 5))

        # Plot mode shapes
        plate.plot_mode_shape(2, 1)

    w1s = np.array(w1s)
    w2s = np.array(w2s)

    for i in range(3):
        print(w1s[:,i].flatten())

    for i in range(3):
        print(w2s[:,i].flatten())