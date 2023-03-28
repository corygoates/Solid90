import numpy as np
import matplotlib.pyplot as plt
from plate import ZeroNinetyZeroPlate


if __name__=="__main__":

    # Loop through boundary conditions
    for BC in ['SSSS', 'SSSF', 'SCSC', 'SCSS']:

        # initialize plate
        plate = ZeroNinetyZeroPlate(a=2.0, b=1.0, E1=155.0e9, E2=12.1e9, v12=0.248, G12=4.4e9, tl=1.5e-3, rho=7860.0, BC=BC)

        if BC == 'SCSC':
            m = 1
            w_trans = plate.get_freq_transition_point(m)
            W = np.logspace(np.log10(w_trans), 3, 1000)
            plt.figure()
            plt.plot(W, abs(plate.high_freq_characteristic(W, m)))
            plt.xscale('log')
            plt.yscale('log')
            plt.show()

        # Get natural frequencies
        print()
        print("--- BC: {0} ---".format(BC))
        w1s = plate.get_natural_freqs(1, 5)
        print("w1s: ", w1s)
        w2s = plate.get_natural_freqs(2, 5)
        print("w2s: ", w2s)

        # Plot mode shapes
        plate.plot_mode_shape(1, 1)