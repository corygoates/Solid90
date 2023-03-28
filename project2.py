import numpy as np
from plate import ZeroNinetyZeroPlate


if __name__=="__main__":

    # Loop through boundary conditions
    for BC in ['SSSS', 'SSSF', 'SCSC', 'SCSS']:

        # initialize plate
        plate = ZeroNinetyZeroPlate(a=2.0, b=1.0, E1=155.0e9, E2=12.1e9, v12=0.248, G12=4.4e9, tl=1.5e-3, rho=7860.0, BC=BC)

        # Get natural frequencies
        print()
        print("--- BC: {0} ---".format(BC))
        w1s = plate.get_natural_freqs(1, 5)
        print("w1s: ", w1s)
        w2s = plate.get_natural_freqs(2, 5)
        print("w2s: ", w2s)

        # Plot mode shapes
        plate.plot_mode_shape(1, 1)