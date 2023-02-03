import matplotlib.pyplot as plt

from plate import Plate
from loads import *


if __name__=="__main__":

    # Initialize plate
    plate = Plate()

    # Initialize load
    load = HydrostaticLoad(plate=plate)

    # Get load at various points
    Nx = 50
    Ny = 50
    X = np.linspace(0.0, plate.a, Nx)
    Y = np.linspace(0.0, plate.b, Ny)
    p = np.zeros((Nx, Ny))
    for i, xi in enumerate(X):
        for j, yj in enumerate(Y):
            p[i,j] = load.p_from_series(xi, yj)

    # Plot load
    plt.figure()
    plt.contourf(X, Y, p, 100)
    plt.show()