import sys

import matplotlib.pyplot as plt

from plate import Plate
from loads import *
from solver import NavierSolver, LevySolver


def load_input(input_file):
    # Loads the input and returns the solver object

    # Get input
    with open(input_file) as input_handle:
        input_lines = input_handle.readlines()

    # Geometry
    a = float(input_lines[0].split()[-1])
    b = float(input_lines[1].split()[-1])
    h = float(input_lines[2].split()[-1])

    print()
    print("---Plate Dimensions---")
    print()
    print("a:", a)
    print("b:", b)
    print("h:", h)

    # Material properties
    E = float(input_lines[3].split()[-1])
    v = float(input_lines[4].split()[-1])
    print()
    print("---Material Properties---")
    print()
    print("E:", E)
    print("v:", v)

    # Initialize plate
    plate = Plate(a=a, b=b, h=h, E=E, v=v)
    print("D:", plate.D)

    # Locations of interest
    x0 = input_lines[5].split('=')[-1].strip()[1:-1]
    x0 = [float(x.strip()) for x in x0.split(',')]
    y0 = input_lines[6].split('=')[-1].strip()[1:-1]
    y0 = [float(y.strip()) for y in y0.split(',')]

    # Load
    load_type = int(input_lines[7].split()[-1])
    p0 = float(input_lines[8].split()[-1])
    print()
    print("---Load---")
    print()
    if load_type == 0:
        load = UniformLoad(p0=p0, plate=plate)
        print("Type: uniform")
    elif load_type == 1:
        load = SinusoidalLoad(p0=p0, plate=plate)
        print("Type: sinusoidal")
    elif load_type == 2:
        load = HydrostaticLoad(p0=p0, plate=plate)
        print("Type: hydrostatic")
    elif load_type == 3:
        c = float(input_lines[9].split()[-1])
        d = float(input_lines[10].split()[-1])
        x = float(input_lines[11].split()[-1])
        y = float(input_lines[12].split()[-1])
        load = PatchLoad(p0=p0, c=c, d=d, x=x, y=y, plate=plate)
        print("Type: patch")
    print("P0: ", p0)

    # Solver parameters
    m_max = int(input_lines[-4].split()[-1])
    n_max = int(input_lines[-3].split()[-1])
    solver = input_lines[-2].split()[-1]
    BC = input_lines[-1].split()[-1]

    print()
    print("---Solver---")
    print()
    print("Type: ", solver.title())
    print("Boundary Conditions: ", BC)
    print("m_max: ", m_max)

    if solver == "navier":
        if BC != 'SSSS':
            raise IOError("Navier solver cannot be used with BCs other than SSSS! Quitting...")
        solver = NavierSolver(plate, load, m_max, n_max)
        print("n_max: ", n_max)
    else:
        if BC[0] != 'S' or BC[2] != 'S':
            raise IOError("Levy solver must have 'S' BCs on x-faces! Quitting...")
        solver = LevySolver(plate, load, BC, m_max)
        print("Symmetric: ", solver.symmetric)

    return solver


if __name__=="__main__":

    # Load input
    input_file = sys.argv[-1]
    solver = load_input(input_file)

    # Print stresses
    print()
    print("---Solution---")
    print()
    print("Bending stresses at center")
    print("    \u03C3_x: {0:.3f}".format(solver.sigma_x(0.5*solver.plate.a, 0.5*solver.plate.b, 0.5*solver.plate.h)))
    print("    \u03C3_y: {0:.3f}".format(solver.sigma_y(0.5*solver.plate.a, 0.5*solver.plate.b, 0.5*solver.plate.h)))
    print()
    print("Shear stresses at corners")
    print("(0,0):")
    print("    \u03C4_xy: {0:.3f}".format(solver.tau_xy(0.0, 0.0, 0.5*solver.plate.h)))
    print("(a,0):")
    print("    \u03C4_xy: {0:.3f}".format(solver.tau_xy(solver.plate.a, 0.0, 0.5*solver.plate.h)))
    print("(a,b):")
    print("    \u03C4_xy: {0:.3f}".format(solver.tau_xy(solver.plate.a, solver.plate.b, 0.5*solver.plate.h)))
    print("(0,b):")
    print("    \u03C4_xy: {0:.3f}".format(solver.tau_xy(0.0, solver.plate.b, 0.5*solver.plate.h)))

    # Get max deflection
    x_w_max, y_w_max, w_max = solver.get_maximum_w()
    print()
    print("Deformations")
    print("    Maximum: {0:.5f} at ({1:.5f}, {2:.5f})".format(w_max, x_w_max, y_w_max))
    print("    At center: {0:.5f} at ({1:.1f}, {2:.1f})".format(solver.w(0.5*solver.plate.a, 0.5*solver.plate.b), 0.5*solver.plate.a, 0.5*solver.plate.b))

    print()
    print("Average corner reactions")
    print("    {0:.2f} per corner".format(solver.corner_reactions()))

    print()
    print("Maximum stresses")
    print("    \u03C3_x: {3:.3f} at ({0:.2f}, {1:.2f}, {2:.2f})".format(*solver.get_maximum_sigma_x()))
    print("    \u03C3_y: {3:.3f} at ({0:.2f}, {1:.2f}, {2:.2f})".format(*solver.get_maximum_sigma_y()))
    print("    \u03C4_xz: {3:.3f} at ({0:.2f}, {1:.2f}, {2:.2f})".format(*solver.get_maximum_tau_xz()))
    print("    \u03C4_yz: {3:.3f} at ({0:.2f}, {1:.2f}, {2:.2f})".format(*solver.get_maximum_tau_yz()))

    #solver.plot_deflection_field()