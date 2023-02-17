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
    print("a:", a)
    print("b:", b)
    print("h:", h)

    # Material properties
    E = float(input_lines[3].split()[-1])
    v = float(input_lines[4].split()[-1])
    print()
    print("---Material Properties---")
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
        load = PatchLoad(p0=p0, c=c, d=d, x=x, y=y)
        print("Type: patch")
    print("P0:", p0)

    # Solver parameters
    m_max = int(input_lines[-4].split()[-1])
    n_max = int(input_lines[-3].split()[-1])
    solver = input_lines[-2].split()[-1]
    BC = input_lines[-1].split()[-1]

    print()
    print("---Solver---")
    print("Type: ", solver.title())
    print("Boundary Conditions: ", BC)
    print("m_max:", m_max)
    print("n_max:", n_max)

    if solver == "navier":
        if BC != 'SSSS':
            raise IOError("Navier solver cannot be used with BCs other than SSSS! Quitting...")
        solver = NavierSolver(plate, load, m_max, n_max)
    else:
        if BC[0] != 'S' or BC[2] != 'S':
            raise IOError("Levy solver must have 'S' BCs on x-faces! Quitting...")
        solver = LevySolver(plate, load, BC, m_max)

    return solver


if __name__=="__main__":

    # Load input
    input_file = sys.argv[-1]
    solver = load_input(input_file)

    # Print stresses
    print()
    print("---Solution---")
    print("Maximum Plane Stresses")
    print("At center:")
    print("    \u03C3x: ", solver.sigma_x(0.5*solver.plate.a, 0.5*solver.plate.b, 0.5*solver.plate.h))
    print("    \u03C3y: ", solver.sigma_y(0.5*solver.plate.a, 0.5*solver.plate.b, 0.5*solver.plate.h))
    print("At (0,0):")
    print("    \u03C4xy: ", solver.tau_xy(0.0, 0.0, 0.5*solver.plate.h))
    print("At (a,0):")
    print("    \u03C4xy: ", solver.tau_xy(solver.plate.a, 0.0, 0.5*solver.plate.h))
    print("At (a,b):")
    print("    \u03C4xy: ", solver.tau_xy(solver.plate.a, solver.plate.b, 0.5*solver.plate.h))
    print("At (0,b):")
    print("    \u03C4xy: ", solver.tau_xy(0.0, solver.plate.b, 0.5*solver.plate.h))

    print("Maximum Transverse Stresses")
    print("At (0,0):")
    print("    \u03C4xz: ", solver.tau_xz(0.0, 0.0, 0.0))
    print("    \u03C4yz: ", solver.tau_yz(0.0, 0.0, 0.0))
    print("At (a/2,0):")
    print("    \u03C4xz: ", solver.tau_xz(0.5*solver.plate.a, 0.0, 0.0))
    print("    \u03C4yz: ", solver.tau_yz(0.5*solver.plate.a, 0.0, 0.0))
    print("At (a,0):")
    print("    \u03C4xz: ", solver.tau_xz(solver.plate.a, 0.0, 0.0))
    print("    \u03C4yz: ", solver.tau_yz(solver.plate.a, 0.0, 0.0))
    print("At (a,b/2):")
    print("    \u03C4xz: ", solver.tau_xz(solver.plate.a, 0.5*solver.plate.b, 0.0))
    print("    \u03C4yz: ", solver.tau_yz(solver.plate.a, 0.5*solver.plate.b, 0.0))
    print("At (a,b):")
    print("    \u03C4xz: ", solver.tau_xz(solver.plate.a, solver.plate.b, 0.0))
    print("    \u03C4yz: ", solver.tau_yz(solver.plate.a, solver.plate.b, 0.0))
    print("At (a/2,b):")
    print("    \u03C4xz: ", solver.tau_xz(0.5*solver.plate.a, solver.plate.b, 0.0))
    print("    \u03C4yz: ", solver.tau_yz(0.5*solver.plate.a, solver.plate.b, 0.0))
    print("At (0,b):")
    print("    \u03C4xz: ", solver.tau_xz(0.0, solver.plate.b, 0.0))
    print("    \u03C4yz: ", solver.tau_yz(0.0, solver.plate.b, 0.0))
    print("At (0,b/2):")
    print("    \u03C4xz: ", solver.tau_xz(0.0, 0.5*solver.plate.b, 0.0))
    print("    \u03C4yz: ", solver.tau_yz(0.0, 0.5*solver.plate.b, 0.0))

    solver.plot_deflection_field()