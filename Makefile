# make for FEM

# Directories
SRC_DIR = ./src

# Compiler
COMPILER = gfortran

# Flags
FLAGS = -O2 -fdefault-real-8 -fbounds-check -fbacktrace

# Program name
PROGRAM = fem.exe

default:
	$(COMPILER) $(FLAGS) -o $(PROGRAM) \
	src/math.f95 \
	src/linalg.f95 \
	src/json.f95 \
	src/json_xtnsn.f95 \
	src/gauss_integration.f95 \
	src/element.f95 \
	src/plate.f95 \
	src/main.f95