CC = g++
F77 = gfort
LIB =  -lc -lm -lgsl -lgslcblas
#CFL = -fopenmp -Wall -W -Wno-unused -Wno-write-strings
#CFL = -g  -openmp -O2 -Wall
CFL = -g -O2 -fopenmp -ffast-math -funroll-loops -fstrict-aliasing -Wall -W -Wno-unused -Wno-write-strings\
      #-march=opteron-sse3 -mtunecore2 -mpc64  

OBJS = main.o dg_monomial_basis.o mesh_utils.o output.o refine.o gauss_pt.o gl_weight.o\
       hyperbolic_solver.o data_functions.o dg_error.o physics.o compute_reference_mesh_error.o\
       parameters.o electron_physics.o limiting_gauss_pts.o positivity_preserve.o

all:  $(OBJS)
	$(CC) $(CFL) -o driver $(OBJS) $(LIB)
.c.o:
	$(CC) $(CFL) -c $< -o $*.o
.f.o:
	$(F77) $(CFL) -c $< -o $*.o
clean:
	rm -f *.o driver
