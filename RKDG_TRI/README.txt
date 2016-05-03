========================================================================
=== 2D RKDG code for hyperbolic equations on unstructured grids  ===
========================================================================
===  Prince Chidyagwai (Loyola University Maryland)  
===  pchidyagwai@loyola.edu
========================================================================

===================================
Description:
- The package includes an implementation of the RKDG method in C for the moment approximations to transport equations, 
working for piecewise constant, linear or quadratic polynomials, yielding an overall first, second or third order convergence rate. 
- The code works with arbitrary triangular meshes in 2 space dimensions. Some sample meshes are provided in /meshes
- Flux functions (closures), boundary and initial conditions can be specified by the user.
- A slope limiter that can be run either primitive or characteristic variables is implemented
- Runs in parallel with OPENMP.

===================================
Prerequesites:
See MAKEFILE for required libraries
For the visualization
- MATLAB
Additionally, to render movies 
- MEncoder ( www.mplayerhq.hu/ )
is being used. It comes with the popular program mplayer. 

===================================
Installation:
Simply type
> make
to compile the software. After changes to the code, one should use
> make clean
before compiling the project again.

===================================
Configuration:

- All parameters are prescribed in the file parameters.c, for detailed descriptions refer to the comments in the file.
- Flux functions, boundary and initial conditions are prescribed in data_functions.c

===================================
Usage:

The program is being started typing
> ./driver

Meshfile describing a triangulation (typically in the folder meshes/) needs to be specified. Regular meshes ranging from 32 to 32768 triangles are included as are 3 irregular meshes and various others.

Visualization is provided in MATLAB. See the reps. documentation in the files solution/plot_dg_solution.m and solution/plot_dg_movie.m .


===================================
Error Computation
There are 2 ways to compute errors using the code:
1. Providing an exact solution at the final time T.
2. Computing a reference solution on a refined grid.
-The code automatically stores the solution coefficients for a given mesh in the directory solution_coefs
-These coefficients can then be used as a reference solution.

===================================
