# PGE 383 Advanced Geomechanics Final Project 

## Wellbore Poroelasticity Simulation

### Overview
In this final project, you will implement a poroelastic simulation of a wellbore in plane strain conditions. You will use the Gridap.jl finite element framework to solve the coupled equations of poroelasticity around a wellbore geometry.

You will simulate the poroelastic response of rock around a wellbore. The wellbore is subjected to a pressure differece from the initial pore pressure in the formation. This creates a pressure gradient that drives fluid flow and causes deformation of the rock matrix.

### Assignment Tasks

1. Complete the implementation of the `wellbore.jl` file using `mandel.jl` as a reference.
2. Use the provided `wellbore.geo` file as the input geometry.
3. Apply the following boundary conditions:
   - At the `top_bottom` boundaries: Fix both $x$ and $y$ displacements ($u_x = 0, u_y = 0$)
   - At the `wellbore` boundary: Apply a constant pressure of $p_b = 30.5$ MPa
   - Initial pressure throughout the domain: $p_0 = 20.0$ MPa
4. Use the material properties and simulation parameters already defined in `wellbore.jl`.
5. Implement the weak form of the poroelastic equations.
6. Set up the solver and time-stepping scheme.
7. Generate output files for visualization.

### Submission Requirements
- Submit your completed `wellbore.jl` file via Canvas.
- Your code should run without errors when executed with the provided `wellbore.geo` file.

### Evaluation Criteria
- Correctness of the implementation
- Proper application of boundary conditions
- Code organization and clarity
- Proper implementation of the weak form
- Successful execution of the simulation

### Note
All submitted code will be executed to verify functionality. Make sure your code runs correctly with the provided geometry file and parameters.

