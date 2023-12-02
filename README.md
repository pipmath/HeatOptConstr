# Heat Equation Optimiziation with Constraint
MATLAB code to test and perform an optimal reconstruction of the Neumann boundary condition from data in a one dimensional Heat Equation. Includes both unconstrained and constrained formulations of the problem, in the latter case with an exact (via retraction) and an approximate (via projection onto the subspace tangent to the constraint manifold) enforcement of the constraint.

# Files
- `Heat_Optimization.m`: Script performs an optimal reconstruction of the Neumann boundary condition from data in a 1D Heat Equation
- `KappaTest.m`: Script of the Kappa Test for validating gradients and normal elements for optimal reconstruction
- `HeatOptFuncs.m`: Functions required to run validation and optimization contained in `Heat_Optimization.m` and `KappaTest.m`

# How to Use
Default spatial and temporal step sizes can be changed in `HeatOptFuncs.m`, and parameters (such as spatial domain, temporal domain, initial condition, initial guess and true Neumann boundary function) for solving heat equation can be changed in subroutine `heat_init`. Run script `KappaTest.m` to run kappa test and validate gradient/ normal element (in the $L^2$ topology) of problem, using the flags *Constr* and *IC* to toggle the enforcement of the constraint and switch the initial condition of the PDE, respectively.

Run the script `Heat_Optimization.m` to perform optimal reconstruction of the Neumann boundary condition. Script also includes flags *Constr* and *IC*, in addition to optimization parameters that can be adjusted. When constrained problem is solved, constraint is exactly enforced via retraction when zero initial condition is used, and approximately enforced via projection onto the subspace tangent to the constraint manifold when non-zero initial condition is used. In the unconstrained versions of the problems we use the Polak-Ribiere version of the conjugate gradient method to accelerate convergence of iterations. Sobolev gradients are employed to ensure regularity of gradients, then, bracketing and line search methods are used to determine optimal gradient step length for each iteration. 
