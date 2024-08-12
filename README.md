# Heat Equation Optimiziation with Constraint
MATLAB code to test and perform an optimal reconstruction of the Neumann boundary condition from data in a one dimensional Heat Equation. Includes both unconstrained and constrained formulations of the problem, in the latter case with an exact (via retraction) and an approximate (via projection onto the subspace tangent to the constraint manifold) enforcement of the constraint.

# Files
- `Heat_Optimization.m`: Script performs an optimal reconstruction of the Neumann boundary condition from data in a 1D Heat Equation
- `KappaTest.m`: Script of the Kappa Test for validating gradients and normal elements for optimal reconstruction
- `HeatOptFuncs.m`: Functions required to run validation and optimization contained in `Heat_Optimization.m` and `KappaTest.m`

# How to Use
Default spatial and temporal step sizes can be changed in `HeatOptFuncs.m`, and parameters (such as spatial domain, temporal domain, initial condition, initial guess and true Neumann boundary function) for solving heat equation can be changed in subroutine `heat_init`. Run script `KappaTest.m` to run kappa test and validate gradient/ normal element (in the $L^2$ topology) of problem, using the flags *Constr* and *IC* to toggle the enforcement of the constraint and switch the initial condition of the PDE, respectively.

Run the script `Heat_Optimization.m` to perform optimal reconstruction of the Neumann boundary condition. Script also includes flags *Constr* and *IC*, in addition to optimization parameters that can be adjusted. When constrained problem is solved, constraint is exactly enforced via retraction when zero initial condition is used, and approximately enforced via projection onto the subspace tangent to the constraint manifold when non-zero initial condition is used. In the unconstrained versions of the problems we use the Polak-Ribiere version of the conjugate gradient method to accelerate convergence of iterations. Sobolev gradients are employed to ensure regularity of gradients, then, bracketing and line search methods are used to determine optimal gradient step length for each iteration. 

# Citing
Work has been published in the Journal of Computational Physics (JCP). The paper can be found [here](https://doi.org/10.1016/j.jcp.2024.113298).

Pritpal Matharu and Bartosz Protas., (2024). *Adjoint-Based Enforcement of State Constraints in PDE Optimization Problems.* Journal of Computational Physics 517, 113298, 2024, https://doi.org/10.1016/j.jcp.2024.113298

Bibtex:
```
  @article{Matharu_JCP517_2024,
  title = {{Adjoint-based enforcement of state constraints in PDE optimization problems}},
  journal = {J. Comput. Phys.},
  fjournal = {Journal of Computational Physics},
  volume = {517},
  pages = {113298},
  year = {2024},
  issn = {0021-9991},
  doi = {https://doi.org/10.1016/j.jcp.2024.113298},
  url = {https://www.sciencedirect.com/science/article/pii/S0021999124005461},
  author = {Matharu, P. and Protas, B.},
  fauthor = {Pritpal Matharu and Bartosz Protas},
  keywords = {PDE optimization, Adjoint analysis, State constraints, Heat transfer, Turbulence},
  abstract = {This study demonstrates how the adjoint-based framework traditionally used to compute gradients in PDE optimization problems can be extended to handle general constraints on the state variables. This is accomplished by constructing a projection of the gradient of the objective functional onto a subspace tangent to the manifold defined by the constraint. This projection is realized by solving an adjoint problem defined in terms of the same adjoint operator as used in the system employed to determine the gradient, but with a different forcing. We focus on the “optimize-then-discretize” paradigm in the infinite-dimensional setting where the required regularity of both the gradient and of the projection is ensured. The proposed approach is illustrated with two examples: a simple test problem describing optimization of heat transfer in one direction and a more involved problem where an optimal closure is found for a turbulent flow described by the Navier-Stokes system in two dimensions, both considered subject to different state constraints. The accuracy of the gradients and projections computed by solving suitable adjoint systems is carefully verified and the presented computational results show that the solutions of the optimization problems obtained with the proposed approach satisfy the state constraints with a good accuracy, although not exactly.}
  }
```
