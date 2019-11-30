# FOCEllipticEquationSolver
Fourth order compact solver for a general elliptic equation 
a(x,y,z)p_{xx} + b(x,y,z)p_{yy} + c(x,y,z)p_{zz} 
+ d(x,y,z)p_{x} + e(x,y,z)p_{y} + f(x,y,z)p_{z} + K(x,y,z) p = g(x,y,z) 
with homogeneous boundary conditions on the domain [0,1]^{3} using multigrid methods. 

By compact, we mean that the implicit finite difference scheme 
uses a 3 x 3 x 3 stiencil when centered at the point (x_i, y_j, z_k) 

We will assume that a, b, c, d, e, f, K, p, g are all smooth functions where a, b, and c are all nonegative (or all nonpositive) 

We will give a short discrption of the modules 
I. TestEPDE.f90 
This is the main subroutine which actually solves the fourth order compact elliptic equation. 
In the preamble we define the settings for the multigrid methods. The following are the options for the 
multigrid method 

A. smoother - the smoother will be either pointwise Gauss-Seidel ('GS') or 
Four Color Gauss-Seidel ('FC-GS'). 
B. cycle Type - the cycle type will be either the V, W, or F cycles ('V', 'W', 'F'). 
C. Test Problem - we fabricated 3 test problems to determine if fourth order convergence is met. 
Choose 1,2, or 3 to choose among the test problems 
D. numGrids - this integer determines the number of auxiliary grids used for multigrid methods. The number of (uniform) subintervals (per axis) are 2^(numGrids) on the finest grid, 2^(numGrids - 1) on the next finest grid, ..., 2 on the coarsest grid. 
E. preSweeps - number of smoothing steps on each grid along the descending portion. We recomend making    preSweeps <= 3 
F. postSweeps - number of smoothing steps on each grid along the ascending portion. We recomend making postSweeps <= 3 
G. totalIterations - the total number of multigrid steps taken. We make the zero vector the initial guess.
H. orderRestriction - The degree of accuracy of the restriction operator. Use 0 for the injection restriction operator or 2 for the Full-Weighting restriction operator 
I. orderInterpolation - The degree of accuracy of the interpolation operator. Use 2 for trilinear interpolation or 4 for tricubic interpolation. 


