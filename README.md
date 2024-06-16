## About
Implementation of Interior Penalty Discontinuous Galerkin methods (IP-DG) for diffusion problems using arbitrary high order Lagrange basis. 

MATLAB R2019 at least is required.

## Quick start
Clone the repository
```[bash]
git clone https://github.com/mobenabd/DG-2D.git
cd DG-2D
```

## File structure

```
DG-2D
|   LICENSE
│   README.md
|   mainPoisson.m                  (solve 2D elliptic equations)
|   mainHeat.m                     (solve 2D parabolic equations with backward Euler scheme in time)
|   mainObstacle.m                 (solve 2D linearized obstacle problem)
|   mainConvexification.m          (solve convexification problems)
|   mainOptimalTransport.m         (Solve 2D L2-MKP problems)
|
└───obstacle_problem               (obstacle problem fixed point algorithm related files)
|   |   ...                     
|
└───optimal_transport              (optimal transport related files)
|   |   ...  
|
└───src_DG                         (IP-DG solver core functions)
|   |   ...
|
└───test                           (various test mains)
|   |   ...
```

