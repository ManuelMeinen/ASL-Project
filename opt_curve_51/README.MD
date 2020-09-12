# Opt Curve radix51
contains mathematical optimizations in point_add and point_double functions. 
- adds instead of mults for small constants
- pre computations
- optimizations for curve parameters (B==1, remove *B)
- prefer add over sub