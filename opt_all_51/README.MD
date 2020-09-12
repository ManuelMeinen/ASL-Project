# Opt all radix51
Combination of various optimizations:
- mathematical optimizations in point_add and point_double functions. (see readme for opt_curve_51 for more information)
- removed call to reduce function in "add_256"
- replaced radix51_t data type with arrays of long
- removed various copy functions (previously necessary to avoid aliasing issues with struct)
- instead of adding an array with lowest limb set to 1, introduced function add_1_256
- scalar replacement in "mul_256" and "square_256"
- vectorized "sub_256"
