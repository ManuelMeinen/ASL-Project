# Opt Curve radix17
contains scalar replacement and vectoritation in 256 arithmetic. 
- radix17 number represented as 2 __256i vectors
- add and sub use 8x vectors 
- mul and square use 4x vectors
- mul and square make use of patterns to group the operations 
- square still introduces lots of permutes compared to the actual computations
- reduce not vectorizable 

