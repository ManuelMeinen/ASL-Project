# Mem Opt radix17
For mem-opt we extended the problem to generating multiple pks for multiple sks. 
The sks are stored as matrix where the first row is the first sk as bitstring.
In this version we iterate column-wise instead of row-wise (use 1st bit of all sks, then 2nd bit of all ..).
The idea behind it is that column-wise access should produce more misses (capacity misses) when the matrix of all sks does not fit in L1 (64KB) 
We profiled the cache with valgrind. 
Result: 0% miss rate in both versions 
