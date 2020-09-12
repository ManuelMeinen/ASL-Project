micro:
	gcc \
	-march=native \
	main_microbenchmarking.c \
	microbenchmarking.c \
	microbenchmark_aux/vec_constants.c \
	microbenchmark_aux/256_17_bignum.c \
	microbenchmark_aux/256_17_opt_all_bignum.c \
	microbenchmark_aux/256_17_vec_bignum.c \
	microbenchmark_aux/256_17_vec_naive_bignum.c \
	microbenchmark_aux/256_51_bignum.c \
	../_common_libs/radix17.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix51.c \
	-o out
	
	
micro-o3:
	gcc -O3 -mavx2 \
	-march=native \
	main_microbenchmarking.c \
	microbenchmarking.c \
	microbenchmark_aux/vec_constants.c \
	microbenchmark_aux/256_17_bignum.c \
	microbenchmark_aux/256_17_opt_all_bignum.c \
	microbenchmark_aux/256_17_vec_bignum.c \
	microbenchmark_aux/256_17_vec_naive_bignum.c \
	microbenchmark_aux/256_51_bignum.c \
	../_common_libs/radix17.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix51.c \
	-o out
	
micro-o2:
	gcc -O2 -mavx2 \
	-march=native \
	main_microbenchmarking.c \
	microbenchmarking.c \
	microbenchmark_aux/vec_constants.c \
	microbenchmark_aux/256_17_bignum.c \
	microbenchmark_aux/256_17_opt_all_bignum.c \
	microbenchmark_aux/256_17_vec_bignum.c \
	microbenchmark_aux/256_17_vec_naive_bignum.c \
	microbenchmark_aux/256_51_bignum.c \
	../_common_libs/radix17.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix51.c \
	-o out

micro-o3-no-vect:
	gcc -O3 -fno-tree-vectorize -mavx2 \
	-march=native \
	main_microbenchmarking.c \
	microbenchmarking.c \
	microbenchmark_aux/vec_constants.c \
	microbenchmark_aux/256_17_bignum.c \
	microbenchmark_aux/256_17_opt_all_bignum.c \
	microbenchmark_aux/256_17_vec_bignum.c \
	microbenchmark_aux/256_17_vec_naive_bignum.c \
	microbenchmark_aux/256_51_bignum.c \
	../_common_libs/radix17.c \
	../_common_libs/radix17_vec.c \
	../_common_libs/radix51.c \
	-o out