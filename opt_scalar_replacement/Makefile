all:
	gcc \
	main_opt_scalar_replacement.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/radix17.c \
	../_common_libs/radix17_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

no-opt:
	gcc \
	main_opt_scalar_replacement.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/radix17.c \
	../_common_libs/radix17_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3:
	gcc -O3 -mavx2 -march=native -mfma \
	main_opt_scalar_replacement.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/radix17.c \
	../_common_libs/radix17_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o2:
	gcc -O2 -mavx2 -march=native -mfma \
	main_opt_scalar_replacement.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/radix17.c \
	../_common_libs/radix17_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3-no-vect:
	gcc -O3 -fno-tree-vectorize -march=native -mfma \
	main_opt_scalar_replacement.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/radix17.c \
	../_common_libs/radix17_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3-avx2:
	gcc -O3 -mavx2 -mfma \
	main_opt_scalar_replacement.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/radix17.c \
	../_common_libs/radix17_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

autoinline:
	gcc -finline-small-functions -finline-functions-called-once -findirect-inlining \
	 -march=native \
	main_opt_scalar_replacement.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_17.c  \
	../_common_libs/radix17.c \
	../_common_libs/radix17_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out
clean:
	rm out
	