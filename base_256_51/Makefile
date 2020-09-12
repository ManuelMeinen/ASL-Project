all:
	gcc \
	main_base_256_51.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_51.c  \
	../_common_libs/radix51.c \
	../_common_libs/radix51_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

no-opt:
	gcc \
	main_base_256_51.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_51.c  \
	../_common_libs/radix51.c \
	../_common_libs/radix51_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3:
	gcc -O3 -mavx2 -march=native -mfma \
	main_base_256_51.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_51.c  \
	../_common_libs/radix51.c \
	../_common_libs/radix51_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o2:
	gcc -O2 -mavx2 -march=native -mfma \
	main_base_256_51.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_51.c  \
	../_common_libs/radix51.c \
	../_common_libs/radix51_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3-no-vect:
	gcc -O3 -fno-tree-vectorize -march=native -mfma \
	main_base_256_51.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_51.c  \
	../_common_libs/radix51.c \
	../_common_libs/radix51_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3-avx2:
	gcc -O3 -mavx2 -mfma \
	main_base_256_51.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_51.c  \
	../_common_libs/radix51.c \
	../_common_libs/radix51_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

clean:
	rm out
	