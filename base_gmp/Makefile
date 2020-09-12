all:
	gcc \
	main_base_gmp.c \
	../_common_libs/common.c \
	../_common_libs/gmp_curveParams.c \
	curve25519/curve25519_base_gmp.c  \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out -lgmp

no-opt:
	gcc \
	main_base_gmp.c \
	../_common_libs/common.c \
	../_common_libs/gmp_curveParams.c \
	curve25519/curve25519_base_gmp.c  \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out -lgmp

o3:
	gcc -O3 -mavx2 -march=native -mfma \
	main_base_gmp.c \
	../_common_libs/common.c \
	../_common_libs/gmp_curveParams.c \
	curve25519/curve25519_base_gmp.c  \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out -lgmp

o2:
	gcc -O2 -mavx2 -march=native -mfma \
	main_base_gmp.c \
	../_common_libs/common.c \
	../_common_libs/gmp_curveParams.c \
	curve25519/curve25519_base_gmp.c  \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out -lgmp

o3-no-vect:
	gcc -O3 -fno-tree-vectorize \
	main_base_gmp.c \
	../_common_libs/common.c \
	../_common_libs/gmp_curveParams.c \
	curve25519/curve25519_base_gmp.c  \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out -lgmp

o3-avx2:
	gcc -O3 -mavx2 -mfma \
	main_base_gmp.c \
	../_common_libs/common.c \
	../_common_libs/gmp_curveParams.c \
	curve25519/curve25519_base_gmp.c  \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out -lgmp

clean:
	rm out