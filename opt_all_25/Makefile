all:
	gcc -march=native \
	main_opt_all_25.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_25.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix25_vec.c \
	../_common_libs/radix25_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

no-opt:
	gcc -march=native -mavx2 \
	main_opt_all_25.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_25.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix25_vec.c \
	../_common_libs/radix25_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3:
	gcc -O3 -mavx2 -march=native -mfma \
	main_opt_all_25.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_25.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix25_vec.c \
	../_common_libs/radix25_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o2:
	gcc -O2 -mavx2 -march=native -mfma \
	main_opt_all_25.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_25.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix25_vec.c \
	../_common_libs/radix25_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3-no-vect:
	gcc -O3 -march=native -fno-tree-vectorize -mfma \
	main_opt_all_25.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_25.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix25_vec.c \
	../_common_libs/radix25_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3-avx2:
	gcc -O3 -mavx2 -mfma \
	main_opt_all_25.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_25.c  \
	../_common_libs/vec_utils.c \
	../_common_libs/radix25_vec.c \
	../_common_libs/radix25_curveParams_vec.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

clean:
	rm out
	