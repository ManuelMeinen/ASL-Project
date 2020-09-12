all:
	gcc \
	main_base_256_25.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_25.c  \
	../_common_libs/radix25.c \
	../_common_libs/radix25_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

no-opt:
	gcc \
	main_base_256_25.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_25.c  \
	../_common_libs/radix25.c \
	../_common_libs/radix25_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

o3:
	gcc -O3 \
	main_base_256_25.c \
	../_common_libs/common.c \
	256arithmetic/bignum.c \
	curve25519/curve25519_base_256_25.c  \
	../_common_libs/radix25.c \
	../_common_libs/radix25_curveParams.c \
	../benchmark_framework/benchmarking.c \
	../benchmark_framework/test_set_reader.c \
	-o out

clean:
	rm out