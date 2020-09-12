all:
	gcc main.c \
	-march=native \
	curve25519_base_256_25.c \
	../256arithmetic/bignum.c \
	../../_common_libs/radix25_vec.c \
	../../_common_libs/radix25_curveParams_vec.c \
	../../_common_libs/vec_utils.c \
	../../_common_libs/common.c \
	-o out
clean:
	rm out