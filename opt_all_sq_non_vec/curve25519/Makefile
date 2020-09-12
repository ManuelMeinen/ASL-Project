all:
	gcc -march=native main.c curve25519_base_256_17.c ../256arithmetic/bignum.c \
	../../_common_libs/radix17_vec.c \
	../../_common_libs/radix17_curveParams_vec.c \
	../../_common_libs/common.c \
	../../_common_libs/vec_utils.c \
	-o out
clean:
	rm out