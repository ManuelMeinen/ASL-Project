all:
	gcc -march=native \
	main.c curve25519_opt_51.c \
	../256arithmetic/bignum.c \
	../../_common_libs/radix51_arr.c \
	../../_common_libs/radix51_curveParams_arr.c \
	../../_common_libs/common.c \
	-o out
clean:
	rm out