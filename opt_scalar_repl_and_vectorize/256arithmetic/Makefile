all: 
	gcc --version
	gcc main.c bignum.c ../../_common_libs/vec_utils.c -o out -march=native -mavx2

test: 
	gcc test.c bignum.c ../../_common_libs/vec_utils.c -o test -lgmp -march=native