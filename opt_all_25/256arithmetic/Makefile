all: 
	gcc -march=native ../../_common_libs/radix25.c ../../_common_libs/radix25_vec.c \
	../../_common_libs/vec_utils.c \
	main.c bignum.c -o out  
test: 
	gcc test.c bignum.c -o test -lgmp 