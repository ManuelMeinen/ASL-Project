all: 
	gcc main.c curve25519_base_gmp.c ../../_common_libs/gmp_curveParams.c -o out -lgmp 
test: 
	gcc test.c curve25519.c -o test -lgmp 