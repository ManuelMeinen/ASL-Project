
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pubkeygen_gmp.h"
#include "../test_set_reader.h"

int main( int argc, char **argv ) {
    // - retrieve input parameter
    if ( argc != 2 ) {
    	printf( "Wrong number of parameters specified. Usage: ./create_ground_truths <test_number>\n" );
    	return -1;
    }
    int test_number = atoi( argv[1] );
        
    int num_test_cases, sk_len;
    // - read values from input files
	if ( ! read_test_set( "..", test_number, 0, &num_test_cases, &sk_len ) ) {
		printf( " -- ERROR: reading of input files failed" );
		return -1;
	}
	
	// - create the public keys
	init_pubkeygen_gmp();
	for ( int i = 0; i < num_test_cases; ++i ) {
		if ( !validate_pubkeygen_gmp( secrets_A[i], secrets_B[i] ) ) {
			printf( " -- ERROR: public key generation is invalid" );
			return -1;
		}
		gmp_public_key_gen( secrets_A[i] );
	}

    return 0;
}
