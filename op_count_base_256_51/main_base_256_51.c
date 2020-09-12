#include "../_common_libs/radix51_curveParams.h"
#include "../_common_libs/common.h"
#include "curve25519/curve25519_base_256_51.h"
#include "../benchmark_framework/benchmarking.h"
#include <string.h>
#include <stdlib.h>

void pubKeyGen_base_256_51( char **sk, int num_test_cases ) {
	radix51_t Px, Py;
	for ( int i = 0; i < num_test_cases; ++i ) {
		// TODO: replace with scalar_mult function
		pub_key_gen(sk[i], &Px, &Py);
	}
}

int main( int argc, char **argv ) {
    // - retrieve input parameter
    if ( argc != 5 ) {
        printUsageHelp();
    	return -1;
    }
    int  test_number   = atoi( argv[1] );
    char run_mode      = argv[2][0];
    int  be_verbose    = atoi( argv[3] );
    char *res_dir_name = argv[4];
    
    // - initialize curve parameters
    init_curve25519_params_rad51();
	
	// - run benchmarking
	run_benchmarking( "../benchmark_framework", test_number, run_mode, be_verbose, res_dir_name, &pubKeyGen_base_256_51, &validate, &checkPk, "radix51", "Base line implementation with 256 arithmetic (radix51)", "base_256_51" );
    return 0;
}

