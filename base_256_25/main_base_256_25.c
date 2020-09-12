#include "../_common_libs/radix25_curveParams.h"
#include "../_common_libs/common.h"
#include "curve25519/curve25519_base_256_25.h"
#include "../benchmark_framework/benchmarking.h"
#include <string.h>
#include <stdlib.h>
#include "../benchmark_framework/tsc_x86.h"

int pubKeyGen_base_256_25( char **sk, int num_test_cases ) {
    myInt64 start, end;
	int cycles = 0;
	
	radix25_t Px, Py;
	for ( int i = 0; i < num_test_cases; ++i ) {
        start = start_tsc();	
		scalar_mult_base_256_25(&Gx_rad25, &Gy_rad25, sk[i], &Px, &Py);
        end = stop_tsc( start );
        cycles += end;
	}
	return cycles;
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
    init_curve25519_params_rad25();
	
	// - run benchmarking
	int ok = run_benchmarking( "../benchmark_framework", test_number, run_mode, be_verbose, res_dir_name, &pubKeyGen_base_256_25, &validate_base_256_25, &checkPk_base_256_25, "radix25", "Base line implementation with 256 arithmetic (radix25.5)", "base_256_25" );
    return !ok;
}

