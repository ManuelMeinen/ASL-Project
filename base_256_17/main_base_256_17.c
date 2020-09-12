#include "../_common_libs/radix17_curveParams.h"
#include "../_common_libs/common.h"
#include "curve25519/curve25519_base_256_17.h"
#include "../benchmark_framework/benchmarking.h"
#include "../benchmark_framework/tsc_x86.h"
#include <string.h>
#include <stdlib.h>

int pubKeyGen_base_256_17( char **sk, int num_test_cases ) {
    myInt64 start, end;
	radix17_t Px, Py;
    int cycles = 0;
	for ( int i = 0; i < num_test_cases; ++i ) {
        start = start_tsc();	
        scalar_mult_base_256_17(&Gx_rad17, &Gy_rad17, sk[i],&Px, &Py);
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
    init_curve25519_params_rad17();
	
	// - run benchmarking
	int ok = run_benchmarking( "../benchmark_framework", test_number, run_mode, be_verbose, res_dir_name, &pubKeyGen_base_256_17, &validate_base_256_17, &checkPk_base_256_17, "radix17", "Base line implementation with 256 arithmetic (radix17)", "base_256_17" );
    return !ok;
}

