#include "../_common_libs/radix51_curveParams_arr.h"
#include "../_common_libs/common.h"
#include "curve25519/curve25519_opt_51.h"
#include "../benchmark_framework/tsc_x86.h"
#include "../benchmark_framework/benchmarking.h"
#include <string.h>
#include <stdlib.h>

int pubKeyGen_opt_all_51( char **sk, int num_test_cases ) {
	long Px[5], Py[5];
     myInt64 start, end;

    int cycles = 0;
	for ( int i = 0; i < num_test_cases; ++i ) {
        start = start_tsc();	
		scalar_mult_opt_all_51(Gx_rad51_arr, Gy_rad51_arr, sk[i], Px, Py);
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
    init_curve25519_params_rad51_arr();
	
	// - run benchmarking
	run_benchmarking( "../benchmark_framework", test_number, run_mode, be_verbose, res_dir_name, &pubKeyGen_opt_all_51, &validate_opt_all_51, &checkPk_opt_all_51, "radix51", "Combined optimizations using radix51 256-bit arithmetic", "opt_all_51" );
    return 0;
}

