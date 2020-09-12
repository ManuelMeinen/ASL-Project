#include "../_common_libs/common.h"
#include "../_common_libs/radix25_curveParams.h"
#include "../_common_libs/radix25_curveParams_vec.h"
#include "curve25519/curve25519_base_256_25.h"
#include "../benchmark_framework/benchmarking.h"
#include "../benchmark_framework/tsc_x86.h"
#include <string.h>
#include <stdlib.h>

int pubKeyGen_opt_all_25( char **sk, int num_test_cases ) {
	__m256i Pxlo, Pxhi, Pylo, Pyhi;
    myInt64 start, end;

    int cycles = 0;
	for ( int i = 0; i < num_test_cases; ++i ) {
        start = start_tsc();	
        scalar_mult_opt_all_25(Gx_rad25lo, Gx_rad25hi, Gy_rad25lo, Gy_rad25hi, sk[i], &Pxlo, &Pxhi, &Pylo, &Pyhi);
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
	run_benchmarking( "../benchmark_framework", test_number, run_mode, be_verbose, res_dir_name, &pubKeyGen_opt_all_25, &validate_opt_all_25, &checkPk_opt_all_25, "radix25", "All optimization using raidx25 255-bit arithmetic", "base_opt_all_25" );
    return 0;
}

