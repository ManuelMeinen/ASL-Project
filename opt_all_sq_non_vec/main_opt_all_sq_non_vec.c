#include "../_common_libs/radix17_curveParams_vec.h"
#include "../_common_libs/common.h"
#include "curve25519/curve25519_base_256_17.h"
#include "../benchmark_framework/benchmarking.h"
#include "../benchmark_framework/tsc_x86.h"
#include <immintrin.h>
#include <string.h>
#include <stdlib.h>

int pubKeyGen_opt_all( char **sk, int num_test_cases ) {
	__m256i Pxlo, Pxhi, Pylo, Pyhi;
    myInt64 start, end;

    int cycles = 0;
	for ( int i = 0; i < num_test_cases; ++i ) {
        start = start_tsc();	
        scalar_mult_base_256_17(Gx_rad17lo, Gx_rad17hi, Gy_rad17lo, Gy_rad17hi, sk[i], &Pxlo, &Pxhi, &Pylo, &Pyhi);
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
	run_benchmarking( "../benchmark_framework", test_number, run_mode, be_verbose, res_dir_name, &pubKeyGen_opt_all, &validate_base_256_17, &checkPk_base_256_17, "radix17", "All optimisations in radix17 256 arithmetics", "opt_all" );
    return 0;
}

