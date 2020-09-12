#include "../_common_libs/gmp_curveParams.h"
#include "../_common_libs/common.h"
#include "curve25519/curve25519_base_gmp.h"
#include "../benchmark_framework/benchmarking.h"
#include "../benchmark_framework/tsc_x86.h"
#include <string.h>
#include <stdlib.h>

int pubKeyGen_base_gmp( char **sk, int num_test_cases ) {
    myInt64 start, end;
	mpz_t Pkx, Pky;
    mpz_init(Pkx);
    mpz_init(Pky);

    int cycles = 0;
	for ( int i = 0; i < num_test_cases; ++i ) {
        start = start_tsc();	
        scalar_mult_base_gmp(Gx_gmp, Gy_gmp, sk[i], Pkx, Pky);
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
    init_curve25519_params_gmp();
	
	// - run benchmarking
	run_benchmarking( "../benchmark_framework", test_number, run_mode, be_verbose, res_dir_name, &pubKeyGen_base_gmp, &validate_base_gmp, &checkPk_base_gmp, "decimal", "Base line implementation using GMP library", "base_gmp" );
    return 0;
}
