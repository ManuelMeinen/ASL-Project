#include "../_common_libs/radix17_curveParams.h"
#include "../_common_libs/common.h"
#include "curve25519/curve25519_base_256_17.h"
#include "../benchmark_framework/benchmarking.h"
#include "256arithmetic/bignum.h"
#include "../_common_libs/radix17.h"
#include <string.h>
#include <stdlib.h>

int public_key_gen_base_256_17_rev( char **sk, int num_test_cases, radix17_t *Px, radix17_t *Py){
        radix17_t *R0x, *R0y, *R1x, *R1y;
    R0x = (radix17_t*)malloc(sizeof(radix17_t)*num_test_cases);
    R0y = (radix17_t*)malloc(sizeof(radix17_t)*num_test_cases);
    R1x = (radix17_t*)malloc(sizeof(radix17_t)*num_test_cases);
    R1y = (radix17_t*)malloc(sizeof(radix17_t)*num_test_cases);
   
    for(int i = 0; i < num_test_cases; i++){
        set_ui_256(&R0x[i], 0);
        set_ui_256(&R0y[i], 0);
        copy_256(&R1x[i], Gx_rad17);
        copy_256(&R1y[i], Gy_rad17);
    }

    for ( int j = 0; j < sk_length; ++j ) {
        radix17_t resX, resY;
        for(int i=0; i<num_test_cases; i++){
            set_ui_256(&resX, 0);
            set_ui_256(&resY, 0);
            if(sk[i][j]=='0'){
                // R1 <- add (R0, R1)
                point_add_base_256_17(&resX, &resY, &R0x[i], &R0y[i], &R1x[i], &R1y[i]);
                copy_256(&R1x[i], resX);
                copy_256(&R1y[i], resY);

                set_ui_256(&resX, 0);
                set_ui_256(&resY, 0);
                // R0 <- double(R0)
                point_double_base_256_17(&resX, &resY, &R0x[i], &R0y[i]);
                copy_256(&R0x[i], resX);
                copy_256(&R0y[i], resY);
            }else{
                // R0 <- add(R0, R1)
                point_add_base_256_17(&resX, &resY, &R0x[i], &R0y[i], &R1x[i], &R1y[i]);
                copy_256(&R0x[i], resX);
                copy_256(&R0y[i], resY);

                set_ui_256(&resX, 0);
                set_ui_256(&resY, 0);
                // R1 <- double(R1)
                point_double_base_256_17(&resX, &resY, &R1x[i], &R1y[i]);
                copy_256(&R1x[i], resX);
                copy_256(&R1y[i], resY);
            }
        }
        
    }

    *Px = *R0x;
    *Py = *R0y;
    return 0;
}

void pubKeyGen_base_256_17_mem_opt( char **sk, int num_test_cases ) {
	radix17_t *Px, *Py;
    Px = (radix17_t*)malloc(sizeof(radix17_t)*num_test_cases);
    Py = (radix17_t*)malloc(sizeof(radix17_t)*num_test_cases);
    public_key_gen_base_256_17_rev(sk, num_test_cases, Px, Py);
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
	int ok = run_benchmarking( "../benchmark_framework", test_number, run_mode, be_verbose, res_dir_name, &pubKeyGen_base_256_17_mem_opt, &validate_base_256_17, &checkPk_base_256_17, "radix17", "Base line implementation with 256 arithmetic mem opt(radix17)", "pubKeyGen_base_256_17_mem_opt" );
    return !ok;
}

