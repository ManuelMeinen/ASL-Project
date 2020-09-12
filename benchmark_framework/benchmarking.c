#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test_set_reader.h"
#include "../_common_libs/common.h"
#include "benchmarking.h"

#include "tsc_x86.h"

#define NUM_CYCLES 15000000

int run_benchmarking( char *rel_path_to_benchmarking_dir, int test_set_no, char run_mode, int be_verbose, char *res_dir_name, pubKeyGen_func *f, validate_func *v, checkPk_func *c, char *pubKeyFormat, char *name, char *id ) {
	myInt64 start, end;
	int is_valid_comp = 0;
	int is_pkCheck_ok = 0;

	int secretK_length, num_test_cases;

	char** pks;

    FILE *fp;
    char res_fName[200];

	// - read input values from test set 'test_set_no'
	if ( ! read_test_set( rel_path_to_benchmarking_dir, test_set_no, 1, &num_test_cases, &secretK_length ) ) {
		return 0;
	}
	set_sk_length( secretK_length );
	
	if ( be_verbose ) { printf( " Running benchmarking of function '%s' on test set %d\n", name, test_set_no ); }
	if ( run_mode == 'v' ) {
		if ( be_verbose ) { printf( "    validation-only mode\n\n" ); }
	}
	else if ( run_mode == 'b' ) {
		if ( be_verbose ) { printf( "    benchmarking-only mode\n\n" ); }
	}
	else if ( run_mode == 'f' ) {
		if ( be_verbose ) { printf( "    running full validation and benchmarking\n\n" ); }
	}
	else {
		printf( "\n -- ERROR: unknown run mode '%c' encountered\n", run_mode );
		return 0;
	}
	
	if ( run_mode == 'v' || run_mode == 'f' ) {
		// - check public key generation results
		if ( be_verbose ) { printf( "  - checking public key generation\n" ); }
		if ( ! strncmp( pubKeyFormat, "decimal", 10 ) ) {
			pks = public_keys_dec;
		}
		else if ( ! strncmp( pubKeyFormat, "radix51", 10 ) ) {
			pks = public_keys_r51;
		}
		else if ( ! strncmp( pubKeyFormat, "radix17", 10 ) ) {
			pks = public_keys_r17;
		}
		else if ( ! strncmp( pubKeyFormat, "radix25", 10 ) ) {
			pks = public_keys_r25;
		}
		else {
			printf( "\n -- ERROR: currently unsupported public key format found '%s'\n", pubKeyFormat );
			return 0;
		}
		if ( pks == NULL ) {
			printf( "\n -- ERROR: No corresponding public key file found for format '%s'\n", pubKeyFormat );
			return 0;
		}
		
		is_pkCheck_ok = (*c) ( secrets_A, pks, num_test_cases );
		if ( ! is_pkCheck_ok ) {
			if ( be_verbose ) {
				printf( "    FALSE\n" );
				printf( "Exiting...\n\n" );
			}
			else {
				printf( "\n Public key generation check failed! \n" );
			}
			return 0;
		}
		else {
			if ( be_verbose ) { printf( "    OK\n" ); }
		}
		if ( be_verbose ) { printf( "\n" ); }
		
		// - validate
		if ( be_verbose ) { printf( "  - validating public key generation with full Diffie-Hellman key exchange\n" ); }
		is_valid_comp = (*v) ( secrets_A, secrets_B, num_test_cases );
		if ( ! is_valid_comp ) {
			if ( be_verbose ) {
				printf( "    FALSE\n" );
				printf( "Exiting...\n\n" );
			}
			else {
				printf( "\n Validation failed! \n" );
			}
			return 0;
		}
		else {
			if ( be_verbose ) { printf( "    OK\n" ); }
		}
		if ( be_verbose ) { printf( "\n" ); }
	}
	
	if ( run_mode == 'b' || run_mode == 'f' ) {
		if ( be_verbose ) { printf( "  - benchmarking\n" ); }
		
		
		double performance;

		int p =0;
		for ( p = 0; p < 50; p++ ) {
			
			performance = (*f) ( secrets_A, 1 );
			
			//performance = (( double ) end); 
			if (performance >= NUM_CYCLES){
				break;
			}
		}
		
		
		
		long num_runs = 16;
		// -- get performance
		// start = start_tsc();		
		for ( int j = 0; j < num_runs; j++ ) {
			performance += (*f) ( secrets_A, num_test_cases );
		}
		// end = stop_tsc( start );
		performance = (( double ) performance ) / num_runs;

		if ( be_verbose ) {
			printf( "    %f cycles\n\n", performance );
		}
		else {
			printf( "%f", performance );
		}
		
		// -- write performance to file
		strcpy( res_fName, res_dir_name );
		strcat( res_fName, id );
		strcat( res_fName, ".csv" );
		fp = fopen( res_fName, "a" );
		if ( fp == NULL ) {
			printf( " -- ERROR: unable to find output file '%s'\n", res_fName );
			return 0;
		}
		fprintf( fp, "%d;%f\n", num_test_cases, performance );
		fclose( fp );
	}

	return 1;
}
