#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <immintrin.h>
#include <time.h>

#include "../_common_libs/common.h"
#include "microbenchmarking.h"

//#include "../_common_libs/vec_utils.h"
#include "microbenchmark_aux/vec_constants.h"

#include "../_common_libs/radix17.h"
#include "../_common_libs/radix51.h"
#include "../_common_libs/radix17_vec.h"

#include "microbenchmark_aux/256_17_bignum.h"
#include "microbenchmark_aux/256_51_bignum.h"
#include "microbenchmark_aux/256_17_vec_bignum.h"
#include "microbenchmark_aux/256_17_vec_naive_bignum.h"
#include "microbenchmark_aux/256_17_opt_all_bignum.h"

#include "tsc_x86.h"

#define NUM_WARMUP_RUNS 50
#define NUM_RUNS 16

// ---------------------------------------------------------------------------------------
//   Initialize and handling of test data
// ---------------------------------------------------------------------------------------

void init_rand() {
	srand( (int) time( NULL ) );
}

int rand_int() {
	return rand() % ( ( 1 << 17 ) - 1 );
}

void print_data( int** test_data, int num_test_cases ) {
	for ( int i = 0; i < num_test_cases; i++ ) {
		for ( int j = 0; j < 16; j++ ) {
			printf( "%d ", test_data[i][j] );
		}
		printf( "\n" );
	}
}

void init_data( int** test_data, int num_test_cases ) {
	for ( int i = 0; i < num_test_cases; i++ ) {
		for ( int j = 0; j < 15; j++ ) {
			test_data[i][j] = rand_int(); //printf( "%d ", test_data[i][j] );
		}
		test_data[i][15] = 0; //printf( "%d\n", test_data[i][15] );
	}
}

// ---------------------------------------------------------------------------------------
//  Output handling (writing to file and to stdOut)
// ---------------------------------------------------------------------------------------

void print_mb_start( char* name ) {
	printf( " Running microbenchmarking for '%s'\n", name );
}

void print_perf( char* op_name, double performance ) {
	printf( "   %s: %f cycles\n", op_name, performance );
}

void write_perf_to_file( FILE *fp, double performance, char* op_name ) {
	fprintf( fp, "%s;%f\n", op_name, performance );
}

// ---------------------------------------------------------------------------------------
//  Performance calculations
// ---------------------------------------------------------------------------------------

typedef void ( single_op_rad17_func )( radix17_t*, radix17_t* );
typedef void ( double_op_rad17_func )( radix17_t*, radix17_t*, radix17_t* );

double get_perf_single_op_rad17( int** test_data, int num_test_cases, single_op_rad17_func *f ) {
	myInt64 start, end;
	double num_cycles = 0;
	
	radix17_t a, r;
	
	for ( int i = 0; i < num_test_cases; i++ ) {
		for ( int j = 0; j < 16; j++ ) { a.d[j] = test_data[i][j]; }
		
		for ( int j = 0; j < NUM_WARMUP_RUNS; j++ ) {
			(*f) ( &r, &a );
		}
		
		for ( int j = 0; j < NUM_RUNS; j++ ) {
			start = start_tsc();
			(*f) ( &r, &a );
			end = stop_tsc( start );
			num_cycles += (( double ) end );
		}
	}
	num_cycles = num_cycles / NUM_RUNS;
	return num_cycles / num_test_cases;
}

double get_perf_double_op_rad17( int** test_data, int num_test_cases, double_op_rad17_func *f ) {
	myInt64 start, end;
	double num_cycles = 0;
	
	radix17_t a, b, r;
	
	for ( int i = 0; i < num_test_cases/2; i++ ) {
		for ( int j = 0; j < 16; j++ ) { a.d[j] = test_data[i  ][j]; }
		for ( int j = 0; j < 16; j++ ) { b.d[j] = test_data[i+1][j]; }
		
		for ( int j = 0; j < NUM_WARMUP_RUNS; j++ ) {
			(*f) ( &r, &a, &b );
		}
		
		for ( int j = 0; j < NUM_RUNS; j++ ) {
			start = start_tsc();
			(*f) ( &r, &a, &b );
			end = stop_tsc( start );
			num_cycles += (( double ) end );
		}
	}
	num_cycles = num_cycles / NUM_RUNS;
	return num_cycles / ( num_test_cases / 2 );
}


double get_perf_radix17_add( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17( test_data, num_test_cases, &add_256_17 );
}

double get_perf_radix17_add_no_reduce( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17( test_data, num_test_cases, &add_256_17_no_reduce );
}

double get_perf_radix17_sub( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17( test_data, num_test_cases, &sub_256_17 );
}

double get_perf_radix17_sub_no_reduce( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17( test_data, num_test_cases, &sub_256_17_no_reduce );
}

double get_perf_radix17_mul( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17( test_data, num_test_cases, &mul_256_17 );
}


double get_perf_radix17_div( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17( test_data, num_test_cases, &div_256_17 );
}

double get_perf_radix17_sqr( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17( test_data, num_test_cases, &square_256_17 );
}

double get_perf_radix17_inv( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17( test_data, num_test_cases, &inverse_256_17 );
}

// ---------------------------------------------------------------------------------------

typedef void ( single_op_rad51_func )( radix51_t*, radix51_t* );
typedef void ( double_op_rad51_func )( radix51_t*, radix51_t*, radix51_t* );

double get_perf_single_op_rad51( int** test_data, int num_test_cases, single_op_rad51_func *f ) {
	myInt64 start, end;
	double num_cycles = 0;
	
	radix51_t a, r;
	
	for ( int i = 0; i < num_test_cases; i++ ) {
		for ( int j = 0; j < 4; j++ ) { a.d[j] = test_data[i][j]; }
		a.d[4] = 0;
		
		for ( int j = 0; j < NUM_WARMUP_RUNS; j++ ) {
			(*f) ( &r, &a );
		}
		
		for ( int j = 0; j < NUM_RUNS; j++ ) {
			start = start_tsc();
			(*f) ( &r, &a );
			end = stop_tsc( start );
			num_cycles += (( double ) end );
		}
	}
	num_cycles = num_cycles / NUM_RUNS;
	return num_cycles / num_test_cases;
}

double get_perf_double_op_rad51( int** test_data, int num_test_cases, double_op_rad51_func *f ) {
	myInt64 start, end;
	double num_cycles = 0;
	
	radix51_t a, b, r;
	
	for ( int i = 0; i < num_test_cases/2; i++ ) {
		for ( int j = 0; j < 4; j++ ) { a.d[j] = test_data[i  ][j]; }
		for ( int j = 0; j < 4; j++ ) { b.d[j] = test_data[i+1][j]; }
		a.d[4] = 0;
		b.d[4] = 0;
		
		for ( int j = 0; j < NUM_WARMUP_RUNS; j++ ) {
			(*f) ( &r, &a, &b );
		}
		
		for ( int j = 0; j < NUM_RUNS; j++ ) {
			start = start_tsc();
			(*f) ( &r, &a, &b );
			end = stop_tsc( start );
			num_cycles += (( double ) end );
		}
	}
	num_cycles = num_cycles / NUM_RUNS;
	return num_cycles / ( num_test_cases / 2 );
}


double get_perf_radix51_add( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad51( test_data, num_test_cases, &add_256_51 );
}

double get_perf_radix51_sub( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad51( test_data, num_test_cases, &sub_256_51 );
}

double get_perf_radix51_mul( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad51( test_data, num_test_cases, &mul_256_51 );
}

double get_perf_radix51_div( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad51( test_data, num_test_cases, &div_256_51 );
}

double get_perf_radix51_sqr( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad51( test_data, num_test_cases, &square_256_51 );
}

double get_perf_radix51_inv( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad51( test_data, num_test_cases, &inverse_256_51 );
}

// ---------------------------------------------------------------------------------------

typedef void ( single_op_rad17_vec_func )( __m256i*, __m256i*, __m256i, __m256i );
typedef void ( double_op_rad17_vec_func )( __m256i*, __m256i*, __m256i, __m256i, __m256i, __m256i );

double get_perf_single_op_rad17_vec( int** test_data, int num_test_cases, single_op_rad17_vec_func *f ) {
	myInt64 start, end;
	double num_cycles = 0;
	
	__m256i rlo, rhi, alo, ahi;
	
	for ( int i = 0; i < num_test_cases; i++ ) {
		alo = _mm256_setr_epi32( test_data[i][0], test_data[i][1], test_data[i][2 ], test_data[i][3 ], test_data[i][4 ], test_data[i][5 ], test_data[i][6 ], test_data[i][7 ] );
		ahi = _mm256_setr_epi32( test_data[i][8], test_data[i][9], test_data[i][10], test_data[i][11], test_data[i][12], test_data[i][13], test_data[i][14], test_data[i][15] );
				
		for ( int j = 0; j < NUM_WARMUP_RUNS; j++ ) {
			(*f) ( &rlo, &rhi, alo, ahi );
		}
		
		for ( int j = 0; j < NUM_RUNS; j++ ) {
			start = start_tsc();
			(*f) ( &rlo, &rhi, alo, ahi );
			end = stop_tsc( start );
			num_cycles += (( double ) end );
		}
	}
	num_cycles = num_cycles / NUM_RUNS;
	return num_cycles / num_test_cases;
}

double get_perf_double_op_rad17_vec( int** test_data, int num_test_cases, double_op_rad17_vec_func *f ) {
	myInt64 start, end;
	double num_cycles = 0;
		
	__m256i rlo, rhi, alo, ahi, blo, bhi;
	
	for ( int i = 0; i < num_test_cases/2; i++ ) {
		alo = _mm256_setr_epi32( test_data[i][0], test_data[i][1], test_data[i][2 ], test_data[i][3 ], test_data[i][4 ], test_data[i][5 ], test_data[i][6 ], test_data[i][7 ] );
		ahi = _mm256_setr_epi32( test_data[i][8], test_data[i][9], test_data[i][10], test_data[i][11], test_data[i][12], test_data[i][13], test_data[i][14], test_data[i][15] );
		blo = _mm256_setr_epi32( test_data[i+1][0], test_data[i+1][1], test_data[i+1][2 ], test_data[i+1][3 ], test_data[i+1][4 ], test_data[i+1][5 ], test_data[i+1][6 ], test_data[i+1][7 ] );
		bhi = _mm256_setr_epi32( test_data[i+1][8], test_data[i+1][9], test_data[i+1][10], test_data[i+1][11], test_data[i+1][12], test_data[i+1][13], test_data[i+1][14], test_data[i+1][15] );
		
		for ( int j = 0; j < NUM_WARMUP_RUNS; j++ ) {
			(*f) ( &rlo, &rhi, alo, ahi, blo, bhi );
		}
		
		for ( int j = 0; j < NUM_RUNS; j++ ) {
			start = start_tsc();
			(*f) ( &rlo, &rhi, alo, ahi, blo, bhi );
			end = stop_tsc( start );
			num_cycles += (( double ) end );
		}
	}
	num_cycles = num_cycles / NUM_RUNS;
	return num_cycles / ( num_test_cases / 2 );
}


double get_perf_radix17_vec_add( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17_vec( test_data, num_test_cases, &add_256_17_vec );
}

double get_perf_radix17_vec_sub( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17_vec( test_data, num_test_cases, &sub_256_17_vec );
}

double get_perf_radix17_vec_mul( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17_vec( test_data, num_test_cases, &mul_256_17_vec );
}

double get_perf_radix17_vec_div( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17_vec( test_data, num_test_cases, &div_256_17_vec );
}

double get_perf_radix17_vec_sqr( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17_vec( test_data, num_test_cases, &square_256_17_vec );
}

double get_perf_radix17_vec_sqr_alt( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17_vec( test_data, num_test_cases, &square_256_17_vec_alt );
}

double get_perf_radix17_vec_inv( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17_vec( test_data, num_test_cases, &inverse_256_17_vec );
}

// ---------------------------------------------------------------------------------------

double get_perf_radix17_vecn_add( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17( test_data, num_test_cases, &add_256_17_vecn );
}

double get_perf_radix17_vecn_sub( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17( test_data, num_test_cases, &sub_256_17_vecn );
}

double get_perf_radix17_vecn_mul( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17( test_data, num_test_cases, &mul_256_17_vecn );
}


double get_perf_radix17_vecn_div( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17( test_data, num_test_cases, &div_256_17_vecn );
}

double get_perf_radix17_vecn_sqr( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17( test_data, num_test_cases, &square_256_17_vecn );
}

double get_perf_radix17_vecn_inv( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17( test_data, num_test_cases, &inverse_256_17_vecn );
}

// ---------------------------------------------------------------------------------------

double get_perf_radix17_opt_vec_add( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17_vec( test_data, num_test_cases, &add_256_17_all );
}

double get_perf_radix17_opt_vec_sub( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17_vec( test_data, num_test_cases, &sub_256_17_all );
}

double get_perf_radix17_opt_vec_mul( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17_vec( test_data, num_test_cases, &mul_256_17_all );
}

double get_perf_radix17_opt_vec_div( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17_vec( test_data, num_test_cases, &div_256_17_all );
}

double get_perf_radix17_opt_vec_div_nvec( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17_vec( test_data, num_test_cases, &div_256_non_vec );
}

double get_perf_radix17_opt_vec_div_alt( int** test_data, int num_test_cases ) {
	return get_perf_double_op_rad17_vec( test_data, num_test_cases, &div_256_non_vec_alt );
}

double get_perf_radix17_opt_vec_sqr( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17_vec( test_data, num_test_cases, &square_256_17_all );
}

double get_perf_radix17_opt_vec_sqr_new( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17_vec( test_data, num_test_cases, &square_256_17_all_new );
}

double get_perf_radix17_opt_vec_sqr_newvec( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17_vec( test_data, num_test_cases, &square_256_opt_vec );
}

double get_perf_radix17_opt_vec_inv( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17_vec( test_data, num_test_cases, &inverse_256_17_all );
}

double get_perf_radix17_opt_vec_inv_nvec( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17_vec( test_data, num_test_cases, &inverse_256_non_vec );
}

double get_perf_radix17_opt_vec_inv_alt( int** test_data, int num_test_cases ) {
	return get_perf_single_op_rad17_vec( test_data, num_test_cases, &inverse_256_non_vec_alt );
}

// ---------------------------------------------------------------------------------------
//  Main microbenchmarking function
// ---------------------------------------------------------------------------------------

int run_microbenchmarking( int be_verbose, char *res_dir_name ) {
    FILE *fp;
    char res_fName[200];

	// -- initialize any required parameters
	init_vec_constants();
	
	// -- initialize test data
	init_rand();
	int num_test_cases = 100;
	int** test_data = (int**)malloc(num_test_cases * sizeof(int*)); for ( int i = 0; i < num_test_cases; ++i ) { test_data[i] = (int*)malloc( (16) * sizeof(int) ); }
	init_data( test_data, num_test_cases );

	// -- run microbenchmarking
	double perf;
	//  - 1. radix_17
	//   -- prepare result file
	strcpy( res_fName, res_dir_name );
	strcat( res_fName, "256_17.csv" );
	fp = fopen( res_fName, "w" );
	if ( fp == NULL ) {
		printf( " -- ERROR: unable to find output file '%s'\n", res_fName );
		return 0;
	}
	//   -- do benchmarking
	if ( be_verbose ) { print_mb_start( "Radix17" ); }
	perf = get_perf_radix17_add          ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "add"           ); print_perf( "add                       ", perf );
	perf = get_perf_radix17_add_no_reduce( test_data, num_test_cases ); write_perf_to_file( fp, perf, "add_no_reduce" ); print_perf( "add no reduction          ", perf );
	perf = get_perf_radix17_sub          ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sub"           ); print_perf( "sub                       ", perf );
	perf = get_perf_radix17_sub_no_reduce( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sub_no_reduce" ); print_perf( "sub no reduction          ", perf );
	perf = get_perf_radix17_mul          ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "mul"           ); print_perf( "mul                       ", perf );
	perf = get_perf_radix17_div          ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "div"           ); print_perf( "div                       ", perf );
	perf = get_perf_radix17_sqr          ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sqr"           ); print_perf( "sqr                       ", perf );
	perf = get_perf_radix17_inv          ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "inv"           ); print_perf( "inv                       ", perf );
	//   -- cleanup
	fclose( fp );
	
	//  - 2. radix_51
	//   -- prepare result file
	strcpy( res_fName, res_dir_name );
	strcat( res_fName, "256_51.csv" );
	fp = fopen( res_fName, "w" );
	if ( fp == NULL ) {
		printf( " -- ERROR: unable to find output file '%s'\n", res_fName );
		return 0;
	}
	//   -- do benchmarking
	if ( be_verbose ) { print_mb_start( "Radix51" ); }
	perf = get_perf_radix51_add( test_data, num_test_cases ); write_perf_to_file( fp, perf, "add" ); print_perf( "add                       ", perf );
	perf = get_perf_radix51_sub( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sub" ); print_perf( "sub                       ", perf );
	perf = get_perf_radix51_mul( test_data, num_test_cases ); write_perf_to_file( fp, perf, "mul" ); print_perf( "mul                       ", perf );
	perf = get_perf_radix51_div( test_data, num_test_cases ); write_perf_to_file( fp, perf, "div" ); print_perf( "div                       ", perf );
	perf = get_perf_radix51_sqr( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sqr" ); print_perf( "sqr                       ", perf );
	perf = get_perf_radix51_inv( test_data, num_test_cases ); write_perf_to_file( fp, perf, "inv" ); print_perf( "inv                       ", perf );
	//   -- cleanup
	fclose( fp );

	//  - 3. radix_17_vec
	//   -- prepare result file
	strcpy( res_fName, res_dir_name );
	strcat( res_fName, "256_17_vec.csv" );
	fp = fopen( res_fName, "w" );
	if ( fp == NULL ) {
		printf( " -- ERROR: unable to find output file '%s'\n", res_fName );
		return 0;
	}
	//   -- do benchmarking
	if ( be_verbose ) { print_mb_start( "Radix17_vec" ); }
	perf = get_perf_radix17_vec_add    ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "add    " ); print_perf( "add                       ", perf );
	perf = get_perf_radix17_vec_sub    ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sub    " ); print_perf( "sub                       ", perf );
	perf = get_perf_radix17_vec_mul    ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "mul    " ); print_perf( "mul                       ", perf );
	perf = get_perf_radix17_vec_div    ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "div    " ); print_perf( "div                       ", perf );
	perf = get_perf_radix17_vec_sqr    ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sqr    " ); print_perf( "sqr                       ", perf );
	perf = get_perf_radix17_vec_sqr_alt( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sqr alt" ); print_perf( "sqr write to arr + non-vec", perf );
	perf = get_perf_radix17_vec_inv    ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "inv    " ); print_perf( "inv                       ", perf );
	//   -- cleanup
	fclose( fp );

	/*
	//  - 4. radix_17_vec_naive
	//   -- prepare result file
	strcpy( res_fName, res_dir_name );
	strcat( res_fName, "256_17_vecn.csv" );
	fp = fopen( res_fName, "w" );
	if ( fp == NULL ) {
		printf( " -- ERROR: unable to find output file '%s'\n", res_fName );
		return 0;
	}
	//   -- do benchmarking
	if ( be_verbose ) { print_mb_start( "Radix17_vec naive" ); }
	perf = get_perf_radix17_vecn_add( test_data, num_test_cases ); write_perf_to_file( fp, perf, "add" ); print_perf( "add                       ", perf );
	perf = get_perf_radix17_vecn_sub( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sub" ); print_perf( "sub                       ", perf );
	perf = get_perf_radix17_vecn_mul( test_data, num_test_cases ); write_perf_to_file( fp, perf, "mul" ); print_perf( "mul                       ", perf );
	perf = get_perf_radix17_vecn_div( test_data, num_test_cases ); write_perf_to_file( fp, perf, "div" ); print_perf( "div                       ", perf );
	perf = get_perf_radix17_vecn_sqr( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sqr" ); print_perf( "sqr                       ", perf );
	perf = get_perf_radix17_vecn_inv( test_data, num_test_cases ); write_perf_to_file( fp, perf, "inv" ); print_perf( "inv                       ", perf );
	//   -- cleanup
	fclose( fp );
	*/

	//  - 5. opt_all_17
	//   -- prepare result file
	strcpy( res_fName, res_dir_name );
	strcat( res_fName, "opt_all.csv" );
	fp = fopen( res_fName, "w" );
	if ( fp == NULL ) {
		printf( " -- ERROR: unable to find output file '%s'\n", res_fName );
		return 0;
	}
	//   -- do benchmarking
	if ( be_verbose ) { print_mb_start( "opt all" ); }
	perf = get_perf_radix17_opt_vec_add       ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "add"         ); print_perf( "add no reduction          ", perf );
	perf = get_perf_radix17_opt_vec_sub       ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sub"         ); print_perf( "sub conditional reduction ", perf );
	perf = get_perf_radix17_opt_vec_mul       ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "mul"         ); print_perf( "mul                       ", perf );
	perf = get_perf_radix17_opt_vec_div       ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "div"         ); print_perf( "div                       ", perf );
	perf = get_perf_radix17_opt_vec_div_nvec  ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "div non-vec" ); print_perf( "div occasional sq. non-vec", perf );
	perf = get_perf_radix17_opt_vec_div_alt   ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "div new"     ); print_perf( "div always sq. non-vec    ", perf ); // uses the square new function everywhere
	perf = get_perf_radix17_opt_vec_sqr       ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sqr"         ); print_perf( "sqr                       ", perf );
	perf = get_perf_radix17_opt_vec_sqr_new   ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sqr new"     ); print_perf( "sqr write to arr + non-vec", perf );
	perf = get_perf_radix17_opt_vec_sqr_newvec( test_data, num_test_cases ); write_perf_to_file( fp, perf, "sqr_new-vec" ); print_perf( "sqr new vec               ", perf );
	perf = get_perf_radix17_opt_vec_inv       ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "inv"         ); print_perf( "inv                       ", perf );
	perf = get_perf_radix17_opt_vec_inv_nvec  ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "inv non-vec" ); print_perf( "inv occasional sq. non-vec", perf );
	perf = get_perf_radix17_opt_vec_inv_alt   ( test_data, num_test_cases ); write_perf_to_file( fp, perf, "inv new"     ); print_perf( "inv always sq. non-vec    ", perf ); // // uses the square new function everywhere
	//   -- cleanup
	fclose( fp );
	return 1;
}


