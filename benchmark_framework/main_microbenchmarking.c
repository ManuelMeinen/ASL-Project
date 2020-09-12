#include "../_common_libs/common.h"
#include "microbenchmarking.h"
#include <string.h>
#include <stdlib.h>


int main( int argc, char **argv ) {
    // - retrieve input parameter
    if ( argc != 3 ) {
		printf( " -- ERROR: Wrong number of parameters specified. Usage: <be_verbose> <res_dir_name>\n" );
		printf( "      Parameter 'be_verbose' may have the following values:\n" );
		printf( "         1 : be verbose (i.e. prints information on all steps)\n" );
		printf( "         0 : silent mode (i.e. prints only the performance or an error if one occurred)\n" );
		printf( "\n" );
    	return -1;
    }
    int  be_verbose    = atoi( argv[1] );
    char *res_dir_name = argv[2];
    
	// - run benchmarking
	run_microbenchmarking( be_verbose, res_dir_name );
    return 0;
}

