#include "common.h"

int sk_length;
int sk_str_length;

void set_sk_length( int len ) {
    sk_length = len;
    sk_str_length = len+1;
}

void printUsageHelp() {
	printf( " -- ERROR: Wrong number of parameters specified. Usage: <test_number> <run_mode> <be_verbose> <res_dir_name>\n" );
	printf( "      Parameter 'run_mode' may have the following values:\n" );
	printf( "         v : validation-only\n" );
	printf( "         b : benchmarking-only\n" );
	printf( "         f : full validation and benchmarking\n" );
	printf( "      Parameter 'be_verbose' may have the following values:\n" );
	printf( "         1 : be verbose (i.e. prints information on all steps)\n" );
	printf( "         0 : silent mode (i.e. prints only the performance or an error if one occurred)\n" );
	printf( "\n" );
}
