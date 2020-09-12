#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "test_set_reader.h"

#define MAX_PUB_KEY_LENGTH 256

// read input values
int get_str_input_value( FILE *fp, char *value, int len, int with_fixed_length ) {
	char c;
	if ( with_fixed_length ) {
		for ( int i = 0; i < len; ++i ) {
			c = getc( fp );
			if ( ( c == '\n' ) || ( c == EOF ) ) {
				return 0;
			}
			value[i] = c;
		}
		// read the trailing new line character to prepare for the retrieval of the next value
		c = getc( fp );
		if ( ( c != '\n' ) && ( c != EOF ) ) {
			return 0;
		}
		// add the string terminator
		value[len] = '\0';
	}
	else {
		int i = 0;
		int cnt = 0;
		int term_found = 0;
		for ( int i = 0; i < len; ++i ) {
			c = getc( fp );
			if ( ( c == '\n' ) || ( c == EOF ) ) {
				term_found = 1;
				break;
			}
			value[i] = c;
			cnt++;
		}
		if ( !term_found ) {
			c = getc( fp );
			// if the maximum length "len" has been read (i.e. no terminating character has been found), but the string has still not terminated, return "false"
			if ( ( c != '\n' ) && ( c != EOF ) ) {
				return 0;
			}
		}
		value[cnt] = '\0';
	}
	// return "true"
	return 1;
}

int read_test_set( char *rel_path_to_benchmarking_dir, int test_set_no, int do_read_pks, int *num_test_cases, int *secret_key_length ) {

    // convert test_number to string
    char test_no_str[5];
    snprintf( test_no_str, 5, "%d", test_set_no );
    
    // - determine test input file directory
	#ifdef _WIN32
		char path_separator[3] = "\\";
	#else
		char path_separator[2] = "/";
	#endif
    //char test_input_dir_name[100] = "../benchmark_framework/benchmark_in_files";
    char test_input_dir_name[100] = "";
    strcat( test_input_dir_name, rel_path_to_benchmarking_dir );
    strcat( test_input_dir_name, path_separator               );
    strcat( test_input_dir_name, "benchmark_in_files"         );
    strcat( test_input_dir_name, path_separator               );
    strcat( test_input_dir_name, "test_set_"                  );
    strcat( test_input_dir_name, test_no_str                  );
    strcat( test_input_dir_name, path_separator               );
        
    // - read values from input files
    // -- retrieve general parameters
    //    (currently only retrieves number of test cases and secret key length, could expand this in the future to include more parameters specific to the chosen test case)
    FILE *fp;
    char fName[200];
    
    strcpy( fName, test_input_dir_name );
    strcat( fName, "general_parameters.txt" );
    fp = fopen( fName, "r" );
    if ( fp == NULL ) {
    	printf( " -- ERROR: unable to find input file 'general_parameters.txt'\n" );
    	return 0;
    }
    if(fscanf( fp, "%d", num_test_cases)!=1){
		printf("scanf had a problem!\n");
	}
	if(fscanf( fp, "%d", secret_key_length)!=1){
		printf("scanf had a problem!\n");
	}
    fclose( fp );
    
    /*
    printf( "secret key length = %d\n", *sk_length );
    printf( "num test cases = %d\n", *num_test_cases );
    */
    
    int input_size = *num_test_cases;
    int sk_length  = *secret_key_length;
    // -- read secret keys
    secrets_A = (char**)malloc(input_size * sizeof(char*)); for ( int i = 0; i < input_size; ++i ) { secrets_A[i] = (char*)malloc( (sk_length+1) * sizeof(char) ); }
    strcpy( fName, test_input_dir_name );
    strcat( fName, "secret_keys.txt" );
    fp = fopen( fName, "r" );
    if ( fp == NULL ) {
    	printf( " -- ERROR: unable to find input file 'secret_keys.txt'\n" );
    	return 0;
    }
    for ( int i = 0; i < input_size ; ++i ) {
    	if ( !get_str_input_value( fp, secrets_A[i], sk_length, 1 ) ) { printf( " -- ERROR: Unable to read secret key from file '%s' (either incorrect or too few secretkeys provided; expected %d secret keys of length %d)\n", fName, input_size, sk_length ); fclose( fp ); return 0; }
    }
    fclose( fp );
	
    // -- read secret keys of B (for validation purposes)
    secrets_B = (char**)malloc(input_size*sizeof(char*)); for ( int i = 0; i < input_size; ++i ) { secrets_B[i] = (char*)malloc( (sk_length+1) * sizeof(char) ); }
    strcpy( fName, test_input_dir_name );
    strcat( fName, "secret_keys_B.txt" );
    fp = fopen( fName, "r" );
    if ( fp == NULL ) {
    	printf( " -- ERROR: unable to find input file 'secret_keys_B.txt'\n" );
    	return 0;
    }
    for ( int i = 0; i < input_size ; ++i ) {
    	if ( !get_str_input_value( fp, secrets_B[i], sk_length, 1 ) ) { printf( " -- ERROR: Unable to read secret key from file '%s' (either incorrect or too few secretkeys provided; expected %d secret keys of length %d)\n", fName, input_size, sk_length ); fclose( fp ); return 0; }
    }
    fclose( fp );
	
	if ( do_read_pks ) {
		// -- read (decimal) public keys (i.e. the expected results)
    	public_keys_dec = (char**)malloc(input_size*sizeof(char*)); for ( int i = 0; i < input_size; ++i ) { public_keys_dec[i] = (char*)malloc( (MAX_PUB_KEY_LENGTH+1) * sizeof(char) ); }
    	strcpy( fName, test_input_dir_name );
    	strcat( fName, "public_keys.txt" );
    	fp = fopen( fName, "r" );
    	if ( fp == NULL ) {
    		printf( " -- ERROR: unable to find input file 'public_keys.txt'\n" );
    		return 0;
    	}
    	for ( int i = 0; i < input_size; ++i ) {
    		if ( !get_str_input_value( fp, public_keys_dec[i], MAX_PUB_KEY_LENGTH, 0 ) ) { printf( " -- ERROR: Unable to read public key from file '%s' (either incorrect or too few public keys provided; expected %d public keys with maximum length %d)\n", fName, input_size, MAX_PUB_KEY_LENGTH ); fclose( fp ); return 0; } 
    	}
    	fclose( fp );
    	
		// -- read (radix51) public keys (i.e. the expected results)
    	strcpy( fName, test_input_dir_name );
    	strcat( fName, "public_keys_radix51.txt" );
    	fp = fopen( fName, "r" );
    	if ( fp == NULL ) {
    		printf( " -- WARNING: unable to find input file 'public_keys_radix51.txt'\n" );
    		public_keys_r51 = NULL;
    	}
    	else {
		    public_keys_r51 = (char**)malloc(input_size*sizeof(char*)); for ( int i = 0; i < input_size; ++i ) { public_keys_r51[i] = (char*)malloc( (MAX_PUB_KEY_LENGTH+1) * sizeof(char) ); }
    		for ( int i = 0; i < input_size; ++i ) {
    			if ( !get_str_input_value( fp, public_keys_r51[i], MAX_PUB_KEY_LENGTH, 0 ) ) { printf( " -- ERROR: Unable to read public key from file '%s' (either incorrect or too few public keys provided; expected %d public keys with maximum length %d)\n", fName, input_size, MAX_PUB_KEY_LENGTH ); fclose( fp ); return 0; } 
    		}
    		fclose( fp );
		}
	
		// -- read (radix17) public keys (i.e. the expected results)
		strcpy( fName, test_input_dir_name );
		strcat( fName, "public_keys_radix17.txt" );
		fp = fopen( fName, "r" );
		if ( fp == NULL ) {
			printf( " -- WARNING: unable to find input file 'public_keys_radix17.txt'\n" );
			public_keys_r17 = NULL;
		}
		else {
		    public_keys_r17 = (char**)malloc(input_size*sizeof(char*)); for ( int i = 0; i < input_size; ++i ) { public_keys_r17[i] = (char*)malloc( (MAX_PUB_KEY_LENGTH+1) * sizeof(char) ); }
			for ( int i = 0; i < input_size; ++i ) {
				if ( !get_str_input_value( fp, public_keys_r17[i], MAX_PUB_KEY_LENGTH, 0 ) ) { printf( " -- ERROR: Unable to read public key from file '%s' (either incorrect or too few public keys provided; expected %d public keys with maximum length %d)\n", fName, input_size, MAX_PUB_KEY_LENGTH ); fclose( fp ); return 0; } 
			}
			fclose( fp );
		}
		
		// -- read (radix25) public keys (i.e. the expected results)
    	strcpy( fName, test_input_dir_name );
    	strcat( fName, "public_keys_radix25.txt" );
    	fp = fopen( fName, "r" );
    	if ( fp == NULL ) {
    		printf( " -- WARNING: unable to find input file 'public_keys_radix25.txt'\n" );
    		public_keys_r25 = NULL;
    	}
    	else {
		    public_keys_r25 = (char**)malloc(input_size*sizeof(char*)); for ( int i = 0; i < input_size; ++i ) { public_keys_r25[i] = (char*)malloc( (MAX_PUB_KEY_LENGTH+1) * sizeof(char) ); }
    		for ( int i = 0; i < input_size; ++i ) {
    			if ( !get_str_input_value( fp, public_keys_r25[i], MAX_PUB_KEY_LENGTH, 0 ) ) { printf( " -- ERROR: Unable to read public key from file '%s' (either incorrect or too few public keys provided; expected %d public keys with maximum length %d)\n", fName, input_size, MAX_PUB_KEY_LENGTH ); fclose( fp ); return 0; } 
    		}
    		fclose( fp );
		}
	}
	else {
		public_keys_dec = NULL;
		public_keys_r51 = NULL;
		public_keys_r17 = NULL;
		public_keys_r25 = NULL;
	}
	
    return 1;
}
