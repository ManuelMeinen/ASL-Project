#include <stdio.h>
#include <string.h>
#include "radix51_arr.h"

void print_radix51_arr( long x[5], char *string ) {
    printf( "%s : %ld, %ld, %ld, %ld, %ld\n", string, x[0], x[1], x[2], x[3], x[4] );
}

void radix51_arr_toString( char* dest, long x[5] ) {
    char tmp[100];
	
	strcpy( dest, "[" );
	for ( int i = 0; i < 5; ++i ) {
		snprintf( tmp, 100, "%ld", x[i] );
		strcat  ( dest, tmp );
		if ( i < 4 ) {
			strcat( dest, ", " );
		}
	}
    strcat( dest, "]" );
}
