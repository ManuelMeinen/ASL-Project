#include <stdio.h>
#include <string.h>
#include "radix17_vec.h"

void print_radix17_vec(__m256i xlo, __m256i xhi, char *string){
	int* i_v = (int*) &xlo;
    printf("%s:\n%d, %d, %d, %d, %d, %d, %d, %d,",string, i_v[0], i_v[1], i_v[2], i_v[3], i_v[4], i_v[5], i_v[6], i_v[7]);
	i_v = (int*) &xhi;
    printf(" %d, %d, %d, %d, %d, %d, %d, %d,\n", i_v[0], i_v[1], i_v[2], i_v[3], i_v[4], i_v[5], i_v[6], i_v[7]);
}


void radix17_vec_toString( char* dest, __m256i xlo, __m256i xhi ) {
    char tmp[100];
	
	strcpy( dest, "[" );
	
    int* i_v0 = (int*) &xlo;
    int* i_v1 = (int*) &xhi;

	for ( int i = 0; i < 8; ++i ) {
		snprintf( tmp, 100, "%d", i_v0[i] );
		strcat  ( dest, tmp );
		strcat  ( dest, ", " );
	}	
	for ( int i = 0; i < 7; ++i ) {
		snprintf( tmp, 100, "%d", i_v1[i] );
		strcat  ( dest, tmp );
		if ( i < 6 ) {
			strcat( dest, ", " );
		}
	}	
    strcat( dest, "]" );
}