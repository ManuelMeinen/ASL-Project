#include <stdio.h>
#include <string.h>
#include "radix51.h"

void print_radix51(radix51_t *x, char *string){
    printf("%s : %lld, %lld, %lld, %lld, %lld\n", string, x->d[0], x->d[1], x->d[2], x->d[3], x->d[4]);
}

void radix51_toString( char* dest, radix51_t x ) {
    char tmp[100];
	
	strcpy( dest, "[" );
	for ( int i = 0; i < 5; ++i ) {
		snprintf( tmp, 100, "%lld", x.d[i] );
		strcat  ( dest, tmp );
		if ( i < 4 ) {
			strcat( dest, ", " );
		}
	}
    strcat( dest, "]" );
}
