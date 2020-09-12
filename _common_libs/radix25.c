#include <stdio.h>
#include <string.h>
#include "radix25.h"

void print_radix25(radix25_t *x, char *string) {
    printf("%s: \n", string);
    for(int i =0;i < 10; i++){
        printf("%d, ",x->d[i]);
    }
    printf("\n");
}

void radix25_toString( char* dest, radix25_t x ) {
    char tmp[100];
	
	strcpy( dest, "[" );
	for ( int i = 0; i < 10; ++i ) {
		snprintf( tmp, 100, "%d", x.d[i] );
		strcat  ( dest, tmp );
		if ( i < 9 ) {
			strcat( dest, ", " );
		}
	}
    strcat( dest, "]" );
}
