#include <stdio.h>
#include <string.h>
#include "radix17.h"

void print_radix17(radix17_t *x, char *string) {
    printf("%s: \n", string);
    for(int i =0;i < 16; i++){
        printf("%d, ",x->d[i]);
    }
    printf("\n");
}

void radix17_toString( char* dest, radix17_t x ) {
    char tmp[100];
	
	strcpy( dest, "[" );
	for ( int i = 0; i < 15; ++i ) {
		snprintf( tmp, 100, "%d", x.d[i] );
		strcat  ( dest, tmp );
		if ( i < 14 ) {
			strcat( dest, ", " );
		}
	}
    strcat( dest, "]" );
}

