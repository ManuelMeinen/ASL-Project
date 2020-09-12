#include <stdlib.h>
#include <stdio.h>

int main( int argc, char **argv ) {
	printf( "integer size         : %lu bytes\n", sizeof( int ) );
	printf( "long size            : %lu bytes\n", sizeof( long ) );
	printf( "long long size       : %lu bytes\n", sizeof( long long ) );
	printf( "signed integer size  : %lu bytes\n", sizeof( signed int ) );
	printf( "signed long size     : %lu bytes\n", sizeof( signed long ) );
	printf( "signed long long size: %lu bytes\n", sizeof( signed long long ) );
	return 0;
}
