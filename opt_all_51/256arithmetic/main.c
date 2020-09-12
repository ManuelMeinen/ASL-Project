#include "bignum.h"

int main(int argc, char const *argv[])
{
	long a[5] = {575160602422936, 1618811318005166, 97863898957929, 1742854989071013, 517724694265016};
    long b[5] = {1704666820865218, 1958688903789407, 999601117524594, 945849795254012, 435221983174332};
    long r[5] = {0, 0, 0, 0, 0};

	long res_add[5] = {28027609602906  , 1325700408109326, 1097465016482524, 436904970639777 , 952946677439349 };
	long res_sub[5] = {1122293595242966, 1911922227901006, 1350062595118582, 797005193817000 , 82502711090684  };
	long res_mul[5] = {1765262878194772, 2160670913896529, 1059500224136984, 2120783277102739, 2078530830153566};
	long res_div[5] = {664433579777367 , 798145316342188 , 1448017419764148, 1540225854329050, 911028651365221 };
	long res_sqr[5] = {1391896342111381, 1351463214626422, 1089264114725229, 2249629422887814, 919810247932394 };
	long res_inv[5] = {793321103500114 , 1755768637347786, 1467874338892072, 104638582436266 , 163445459748874 };

	add_256( r, a, b ); reduce_2_to_2( r ); // reduce done separately as no longer done in add_256 for optimization purposes (only yields correct result if used as an input to a mult)
	if ( cmp_coord_256( r, res_add ) ) { printf( "Add   : OK\n" ); } else { printf( " -- Add   : failed\n" ); }

	sub_256( r, a, b );
	if ( cmp_coord_256( r, res_sub ) ) { printf( "Sub   : OK\n" ); } else { printf( " -- Sub   : failed\n" ); }

	mul_256( r, a, b );
	if ( cmp_coord_256( r, res_mul ) ) { printf( "Mul   : OK\n" ); } else { printf( " -- Mul   : failed\n" ); }

	div_256( r, a, b );
	if ( cmp_coord_256( r, res_div ) ) { printf( "Div   : OK\n" ); } else { printf( " -- Div   : failed\n" ); }

	square_256( r, a );
	if ( cmp_coord_256( r, res_sqr ) ) { printf( "Square: OK\n" ); } else { printf( " -- Square: failed\n" ); }

	inverse_256( r, a );
	if ( cmp_coord_256( r, res_inv ) ) { printf( "Inv   : OK\n" ); } else { printf( " -- Inv   : failed\n" ); }

    return 0;
}
