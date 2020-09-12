#include "curve25519_base_256_17.h"
#include "../256arithmetic/bignum.h"
#include "../../_common_libs/radix17.h"
#include "../../_common_libs/radix17_curveParams.h"
#include "../../_common_libs/common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>

void point_double_base_256_17(radix17_t *resX, radix17_t *resY, radix17_t *x, radix17_t *y)
{
    /*
    res = 2*(x,y)
    */
    if(iszero_256(x) && iszero_256(y)){  
        copy_256(resX, *x);
        copy_256(resY, *y);
        return;
    }

    radix17_t upperfrac, lowerfrac, tmp1, tmp2, tmp3, l;
    set_ui_256(&upperfrac, 0);
    set_ui_256(&lowerfrac, 0);
    set_ui_256(&tmp1, 0);
    set_ui_256(&tmp2, 0);
    set_ui_256(&tmp3, 0);
    set_ui_256(&l, 0);
    

    //x - coordinate
    mul_256(&tmp1, &A_rad17, x);               // tmp1 = A*x
    add_256(&tmp1, &tmp1, &tmp1);          // tmp1 = 2*tmp1
    add_256(&tmp1, &tmp1, &const_1_rad17);       // tmp1 = tmp1 + 1 

    square_256(&tmp2, x);                // tmp2 = x^2
    mul_256(&tmp2, &tmp2, &const_3_rad17);       // tmp2 = tmp2*3

    add_256(&upperfrac, &tmp2, &tmp1);     // upperfrac = tmp2+tmp1

    mul_256(&lowerfrac, &B_rad17, y);          // lowerfrac = B*y
    add_256(&lowerfrac, &lowerfrac, &lowerfrac); // lowwerfrac = 2*lowerfrac

    div_256(&l, &upperfrac, &lowerfrac);   // l = uperfrac/lowerfrac

    square_256(&l, &l);                    // l=l^2

    mul_256(resX, &B_rad17, &l);           // resX = B*l
    sub_256(resX, resX, &A_rad17);    // resX = resX - A
    sub_256(resX, resX, x);     // resX = resX - x
    sub_256(resX, resX, x);     // resX = resX - x

    //y - coordinate
    // l, upperfrac, lowerfrac, tmp1, tmp2 don't need to be reinitialized
    // l = 3x^2 + 2 Ax + 1 / 2By   
    mul_256(&tmp1, &A_rad17, x);               // tmp1 = A*x
    add_256(&tmp1, &tmp1, &tmp1);          // tmp1 = 2*tmp1
    add_256(&tmp1, &tmp1, &const_1_rad17);       // tmp1 = tmp1 + 1

    square_256(&tmp2, x);                // tmp2 = x^2
    mul_256(&tmp2, &tmp2, &const_3_rad17);       // tmp2 = tmp2*3

    add_256(&upperfrac, &tmp2, &tmp1);     // upperfrac = tmp2/tmp1

    mul_256(&lowerfrac, &B_rad17, y);          // lowerfrac = B*y
    add_256(&lowerfrac, &lowerfrac, &lowerfrac);   // lowerfrac = 2*lowerfrac

    div_256(&l, &upperfrac, &lowerfrac);   // l = upperfrac/lowerfrac

    mul_256(&tmp1, x, &const_3_rad17);           // tmp1 = 3*x
    add_256(&tmp1, &tmp1, &A_rad17);                // tmp1 = tmp1+A

    mul_256(&tmp1, &tmp1, &l);             // tmp1 = tmp1*l

    square_256(&tmp3,&l);                 // tmp3=l^2
    mul_256(&l, &tmp3, &l);                // l = l^3

    mul_256(&l, &l, &B_rad17);                  // l = l*B

    sub_256(resY, &tmp1, &l);         // resY = tmp1 - l

    sub_256(resY, resY, y);     // resY = resY - y
}

void point_add_base_256_17(radix17_t *resX, radix17_t *resY, radix17_t *op1x, radix17_t *op1y, radix17_t *op2x, radix17_t *op2y)
{    
    // resX <- op1x + op2x
    // resY <- op1y + op2y
    // https://en.wikipedia.org/wiki/Montgomery_curve#Addition
    
    // return (0,0) if ( op1x == op2x AND op1y != op2y )
    if ( cmp_coord_256( op1x, op2x ) && !cmp_coord_256(op1y, op2y)) {
    	set_ui_256(resX, 0 );
    	set_ui_256(resY, 0 );
    	return;
    }
    
    // return op2 if the point op1 is (0,0)
    if ( iszero_256( op1x ) && iszero_256( op1y ) ) {  
    	// probably not side channel safe anymore ?
    	copy_256( resX, *op2x );
    	copy_256( resY, *op2y );
    	return;
    }
    
    // return op1 if the point op2 is (0,0)
    if ( iszero_256( op2x ) && iszero_256( op2y ) ) {  
    	// probably not side channel safe anymore ?
    	copy_256( resX, *op1x );
    	copy_256( resY, *op1y );
    	return;
    }

    if ( cmp_point_256( op1x, op1y, op2x, op2y ) ) { 
    	// point_add_base_256_17(p1, p1) = double(p1), otherwise exception div by 0
    	point_double_base_256_17( resX, resY, op1x, op1y );
    	return;
    }

    radix17_t diffx, diffy, tmp1, tmp2;
    set_ui_256( &diffx   , 0 ); 
    set_ui_256( &diffy   , 0 );
    set_ui_256( &tmp1    , 0 );
    set_ui_256( &tmp2    , 0 );
    

    // x - coordinate
    sub_256( &diffy, op2y, op1y );			// diffy = y2 - y1
    sub_256( &diffx, op2x, op1x );			// diffx = x2 - x1
    
    square_256( &diffy, &diffy );				// diffy = diffy^2
    square_256( &diffx, &diffx );				// diffx = diffx^2
    
    mul_256( &diffy, &diffy, &B_rad17 );			// diffy = B * diffy
    div_256( resX, &diffy, &diffx );		// resX = diffy / diffx
    
    sub_256( resX, resX, &A_rad17 );		// resX -= A
    sub_256( resX, resX, op1x );	// resX -= x1
    sub_256( resX, resX, op2x );	// resX -= x2

    // y - coordinate
	// Note: no need to re-initialize the local variables diffx and diffy as they are always overwritten
	
	add_256( &tmp1, op1x, op2x );			// tmp1 = x1 + x2
	add_256( &tmp1, &tmp1, op1x );			// tmp1 = tmp1 + x1
	add_256( &tmp1, &tmp1, &A_rad17   );			// tmp1 = tmp1 + A, i.e. tmp1 = 2*x1 + x2 + A
	
	sub_256( &diffy, op2y, op1y );			// diffy = y2 - y1
	mul_256( &tmp1, &tmp1, &diffy );			// tmp1 = tmp1 * diffy
	
	sub_256( &diffx, op2x, op1x );			// diffx = x2 - x1
	div_256( &tmp1, &tmp1, &diffx );			// tmp1 = tmp1 / diffx
	
	square_256( &tmp2, &diffx );
	mul_256   ( &diffx, &tmp2, &diffx );		// diffx = diffx ^ 3
	square_256( &tmp2, &diffy );
	mul_256   ( &diffy, &tmp2, &diffy );		// diffy = diffy ^ 3
	mul_256   ( &diffy, &diffy, &B_rad17 );			// diffy = diffy * B
	
	div_256( &tmp2, &diffy, &diffx );		// tmp2 = diffy / diffx
	sub_256( resY, &tmp1, &tmp2 );		// resY = tmp1 - tmp2
	sub_256( resY, resY, op1y );	    // resY -= y1
}

void scalar_mult_base_256_17(radix17_t *base_point_x, radix17_t *base_point_y, char *sk, radix17_t *Pkx, radix17_t *Pky)
{
    /*
    Pk = sk*base_point
    */
    int n = sk_length;
    radix17_t R0x, R0y, R1x, R1y, resX, resY;
    set_ui_256(&R0x, 0);
    set_ui_256(&R0y, 0);


    copy_256(&R1x, *base_point_x);
    copy_256(&R1y, *base_point_y);

    for(int i=0; i<n; i++){

        set_ui_256(&resX, 0);
        set_ui_256(&resY, 0);


        if(sk[i]=='0'){
            // R1 <- add (R0, R1)
            point_add_base_256_17(&resX, &resY, &R0x, &R0y, &R1x, &R1y);
            copy_256(&R1x, resX);
            copy_256(&R1y, resY);

            set_ui_256(&resX, 0);
            set_ui_256(&resY, 0);
            // R0 <- double(R0)
            point_double_base_256_17(&resX, &resY, &R0x, &R0y);
            copy_256(&R0x, resX);
            copy_256(&R0y, resY);
        }else{
            // R0 <- add(R0, R1)
            point_add_base_256_17(&resX, &resY, &R0x, &R0y, &R1x, &R1y);
            copy_256(&R0x, resX);
            copy_256(&R0y, resY);

            set_ui_256(&resX, 0);
            set_ui_256(&resY, 0);
            // R1 <- double(R1)
            point_double_base_256_17(&resX, &resY, &R1x, &R1y);
            copy_256(&R1x, resX);
            copy_256(&R1y, resY);
        }
    }
    // Assign Pk
    
    copy_256(Pkx, R0x);
    copy_256(Pky, R0y);
}

void public_key_gen_base_256_17(char *sk, radix17_t *x, radix17_t *y)
{
    /*
    (x,y) = sk*G
    */
    scalar_mult_base_256_17(&Gx_rad17, &Gy_rad17, sk, x, y);
    /*
    print_radix17(x, "Pk(x)");
    print_radix17(y, "Pk(y)");
    */
}


int validate_base_256_17(char **sk_A, char **sk_B, int num_test_cases) {
    /*
    validate by running an instance of DH
    */
	radix17_t Pkx_A, Pky_A, Pkx_B, Pky_B, ssx_A, ssy_A, ssx_B, ssy_B;
	set_ui_256( &Pkx_A, 0 );
	set_ui_256( &Pky_A, 0 );
	set_ui_256( &Pkx_B, 0 );
	set_ui_256( &Pky_B, 0 );
	set_ui_256( &ssx_A, 0 );
	set_ui_256( &ssy_A, 0 );
	set_ui_256( &ssx_B, 0 );
	set_ui_256( &ssy_B, 0 );
	
	int ok = 1;
	
	for ( int i = 0; i < num_test_cases; ++i ) {
		// 1. generate Alice's public key
		public_key_gen_base_256_17(sk_A[i], &Pkx_A, &Pky_A);
    	// 2. generate Bob's public key
		public_key_gen_base_256_17(sk_B[i], &Pkx_B, &Pky_B);
		// 3. Generate the shared secret
		scalar_mult_base_256_17(&Pkx_B, &Pky_B, sk_A[i], &ssx_A, &ssy_A);
		scalar_mult_base_256_17(&Pkx_A, &Pky_A, sk_B[i], &ssx_B, &ssy_B);
		if ( ! cmp_coord_256( &ssx_A, &ssx_B ) ) {
			printf( "     Validation failed for key '%s' found in row '%d'\n", sk_A[i], i+1 );
			ok = 0;
		}
	}
	return ok;
}

int checkPk_base_256_17(char **sk, char **pk, int num_test_cases) {
    /*
    Generate PKs and compare them
    */
	radix17_t Pkx, Pky;
	set_ui_256( &Pkx, 0 );
	set_ui_256( &Pky, 0 );
	
	int ok = 1;
	
	char pkAsString[513];
	
	for ( int i = 0; i < num_test_cases; ++i ) {
		public_key_gen_base_256_17(sk[i], &Pkx, &Pky);
		radix17_toString( pkAsString, Pkx );
		if ( strncmp( pk[i], pkAsString, 513 ) ) {
			printf( "     Public key generation failed for key found in row %d (computed key: '%s', expected key: '%s')\n", i+1, pkAsString, pk[i] );
			ok = 0;
		}
	}
	return ok;
}
