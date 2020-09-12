#include "curve25519_base_256_17.h"
#include "../256arithmetic/bignum.h"
#include "../../_common_libs/vec_utils.h"
#include "../../_common_libs/radix17_vec.h"
#include "../../_common_libs/radix17_curveParams_vec.h"
#include "../../_common_libs/common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>

void point_double_base_256_17(__m256i *resXlo, __m256i *resXhi, __m256i *resYlo,  __m256i *resYhi, __m256i xlo, __m256i xhi, __m256i ylo, __m256i yhi)
{
    /*
    res = 2*(x,y)
    */
    if(iszero_256(xlo, xhi) && iszero_256(ylo, yhi)){  
        copy_256(resXlo, resXhi, xlo, xhi);
        copy_256(resYlo, resYhi, ylo, yhi);
        return;
    }

    __m256i upperfraclo, upperfrachi, lowerfraclo, lowerfrachi, tmp1lo, tmp1hi, tmp2lo, tmp2hi, tmp3lo, tmp3hi, llo, lhi;
    set_zero(&upperfraclo, &upperfrachi);
    set_zero(&lowerfraclo, &lowerfrachi);
    set_zero(&tmp1lo, &tmp1hi);
    set_zero(&tmp2lo, &tmp2hi);
    set_zero(&tmp3lo, &tmp3hi);
    set_zero(&llo, &lhi);    

    //x - coordinate
    mul_256(&tmp1lo, &tmp1hi, A_rad17lo, A_rad17hi, xlo, xhi);               // tmp1 = A*x
    add_256(&tmp1lo, &tmp1hi, tmp1lo, tmp1hi, tmp1lo, tmp1hi);          // tmp1 = 2*tmp1
    add_256(&tmp1lo, &tmp1hi, tmp1lo, tmp1hi, const_1_rad17lo, const_1_rad17hi);       // tmp1 = tmp1 + 1 

    square_256(&tmp2lo, &tmp2hi, xlo, xhi);                // tmp2 = x^2
    mul_256(&tmp2lo, &tmp2hi, tmp2lo, tmp2hi, const_3_rad17lo, const_3_rad17hi);       // tmp2 = tmp2*3

    add_256(&upperfraclo, &upperfrachi, tmp2lo, tmp2hi, tmp1lo, tmp1hi);     // upperfrac = tmp2+tmp1

    mul_256(&lowerfraclo, &lowerfrachi, B_rad17lo, B_rad17hi, ylo, yhi);          // lowerfrac = B*y
    add_256(&lowerfraclo, &lowerfrachi, lowerfraclo, lowerfrachi, lowerfraclo, lowerfrachi); // lowwerfrac = 2*lowerfrac

    div_256(&llo, &lhi, upperfraclo, upperfrachi, lowerfraclo, lowerfrachi);   // l = uperfrac/lowerfrac

    square_256(&llo, &lhi, llo, lhi);                    // l=l^2

    mul_256(resXlo, resXhi, B_rad17lo, B_rad17hi, llo, lhi);           // resX = B*l
    sub_256(resXlo, resXhi, *resXlo, *resXhi, A_rad17lo, A_rad17hi);    // resX = resX - A
    sub_256(resXlo, resXhi, *resXlo, *resXhi, xlo, xhi);     // resX = resX - x
    sub_256(resXlo, resXhi, *resXlo, *resXhi, xlo, xhi);     // resX = resX - x

    //y - coordinate
    // l, upperfrac, lowerfrac, tmp1, tmp2 don't need to be reinitialized
    // l = 3x^2 + 2 Ax + 1 / 2By   
    mul_256(&tmp1lo, &tmp1hi, A_rad17lo, A_rad17hi, xlo, xhi);               // tmp1 = A*x
    add_256(&tmp1lo, &tmp1hi, tmp1lo, tmp1hi, tmp1lo, tmp1hi);          // tmp1 = 2*tmp1
    add_256(&tmp1lo, &tmp1hi, tmp1lo, tmp1hi, const_1_rad17lo, const_1_rad17hi);       // tmp1 = tmp1 + 1

    square_256(&tmp2lo, &tmp2hi, xlo, xhi);                // tmp2 = x^2
    mul_256(&tmp2lo, &tmp2hi, tmp2lo, tmp2hi, const_3_rad17lo, const_3_rad17hi);       // tmp2 = tmp2*3

    add_256(&upperfraclo, &upperfrachi, tmp2lo, tmp2hi, tmp1lo, tmp1hi);     // upperfrac = tmp2/tmp1

    mul_256(&lowerfraclo, &lowerfrachi, B_rad17lo, B_rad17hi, ylo, yhi);          // lowerfrac = B*y
    add_256(&lowerfraclo, &lowerfrachi, lowerfraclo, lowerfrachi, lowerfraclo, lowerfrachi);   // lowerfrac = 2*lowerfrac

    div_256(&llo, &lhi, upperfraclo, upperfrachi, lowerfraclo, lowerfrachi);   // l = upperfrac/lowerfrac

    mul_256(&tmp1lo, &tmp1hi, xlo, xhi, const_3_rad17lo, const_3_rad17hi);           // tmp1 = 3*x
    add_256(&tmp1lo, &tmp1hi, tmp1lo, tmp1hi, A_rad17lo, A_rad17hi);                // tmp1 = tmp1+A

    mul_256(&tmp1lo, &tmp1hi, tmp1lo, tmp1hi, llo, lhi);             // tmp1 = tmp1*l

    square_256(&tmp3lo, &tmp3hi,llo, lhi);                 // tmp3=l^2
    mul_256(&llo, &lhi, tmp3lo, tmp3hi, llo, lhi);                // l = l^3

    mul_256(&llo, &lhi, llo, lhi, B_rad17lo, B_rad17hi);                  // l = l*B

    sub_256(resYlo, resYhi, tmp1lo, tmp1hi, llo, lhi);         // resY = tmp1 - l

    sub_256(resYlo, resYhi, *resYlo, *resYhi, ylo, yhi);     // resY = resY - y

    
}

void point_add_base_256_17(__m256i *resXlo, __m256i *resXhi, __m256i *resYlo, __m256i *resYhi, __m256i op1xlo, __m256i op1xhi, __m256i op1ylo, __m256i op1yhi, __m256i op2xlo, __m256i op2xhi, __m256i op2ylo, __m256i op2yhi)
{    
    // resX <- op1x + op2x
    // resY <- op1y + op2y
    // https://en.wikipedia.org/wiki/Montgomery_curve#Addition
    
    // return (0,0) if ( op1x == op2x AND op1y != op2y )
    if ( cmp_coord_256( op1xlo, op1xhi, op2xlo, op2xhi ) && !cmp_coord_256(op1ylo, op1yhi, op2ylo, op2yhi)) {
        set_zero(resXlo, resXhi);
        set_zero(resYlo, resYhi);
    	return;
    }
    
    // return op2 if the point op1 is (0,0)
    if ( iszero_256(op1xlo, op1xhi) && iszero_256( op1ylo, op1yhi ) ) {  
    	// probably not side channel safe anymore ?
    	copy_256( resXlo, resXhi, op2xlo, op2xhi );
    	copy_256( resYlo, resYhi, op2ylo, op2yhi );
    	return;
    }
    
    // return op1 if the point op2 is (0,0)
    if ( iszero_256( op2xlo, op2xhi ) && iszero_256( op2ylo, op2yhi ) ) {  
    	// probably not side channel safe anymore ?
    	copy_256( resXlo, resXhi, op1xlo, op1xhi );
    	copy_256( resYlo, resYhi, op1ylo, op1yhi );
    	return;
    }

    if ( cmp_point_256( op1xlo, op1xhi, op1ylo, op1yhi, op2xlo, op2xhi, op2ylo, op2yhi ) ) { 
    	// point_add_base_256_17(p1, p1) = double(p1), otherwise exception div by 0
    	point_double_base_256_17( resXlo, resXhi, resYlo, resYhi, op1xlo, op1xhi, op1ylo, op1yhi );
    	return;
    }

    __m256i diffxlo, diffxhi, diffylo, diffyhi, tmp1lo, tmp1hi, tmp2lo, tmp2hi;
    set_zero(&diffxlo, &diffxhi);
    set_zero(&diffylo, &diffyhi);
    set_zero(&tmp1lo, &tmp1hi);
    set_zero(&tmp2lo, &tmp2hi);    

    // x - coordinate
    sub_256( &diffylo, &diffyhi, op2ylo, op2yhi, op1ylo, op1yhi );			// diffy = y2 - y1
    sub_256( &diffxlo, &diffxhi, op2xlo, op2xhi, op1xlo, op1xhi );			// diffx = x2 - x1
    
    square_256( &diffylo, &diffyhi, diffylo, diffyhi );				// diffy = diffy^2
    square_256( &diffxlo, &diffxhi, diffxlo, diffxhi );				// diffx = diffx^2
    
    mul_256( &diffylo, &diffyhi, diffylo, diffyhi, B_rad17lo, B_rad17hi );			// diffy = B * diffy
    div_256( resXlo, resXhi, diffylo, diffyhi, diffxlo, diffxhi );		// resX = diffy / diffx
    
    sub_256( resXlo, resXhi, *resXlo, *resXhi, A_rad17lo, A_rad17hi );		// resX -= A
    sub_256( resXlo, resXhi, *resXlo, *resXhi, op1xlo, op1xhi );	// resX -= x1
    sub_256( resXlo, resXhi, *resXlo, *resXhi, op2xlo, op2xhi );	// resX -= x2

    // y - coordinate
	// Note: no need to re-initialize the local variables diffx and diffy as they are always overwritten
	
	add_256( &tmp1lo, &tmp1hi, op1xlo, op1xhi, op2xlo, op2xhi );			// tmp1 = x1 + x2
	add_256( &tmp1lo, &tmp1hi, tmp1lo, tmp1hi, op1xlo, op1xhi );			// tmp1 = tmp1 + x1
	add_256( &tmp1lo, &tmp1hi, tmp1lo, tmp1hi, A_rad17lo, A_rad17hi   );			// tmp1 = tmp1 + A, i.e. tmp1 = 2*x1 + x2 + A
	
	sub_256( &diffylo, &diffyhi, op2ylo, op2yhi, op1ylo, op1yhi );			// diffy = y2 - y1
	mul_256( &tmp1lo, &tmp1hi, tmp1lo, tmp1hi, diffylo, diffyhi );			// tmp1 = tmp1 * diffy
	
	sub_256( &diffxlo, &diffxhi, op2xlo, op2xhi, op1xlo, op1xhi );			// diffx = x2 - x1
	div_256( &tmp1lo, &tmp1hi, tmp1lo, tmp1hi, diffxlo, diffxhi );			// tmp1 = tmp1 / diffx
	
	square_256( &tmp2lo, &tmp2hi, diffxlo, diffxhi );
	mul_256   ( &diffxlo, &diffxhi, tmp2lo, tmp2hi, diffxlo, diffxhi );		// diffx = diffx ^ 3
	square_256( &tmp2lo, &tmp2hi, diffylo, diffyhi );
	mul_256   ( &diffylo, &diffyhi, tmp2lo, tmp2hi, diffylo, diffyhi );		// diffy = diffy ^ 3
	mul_256   ( &diffylo, &diffyhi, diffylo, diffyhi, B_rad17lo, B_rad17hi );			// diffy = diffy * B
	
	div_256( &tmp2lo, &tmp2hi, diffylo, diffyhi, diffxlo, diffxhi );		// tmp2 = diffy / diffx
	sub_256( resYlo, resYhi, tmp1lo, tmp1hi, tmp2lo, tmp2hi );		// resY = tmp1 - tmp2
	sub_256( resYlo, resYhi, *resYlo, *resYhi, op1ylo, op1yhi );	    // resY -= y1
	
	
}

void scalar_mult_base_256_17(__m256i base_point_xlo, __m256i base_point_xhi, __m256i base_point_ylo, __m256i base_point_yhi, char *sk, __m256i *Pkxlo, __m256i *Pkxhi, __m256i *Pkylo, __m256i *Pkyhi)
{
    /*
    Pk = sk*base_point
    */
    int n = sk_length;
    __m256i R0xlo, R0xhi, R0ylo, R0yhi, R1xlo, R1xhi, R1ylo, R1yhi, resXlo, resXhi, resYlo, resYhi;
    set_zero(&R0xlo, &R0xhi);
    set_zero(&R0ylo, &R0yhi);

    copy_256(&R1xlo, &R1xhi, base_point_xlo, base_point_xhi);
    copy_256(&R1ylo, &R1yhi, base_point_ylo, base_point_yhi);

    for(int i=0; i<n; i++){

        set_zero(&resXlo, &resXhi);
        set_zero(&resYlo, &resYhi);


        if(sk[i]=='0'){
            // R1 <- add (R0, R1)
            point_add_base_256_17(&resXlo, &resXhi, &resYlo, &resYhi, R0xlo, R0xhi, R0ylo, R0yhi, R1xlo, R1xhi, R1ylo, R1yhi);
            copy_256(&R1xlo, &R1xhi, resXlo, resXhi);
            copy_256(&R1ylo, &R1yhi, resYlo, resYhi);

            set_zero(&resXlo, &resXhi);
            set_zero(&resYlo, &resYhi);
            // R0 <- double(R0)
            point_double_base_256_17(&resXlo, &resXhi, &resYlo, &resYhi, R0xlo, R0xhi, R0ylo, R0yhi);
            copy_256(&R0xlo, &R0xhi, resXlo, resXhi);
            copy_256(&R0ylo, &R0yhi, resYlo, resYhi);
                
        }else{
            // R0 <- add(R0, R1)
            point_add_base_256_17(&resXlo, &resXhi, &resYlo, &resYhi, R0xlo, R0xhi, R0ylo, R0yhi, R1xlo, R1xhi, R1ylo, R1yhi);
            copy_256(&R0xlo, &R0xhi, resXlo, resXhi);
            copy_256(&R0ylo, &R0yhi, resYlo, resYhi);

            set_zero(&resXlo, &resXhi);
            set_zero(&resYlo, &resYhi);
            // R1 <- double(R1)
            point_double_base_256_17(&resXlo, &resXhi, &resYlo, &resYhi, R1xlo, R1xhi, R1ylo, R1yhi);
            copy_256(&R1xlo, &R1xhi, resXlo, resXhi);
            copy_256(&R1ylo, &R1yhi, resYlo, resYhi);
                
        }
    }
    // Assign Pk 
    copy_256(Pkxlo, Pkxhi, R0xlo, R0xhi);
    copy_256(Pkylo, Pkyhi, R0ylo, R0yhi);
}

void public_key_gen_base_256_17(char *sk, __m256i *xlo, __m256i *xhi, __m256i *ylo, __m256i *yhi)
{
    /*
    (x,y) = sk*G
    */
    scalar_mult_base_256_17(Gx_rad17lo, Gx_rad17hi, Gy_rad17lo, Gy_rad17hi, sk, xlo, xhi, ylo, yhi);
}


int validate_base_256_17(char **sk_A, char **sk_B, int num_test_cases) {
    /*
    validate by running an instance of DH
    */
	__m256i Pkx_Alo, Pkx_Ahi, Pky_Alo, Pky_Ahi, Pkx_Blo, Pkx_Bhi, Pky_Blo, Pky_Bhi, ssx_Alo, ssx_Ahi, ssy_Alo, ssy_Ahi, ssx_Blo, ssx_Bhi, ssy_Blo, ssy_Bhi;
    set_zero(&Pkx_Alo, &Pkx_Ahi);
    set_zero(&Pky_Alo, &Pky_Ahi);
    set_zero(&Pkx_Blo, &Pkx_Bhi);
    set_zero(&Pky_Blo, &Pky_Bhi);
    set_zero(&ssx_Alo, &ssx_Ahi);
    set_zero(&ssy_Alo, &ssy_Ahi);
    set_zero(&ssx_Blo, &ssx_Bhi);
    set_zero(&ssy_Blo, &ssy_Bhi);
	
	int ok = 1;
	
	for ( int i = 0; i < num_test_cases; ++i ) {
		// 1. generate Alice's public key
		public_key_gen_base_256_17(sk_A[i], &Pkx_Alo, &Pkx_Ahi, &Pky_Alo, &Pky_Ahi);
    	// 2. generate Bob's public key
		public_key_gen_base_256_17(sk_B[i], &Pkx_Blo, &Pkx_Bhi, &Pky_Blo, &Pky_Bhi);
		// 3. Generate the shared secret
		scalar_mult_base_256_17(Pkx_Blo, Pkx_Bhi, Pky_Blo, Pky_Bhi, sk_A[i], &ssx_Alo, &ssx_Ahi, &ssy_Alo, &ssy_Ahi);
		scalar_mult_base_256_17(Pkx_Alo, Pkx_Ahi, Pky_Alo, Pky_Ahi, sk_B[i], &ssx_Blo, &ssx_Bhi, &ssy_Blo, &ssy_Bhi);
		if ( ! cmp_coord_256( ssx_Alo, ssx_Ahi, ssx_Blo, ssx_Bhi ) ) {
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
	__m256i Pkxlo, Pkxhi, Pkylo, Pkyhi;
	set_zero( &Pkxlo, &Pkxhi);
	set_zero( &Pkylo, &Pkyhi);
	
	int ok = 1;
	
	char pkAsString[513];
	
	for ( int i = 0; i < num_test_cases; ++i ) {
		public_key_gen_base_256_17(sk[i], &Pkxlo, &Pkxhi, &Pkylo, &Pkyhi);
		radix17_vec_toString( pkAsString, Pkxlo, Pkxhi );
		if ( strncmp( pk[i], pkAsString, 513 ) ) {
			printf( "     Public key generation failed for key found in row %d (computed key: '%s', expected key: '%s')\n", i+1, pkAsString, pk[i] );
			ok = 0;
		}
	}
	return ok;
}
