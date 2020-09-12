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
    if (iszero_256(xlo, xhi) && iszero_256(ylo, yhi))
    {
        copy_256(resXlo, resXhi, xlo, xhi);
        copy_256(resYlo, resYhi, ylo, yhi);
        return;
    }

    __m256i _2xlo, Axxlo, x2lo, _2x2lo, _3x2lo, _2ylo, llo, l2lo, l3lo, tmp1lo, tmp2lo, onelo,
            _2xhi, Axxhi, x2hi, _2x2hi, _3x2hi, _2yhi, lhi, l2hi, l3hi, tmp1hi, tmp2hi, onehi;
    // 2x = x+x
    add_256(&_2xlo, &_2xhi, xlo, xhi, xlo, xhi);
    // Axx = A + 2x
    add_256(&Axxlo, &Axxhi, A_rad17lo, A_rad17hi, _2xlo, _2xhi);
    // x2 = x^2
    square_256(&x2lo, &x2hi, xlo, xhi);
    // 2x^2 = x^2 + x^2
    add_256(&_2x2lo, &_2x2hi, x2lo, x2hi, x2lo, x2hi);
    // 3x^2 = 2x^2 + x^2
    add_256(&_3x2lo, &_3x2hi, _2x2lo, _2x2hi, x2lo, x2hi);
    // 2y = y + y
    add_256(&_2ylo, &_2yhi, ylo, yhi, ylo, yhi);

    // l = (_3x2 + A*_2x + 1)/(_2y)
    mul_256(&tmp1lo, &tmp1hi, A_rad17lo, A_rad17hi, _2xlo, _2xhi);
    add_256(&tmp1lo, &tmp1hi, tmp1lo, tmp1hi, _3x2lo, _3x2hi);
    add_256(&tmp1lo, &tmp1hi, tmp1lo, tmp1hi, const_1_rad17lo, const_1_rad17hi);
    div_256(&llo, &lhi, tmp1lo, tmp1hi, _2ylo, _2yhi);
    // l2 = l^2
    square_256(&l2lo, &l2hi, llo, lhi);
    // l3 = l2*l
    mul_256(&l3lo, &l3hi, l2lo, l2hi, llo, lhi);

    // x - coordinate
    // resX = l2 - Axx
    sub_256(resXlo, resXhi, l2lo, l2hi, Axxlo, Axxhi);

    // y - coordinate
    // resY = (Axx + x)*l - (l3 + y)
    add_256(&tmp1lo, &tmp1hi, Axxlo, Axxhi, xlo, xhi);
    mul_256(&tmp1lo, &tmp1hi, tmp1lo, tmp1hi, llo, lhi);
    add_256(&tmp2lo, &tmp2hi, l3lo, l3hi, ylo, yhi);
    sub_256(resYlo, resYhi, tmp1lo, tmp1hi, tmp2lo, tmp2hi);
}

void point_add_base_256_17(__m256i *resXlo, __m256i *resXhi, __m256i *resYlo, __m256i *resYhi, __m256i op1xlo, __m256i op1xhi, __m256i op1ylo, __m256i op1yhi, __m256i op2xlo, __m256i op2xhi, __m256i op2ylo, __m256i op2yhi)
{
    // resX <- op1x + op2x
    // resY <- op1y + op2y
    // https://en.wikipedia.org/wiki/Montgomery_curve#Addition

    // return (0,0) if ( op1x == op2x AND op1y != op2y )
    if (cmp_coord_256(op1xlo, op1xhi, op2xlo, op2xhi) && !cmp_coord_256(op1ylo, op1yhi, op2ylo, op2yhi))
    {
        set_ui_256(resXlo, resXhi, 0);
        set_ui_256(resYlo, resYhi, 0);
        return;
    }

    // return op2 if the point op1 is (0,0)
    if (iszero_256(op1xlo, op1xhi) && iszero_256(op1ylo, op1yhi))
    {
        // probably not side channel safe anymore ?
        copy_256(resXlo, resXhi, op2xlo, op2xhi);
        copy_256(resYlo, resYhi, op2ylo, op2yhi);
        return;
    }

    // return op1 if the point op2 is (0,0)
    if (iszero_256(op2xlo, op2xhi) && iszero_256(op2ylo, op2yhi))
    {
        // probably not side channel safe anymore ?
        copy_256(resXlo, resXhi, op1xlo, op1xhi);
        copy_256(resYlo, resYhi, op1ylo, op1yhi);
        return;
    }

    if (cmp_point_256(op1xlo, op1xhi, op1ylo, op1yhi, op2xlo, op2xhi, op2ylo, op2yhi))
    {
        // point_add_base_256_17(p1, p1) = double(p1), otherwise exception div by 0
        point_double_base_256_17(resXlo, resXhi, resYlo, resYhi, op1xlo, op1xhi, op1ylo, op1yhi);
        return;
    }

    __m256i   diffxlo, diffylo, fraclo, frac2lo, frac3lo, Axxlo, tmp1lo, tmp2lo,
              diffxhi, diffyhi, frachi, frac2hi, frac3hi, Axxhi, tmp1hi, tmp2hi;
    // diffx = x2-x1
    sub_256(&diffxlo, &diffxhi, op2xlo, op2xhi, op1xlo, op1xhi);
    // diffy = y2-y1
    sub_256(&diffylo, &diffyhi, op2ylo, op2yhi, op1ylo, op1yhi);
    // frac = diffy/diffx
    div_256(&fraclo, &frachi, diffylo, diffyhi, diffxlo, diffxhi);
    // frac2 = frac^2
    square_256(&frac2lo, &frac2hi, fraclo, frachi);
    // frac3 = frac*frac^2
    mul_256(&frac3lo, &frac3hi, fraclo, frachi, frac2lo, frac2hi);
    // Axx = A + x1 + x2
    add_256(&tmp1lo, &tmp1hi, A_rad17lo, A_rad17hi, op1xlo, op1xhi);
    add_256(&Axxlo, &Axxhi, tmp1lo, tmp1hi, op2xlo, op2xhi);
    // x - coordinate
    // x3 = frac2-Axx
    sub_256(resXlo, resXhi, frac2lo, frac2hi, Axxlo, Axxhi);
    // y - coordinate
    // y3 = (x1 + Axx)*frac - (frac3 + y1)
    add_256(&tmp1lo, &tmp1hi, op1xlo, op1xhi, Axxlo, Axxhi);
    mul_256(&tmp1lo, &tmp1hi, tmp1lo, tmp1hi, fraclo, frachi);
    add_256(&tmp2lo, &tmp2hi, frac3lo, frac3hi, op1ylo, op1yhi);
    sub_256(resYlo, resYhi, tmp1lo, tmp1hi, tmp2lo, tmp2hi);
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
            point_add_base_256_17(&R1xlo, &R1xhi, &R1ylo, &R1yhi, R0xlo, R0xhi, R0ylo, R0yhi, R1xlo, R1xhi, R1ylo, R1yhi);
            //copy_256(&R1xlo, &R1xhi, resXlo, resXhi);
            //copy_256(&R1ylo, &R1yhi, resYlo, resYhi);
            // R0 <- double(R0)
            point_double_base_256_17(&R0xlo, &R0xhi, &R0ylo, &R0yhi, R0xlo, R0xhi, R0ylo, R0yhi);
            // copy_256(&R0xlo, &R0xhi, resXlo, resXhi);
            // copy_256(&R0ylo, &R0yhi, resYlo, resYhi);
                
        }else{
            // R0 <- add(R0, R1)
            point_add_base_256_17(&R0xlo, &R0xhi, &R0ylo, &R0yhi, R0xlo, R0xhi, R0ylo, R0yhi, R1xlo, R1xhi, R1ylo, R1yhi);
            // copy_256(&R0xlo, &R0xhi, resXlo, resXhi);
            // copy_256(&R0ylo, &R0yhi, resYlo, resYhi);
            // R1 <- double(R1)
            point_double_base_256_17(&R1xlo, &R1xhi, &R1ylo, &R1yhi, R1xlo, R1xhi, R1ylo, R1yhi);
            // copy_256(&R1xlo, &R1xhi, resXlo, resXhi);
            // copy_256(&R1ylo, &R1yhi, resYlo, resYhi);
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
