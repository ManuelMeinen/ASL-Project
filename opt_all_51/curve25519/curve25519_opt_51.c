#include "curve25519_opt_51.h"
#include "../256arithmetic/bignum.h"
#include "../../_common_libs/radix51_arr.h"
#include "../../_common_libs/radix51_curveParams_arr.h"
#include "../../_common_libs/common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>

void point_double_opt_all_51(long resX[5], long resY[5], long x[5], long y[5]) {
    /*
    res = 2*(x,y)
    */
    if (iszero_256(x) && iszero_256(y))
    {
        copy_256(resX, x);
        copy_256(resY, y);
        return;
    }

    long _2x[5], Axx[5], x2[5], _2x2[5], _3x2[5], _2y[5], l[5], l2[5], l3[5], tmp1[5], tmp2[5], one[5];
    // 2x = x+x
    add_256(_2x, x, x);
    // Axx = A + 2x
    add_256(Axx, A_rad51_arr, _2x);
    // x2 = x^2
    square_256(x2, x);
    // 2x^2 = x^2 + x^2
    add_256(_2x2, x2, x2);
    // 3x^2 = 2x^2 + x^2
    add_256(_3x2, _2x2, x2);
    // 2y = y + y
    add_256(_2y, y, y);

    // l = (_3x2 + A*_2x + 1)/(_2y)
    mul_256(tmp1, A_rad51_arr, _2x);
    add_256(tmp1, tmp1, _3x2);
    add_1_256(tmp1);
    div_256(l, tmp1, _2y);
    // l2 = l^2
    square_256(l2, l);
    // l3 = l2*l
    mul_256(l3, l2, l);

	// save the value of x
	long xcopy[5];
	copy_256( xcopy, x );
	
    // x - coordinate
    // resX = l2 - Axx
    sub_256(resX, l2, Axx);

    // y - coordinate
    // resY = (Axx + x)*l - (l3 + y)
    add_256(tmp1, Axx, xcopy);
    mul_256(tmp1, tmp1, l);
    add_256(tmp2, l3, y);
    sub_256(resY, tmp1, tmp2);
}

void point_add_opt_all_51(long resX[5], long resY[5], long op1x[5], long op1y[5], long op2x[5], long op2y[5]) {
    // resX <- op1x + op2x
    // resY <- op1y + op2y
    // https://en.wikipedia.org/wiki/Montgomery_curve#Addition

    // return (0,0) if ( op1x == op2x AND op1y != op2y )
    if (cmp_coord_256(op1x, op2x) && !cmp_coord_256(op1y, op2y))
    {
        set_ui_256(resX, 0);
        set_ui_256(resY, 0);
        return;
    }

    // return op2 if the point op1 is (0,0)
    if (iszero_256(op1x) && iszero_256(op1y))
    {
        // probably not side channel safe anymore ?
        copy_256(resX, op2x);
        copy_256(resY, op2y);
        return;
    }

    // return op1 if the point op2 is (0,0)
    if (iszero_256(op2x) && iszero_256(op2y))
    {
        // probably not side channel safe anymore ?
        copy_256(resX, op1x);
        copy_256(resY, op1y);
        return;
    }

    if (cmp_point_256(op1x, op1y, op2x, op2y))
    {
        // point_add_opt_all_51(p1, p1) = double(p1), otherwise exception div by 0
        point_double_opt_all_51(resX, resY, op1x, op1y);
        return;
    }

    long diffx[5], diffy[5], frac[5], frac2[5], frac3[5], Axx[5], tmp1[5], tmp2[5];
    // diffx = x2-x1
    sub_256(diffx, op2x, op1x);
    // diffy = y2-y1
    sub_256(diffy, op2y, op1y);
    // frac = diffy/diffx
    div_256(frac, diffy, diffx);
    // frac2 = frac^2
    square_256(frac2, frac);
    // frac3 = frac*frac^2
    mul_256(frac3, frac, frac2);
    // Axx = A + x1 + x2
    add_256(tmp1, A_rad51_arr, op1x);
    add_256(Axx, tmp1, op2x);
    // x - coordinate
    // x3 = frac2-Axx
    sub_256(resX, frac2, Axx);
    // y - coordinate
    // y3 = (x1 + Axx)*frac - (frac3 + y1)
    add_256(tmp1, op1x, Axx);
    mul_256(tmp1, tmp1, frac);
    add_256(tmp2, frac3, op1y);
    sub_256(resY, tmp1, tmp2);
}

void scalar_mult_opt_all_51(long base_point_x[5], long base_point_y[5], char *sk, long Pkx[5], long Pky[5]) {
    /*
    Pk = sk*base_point
    */
    int n = sk_length;
    long R0x[5], R0y[5], R1x[5], R1y[5], resX[5], resY[5];
    
    // initialize R0x, R0y, R1x, R1y
    set_ui_256(R0x, 0);
    set_ui_256(R0y, 0);
    copy_256(R1x, base_point_x);
    copy_256(R1y, base_point_y);

    for (int i = 0; i < n; i++)
    {
        if (sk[i] == '0')
        {
            // R1 <- add (R0, R1)
            point_add_opt_all_51(R1x, R1y, R0x, R0y, R1x, R1y);
            // R0 <- double(R0)
            point_double_opt_all_51(R0x, R0y, R0x, R0y);
        }
        else
        {
            // R0 <- add(R0, R1)
            // Note: R1x and R1y must be the first operand, otherwise aliasing issues will arise due to R0x and R0y being overwritten
            point_add_opt_all_51(R0x, R0y, R1x, R1y, R0x, R0y);
            // R1 <- double(R1)
            point_double_opt_all_51(R1x, R1y, R1x, R1y);
            /*
            point_double_opt_all_51(resX, resY, R1x, R1y);
            copy_256(R1x, resX);
            copy_256(R1y, resY);
            */
        }
    }
    // Assign Pk
    copy_256(Pkx, R0x);
    copy_256(Pky, R0y);
}

void public_key_gen_opt_all_51( char *sk, long Pkx[5], long Pky[5] ) {
	scalar_mult_opt_all_51( Gx_rad51_arr, Gy_rad51_arr, sk, Pkx, Pky );
}

int validate_opt_all_51(char **sk_A, char **sk_B, int num_test_cases)
{
    /*
    validate by running an instance of DH
    */
    long Pkx_A[5], Pky_A[5], Pkx_B[5], Pky_B[5], ssx_A[5], ssy_A[5], ssx_B[5], ssy_B[5];
    set_ui_256(Pkx_A, 0);
    set_ui_256(Pky_A, 0);
    set_ui_256(Pkx_B, 0);
    set_ui_256(Pky_B, 0);
    set_ui_256(ssx_A, 0);
    set_ui_256(ssy_A, 0);
    set_ui_256(ssx_B, 0);
    set_ui_256(ssy_B, 0);

    int ok = 1;

    for (int i = 0; i < num_test_cases; ++i)
    {
        // 1. generate Alice's public key
        public_key_gen_opt_all_51(sk_A[i], Pkx_A, Pky_A);
        // 2. generate Bob's public key
        public_key_gen_opt_all_51(sk_B[i], Pkx_B, Pky_B);
        // 3. Generate the shared secret
        scalar_mult_opt_all_51(Pkx_B, Pky_B, sk_A[i], ssx_A, ssy_A);
        scalar_mult_opt_all_51(Pkx_A, Pky_A, sk_B[i], ssx_B, ssy_B);
        if (!cmp_coord_256(ssx_A, ssx_B))
        {
            printf("     Validation failed for key '%s' found in row '%d'\n", sk_A[i], i + 1);
            ok = 0;
        }
    }
    return ok;
}

int checkPk_opt_all_51(char **sk, char **pk, int num_test_cases)
{
    /*
    Generate PKs and compare them
    */
    long Pkx[5], Pky[5];
    set_ui_256(Pkx, 0);
    set_ui_256(Pky, 0);

    int ok = 1;

    char pkAsString[513];

    for (int i = 0; i < num_test_cases; ++i)
    {
        public_key_gen_opt_all_51(sk[i], Pkx, Pky);
        radix51_arr_toString(pkAsString, Pkx);
        if (strncmp(pk[i], pkAsString, 513))
        {
            printf("     Public key generation failed for key found in row %d (computed key: '%s', expected key: '%s')\n", i + 1, pkAsString, pk[i]);
            ok = 0;
        }
    }
    return ok;
}
