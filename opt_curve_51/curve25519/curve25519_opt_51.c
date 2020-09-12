#include "curve25519_opt_51.h"
#include "../../base_256_51/256arithmetic/bignum.h"
#include "../../_common_libs/radix51.h"
#include "../../_common_libs/radix51_curveParams.h"
#include "../../_common_libs/common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>

void point_double_base_256_51(radix51_t *resX, radix51_t *resY, radix51_t *x, radix51_t *y)
{
    /*
    res = 2*(x,y)
    */
    if (iszero_256(x) && iszero_256(y))
    {
        copy_256(resX, *x);
        copy_256(resY, *y);
        return;
    }

    radix51_t _2x, Axx, x2, _2x2, _3x2, _2y, l, l2, l3, tmp1, tmp2, one;
    // 2x = x+x
    add_256(&_2x, x, x);
    // Axx = A + 2x
    add_256(&Axx, &A_rad51, &_2x);
    // x2 = x^2
    square_256(&x2, x);
    // 2x^2 = x^2 + x^2
    add_256(&_2x2, &x2, &x2);
    // 3x^2 = 2x^2 + x^2
    add_256(&_3x2, &_2x2, &x2);
    // 2y = y + y
    add_256(&_2y, y, y);

    // l = (_3x2 + A*_2x + 1)/(_2y)
    mul_256(&tmp1, &A_rad51, &_2x);
    add_256(&tmp1, &tmp1, &_3x2);
    add_256(&tmp1, &tmp1, &const_1_rad51);
    div_256(&l, &tmp1, &_2y);
    // l2 = l^2
    square_256(&l2, &l);
    // l3 = l2*l
    mul_256(&l3, &l2, &l);

    // x - coordinate
    // resX = l2 - Axx
    sub_256(resX, &l2, &Axx);

    // y - coordinate
    // resY = (Axx + x)*l - (l3 + y)
    add_256(&tmp1, &Axx, x);
    mul_256(&tmp1, &tmp1, &l);
    add_256(&tmp2, &l3, y);
    sub_256(resY, &tmp1, &tmp2);
}

void point_add_base_256_51(radix51_t *resX, radix51_t *resY, radix51_t *op1x, radix51_t *op1y, radix51_t *op2x, radix51_t *op2y)
{
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
        copy_256(resX, *op2x);
        copy_256(resY, *op2y);
        return;
    }

    // return op1 if the point op2 is (0,0)
    if (iszero_256(op2x) && iszero_256(op2y))
    {
        // probably not side channel safe anymore ?
        copy_256(resX, *op1x);
        copy_256(resY, *op1y);
        return;
    }

    if (cmp_point_256(op1x, op1y, op2x, op2y))
    {
        // point_add_base_256_51(p1, p1) = double(p1), otherwise exception div by 0
        point_double_base_256_51(resX, resY, op1x, op1y);
        return;
    }

    radix51_t diffx, diffy, frac, frac2, frac3, Axx, tmp1, tmp2;
    // diffx = x2-x1
    sub_256(&diffx, op2x, op1x);
    // diffy = y2-y1
    sub_256(&diffy, op2y, op1y);
    // frac = diffy/diffx
    div_256(&frac, &diffy, &diffx);
    // frac2 = frac^2
    square_256(&frac2, &frac);
    // frac3 = frac*frac^2
    mul_256(&frac3, &frac, &frac2);
    // Axx = A + x1 + x2
    add_256(&tmp1, &A_rad51, op1x);
    add_256(&Axx, &tmp1, op2x);
    // x - coordinate
    // x3 = frac2-Axx
    sub_256(resX, &frac2, &Axx);
    // y - coordinate
    // y3 = (x1 + Axx)*frac - (frac3 + y1)
    add_256(&tmp1, op1x, &Axx);
    mul_256(&tmp1, &tmp1, &frac);
    add_256(&tmp2, &frac3, op1y);
    sub_256(resY, &tmp1, &tmp2);
}

void scalar_mult_base_256_51(radix51_t *base_point_x, radix51_t *base_point_y, char *sk, radix51_t *Pkx, radix51_t *Pky)
{
    /*
    Pk = sk*base_point
    */
    int n = sk_length;
    radix51_t R0x, R0y, R1x, R1y, resX, resY;
    set_ui_256(&R0x, 0);
    set_ui_256(&R0y, 0);

    copy_256(&R1x, *base_point_x);
    copy_256(&R1y, *base_point_y);

    for (int i = 0; i < n; i++)
    {
        if (sk[i] == '0')
        {
            // R1 <- add (R0, R1)
            point_add_base_256_51(&resX, &resY, &R0x, &R0y, &R1x, &R1y);
            copy_256(&R1x, resX);
            copy_256(&R1y, resY);
            // R0 <- double(R0)
            point_double_base_256_51(&resX, &resY, &R0x, &R0y);
            copy_256(&R0x, resX);
            copy_256(&R0y, resY);
        }
        else
        {
            // R0 <- add(R0, R1)
            point_add_base_256_51(&resX, &resY, &R0x, &R0y, &R1x, &R1y);
            copy_256(&R0x, resX);
            copy_256(&R0y, resY);
            // R1 <- double(R1)
            point_double_base_256_51(&resX, &resY, &R1x, &R1y);
            copy_256(&R1x, resX);
            copy_256(&R1y, resY);
        }
    }
    // Assign Pk
    copy_256(Pkx, R0x);
    copy_256(Pky, R0y);
}

void public_key_gen_base_256_51(char *sk, radix51_t *x, radix51_t *y)
{
    /*
    (x,y) = sk*G
    */
    scalar_mult_base_256_51(&Gx_rad51, &Gy_rad51, sk, x, y);
}

int validate_base_256_51(char **sk_A, char **sk_B, int num_test_cases)
{
    /*
    validate by running an instance of DH
    */
    radix51_t Pkx_A, Pky_A, Pkx_B, Pky_B, ssx_A, ssy_A, ssx_B, ssy_B;
    set_ui_256(&Pkx_A, 0);
    set_ui_256(&Pky_A, 0);
    set_ui_256(&Pkx_B, 0);
    set_ui_256(&Pky_B, 0);
    set_ui_256(&ssx_A, 0);
    set_ui_256(&ssy_A, 0);
    set_ui_256(&ssx_B, 0);
    set_ui_256(&ssy_B, 0);

    int ok = 1;

    for (int i = 0; i < num_test_cases; ++i)
    {
        // 1. generate Alice's public key
        public_key_gen_base_256_51(sk_A[i], &Pkx_A, &Pky_A);
        // 2. generate Bob's public key
        public_key_gen_base_256_51(sk_B[i], &Pkx_B, &Pky_B);
        // 3. Generate the shared secret
        scalar_mult_base_256_51(&Pkx_B, &Pky_B, sk_A[i], &ssx_A, &ssy_A);
        scalar_mult_base_256_51(&Pkx_A, &Pky_A, sk_B[i], &ssx_B, &ssy_B);
        if (!cmp_coord_256(&ssx_A, &ssx_B))
        {
            printf("     Validation failed for key '%s' found in row '%d'\n", sk_A[i], i + 1);
            ok = 0;
        }
    }
    return ok;
}

int checkPk_base_256_51(char **sk, char **pk, int num_test_cases)
{
    /*
    Generate PKs and compare them
    */
    radix51_t Pkx, Pky;
    set_ui_256(&Pkx, 0);
    set_ui_256(&Pky, 0);

    int ok = 1;

    char pkAsString[513];

    for (int i = 0; i < num_test_cases; ++i)
    {
        public_key_gen_base_256_51(sk[i], &Pkx, &Pky);
        radix51_toString(pkAsString, Pkx);
        if (strncmp(pk[i], pkAsString, 513))
        {
            printf("     Public key generation failed for key found in row %d (computed key: '%s', expected key: '%s')\n", i + 1, pkAsString, pk[i]);
            ok = 0;
        }
    }
    return ok;
}
