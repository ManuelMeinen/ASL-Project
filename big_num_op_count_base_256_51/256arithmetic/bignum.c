#include <stdio.h>
#include "../../_common_libs/radix51.h"

int bignum_ops = 0;
void copy_256(radix51_t *dest, radix51_t src){
    /*
    set dest to src
    */

    for (size_t i = 0; i < 5; i++)
    {
        dest->d[i] = src.d[i];
    }
    
}
void reduce_2_to_2(radix51_t *x){
    /*
    reduce x
    */
    signed long long carry;
    for (int i = 0; i < 4; i++)
    {
        carry = x->d[i] >> 51;
        x->d[i+1] += carry;
        carry <<= 51;
        x->d[i] -= carry;
    }
    
    carry = x->d[4] >> 51;
    x->d[0] += 19*carry;
    carry <<= 51;
    x->d[4] -= carry;  
}

void reduce_4_to_2(__int128 r[9]){
    /*
    reduce r
    */
    signed __int128 carry;
    for (int i = 0; i < 4; i++)
    {
        carry = r[i] >> 51;
        r[i+1] += carry;
        carry <<= 51;
        r[i] -= carry;
    }
    
    carry = r[4] >> 51;
    r[0] += 19*carry;
    carry <<= 51;
    r[4] -= carry;  
}

void add_256(radix51_t *res, radix51_t *a, radix51_t *b){
    /*
    res = a + b
    */
   bignum_ops += 1;
    radix51_t tmp_a, tmp_b; 
    copy_256(&tmp_a, *a);
    copy_256(&tmp_b, *b);
    
    for (int i = 0; i < 5; i++)
    {
        res->d[i] = tmp_a.d[i] + tmp_b.d[i];
    }
    reduce_2_to_2(res);
}


void sub_256(radix51_t *res, radix51_t *a, radix51_t *b){
    /*
    res = a - b
    */
   bignum_ops += 1;
    radix51_t tmp_a, tmp_b; 
    copy_256(&tmp_a, *a);
    copy_256(&tmp_b, *b);

    for (int i = 0; i < 5; i++)
    {
        res->d[i] = tmp_a.d[i] - tmp_b.d[i];
    }
    reduce_2_to_2(res);

}

void reduce_mul_512(__int128 r[9]){
    /*
    reduce r
    */
    r[0] += 19*r[5];
    r[1] += 19*r[6];
    r[2] += 19*r[7];
    r[3] += 19*r[8];

}

void mul_256(radix51_t *res, radix51_t *x, radix51_t *y){
    /*
    res = x * y
    */
    bignum_ops += 1;
    unsigned __int128 r[9];

    radix51_t tmp_x, tmp_y;
    copy_256(&tmp_x, *x);
    copy_256(&tmp_y, *y);

    const signed long long *a = tmp_x.d;
    const signed long long *b = tmp_y.d;

    
    r[0] = (__int128_t) a[0]*b[0];
    r[1] = (__int128_t) a[0]*b[1] + (__int128_t) a[1]*b[0];
    r[2] = (__int128_t) a[0]*b[2] + (__int128_t) a[1]*b[1] + (__int128_t) a[2]*b[0];
    r[3] = (__int128_t) a[0]*b[3] + (__int128_t) a[1]*b[2] + (__int128_t) a[2]*b[1] + (__int128_t) a[3]*b[0];
    r[4] = (__int128_t) a[0]*b[4] + (__int128_t) a[1]*b[3] + (__int128_t) a[2]*b[2] + (__int128_t) a[3]*b[1] + (__int128_t) a[4]*b[0];
    r[5] = (__int128_t) a[1]*b[4] + (__int128_t) a[2]*b[3] + (__int128_t) a[3]*b[2] + (__int128_t) a[4]*b[1];
    r[6] = (__int128_t) a[2]*b[4] + (__int128_t) a[3]*b[3] + (__int128_t) a[4]*b[2];
    r[7] = (__int128_t) a[3]*b[4] + (__int128_t) a[4]*b[3];
    r[8] = (__int128_t) a[4]*b[4];

    reduce_mul_512(r);

    reduce_4_to_2(r);
    
    res->d[0] = r[0];
    res->d[1] = r[1];
    res->d[2] = r[2];
    res->d[3] = r[3];
    res->d[4] = r[4];

    
    reduce_2_to_2(res);
}



int cmp_coord_256(radix51_t *x1, radix51_t *x2){
    /*
    return 1 if x1 == x2
    */ 
    return (x1->d[0] == x2->d[0]) && (x1->d[1] == x2->d[1]) && (x1->d[2] == x2->d[2]) && (x1->d[3] == x2->d[3]) && (x1->d[4] == x2->d[4]);
}

int cmp_point_256(radix51_t *x1, radix51_t *y1, radix51_t *x2, radix51_t *y2){
    /*
    return 1 if x1 == y1 and x2 == y2
    */
    return cmp_coord_256(x1, y1) && cmp_coord_256(x2, y2);
}

int iszero_256(radix51_t *x){
    /*
    return 1 if x == y == 0
    */
    return (x->d[0] == 0) && (x->d[1] == 0) && (x-> d[2] == 0) && (x->d[3] == 0) && (x->d[4] == 0);
}

void set_ui_256(radix51_t *x, signed int ui){
    /*
    set radix51_t to a unsigned int ui (lowest limb only, others 0)
    maximum 51 bit numbers !!!
    */
    x->d[0] = (signed long long) ui;
    x->d[1] = (signed long long) 0;
    x->d[2] = (signed long long) 0;
    x->d[3] = (signed long long) 0;
    x->d[4] = (signed long long) 0;
}

void square_256(radix51_t *res, radix51_t *x){
    bignum_ops += 1;
    /*
    res = x^2
    */
    __int128_t r[9];

    radix51_t tmp_x;
    copy_256(&tmp_x, *x);

    const signed long long *a = tmp_x.d;
    signed long long _2a[4]; 

    _2a[0] = a[0] << 1; 
    _2a[1] = a[1] << 1; 
    _2a[2] = a[2] << 1; 
    _2a[3] = a[3] << 1;

    r[0] = (__int128_t) a[0]*a[0];
    r[1] = (__int128_t) _2a[0]*a[1];
    r[2] = (__int128_t) _2a[0]*a[2] + (__int128_t) a[1]*a[1];
    r[3] = (__int128_t) _2a[0]*a[3] + (__int128_t) _2a[1]*a[2];
    r[4] = (__int128_t) _2a[0]*a[4] + (__int128_t) _2a[1]*a[3] + (__int128_t) a[2]*a[2]; 
    r[5] = (__int128_t) _2a[1]*a[4] + (__int128_t) _2a[2]*a[3];
    r[6] = (__int128_t) _2a[2]*a[4] + (__int128_t) a[3]*a[3];
    r[7] = (__int128_t) _2a[3]*a[4];
    r[8] = (__int128_t) a[4]*a[4];

    reduce_mul_512(r);

    reduce_4_to_2(r);
    
    res->d[0] = r[0];
    res->d[1] = r[1];
    res->d[2] = r[2];
    res->d[3] = r[3];
    res->d[4] = r[4];

    reduce_2_to_2(res);
}

void inverse_256(radix51_t *res, radix51_t *x){
    /*
    res = x^-1 (multiplicative inverse modulo p)
    */
    radix51_t z2, z9, z11, z2_5_0, z2_10_0, z2_20_0, z2_50_0, z2_100_0, t;
    int i;
    /* 2 */ square_256(&z2,x);
    /* 4 */ square_256(&t,&z2);
    /* 8 */ square_256(&t,&t);
    /* 9 */ mul_256(&z9,&t,x);
    /* 11 */ mul_256(&z11,&z9,&z2);
    /* 22 */ square_256(&t,&z11);
    /* 2^5 - 2^0 = 31 */ mul_256(&z2_5_0,&t,&z9);
    /* 2^6 - 2^1 */ square_256(&t,&z2_5_0);
    /* 2^10 - 2^5 */ for (i = 1;i < 5;i++) { square_256(&t,&t); }
    /* 2^10 - 2^0 */ mul_256(&z2_10_0,&t,&z2_5_0);
    /* 2^11 - 2^1 */ square_256(&t,&z2_10_0);
    /* 2^20 - 2^10 */ for (i = 1;i < 10;i++) { square_256(&t,&t); }
    /* 2^20 - 2^0 */ mul_256(&z2_20_0,&t,&z2_10_0);
    /* 2^21 - 2^1 */ square_256(&t,&z2_20_0);
    /* 2^40 - 2^20 */ for (i = 1;i < 20;i++) { square_256(&t,&t); }
    /* 2^40 - 2^0 */ mul_256(&t,&t,&z2_20_0);
    /* 2^41 - 2^1 */ square_256(&t,&t);
    /* 2^50 - 2^10 */ for (i = 1;i < 10;i++) { square_256(&t,&t); }
    /* 2^50 - 2^0 */ mul_256(&z2_50_0,&t,&z2_10_0);
    /* 2^51 - 2^1 */ square_256(&t,&z2_50_0);
    /* 2^100 - 2^50 */ for (i = 1;i < 50;i++) { square_256(&t,&t); }
    /* 2^100 - 2^0 */ mul_256(&z2_100_0,&t,&z2_50_0);
    /* 2^101 - 2^1 */ square_256(&t,&z2_100_0);
    /* 2^200 - 2^100 */ for (i = 1;i < 100;i++) { square_256(&t,&t); }
    /* 2^200 - 2^0 */ mul_256(&t,&t,&z2_100_0);
    /* 2^201 - 2^1 */ square_256(&t,&t);
    /* 2^250 - 2^50 */ for (i = 1;i < 50;i++) { square_256(&t,&t); }
    /* 2^250 - 2^0 */ mul_256(&t,&t,&z2_50_0);
    /* 2^251 - 2^1 */ square_256(&t,&t);
    /* 2^252 - 2^2 */ square_256(&t,&t);
    /* 2^253 - 2^3 */ square_256(&t,&t);
    /* 2^254 - 2^4 */ square_256(&t,&t);
    /* 2^255 - 2^5 */ square_256(&t,&t);
    /* 2^255 - 21 */ mul_256(res,&t,&z11);
}


int get_bignum_ops(){
    return bignum_ops;
}

void div_256(radix51_t *res, radix51_t *a, radix51_t *b){
    /*
    res = a/b
    */
    // inverse and mult
    radix51_t inv;
    inverse_256(&inv, b);
    mul_256(res, a, &inv);
}