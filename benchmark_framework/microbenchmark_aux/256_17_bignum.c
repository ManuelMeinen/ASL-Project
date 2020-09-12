#include <stdio.h>
#include "256_17_bignum.h"
#include "../../_common_libs/radix17.h"

 const int NO_LIMBS = 15;
 const int LIMB_SIZE = 17;


void copy_256_17(radix17_t *dest, radix17_t src){
    //set dest to src

    for (size_t i = 0; i < NO_LIMBS; i++)
    {
        dest->d[i] = src.d[i];
    }
    
}
void reduce_256_17(radix17_t *x){

    signed long carry;
    for (int i = 0; i < NO_LIMBS-1; i++)
    {
        carry = x->d[i] >> LIMB_SIZE;
        x->d[i+1] += carry;
        carry <<= LIMB_SIZE;
        x->d[i] -= carry;
    }
    
    carry = x->d[NO_LIMBS - 1] >> LIMB_SIZE;
    x->d[0] += 19*carry;
    carry <<= LIMB_SIZE;
    x->d[NO_LIMBS - 1] -= carry;  

}

void reduce_mul_256_17(unsigned long r[29]){

    signed long carry;
    for (int i = 0; i < NO_LIMBS-1; i++)
    {
        carry = r[i] >> LIMB_SIZE;
        r[i+1] += carry;
        carry <<= LIMB_SIZE;
        r[i] -= carry;
    }
    
    carry = r[NO_LIMBS - 1] >> LIMB_SIZE;
    r[0] += 19*carry;
    carry <<= LIMB_SIZE;
    r[NO_LIMBS - 1] -= carry;  

}

void add_256_17(radix17_t *res, radix17_t *a, radix17_t *b){
	
    radix17_t tmp_a, tmp_b; 
    copy_256_17(&tmp_a, *a);
    copy_256_17(&tmp_b, *b);
    
    for (int i = 0; i < NO_LIMBS; i++)
    {
        res->d[i] = tmp_a.d[i] + tmp_b.d[i];
    }
    reduce_256_17(res);
}

void add_256_17_no_reduce(radix17_t *res, radix17_t *a, radix17_t *b){
	
    radix17_t tmp_a, tmp_b; 
    copy_256_17(&tmp_a, *a);
    copy_256_17(&tmp_b, *b);
    
    for (int i = 0; i < NO_LIMBS; i++)
    {
        res->d[i] = tmp_a.d[i] + tmp_b.d[i];
    }
}


void sub_256_17(radix17_t *res, radix17_t *a, radix17_t *b){

    radix17_t tmp_a, tmp_b; 
    copy_256_17(&tmp_a, *a);
    copy_256_17(&tmp_b, *b);

    for (int i = 0; i < NO_LIMBS; i++)
    {
        res->d[i] = tmp_a.d[i] - tmp_b.d[i];
    }
    reduce_256_17(res);
    reduce_256_17(res);
}

void sub_256_17_no_reduce(radix17_t *res, radix17_t *a, radix17_t *b){

    radix17_t tmp_a, tmp_b; 
    copy_256_17(&tmp_a, *a);
    copy_256_17(&tmp_b, *b);

    for (int i = 0; i < NO_LIMBS; i++)
    {
        res->d[i] = tmp_a.d[i] - tmp_b.d[i];
    }
}

void reduce_mul_512_17(unsigned long r[29]){

    for(int i = 0; i<NO_LIMBS-1; i++){
        r[i] += 19*r[i+NO_LIMBS];
    }
}

void mul_256_17(radix17_t *res, radix17_t *x, radix17_t *y){

    unsigned long r[NO_LIMBS*2 - 1];

    radix17_t tmp_x, tmp_y;
    copy_256_17(&tmp_x, *x);
    copy_256_17(&tmp_y, *y);

    const signed int *a = tmp_x.d;
    const signed int *b = tmp_y.d;

    
    r[0]  =  (unsigned long) a[0] * b[0];
    r[1]  =  (unsigned long) a[0] * b[1]   + (unsigned long) a[1] * b[0];
    r[2]  =  (unsigned long) a[0] * b[2]   + (unsigned long) a[1] * b[1]   + (unsigned long) a[2] * b[0];
    r[3]  =  (unsigned long) a[0] * b[3]   + (unsigned long) a[1] * b[2]   + (unsigned long) a[2] * b[1]   + (unsigned long) a[3] * b[0];
    r[4]  =  (unsigned long) a[0] * b[4]   + (unsigned long) a[1] * b[3]   + (unsigned long) a[2] * b[2]   + (unsigned long) a[3] * b[1]   + (unsigned long) a[4] * b[0];
    r[5]  =  (unsigned long) a[0] * b[5]   + (unsigned long) a[1] * b[4]   + (unsigned long) a[2] * b[3]   + (unsigned long) a[3] * b[2]   + (unsigned long) a[4] * b[1]   + (unsigned long) a[5] * b[0];
    r[6]  =  (unsigned long) a[0] * b[6]   + (unsigned long) a[1] * b[5]   + (unsigned long) a[2] * b[4]   + (unsigned long) a[3] * b[3]   + (unsigned long) a[4] * b[2]   + (unsigned long) a[5] * b[1]  + (unsigned long) a[6] * b[0];
    r[7]  =  (unsigned long) a[0] * b[7]   + (unsigned long) a[1] * b[6]   + (unsigned long) a[2] * b[5]   + (unsigned long) a[3] * b[4]   + (unsigned long) a[4] * b[3]   + (unsigned long) a[5] * b[2]  + (unsigned long) a[6] * b[1]  + (unsigned long) a[7] * b[0];
    r[8]  =  (unsigned long) a[0] * b[8]   + (unsigned long) a[1] * b[7]   + (unsigned long) a[2] * b[6]   + (unsigned long) a[3] * b[5]   + (unsigned long) a[4] * b[4]   + (unsigned long) a[5] * b[3]  + (unsigned long) a[6] * b[2]  + (unsigned long) a[7] * b[1]  + (unsigned long) a[8] * b[0];
    r[9]  =  (unsigned long) a[0] * b[9]   + (unsigned long) a[1] * b[8]   + (unsigned long) a[2] * b[7]   + (unsigned long) a[3] * b[6]   + (unsigned long) a[4] * b[5]   + (unsigned long) a[5] * b[4]  + (unsigned long) a[6] * b[3]  + (unsigned long) a[7] * b[2]  + (unsigned long) a[8] * b[1]  + (unsigned long) a[9] * b[0];
    r[10] =  (unsigned long) a[0] * b[10]  + (unsigned long) a[1] * b[9]   + (unsigned long) a[2] * b[8]   + (unsigned long) a[3] * b[7]   + (unsigned long) a[4] * b[6]   + (unsigned long) a[5] * b[5]  + (unsigned long) a[6] * b[4]  + (unsigned long) a[7] * b[3]  + (unsigned long) a[8] * b[2]  + (unsigned long) a[9] * b[1]  + (unsigned long) a[10] * b[0];
    r[11] =  (unsigned long) a[0] * b[11]  + (unsigned long) a[1] * b[10]  + (unsigned long) a[2] * b[9]   + (unsigned long) a[3] * b[8]   + (unsigned long) a[4] * b[7]   + (unsigned long) a[5] * b[6]  + (unsigned long) a[6] * b[5]  + (unsigned long) a[7] * b[4]  + (unsigned long) a[8] * b[3]  + (unsigned long) a[9] * b[2]  + (unsigned long) a[10] * b[1] + (unsigned long) a[11] * b[0];
    r[12] =  (unsigned long) a[0] * b[12]  + (unsigned long) a[1] * b[11]  + (unsigned long) a[2] * b[10]  + (unsigned long) a[3] * b[9]   + (unsigned long) a[4] * b[8]   + (unsigned long) a[5] * b[7]  + (unsigned long) a[6] * b[6]  + (unsigned long) a[7] * b[5]  + (unsigned long) a[8] * b[4]  + (unsigned long) a[9] * b[3]  + (unsigned long) a[10] * b[2] + (unsigned long) a[11] * b[1] + (unsigned long) a[12] * b[0];
    r[13] =  (unsigned long) a[0] * b[13]  + (unsigned long) a[1] * b[12]  + (unsigned long) a[2] * b[11]  + (unsigned long) a[3] * b[10]  + (unsigned long) a[4] * b[9]   + (unsigned long) a[5] * b[8]  + (unsigned long) a[6] * b[7]  + (unsigned long) a[7] * b[6]  + (unsigned long) a[8] * b[5]  + (unsigned long) a[9] * b[4]  + (unsigned long) a[10] * b[3] + (unsigned long) a[11] * b[2] + (unsigned long) a[12] * b[1] + (unsigned long) a[13] * b[0];
    r[14] =  (unsigned long) a[0] * b[14]  + (unsigned long) a[1] * b[13]  + (unsigned long) a[2] * b[12]  + (unsigned long) a[3] * b[11]  + (unsigned long) a[4] * b[10]  + (unsigned long) a[5] * b[9]  + (unsigned long) a[6] * b[8]  + (unsigned long) a[7] * b[7]  + (unsigned long) a[8] * b[6]  + (unsigned long) a[9] * b[5]  + (unsigned long) a[10] * b[4] + (unsigned long) a[11] * b[3] + (unsigned long) a[12] * b[2] + (unsigned long) a[13] * b[1] + (unsigned long) a[14] * b[0];
    r[15] =  (unsigned long) a[1] * b[14]  + (unsigned long) a[2] * b[13]  + (unsigned long) a[3] * b[12]  + (unsigned long) a[4] * b[11]  + (unsigned long) a[5] * b[10]  + (unsigned long) a[6] * b[9]  + (unsigned long) a[7] * b[8]  + (unsigned long) a[8] * b[7]  + (unsigned long) a[9] * b[6]  + (unsigned long) a[10] * b[5] + (unsigned long) a[11] * b[4] + (unsigned long) a[12] * b[3] + (unsigned long) a[13] * b[2] + (unsigned long) a[14] * b[1];
    r[16] =  (unsigned long) a[2] * b[14]  + (unsigned long) a[3] * b[13]  + (unsigned long) a[4] * b[12]  + (unsigned long) a[5] * b[11]  + (unsigned long) a[6] * b[10]  + (unsigned long) a[7] * b[9]  + (unsigned long) a[8] * b[8]  + (unsigned long) a[9] * b[7]  + (unsigned long) a[10] * b[6] + (unsigned long) a[11] * b[5] + (unsigned long) a[12] * b[4] + (unsigned long) a[13] * b[3] + (unsigned long) a[14] * b[2];
    r[17] =  (unsigned long) a[3] * b[14]  + (unsigned long) a[4] * b[13]  + (unsigned long) a[5] * b[12]  + (unsigned long) a[6] * b[11]  + (unsigned long) a[7] * b[10]  + (unsigned long) a[8] * b[9]  + (unsigned long) a[9] * b[8]  + (unsigned long) a[10] * b[7] + (unsigned long) a[11] * b[6] + (unsigned long) a[12] * b[5] + (unsigned long) a[13] * b[4] + (unsigned long) a[14] * b[3];
    r[18] =  (unsigned long) a[4] * b[14]  + (unsigned long) a[5] * b[13]  + (unsigned long) a[6] * b[12]  + (unsigned long) a[7] * b[11]  + (unsigned long) a[8] * b[10]  + (unsigned long) a[9] * b[9]  + (unsigned long) a[10] * b[8] + (unsigned long) a[11] * b[7] + (unsigned long) a[12] * b[6] + (unsigned long) a[13] * b[5] + (unsigned long) a[14] * b[4];
    r[19] =  (unsigned long) a[5] * b[14]  + (unsigned long) a[6] * b[13]  + (unsigned long) a[7] * b[12]  + (unsigned long) a[8] * b[11]  + (unsigned long) a[9] * b[10]  + (unsigned long) a[10] * b[9] + (unsigned long) a[11] * b[8] + (unsigned long) a[12] * b[7] + (unsigned long) a[13] * b[6] + (unsigned long) a[14] * b[5];
    r[20] =  (unsigned long) a[6] * b[14]  + (unsigned long) a[7] * b[13]  + (unsigned long) a[8] * b[12]  + (unsigned long) a[9] * b[11]  + (unsigned long) a[10] * b[10] + (unsigned long) a[11] * b[9] + (unsigned long) a[12] * b[8] + (unsigned long) a[13] * b[7] + (unsigned long) a[14] * b[6];
    r[21] =  (unsigned long) a[7] * b[14]  + (unsigned long) a[8] * b[13]  + (unsigned long) a[9] * b[12]  + (unsigned long) a[10] * b[11] + (unsigned long) a[11] * b[10] + (unsigned long) a[12] * b[9] + (unsigned long) a[13] * b[8] + (unsigned long) a[14] * b[7];
    r[22] =  (unsigned long) a[8] * b[14]  + (unsigned long) a[9] * b[13]  + (unsigned long) a[10] * b[12] + (unsigned long) a[11] * b[11] + (unsigned long) a[12] * b[10] + (unsigned long) a[13] * b[9] + (unsigned long) a[14] * b[8];
    r[23] =  (unsigned long) a[9] * b[14]  + (unsigned long) a[10] * b[13] + (unsigned long) a[11] * b[12] + (unsigned long) a[12] * b[11] + (unsigned long) a[13] * b[10] + (unsigned long) a[14] * b[9];
    r[24] =  (unsigned long) a[10] * b[14] + (unsigned long) a[11] * b[13] + (unsigned long) a[12] * b[12] + (unsigned long) a[13] * b[11] + (unsigned long) a[14] * b[10];
    r[25] =  (unsigned long) a[11] * b[14] + (unsigned long) a[12] * b[13] + (unsigned long) a[13] * b[12] + (unsigned long) a[14] * b[11];
    r[26] =  (unsigned long) a[12] * b[14] + (unsigned long) a[13] * b[13] + (unsigned long) a[14] * b[12];
    r[27] =  (unsigned long) a[13] * b[14] + (unsigned long) a[14] * b[13];
    r[28] =  (unsigned long) a[14] * b[14];

    reduce_mul_512_17(r);

    reduce_mul_256_17(r);
    
    for(int i = 0; i<NO_LIMBS; i++){
        res->d[i] = r[i];
    }
    reduce_256_17(res);
}



int cmp_coord_256_17(radix17_t *x1, radix17_t *x2){
    //return 1 if x1 == x2 

    int bool;
    bool = (x1->d[0] == x2->d[0]);

    for (size_t i = 1; i < NO_LIMBS; i++)
    {
        bool = bool && (x1->d[i] == x2->d[i]);
    }
    return bool;

}
int cmp_point_256_17(radix17_t *x1, radix17_t *y1, radix17_t *x2, radix17_t *y2){
    //return 1 if x1 == y1 and x2 == y2

    return cmp_coord_256_17(x1, y1) && cmp_coord_256_17(x2, y2);

}
int iszero_256_17(radix17_t *x){
    //return 1 if x == y == 0

    int bool;
    bool = (x->d[0] == 0);

    for (size_t i = 1; i < NO_LIMBS; i++)
    {
        bool = bool && (x->d[i] == 0);
    }
    return bool;
    
}
void set_ui_256_17(radix17_t *x, signed int ui){
    //set radix17_t to a unsigned int ui (lowest limb only, others 0)
    // maximum 17 bit numbers !!!

    x->d[0] = (signed int) ui;

    for (size_t i = 1; i < NO_LIMBS; i++)
    {
        x->d[i] = (signed int) 0;
    }
    
}


void square_256_17(radix17_t *res, radix17_t *x){
	
    unsigned long r[NO_LIMBS*2 - 1];

    radix17_t tmp_x;
    copy_256_17(&tmp_x, *x);

    const signed int *a = tmp_x.d;
    signed long _2a[NO_LIMBS*2 - 1]; 

    _2a[0] = a[0] << 1;
    _2a[1] = a[1] << 1;
    _2a[2] = a[2] << 1;
    _2a[3] = a[3] << 1;
    _2a[4] = a[4] << 1;
    _2a[5] = a[5] << 1;
    _2a[6] = a[6] << 1;
    _2a[7] = a[7] << 1;
    _2a[8] = a[8] << 1;
    _2a[9] = a[9] << 1;
    _2a[10] = a[10] << 1;
    _2a[11] = a[11] << 1;
    _2a[12] = a[12] << 1;
    _2a[13] = a[13] << 1;

    r[0]  = (unsigned long)  a[0]*a[0];
    r[1]  = (unsigned long)  _2a[0]*a[1];
    r[2]  = (unsigned long) _2a[0]*a[2]   + (unsigned long) a[1]*a[1];
    r[3]  = (unsigned long) _2a[0]*a[3]   + (unsigned long) _2a[1]*a[2];
    r[4]  = (unsigned long) _2a[0]*a[4]   + (unsigned long) _2a[1]*a[3]   + (unsigned long) a[2]*a[2];
    r[5]  = (unsigned long) _2a[0]*a[5]   + (unsigned long) _2a[1]*a[4]   + (unsigned long) _2a[2]*a[3];
    r[6]  = (unsigned long) _2a[0]*a[6]   + (unsigned long) _2a[1]*a[5]   + (unsigned long) _2a[2]*a[4]   + (unsigned long) a[3]*a[3];
    r[7]  = (unsigned long) _2a[0]*a[7]   + (unsigned long) _2a[1]*a[6]   + (unsigned long) _2a[2]*a[5]   + (unsigned long) _2a[3]*a[4];
    r[8]  = (unsigned long) _2a[0]*a[8]   + (unsigned long) _2a[1]*a[7]   + (unsigned long) _2a[2]*a[6]   + (unsigned long) _2a[3]*a[5]  + (unsigned long) a[4]*a[4];
    r[9]  = (unsigned long) _2a[0]*a[9]   + (unsigned long) _2a[1]*a[8]   + (unsigned long) _2a[2]*a[7]   + (unsigned long) _2a[3]*a[6]  + (unsigned long) _2a[4]*a[5];
    r[10] = (unsigned long) _2a[0]*a[10]  + (unsigned long) _2a[1]*a[9]   + (unsigned long) _2a[2]*a[8]   + (unsigned long) _2a[3]*a[7]  + (unsigned long) _2a[4]*a[6]  + (unsigned long) a[5]*a[5];
    r[11] = (unsigned long) _2a[0]*a[11]  + (unsigned long) _2a[1]*a[10]  + (unsigned long) _2a[2]*a[9]   + (unsigned long) _2a[3]*a[8]  + (unsigned long) _2a[4]*a[7]  + (unsigned long) _2a[5]*a[6];
    r[12] = (unsigned long) _2a[0]*a[12]  + (unsigned long) _2a[1]*a[11]  + (unsigned long) _2a[2]*a[10]  + (unsigned long) _2a[3]*a[9]  + (unsigned long) _2a[4]*a[8]  + (unsigned long) _2a[5]*a[7] + (unsigned long) a[6]*a[6];
    r[13] = (unsigned long) _2a[0]*a[13]  + (unsigned long) _2a[1]*a[12]  + (unsigned long) _2a[2]*a[11]  + (unsigned long) _2a[3]*a[10] + (unsigned long) _2a[4]*a[9]  + (unsigned long) _2a[5]*a[8] + (unsigned long) _2a[6]*a[7];
    r[14] = (unsigned long) _2a[0]*a[14]  + (unsigned long) _2a[1]*a[13]  + (unsigned long) _2a[2]*a[12]  + (unsigned long) _2a[3]*a[11] + (unsigned long) _2a[4]*a[10] + (unsigned long) _2a[5]*a[9] + (unsigned long) _2a[6]*a[8] + (unsigned long) a[7]*a[7];
    r[15] = (unsigned long) _2a[1]*a[14]  + (unsigned long) _2a[2]*a[13]  + (unsigned long) _2a[3]*a[12]  + (unsigned long) _2a[4]*a[11] + (unsigned long) _2a[5]*a[10] + (unsigned long) _2a[6]*a[9] + (unsigned long) _2a[7]*a[8];
    r[16] = (unsigned long) _2a[2]*a[14]  + (unsigned long) _2a[3]*a[13]  + (unsigned long) _2a[4]*a[12]  + (unsigned long) _2a[5]*a[11] + (unsigned long) _2a[6]*a[10] + (unsigned long) _2a[7]*a[9] + (unsigned long) a[8]*a[8];
    r[17] = (unsigned long) _2a[3]*a[14]  + (unsigned long) _2a[4]*a[13]  + (unsigned long) _2a[5]*a[12]  + (unsigned long) _2a[6]*a[11] + (unsigned long) _2a[7]*a[10] + (unsigned long) _2a[8]*a[9];
    r[18] = (unsigned long) _2a[4]*a[14]  + (unsigned long) _2a[5]*a[13]  + (unsigned long) _2a[6]*a[12]  + (unsigned long) _2a[7]*a[11] + (unsigned long) _2a[8]*a[10] + (unsigned long) a[9]*a[9];
    r[19] = (unsigned long) _2a[5]*a[14]  + (unsigned long) _2a[6]*a[13]  + (unsigned long) _2a[7]*a[12]  + (unsigned long) _2a[8]*a[11] + (unsigned long) _2a[9]*a[10];
    r[20] = (unsigned long) _2a[6]*a[14]  + (unsigned long) _2a[7]*a[13]  + (unsigned long) _2a[8]*a[12]  + (unsigned long) _2a[9]*a[11] + (unsigned long) a[10]*a[10];
    r[21] = (unsigned long) _2a[7]*a[14]  + (unsigned long) _2a[8]*a[13]  + (unsigned long) _2a[9]*a[12]  + (unsigned long) _2a[10]*a[11];
    r[22] = (unsigned long) _2a[8]*a[14]  + (unsigned long) _2a[9]*a[13]  + (unsigned long) _2a[10]*a[12] + (unsigned long) a[11]*a[11];
    r[23] = (unsigned long) _2a[9]*a[14]  + (unsigned long) _2a[10]*a[13] + (unsigned long) _2a[11]*a[12];
    r[24] = (unsigned long) _2a[10]*a[14] + (unsigned long) _2a[11]*a[13] + (unsigned long) a[12]*a[12];
    r[25] = (unsigned long) _2a[11]*a[14] + (unsigned long) _2a[12]*a[13];
    r[26] = (unsigned long) _2a[12]*a[14] + (unsigned long) a[13]*a[13];
    r[27] = (unsigned long) _2a[13]*a[14];
    r[28] = (unsigned long) a[14]*a[14];

    reduce_mul_512_17(r);

    reduce_mul_256_17(r);
    
     for(int i = 0; i<NO_LIMBS; i++){
        res->d[i] = r[i];
    }

    reduce_256_17(res);
}
void inverse_256_17(radix17_t *res, radix17_t *x){
    radix17_t z2, z9, z11, z2_5_0, z2_10_0, z2_20_0, z2_50_0, z2_100_0, t;
    int i;
    /* 2 */ square_256_17(&z2,x);
    /* 4 */ square_256_17(&t,&z2);
    /* 8 */ square_256_17(&t,&t);
    /* 9 */ mul_256_17(&z9,&t,x);
    /* 11 */ mul_256_17(&z11,&z9,&z2);
    /* 22 */ square_256_17(&t,&z11);
    /* 2^5 - 2^0 = 31 */ mul_256_17(&z2_5_0,&t,&z9);
    /* 2^6 - 2^1 */ square_256_17(&t,&z2_5_0);
    /* 2^10 - 2^5 */ for (i = 1;i < 5;i++) { square_256_17(&t,&t); }
    /* 2^10 - 2^0 */ mul_256_17(&z2_10_0,&t,&z2_5_0);
    /* 2^11 - 2^1 */ square_256_17(&t,&z2_10_0);
    /* 2^20 - 2^10 */ for (i = 1;i < 10;i++) { square_256_17(&t,&t); }
    /* 2^20 - 2^0 */ mul_256_17(&z2_20_0,&t,&z2_10_0);
    /* 2^21 - 2^1 */ square_256_17(&t,&z2_20_0);
    /* 2^40 - 2^20 */ for (i = 1;i < 20;i++) { square_256_17(&t,&t); }
    /* 2^40 - 2^0 */ mul_256_17(&t,&t,&z2_20_0);
    /* 2^41 - 2^1 */ square_256_17(&t,&t);
    /* 2^50 - 2^10 */ for (i = 1;i < 10;i++) { square_256_17(&t,&t); }
    /* 2^50 - 2^0 */ mul_256_17(&z2_50_0,&t,&z2_10_0);
    /* 2^51 - 2^1 */ square_256_17(&t,&z2_50_0);
    /* 2^100 - 2^50 */ for (i = 1;i < 50;i++) { square_256_17(&t,&t); }
    /* 2^100 - 2^0 */ mul_256_17(&z2_100_0,&t,&z2_50_0);
    /* 2^101 - 2^1 */ square_256_17(&t,&z2_100_0);
    /* 2^200 - 2^100 */ for (i = 1;i < 100;i++) { square_256_17(&t,&t); }
    /* 2^200 - 2^0 */ mul_256_17(&t,&t,&z2_100_0);
    /* 2^201 - 2^1 */ square_256_17(&t,&t);
    /* 2^250 - 2^50 */ for (i = 1;i < 50;i++) { square_256_17(&t,&t); }
    /* 2^250 - 2^0 */ mul_256_17(&t,&t,&z2_50_0);
    /* 2^251 - 2^1 */ square_256_17(&t,&t);
    /* 2^252 - 2^2 */ square_256_17(&t,&t);
    /* 2^253 - 2^3 */ square_256_17(&t,&t);
    /* 2^254 - 2^4 */ square_256_17(&t,&t);
    /* 2^255 - 2^5 */ square_256_17(&t,&t);
    /* 2^255 - 21 */ mul_256_17(res,&t,&z11);
    //print_radix17( res, "inv" );
}

void div_256_17(radix17_t *res, radix17_t *a, radix17_t *b){
    // inverse and mult
    radix17_t inv;
    inverse_256_17(&inv, b);
    mul_256_17(res, a, &inv);
}