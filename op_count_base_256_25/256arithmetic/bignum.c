#include <stdio.h>
#include "bignum.h"
#include "../../_common_libs/radix25.h"

 const int NO_LIMBS = 10;
 const int LIMB_SIZE = 25;

int total_ops = 0;
void copy_256(radix25_t *dest, radix25_t src){
    //set dest to src

    for (size_t i = 0; i < NO_LIMBS; i++)
    {
        dest->d[i] = src.d[i];
    }
    
}
void reduce_256(radix25_t *x){

    signed long carry;
    for (int i = 0; i < NO_LIMBS-2; i+=2)
    {
        carry = x->d[i] >> LIMB_SIZE+1;  // 26
        total_ops += 1;
        x->d[i+1] += carry;
        total_ops += 1;
        carry <<= LIMB_SIZE+1;
        total_ops += 1;
        x->d[i] -= carry;
        total_ops += 1;

        carry = x->d[i+1] >> LIMB_SIZE; // 25
        total_ops += 1;
        x->d[i+2] += carry;
        total_ops += 1;
        carry <<= LIMB_SIZE;
        total_ops += 1;
        x->d[i+1] -= carry;
        total_ops += 1;
    }
    
    carry = x->d[8] >> LIMB_SIZE+1;  // 26
    total_ops += 1;
    x->d[9] += carry;
    total_ops += 1;
    carry <<= LIMB_SIZE+1;
    total_ops += 1;
    x->d[8] -= carry;
    total_ops += 1;

    carry = x->d[NO_LIMBS - 1] >> LIMB_SIZE; //25
    total_ops += 1;
    x->d[0] += 19*carry;
    total_ops += 1;
    carry <<= LIMB_SIZE;
    total_ops += 1;
    x->d[NO_LIMBS - 1] -= carry;
    total_ops += 1;  

}

void reduce_mul_256(unsigned long r[29]){

    signed long carry;
    for (int i = 0; i < NO_LIMBS-2; i+=2)
    {
        carry = r[i] >> LIMB_SIZE+1; //26
        total_ops += 1;
        r[i+1] += carry;
        total_ops += 1;
        carry <<= LIMB_SIZE+1;
        total_ops += 1;
        r[i] -= carry;
        total_ops += 1;

        carry = r[i+1] >> LIMB_SIZE; //25
        total_ops += 1;
        r[i+2] += carry;
        total_ops += 1;
        carry <<= LIMB_SIZE;
        total_ops += 1;
        r[i+1] -= carry;
        total_ops += 1;
    }
    
    carry = r[8] >> LIMB_SIZE+1; //26
    total_ops += 1;
    r[9] += carry;
    total_ops += 1;
    carry <<= LIMB_SIZE+1;
    total_ops += 1;
    r[8] -= carry;
    total_ops += 1;

    carry = r[NO_LIMBS - 1] >> LIMB_SIZE; //25
    total_ops += 1;
    r[0] += 19*carry;
    total_ops += 1;
    carry <<= LIMB_SIZE;
    total_ops += 1;
    r[NO_LIMBS - 1] -= carry;
    total_ops += 1;  

}

void add_256(radix25_t *res, radix25_t *a, radix25_t *b){

    radix25_t tmp_a, tmp_b; 
    copy_256(&tmp_a, *a);
    copy_256(&tmp_b, *b);
    
    for (int i = 0; i < NO_LIMBS; i++)
    {
        res->d[i] = tmp_a.d[i] + tmp_b.d[i];
        total_ops += 1;
    }
    reduce_256(res);
    
}


void sub_256(radix25_t *res, radix25_t *a, radix25_t *b){

    radix25_t tmp_a, tmp_b; 
    copy_256(&tmp_a, *a);
    copy_256(&tmp_b, *b);

    for (int i = 0; i < NO_LIMBS; i++)
    {
        res->d[i] = tmp_a.d[i] - tmp_b.d[i];
        total_ops += 1;
    }

    reduce_256(res);
    //reduce_256(res);
    
}
void reduce_mul_512(unsigned long r[29]){

    for(int i = 0; i<NO_LIMBS-1; i++){
        r[i] += 19*r[i+NO_LIMBS];
        total_ops += 1;
    }
}

void mul_256(radix25_t *res, radix25_t *x, radix25_t *y){

    unsigned long r[NO_LIMBS*2 - 1];

    radix25_t tmp_x, tmp_y;
    copy_256(&tmp_x, *x);
    copy_256(&tmp_y, *y);

    const signed int *a = tmp_x.d;
    const signed int *b = tmp_y.d;

    // long m = 38*a[5];
    // printf("%d * %d = %ld \n", 38, a[5],  m);
    // printf("%ld * %ld = %ld \n",m , b[5], m*b[5] );
    r[0]  =  (unsigned long) a[0] * b[0] +   (unsigned long) 38*a[1]*b[9] + (unsigned long) 19*a[2]*b[8] +  (unsigned long) 38*a[3]*b[7] + (unsigned long) 19*a[4]*b[6] + (unsigned long) 38*a[5]*b[5] + (unsigned long) 19*a[6]*b[4] + (unsigned long) 38 * a[7] * b[3] + (unsigned long) 19 * a[8] * b[2] + (unsigned long) 38 * a[9] *b[1];
    r[1]  =  (unsigned long) a[0] * b[1]   + (unsigned long) a[1] * b[0] +  (unsigned long) 19 * a[2] * b[9] + (unsigned long) 19 * a[3] * b[8] + (unsigned long) 19 * a[4] * b[7] + (unsigned long) 19 * a[5] * b[6] +  (unsigned long) 19 * a[6] * b[5] + (unsigned long) 19 * a[7] * b[4] + (unsigned long) 19 * a[8] * b[3] + (unsigned long) 19 * a[9] * b[2];
    r[2]  =  (unsigned long) a[0] * b[2]   + (unsigned long) 2 * a[1] * b[1]   + (unsigned long) a[2] * b[0] + (unsigned long) 38 * a[3] * b[9] + (unsigned long) 19 * a[4] * b[8] + (unsigned long) 38 * a[5] * b[7] + (unsigned long) 19 * a[6] * b[6] + (unsigned long) 38 * a[7] * b[5] + (unsigned long) 19 * a[8] * b[4] + (unsigned long) 38 * a[9] * b[3];
    r[3]  =  (unsigned long) a[0] * b[3]   + (unsigned long) a[1] * b[2]   + (unsigned long) a[2] * b[1]   + (unsigned long) a[3] * b[0] + (unsigned long) 19 * a[4] * b[9] + (unsigned long) 19 * a[5] * b[8] + (unsigned long) 19 * a[6] * b[7] + (unsigned long) 19 * a[7] * b[6] + (unsigned long) 19 * a[8] * b[5] + 19 * (unsigned long) a[9] * b[4];
    r[4]  =  (unsigned long) a[0] * b[4]   + (unsigned long) 2 * a[1] * b[3]   + (unsigned long) a[2] * b[2]   + (unsigned long) 2 * a[3] * b[1]   + (unsigned long) a[4] * b[0] + (unsigned long) 38 * a[5] *b[9] + (unsigned long) 19 * a[6] * b[8] + (unsigned long) 38 * a[7] * b[7] + (unsigned long) 19 * a[8] * b[6] + (unsigned long) 38 * a[9] * b[5];
    r[5]  =  (unsigned long) a[0] * b[5]   + (unsigned long) a[1] * b[4]   + (unsigned long) a[2] * b[3]   + (unsigned long) a[3] * b[2]   + (unsigned long) a[4] * b[1]   + (unsigned long) a[5] * b[0] + (unsigned long) 19 * a[6] * b[9] + (unsigned long) 19 * a[7] * b[8] + (unsigned long) 19 * a[8] * b[7] + (unsigned long) 19 * a[9] * b[6];
    r[6]  =  (unsigned long) a[0] * b[6]   + (unsigned long) 2 *  a[1] * b[5]   + (unsigned long) a[2] * b[4]   + (unsigned long) 2 * a[3] * b[3]   + (unsigned long) a[4] * b[2]   + (unsigned long) 2 * a[5] * b[1]  + (unsigned long) a[6] * b[0] + (unsigned long) 38 * a[7] * b[9] + (unsigned long) 19 * a[8] * b[8] + (unsigned long) 38 * a[9] * b[7];
    r[7]  =  (unsigned long) a[0] * b[7]   + (unsigned long) a[1] * b[6]   + (unsigned long) a[2] * b[5]   + (unsigned long) a[3] * b[4]   + (unsigned long) a[4] * b[3]   + (unsigned long) a[5] * b[2]  + (unsigned long) a[6] * b[1]  + (unsigned long) a[7] * b[0] + (unsigned long) 19 * a[8] * b[9] + (unsigned long) 19 * a[9] * b[8];
    r[8]  =  (unsigned long) a[0] * b[8]   + (unsigned long) 2 * a[1] * b[7]   + (unsigned long) a[2] * b[6]   + (unsigned long) 2 * a[3] * b[5]   + (unsigned long) a[4] * b[4]   + (unsigned long) 2 * a[5] * b[3]  + (unsigned long) a[6] * b[2]  + (unsigned long) 2 * a[7] * b[1]  + (unsigned long) a[8] * b[0] + (unsigned long) 38 * a[9] * b[9];
    r[9]  =  (unsigned long) a[0] * b[9]   + (unsigned long) a[1] * b[8]   + (unsigned long) a[2] * b[7]   + (unsigned long) a[3] * b[6]   + (unsigned long) a[4] * b[5]   + (unsigned long) a[5] * b[4]  + (unsigned long) a[6] * b[3]  + (unsigned long) a[7] * b[2]  + (unsigned long) a[8] * b[1]  + (unsigned long) a[9] * b[0];
    total_ops += 245;

    //reduce_mul_512(r);

    reduce_mul_256(r);
    reduce_mul_256(r);
    
    for(int i = 0; i<NO_LIMBS; i++){
        res->d[i] = r[i];
    }
    // res->d[0] = r[0];
    // res->d[1] = r[1];
    // res->d[2] = r[2];
    // res->d[3] = r[3];
    // res->d[4] = r[4];

    
    reduce_256(res);
    //print_radix17( res, "mul" );
}



int cmp_coord_256(radix25_t *x1, radix25_t *x2){
    //return 1 if x1 == x2 

    int bool;
    bool = (x1->d[0] == x2->d[0]);
    total_ops += 1;
    for (size_t i = 1; i < NO_LIMBS; i++)
    {
        bool = bool && (x1->d[i] == x2->d[i]);
        total_ops += 2;
    }
    return bool;

}
int cmp_point_256(radix25_t *x1, radix25_t *y1, radix25_t *x2, radix25_t *y2){
    //return 1 if x1 == y1 and x2 == y2

    return cmp_coord_256(x1, y1) && cmp_coord_256(x2, y2);

}
int iszero_256(radix25_t *x){
    //return 1 if x == y == 0

    int bool;
    bool = (x->d[0] == 0);
     total_ops += 1;
    for (size_t i = 1; i < NO_LIMBS; i++)
    {
        bool = bool && (x->d[i] == 0);
         total_ops += 2;
    }
    return bool;
    
}
void set_ui_256(radix25_t *x, signed int ui){
    //set radix25_t to a unsigned int ui (lowest limb only, others 0)
    // maximum 17 bit numbers !!!

    x->d[0] = (signed int) ui;

    for (size_t i = 1; i < NO_LIMBS; i++)
    {
        x->d[i] = (signed int) 0;
    }
    
}


void square_256(radix25_t *res, radix25_t *x){
	unsigned long r[NO_LIMBS*2 - 1];
    // mul_256(res, x, x);
    // return;
    const signed int *a = (*x).d;

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
     total_ops += 9;
    r[0]  =   (unsigned long)a[0]*a[0]  + (unsigned long)   38*_2a[1]*a[9]+ (unsigned long)   19*_2a[2]*a[8]+ (unsigned long)   38*_2a[3]*a[7]+ (unsigned long)   19*_2a[4]*a[6]+ (unsigned long)   38*a[5]*a[5];
    r[1]  =   (unsigned long)_2a[0]*a[1]+ (unsigned long)   19*_2a[2]*a[9]+ (unsigned long)   19*_2a[3]*a[8]+ (unsigned long)   19*_2a[4]*a[7]+ (unsigned long)   19*_2a[5]*a[6];
    r[2]  =   (unsigned long)_2a[0]*a[2]+ (unsigned long)   _2a[1]*a[1]   + (unsigned long)   38*_2a[3]*a[9]+ (unsigned long)   19*_2a[4]*a[8]+ (unsigned long)   38*_2a[5]*a[7]+ (unsigned long)   19*a[6]*a[6];
    r[3]  =   (unsigned long)_2a[0]*a[3]+ (unsigned long)   _2a[1]*a[2]   + (unsigned long)   19*_2a[4]*a[9]+ (unsigned long)   19*_2a[5]*a[8]+ (unsigned long)   19*_2a[6]*a[7];
    r[4]  =   (unsigned long)_2a[0]*a[4]+ (unsigned long)   2*_2a[1]*a[3] + (unsigned long)   a[2]*a[2]     + (unsigned long)   38*_2a[5]*a[9]+ (unsigned long)   19*_2a[6]*a[8]+ (unsigned long)   38*a[7]*a[7];
    r[5]  =   (unsigned long)_2a[0]*a[5]+ (unsigned long)   _2a[1]*a[4]   + (unsigned long)   _2a[2]*a[3]   + (unsigned long)   19*_2a[6]*a[9]+ (unsigned long)   19*_2a[7]*a[8];
    r[6]  =   (unsigned long)_2a[0]*a[6]+ (unsigned long)   2*_2a[1]*a[5] + (unsigned long)   _2a[2]*a[4]   + (unsigned long)   _2a[3]*a[3]   + (unsigned long)   38*_2a[7]*a[9]+ (unsigned long)   19*a[8]*a[8];
    r[7]  =   (unsigned long)_2a[0]*a[7]+ (unsigned long)   _2a[1]*a[6]   + (unsigned long)   _2a[2]*a[5]   + (unsigned long)   _2a[3]*a[4]   + (unsigned long)   19*_2a[8]*a[9];
    r[8]  =   (unsigned long)_2a[0]*a[8]+ (unsigned long)   2*_2a[1]*a[7] + (unsigned long)   _2a[2]*a[6]   + (unsigned long)   2*_2a[3]*a[5] + (unsigned long)   a[4]*a[4]     + (unsigned long)   38*a[9]*a[9];
    r[9]  =   (unsigned long)_2a[0]*a[9]+ (unsigned long)   _2a[1]*a[8]   + (unsigned long)   _2a[2]*a[7]   + (unsigned long)   _2a[3]*a[6]   + (unsigned long)   _2a[4]*a[5];
    total_ops += 129;
    reduce_mul_256(r);
    reduce_mul_256(r);
    
     for(int i = 0; i<NO_LIMBS; i++){
        res->d[i] = r[i];
    }

    // res->d[0] = r[0];
    // res->d[1] = r[1];
    // res->d[2] = r[2];
    // res->d[3] = r[3];
    // res->d[4] = r[4];

    
    reduce_256(res);
	//print_radix17( res, "square" );
}
void inverse_256(radix25_t *res, radix25_t *x){
    radix25_t z2, z9, z11, z2_5_0, z2_10_0, z2_20_0, z2_50_0, z2_100_0, t;
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
    //print_radix17( res, "inv" );
}

void div_256(radix25_t *res, radix25_t *a, radix25_t *b){
    // inverse and mult
    radix25_t inv;
    inverse_256(&inv, b);
    mul_256(res, a, &inv);
}

int get_total_ops(){
    return total_ops;
}