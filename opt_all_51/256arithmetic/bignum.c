#include <stdio.h>
#include <immintrin.h>

void copy_256(long dest[5], long src[5]){
    /*
    set dest to src
    */
    for (int i = 0; i < 5; i++) {
        dest[i] = src[i];
    } 
}

void reduce_2_to_2(long x[5]){
    /*
    reduce x
    */
    signed long carry;
    for (int i = 0; i < 4; i++) {
        carry = x[i] >> 51;
        x[i+1] += carry;
        carry <<= 51;
        x[i] -= carry;
    }
    
    carry = x[4] >> 51;
    x[0] += 19*carry;
    carry <<= 51;
    x[4] -= carry;  
}

void reduce_4_to_2(__int128 r[9]){
    /*
    reduce r
    */
    signed __int128 carry;
    for (int i = 0; i < 4; i++) {
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

void add_256(long res[5], long a[5], long b[5]) {
    /*
    res = a + b
    */
    /*
    for (int i = 0; i < 5; i++) {
        res[i] = a[i] + b[i];
    }
    */
    __m256i aVec = _mm256_loadu_si256( (const __m256i*) a );
    __m256i bVec = _mm256_loadu_si256( (const __m256i*) b );
    aVec = _mm256_add_epi64( aVec, bVec );
    _mm256_storeu_si256( (__m256i *) res, aVec );
    res[4] = a[4] + b[4];
}

void add_1_256( long a[5] ) {
	/*
	res = a + 1
	*/
	a[0] = a[0] + 1;
}

void sub_256(long res[5], long a[5], long b[5]) {
    /*
    res = a - b
    */
    for (int i = 0; i < 5; i++) {
        res[i] = a[i] - b[i];
    }
    reduce_2_to_2(res);
}

void reduce_mul_512(__int128 r[9]) {
    /*
    reduce r
    */
    r[0] += 19*r[5];
    r[1] += 19*r[6];
    r[2] += 19*r[7];
    r[3] += 19*r[8];
}

void mul_256(long res[5], long x[5], long y[5]) {
    /*
    res = x * y
    */
    unsigned __int128 r[9];
    
    long x0 = x[0];
    long x1 = x[1];
    long x2 = x[2];
    long x3 = x[3];
    long x4 = x[4];
    
    long y0 = y[0];
    long y1 = y[1];
    long y2 = y[2];
    long y3 = y[3];
    long y4 = y[4];
    
    r[0] = (__int128_t) x0*y0;
    r[1] = (__int128_t) x0*y1 + (__int128_t) x1*y0;
    r[2] = (__int128_t) x0*y2 + (__int128_t) x1*y1 + (__int128_t) x2*y0;
    r[3] = (__int128_t) x0*y3 + (__int128_t) x1*y2 + (__int128_t) x2*y1 + (__int128_t) x3*y0;
    r[4] = (__int128_t) x0*y4 + (__int128_t) x1*y3 + (__int128_t) x2*y2 + (__int128_t) x3*y1 + (__int128_t) x4*y0;
    r[5] = (__int128_t) x1*y4 + (__int128_t) x2*y3 + (__int128_t) x3*y2 + (__int128_t) x4*y1;
    r[6] = (__int128_t) x2*y4 + (__int128_t) x3*y3 + (__int128_t) x4*y2;
    r[7] = (__int128_t) x3*y4 + (__int128_t) x4*y3;
    r[8] = (__int128_t) x4*y4;

    reduce_mul_512(r);
	
	reduce_4_to_2(r);
	    
    res[0] = r[0];
    res[1] = r[1];
    res[2] = r[2];
    res[3] = r[3];
    res[4] = r[4];
    
    reduce_2_to_2(res);
}



int cmp_coord_256(long x1[5], long x2[5]) {
    /*
    return 1 if x1 == x2
    */ 
    return (x1[0] == x2[0]) && (x1[1] == x2[1]) && (x1[2] == x2[2]) && (x1[3] == x2[3]) && (x1[4] == x2[4]);
}

int cmp_point_256(long x1[5], long y1[5], long x2[5], long y2[5]) {
    /*
    return 1 if x1 == y1 and x2 == y2
    */
    return cmp_coord_256(x1, y1) && cmp_coord_256(x2, y2);
}

int iszero_256(long x[5]) {
    /*
    return 1 if x == y == 0
    */
    return (x[0] == 0) && (x[1] == 0) && (x[2] == 0) && (x[3] == 0) && (x[4] == 0);
}

void set_ui_256(long x[5], signed int ui){
    /*
    set lowest limb to ui and other limbs to 0
    maximum 51 bit numbers !!!
    */
    x[0] = (signed long long) ui;
    x[1] = (signed long long) 0;
    x[2] = (signed long long) 0;
    x[3] = (signed long long) 0;
    x[4] = (signed long long) 0;
}

void square_256(long res[5], long a[5]){
    /*
    res = x^2
    */
    __int128_t r[9];

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
    
    res[0] = r[0];
    res[1] = r[1];
    res[2] = r[2];
    res[3] = r[3];
    res[4] = r[4];

    reduce_2_to_2(res);
}

void inverse_256(long res[5], long x[5]) {
    /*
    res = x^-1 (multiplicative inverse modulo p)
    */
    long z2[5], z9[5], z11[5], z2_5_0[5], z2_10_0[5], z2_20_0[5], z2_50_0[5], z2_100_0[5], t[5];
    int i;
    /* 2 */              square_256(z2,x );
    /* 4 */              square_256(t ,z2);
    /* 8 */              square_256(t ,t );
    /* 9 */              mul_256   (z9 ,t ,x );
    /* 11 */             mul_256   (z11,z9,z2);
    /* 22 */             square_256(t  ,z11);
    /* 2^5 - 2^0 = 31 */ mul_256   (z2_5_0,t,z9);
    /* 2^6 - 2^1 */      square_256(t     ,z2_5_0);
    /* 2^10 - 2^5 */     for (i = 1;i < 5;i++) { square_256(t,t); }
    /* 2^10 - 2^0 */     mul_256   (z2_10_0,t,z2_5_0);
    /* 2^11 - 2^1 */     square_256(t,z2_10_0);
    /* 2^20 - 2^10 */    for (i = 1;i < 10;i++) { square_256(t,t); }
    /* 2^20 - 2^0 */     mul_256   (z2_20_0,t,z2_10_0);
    /* 2^21 - 2^1 */     square_256(t      ,z2_20_0);
    /* 2^40 - 2^20 */    for (i = 1;i < 20;i++) { square_256(t,t); }
    /* 2^40 - 2^0 */     mul_256   (t,t,z2_20_0);
    /* 2^41 - 2^1 */     square_256(t,t);
    /* 2^50 - 2^10 */    for (i = 1;i < 10;i++) { square_256(t,t); }
    /* 2^50 - 2^0 */     mul_256   (z2_50_0,t,z2_10_0);
    /* 2^51 - 2^1 */     square_256(t      ,z2_50_0);
    /* 2^100 - 2^50 */   for (i = 1;i < 50;i++) { square_256(t,t); }
    /* 2^100 - 2^0 */    mul_256   (z2_100_0,t,z2_50_0);
    /* 2^101 - 2^1 */    square_256(t       ,z2_100_0);
    /* 2^200 - 2^100 */  for (i = 1;i < 100;i++) { square_256(t,t); }
    /* 2^200 - 2^0 */    mul_256   (t,t,z2_100_0);
    /* 2^201 - 2^1 */    square_256(t,t);
    /* 2^250 - 2^50 */   for (i = 1;i < 50;i++) { square_256(t,t); }
    /* 2^250 - 2^0 */    mul_256   (t,t,z2_50_0);
    /* 2^251 - 2^1 */    square_256(t,t);
    /* 2^252 - 2^2 */    square_256(t,t);
    /* 2^253 - 2^3 */    square_256(t,t);
    /* 2^254 - 2^4 */    square_256(t,t);
    /* 2^255 - 2^5 */    square_256(t,t);
    /* 2^255 - 21 */     mul_256   (res,t,z11);
}


void div_256(long res[5], long a[5], long b[5]){
    /*
    res = a/b
    */
    // inverse and mult
    long inv[5];
    inverse_256(inv, b);
    mul_256(res, a, inv);
}