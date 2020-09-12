#include <stdio.h>
#include "bignum.h"
#include <immintrin.h>
#include <emmintrin.h>
#include "../../_common_libs/radix17.h"

const int NO_LIMBS = 15;
const int LIMB_SIZE = 17;

unsigned long times19(unsigned long a)
{
  return (a << 4) + (a << 1) + a;
}

void print_vec_int(__m256i vec, char* str)
{
    /*
    Print the content of a __m256i vector (8xint)
    */
    int* i_v = (int*) &vec;
    printf("%s: %d %d %d %d %d %d %d %d\n",str, i_v[0], i_v[1], i_v[2], i_v[3], i_v[4], i_v[5], i_v[6], i_v[7]);
}

void print_vec_long(__m256i vec, char* str)
{
    /*
    Print the content of a __m256i vector (4xlong)
    */
    long int* i_v = (long int*) &vec;
    printf("%s: %lu %lu %lu %lu\n",str, i_v[0], i_v[1], i_v[2], i_v[3]);
}

void check_vec_long(__m256i vec, long v0, long v1, long v2, long v3){
    /*
    check if vec (4xlong) is equal to (v0, v1, v2, v3)
    */
    long int* i_v = (long int*) &vec;
    if(i_v[0] == v0 && i_v[1] == v1 && i_v[2] == v2 && i_v[3] == v3){
        printf("OK\n");
    }else{
        printf("Not OK!\n");
        print_vec_long(vec, "wrong vector");
        printf("right values: %lu %lu %lu %lu\n", v0, v1, v2, v3);
    }
}


void copy_256(radix17_t *dest, radix17_t src)
{
    //set dest to src

    for (size_t i = 0; i < NO_LIMBS; i++)
    {
        dest->d[i] = src.d[i];
    }
}
void reduce_256(radix17_t *x)
{

    signed long carry;
    signed int d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14;

    d0 = x->d[0];
    d1 = x->d[1];
    d2 = x->d[2];
    d3 = x->d[3];
    d4 = x->d[4];
    d5 = x->d[5];
    d6 = x->d[6];
    d7 = x->d[7];
    d8 = x->d[8];
    d9 = x->d[9];
    d10 = x->d[10];
    d11 = x->d[11];
    d12 = x->d[12];
    d13 = x->d[13];
    d14 = x->d[14];
    
    carry = d0 >> LIMB_SIZE;
    d1 += carry;
    carry <<= LIMB_SIZE;
    d0 -= carry;

    carry = d1 >> LIMB_SIZE;
    d2 += carry;
    carry <<= LIMB_SIZE;
    d1 -= carry;

    carry = d2 >> LIMB_SIZE;
    d3 += carry;
    carry <<= LIMB_SIZE;
    d2 -= carry;

    carry = d3 >> LIMB_SIZE;
    d4 += carry;
    carry <<= LIMB_SIZE;
    d3 -= carry;

    carry = d4 >> LIMB_SIZE;
    d5 += carry;
    carry <<= LIMB_SIZE;
    d4 -= carry;

    carry = d5 >> LIMB_SIZE;
    d6 += carry;
    carry <<= LIMB_SIZE;
    d5 -= carry;

    carry = d6 >> LIMB_SIZE;
    d7 += carry;
    carry <<= LIMB_SIZE;
    d6 -= carry;

    carry = d7 >> LIMB_SIZE;
    d8 += carry;
    carry <<= LIMB_SIZE;
    d7 -= carry;

    carry = d8 >> LIMB_SIZE;
    d9 += carry;
    carry <<= LIMB_SIZE;
    d8 -= carry;

    carry = d9 >> LIMB_SIZE;
    d10 += carry;
    carry <<= LIMB_SIZE;
    d9 -= carry;

    carry = d10 >> LIMB_SIZE;
    d11 += carry;
    carry <<= LIMB_SIZE;
    d10 -= carry;

    carry = d11 >> LIMB_SIZE;
    d12 += carry;
    carry <<= LIMB_SIZE;
    d11 -= carry;

    carry = d12 >> LIMB_SIZE;
    d13 += carry;
    carry <<= LIMB_SIZE;
    d12 -= carry;

    carry = d13 >> LIMB_SIZE;
    d14 += carry;
    carry <<= LIMB_SIZE;
    d13 -= carry;

    carry = d14 >> LIMB_SIZE;
    d0 += 19 * carry;
    carry <<= LIMB_SIZE;
    d14 -= carry;

    x->d[0] = d0;
    x->d[1] = d1;
    x->d[2] = d2;
    x->d[3] = d3;
    x->d[4] = d4;
    x->d[5] = d5;
    x->d[6] = d6;
    x->d[7] = d7;
    x->d[8] = d8;
    x->d[9] = d9;
    x->d[10] = d10;
    x->d[11] = d11;
    x->d[12] = d12;
    x->d[13] = d13;
    x->d[14] = d14;
}

void reduce_mul_256(unsigned long r[30])
{

    signed long carry;
    unsigned long r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14;

    r0 = r[0];
    r1 = r[1];
    r2 = r[2];
    r3 = r[3];
    r4 = r[4];
    r5 = r[5];
    r6 = r[6];
    r7 = r[7];
    r8 = r[8];
    r9 = r[9];
    r10 = r[10];
    r11 = r[11];
    r12 = r[12];
    r13 = r[13];
    r14 = r[14];
    
    carry = r0 >> LIMB_SIZE;
    r1 += carry;
    carry <<= LIMB_SIZE;
    r0 -= carry;

    carry = r1 >> LIMB_SIZE;
    r2 += carry;
    carry <<= LIMB_SIZE;
    r1 -= carry;

    carry = r2 >> LIMB_SIZE;
    r3 += carry;
    carry <<= LIMB_SIZE;
    r2 -= carry;

    carry = r3 >> LIMB_SIZE;
    r4 += carry;
    carry <<= LIMB_SIZE;
    r3 -= carry;

    carry = r4 >> LIMB_SIZE;
    r5 += carry;
    carry <<= LIMB_SIZE;
    r4 -= carry;

    carry = r5 >> LIMB_SIZE;
    r6 += carry;
    carry <<= LIMB_SIZE;
    r5 -= carry;

    carry = r6 >> LIMB_SIZE;
    r7 += carry;
    carry <<= LIMB_SIZE;
    r6 -= carry;

    carry = r7 >> LIMB_SIZE;
    r8 += carry;
    carry <<= LIMB_SIZE;
    r7 -= carry;

    carry = r8 >> LIMB_SIZE;
    r9 += carry;
    carry <<= LIMB_SIZE;
    r8 -= carry;

    carry = r9 >> LIMB_SIZE;
    r10 += carry;
    carry <<= LIMB_SIZE;
    r9 -= carry;

    carry = r10 >> LIMB_SIZE;
    r11 += carry;
    carry <<= LIMB_SIZE;
    r10 -= carry;

    carry = r11 >> LIMB_SIZE;
    r12 += carry;
    carry <<= LIMB_SIZE;
    r11 -= carry;

    carry = r12 >> LIMB_SIZE;
    r13 += carry;
    carry <<= LIMB_SIZE;
    r12 -= carry;

    carry = r13 >> LIMB_SIZE;
    r14 += carry;
    carry <<= LIMB_SIZE;
    r13 -= carry;

    carry = r14 >> LIMB_SIZE;
    r0 += 19 * carry;
    carry <<= LIMB_SIZE;
    r14 -= carry;

    r[0] = r0;
    r[1] = r1;
    r[2] = r2;
    r[3] = r3;
    r[4] = r4;
    r[5] = r5;
    r[6] = r6;
    r[7] = r7;
    r[8] = r8;
    r[9] = r9;
    r[10] = r10;
    r[11] = r11;
    r[12] = r12;
    r[13] = r13;
    r[14] = r14; 
}

void add_256(radix17_t *res, radix17_t *a, radix17_t *b)
{
    __m256i a1, a2, b1, b2, res1, res2;
    // Load a and b into vectors
    a1 = _mm256_loadu_si256((const __m256i*) a->d);
    a2 = _mm256_loadu_si256((const __m256i*) ((a->d)+8));
    b1 = _mm256_loadu_si256((const __m256i*) b->d);
    b2 = _mm256_loadu_si256((const __m256i*) ((b->d)+8));
    // a + b
    res1 = _mm256_add_epi32(a1,b1);
    res2 = _mm256_add_epi32(a2,b2);
    // store the result back into res
    _mm256_storeu_si256 ((__m256i *) res->d, res1);
    _mm256_storeu_si256 ((__m256i *) ((res->d)+8), res2);
    // reduce res
    reduce_256(res);
}

void sub_256(radix17_t *res, radix17_t *a, radix17_t *b)
{

     __m256i a1, a2, b1, b2, res1, res2;
    // Load a and b into vectors
    a1 = _mm256_loadu_si256((const __m256i*) a->d);
    a2 = _mm256_loadu_si256((const __m256i*) ((a->d)+8));
    b1 = _mm256_loadu_si256((const __m256i*) b->d);
    b2 = _mm256_loadu_si256((const __m256i*) ((b->d)+8));
    // a - b
    res1 = _mm256_sub_epi32(a1,b1);
    res2 = _mm256_sub_epi32(a2,b2);
    // store the result back into res
    _mm256_storeu_si256 ((__m256i *) res->d, res1);
    _mm256_storeu_si256 ((__m256i *) ((res->d)+8), res2);
    // reduce res
    reduce_256(res);
    reduce_256(res);
}
void reduce_mul_512(unsigned long r[30])
{
    // backup r[14]
    long r14 = r[14];
    __m256i r0_3, r4_7,r8_11, r12_15, 
    r15_18, r19_22, r23_16, r27_30,
    r15_18x16, r19_22x16, r23_16x16, r27_30x16,
    r15_18x2, r19_22x2, r23_16x2, r27_30x2;
    r0_3 = _mm256_loadu_si256((const __m256i*) &(r[0]));
    r4_7 = _mm256_loadu_si256((const __m256i*) &(r[4]));
    r8_11 = _mm256_loadu_si256((const __m256i*) &(r[8]));
    r12_15 = _mm256_loadu_si256((const __m256i*) &(r[12]));
    // The following needs to be multiplyed with 19 --> 19x = 16x + 2x + x --> x<<4 + x<<1 +x
    r15_18 = _mm256_loadu_si256((const __m256i*) &(r[15]));
    r19_22 = _mm256_loadu_si256((const __m256i*) &(r[19]));
    r23_16 = _mm256_loadu_si256((const __m256i*) &(r[23]));
    r27_30 = _mm256_loadu_si256((const __m256i*) &(r[27]));
    // 16x
    r15_18x16 = _mm256_slli_epi64(r15_18, 4);
    r19_22x16 = _mm256_slli_epi64(r19_22, 4);
    r23_16x16 = _mm256_slli_epi64(r23_16, 4);
    r27_30x16 = _mm256_slli_epi64(r27_30, 4);
    // 2x
    r15_18x2 = _mm256_slli_epi64(r15_18, 1);
    r19_22x2 = _mm256_slli_epi64(r19_22, 1);
    r23_16x2 = _mm256_slli_epi64(r23_16, 1);
    r27_30x2 = _mm256_slli_epi64(r27_30, 1);
    // x = x + 16x
    r15_18 = _mm256_add_epi64(r15_18, r15_18x16);
    r19_22 = _mm256_add_epi64(r19_22, r19_22x16);
    r23_16 = _mm256_add_epi64(r23_16, r23_16x16);
    r27_30 = _mm256_add_epi64(r27_30, r27_30x16);
    // x = x + 2x
    r15_18 = _mm256_add_epi64(r15_18, r15_18x2);
    r19_22 = _mm256_add_epi64(r19_22, r19_22x2);
    r23_16 = _mm256_add_epi64(r23_16, r23_16x2);
    r27_30 = _mm256_add_epi64(r27_30, r27_30x2);
    // add upper limbs to lower ones
    r0_3 = _mm256_add_epi64(r15_18, r0_3);
    r4_7 = _mm256_add_epi64(r19_22, r4_7);
    r8_11 = _mm256_add_epi64(r23_16, r8_11);
    r12_15 = _mm256_add_epi64(r27_30, r12_15);
    // store back into r
    _mm256_storeu_si256 ((__m256i *) &(r[0]), r0_3);
    _mm256_storeu_si256 ((__m256i *) &(r[4]), r4_7);
    _mm256_storeu_si256 ((__m256i *) &(r[8]), r8_11);
    _mm256_storeu_si256 ((__m256i *) &(r[12]), r12_15);
    // restore original r[14]
    r[14] = r14;
}


void mul_256(radix17_t *res, radix17_t *x, radix17_t *y)
{
	
    unsigned long r[NO_LIMBS*2];
	/*
	//non-vectorized code as reference
	
    //unsigned long r[NO_LIMBS*2 - 1];

    signed long a0 = x->d[0];
    signed long a1 = x->d[1];
    signed long a2 = x->d[2];
    signed long a3 = x->d[3];
    signed long a4 = x->d[4];
    signed long a5 = x->d[5];
    signed long a6 = x->d[6];
    signed long a7 = x->d[7];
    signed long a8 = x->d[8];
    signed long a9 = x->d[9];
    signed long a10 = x->d[10];
    signed long a11 = x->d[11];
    signed long a12 = x->d[12];
    signed long a13 = x->d[13];
    signed long a14 = x->d[14];

    signed long b0 = y->d[0];
    signed long b1 = y->d[1];
    signed long b2 = y->d[2];
    signed long b3 = y->d[3];
    signed long b4 = y->d[4];
    signed long b5 = y->d[5];
    signed long b6 = y->d[6];
    signed long b7 = y->d[7];
    signed long b8 = y->d[8];
    signed long b9 = y->d[9];
    signed long b10 = y->d[10];
    signed long b11 = y->d[11];
    signed long b12 = y->d[12];
    signed long b13 = y->d[13];
    signed long b14 = y->d[14];
	
	
	//       column1          column2          column3          column4          column5         column6       column7         column8        column9         column10       column11        column12       column13       column14      column15
    r[0]  =  a0  * b0;
    r[1]  =  a0  * b1     +   a1  * b0;
    r[2]  =  a0  * b2     +   a1  * b1     +   a2  * b0;
    r[3]  =  a0  * b3     +   a1  * b2     +   a2  * b1     +   a3  * b0;
    r[4]  =  a0  * b4     +   a1  * b3     +   a2  * b2     +   a3  * b1     +   a4  * b0;
    r[5]  =  a0  * b5     +   a1  * b4     +   a2  * b3     +   a3  * b2     +   a4  * b1    +   a5  * b0;
    r[6]  =  a0  * b6     +   a1  * b5     +   a2  * b4     +   a3  * b3     +   a4  * b2    +   a5  * b1   +   a6  * b0;
    r[7]  =  a0  * b7     +   a1  * b6     +   a2  * b5     +   a3  * b4     +   a4  * b3    +   a5  * b2   +   a6  * b1   +   a7  * b0;
    r[8]  =  a0  * b8     +   a1  * b7     +   a2  * b6     +   a3  * b5     +   a4  * b4    +   a5  * b3   +   a6  * b2   +   a7  * b1   +   a8  * b0;
    r[9]  =  a0  * b9     +   a1  * b8     +   a2  * b7     +   a3  * b6     +   a4  * b5    +   a5  * b4   +   a6  * b3   +   a7  * b2   +   a8  * b1   +    a9  * b0;
    r[10] =  a0  * b10    +   a1  * b9     +   a2  * b8     +   a3  * b7     +   a4  * b6    +   a5  * b5   +   a6  * b4   +   a7  * b3   +   a8  * b2   +    a9  * b1   +   a10 * b0;
    r[11] =  a0  * b11    +   a1  * b10    +   a2  * b9     +   a3  * b8     +   a4  * b7    +   a5  * b6   +   a6  * b5   +   a7  * b4   +   a8  * b3   +    a9  * b2   +   a10 * b1   +   a11 * b0;
    r[12] =  a0  * b12    +   a1  * b11    +   a2  * b10    +   a3  * b9     +   a4  * b8    +   a5  * b7   +   a6  * b6   +   a7  * b5   +   a8  * b4   +    a9  * b3   +   a10 * b2   +   a11 * b1   +   a12 * b0;
    r[13] =  a0  * b13    +   a1  * b12    +   a2  * b11    +   a3  * b10    +   a4  * b9    +   a5  * b8   +   a6  * b7   +   a7  * b6   +   a8  * b5   +    a9  * b4   +   a10 * b3   +   a11 * b2   +   a12 * b1   +   a13 * b0;
    r[14] =  a0  * b14    +   a1  * b13    +   a2  * b12    +   a3  * b11    +   a4  * b10   +   a5  * b9   +   a6  * b8   +   a7  * b7   +   a8  * b6   +    a9  * b5   +   a10 * b4   +   a11 * b3   +   a12 * b2   +   a13 * b1   +   a14 * b0;
    
    r[15] =  a1  * b14    +   a2  * b13    +   a3  * b12    +   a4  * b11    +   a5  * b10   +   a6  * b9   +   a7  * b8   +   a8  * b7   +   a9  * b6   +    a10 * b5   +   a11 * b4   +   a12 * b3   +   a13 * b2   +   a14 * b1;
    r[16] =  a2  * b14    +   a3  * b13    +   a4  * b12    +   a5  * b11    +   a6  * b10   +   a7  * b9   +   a8  * b8   +   a9  * b7   +   a10 * b6   +    a11 * b5   +   a12 * b4   +   a13 * b3   +   a14 * b2;
    r[17] =  a3  * b14    +   a4  * b13    +   a5  * b12    +   a6  * b11    +   a7  * b10   +   a8  * b9   +   a9  * b8   +   a10 * b7   +   a11 * b6   +    a12 * b5   +   a13 * b4   +   a14 * b3;
    r[18] =  a4  * b14    +   a5  * b13    +   a6  * b12    +   a7  * b11    +   a8  * b10   +   a9  * b9   +   a10 * b8   +   a11 * b7   +   a12 * b6   +    a13 * b5   +   a14 * b4;
    r[19] =  a5  * b14    +   a6  * b13    +   a7  * b12    +   a8  * b11    +   a9  * b10   +   a10 * b9   +   a11 * b8   +   a12 * b7   +   a13 * b6   +    a14 * b5;
    r[20] =  a6  * b14    +   a7  * b13    +   a8  * b12    +   a9  * b11    +   a10 * b10   +   a11 * b9   +   a12 * b8   +   a13 * b7   +   a14 * b6;
    r[21] =  a7  * b14    +   a8  * b13    +   a9  * b12    +   a10 * b11    +   a11 * b10   +   a12 * b9   +   a13 * b8   +   a14 * b7;
    r[22] =  a8  * b14    +   a9  * b13    +   a10 * b12    +   a11 * b11    +   a12 * b10   +   a13 * b9   +   a14 * b8;
    r[23] =  a9  * b14    +   a10 * b13    +   a11 * b12    +   a12 * b11    +   a13 * b10   +   a14 * b9;
    r[24] =  a10 * b14    +   a11 * b13    +   a12 * b12    +   a13 * b11    +   a14 * b10;
    r[25] =  a11 * b14    +   a12 * b13    +   a13 * b12    +   a14 * b11;
    r[26] =  a12 * b14    +   a13 * b13    +   a14 * b12;
    r[27] =  a13 * b14    +   a14 * b13;
    r[28] =  a14 * b14;
	*/

    signed int a0 = x->d[0];
    signed int a1 = x->d[1];
    signed int a2 = x->d[2];
    signed int a3 = x->d[3];
    signed int a4 = x->d[4];
    signed int a5 = x->d[5];
    signed int a6 = x->d[6];
    signed int a7 = x->d[7];
    signed int a8 = x->d[8];
    signed int a9 = x->d[9];
    signed int a10 = x->d[10];
    signed int a11 = x->d[11];
    signed int a12 = x->d[12];
    signed int a13 = x->d[13];
    signed int a14 = x->d[14];

    signed int b0 = y->d[0];
    signed int b1 = y->d[1];
    signed int b2 = y->d[2];
    signed int b3 = y->d[3];
    signed int b4 = y->d[4];
    signed int b5 = y->d[5];
    signed int b6 = y->d[6];
    signed int b7 = y->d[7];
    signed int b8 = y->d[8];
    signed int b9 = y->d[9];
    signed int b10 = y->d[10];
    signed int b11 = y->d[11];
    signed int b12 = y->d[12];
    signed int b13 = y->d[13];
    signed int b14 = y->d[14];

	// prepare vector data
	// note: these vectors must be reset every time in order to avoid costly permutes
	//       in order to shuffle between 128-bit lanes
	__m256i a_0_3, a_4_7, a_8_11, a_12_15;
	__m256i b_0_3, b_4_7, b_8_11, b_12_15;
	
	__m256i all_vals;	
	__m256i zeroes  = _mm256_set1_epi32( 0 );
	
	// temporary variables
	__m256i tmp_0_3, tmp_4_7, tmp_8_11, tmp_12_15;
	
	// result variables
	__m256i r_0_3, r_4_7, r_8_11, r_12_15, r_16_19, r_20_23, r_24_27, r_28_31;
	
	// multiply
	// TOP HALF
	// column1:
	//   for i=0 to 14: r[i] = a[0] * b[i]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a0, 0, a0, 0, a0, 0, a0, 0 );
	//   - set b
	b_0_3   = _mm256_setr_epi32( b0 , 0, b1 , 0, b2 , 0, b3 , 0 );
	b_4_7   = _mm256_setr_epi32( b4 , 0, b5 , 0, b6 , 0, b7 , 0 );
	b_8_11  = _mm256_setr_epi32( b8 , 0, b9 , 0, b10, 0, b11, 0 );
	b_12_15 = _mm256_setr_epi32( b12, 0, b13, 0, b14, 0, 0  , 0 );
	//   - multiply
	r_0_3   = _mm256_mul_epi32( all_vals, b_0_3   );
	r_4_7   = _mm256_mul_epi32( all_vals, b_4_7   );
	r_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	r_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	// column2:
	//   for i=1 to 14: r[i] += a[1] * b[i-1]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a1, 0, a1, 0, a1, 0, a1, 0 );
	//   - set b
	b_0_3   = _mm256_setr_epi32( 0  , 0, b0 , 0, b1 , 0, b2 , 0 );
	b_4_7   = _mm256_setr_epi32( b3 , 0, b4 , 0, b5 , 0, b6 , 0 );
	b_8_11  = _mm256_setr_epi32( b7 , 0, b8 , 0, b9 , 0, b10, 0 );
	b_12_15 = _mm256_setr_epi32( b11, 0, b12, 0, b13, 0, 0  , 0 );
	//   - multiply
	tmp_0_3   = _mm256_mul_epi32( all_vals, b_0_3   );
	tmp_4_7   = _mm256_mul_epi32( all_vals, b_4_7   );
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_0_3   = _mm256_add_epi64( r_0_3  , tmp_0_3   );
	r_4_7   = _mm256_add_epi64( r_4_7  , tmp_4_7   );
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column3:
	//   for i=2 to 14: r[i] += a[2] * b[i-2]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a2, 0, a2, 0, a2, 0, a2, 0 );
	//   - set b
	b_0_3   = _mm256_setr_epi32( 0  , 0, 0  , 0, b0 , 0, b1 , 0 );
	b_4_7   = _mm256_setr_epi32( b2 , 0, b3 , 0, b4 , 0, b5 , 0 );
	b_8_11  = _mm256_setr_epi32( b6 , 0, b7 , 0, b8 , 0, b9 , 0 );
	b_12_15 = _mm256_setr_epi32( b10, 0, b11, 0, b12, 0, 0  , 0 );
	//   - multiply
	tmp_0_3   = _mm256_mul_epi32( all_vals, b_0_3   );
	tmp_4_7   = _mm256_mul_epi32( all_vals, b_4_7   );
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_0_3   = _mm256_add_epi64( r_0_3  , tmp_0_3   );
	r_4_7   = _mm256_add_epi64( r_4_7  , tmp_4_7   );
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column4:
	//   for i=3 to 14: r[i] += a[3] * b[i-3]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a3, 0, a3, 0, a3, 0, a3, 0 );
	//   - set b
	b_0_3   = _mm256_setr_epi32( 0  , 0, 0  , 0, 0  , 0, b0 , 0 );
	b_4_7   = _mm256_setr_epi32( b1 , 0, b2 , 0, b3 , 0, b4 , 0 );
	b_8_11  = _mm256_setr_epi32( b5 , 0, b6 , 0, b7 , 0, b8 , 0 );
	b_12_15 = _mm256_setr_epi32( b9 , 0, b10, 0, b11, 0, 0  , 0 );
	//   - multiply
	tmp_0_3   = _mm256_mul_epi32( all_vals, b_0_3   );
	tmp_4_7   = _mm256_mul_epi32( all_vals, b_4_7   );
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_0_3   = _mm256_add_epi64( r_0_3  , tmp_0_3   );
	r_4_7   = _mm256_add_epi64( r_4_7  , tmp_4_7   );
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column5:
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a4, 0, a4, 0, a4, 0, a4, 0 );
	//   for i=4 to 14: r[i] += a[4] * b[i-4]
	b_4_7   = _mm256_setr_epi32( b0 , 0, b1 , 0, b2 , 0, b3 , 0 );
	b_8_11  = _mm256_setr_epi32( b4 , 0, b5 , 0, b6 , 0, b7 , 0 );
	b_12_15 = _mm256_setr_epi32( b8 , 0, b9 , 0, b10, 0, 0  , 0 );
	//   - multiply
	tmp_4_7   = _mm256_mul_epi32( all_vals, b_4_7   );
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_4_7   = _mm256_add_epi64( r_4_7  , tmp_4_7   );
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column6:
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a5, 0, a5, 0, a5, 0, a5, 0 );
	//   for i=5 to 14: r[i] += a[5] * b[i-5]
	b_4_7   = _mm256_setr_epi32( 0  , 0, b0 , 0, b1 , 0, b2 , 0 );
	b_8_11  = _mm256_setr_epi32( b3 , 0, b4 , 0, b5 , 0, b6 , 0 );
	b_12_15 = _mm256_setr_epi32( b7 , 0, b8 , 0, b9 , 0, 0  , 0 );
	//   - multiply
	tmp_4_7   = _mm256_mul_epi32( all_vals, b_4_7   );
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_4_7   = _mm256_add_epi64( r_4_7  , tmp_4_7   );
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column7:
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a6, 0, a6, 0, a6, 0, a6, 0 );
	//   for i=6 to 14: r[i] += a[6] * b[i-6]
	b_4_7   = _mm256_setr_epi32( 0  , 0, 0  , 0, b0 , 0, b1 , 0 );
	b_8_11  = _mm256_setr_epi32( b2 , 0, b3 , 0, b4 , 0, b5 , 0 );
	b_12_15 = _mm256_setr_epi32( b6 , 0, b7 , 0, b8 , 0, 0  , 0 );
	//   - multiply
	tmp_4_7   = _mm256_mul_epi32( all_vals, b_4_7   );
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_4_7   = _mm256_add_epi64( r_4_7  , tmp_4_7   );
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column8:
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a7, 0, a7, 0, a7, 0, a7, 0 );
	//   for i=7 to 14: r[i] += a[7] * b[i-7]
	b_4_7   = _mm256_setr_epi32( 0  , 0, 0  , 0, 0  , 0, b0 , 0 );
	b_8_11  = _mm256_setr_epi32( b1 , 0, b2 , 0, b3 , 0, b4 , 0 );
	b_12_15 = _mm256_setr_epi32( b5 , 0, b6 , 0, b7 , 0, 0  , 0 );
	//   - multiply
	tmp_4_7   = _mm256_mul_epi32( all_vals, b_4_7   );
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_4_7   = _mm256_add_epi64( r_4_7  , tmp_4_7   );
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column9:
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a8, 0, a8, 0, a8, 0, a8, 0 );
	//   for i=8 to 14: r[i] += a[8] * b[i-8]
	b_8_11  = _mm256_setr_epi32( b0 , 0, b1 , 0, b2 , 0, b3 , 0 );
	b_12_15 = _mm256_setr_epi32( b4 , 0, b5 , 0, b6 , 0, 0  , 0 );
	//   - multiply
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column10:
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a9, 0, a9, 0, a9, 0, a9, 0 );
	//   for i=9 to 14: r[i] += a[9] * b[i-9]
	b_8_11  = _mm256_setr_epi32( 0  , 0, b0 , 0, b1 , 0, b2 , 0 );
	b_12_15 = _mm256_setr_epi32( b3 , 0, b4 , 0, b5 , 0, 0  , 0 );
	//   - multiply
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column11:
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a10, 0, a10, 0, a10, 0, a10, 0 );
	//   for i=10 to 14: r[i] += a[10] * b[i-10]
	b_8_11  = _mm256_setr_epi32( 0  , 0, 0  , 0, b0 , 0, b1 , 0 );
	b_12_15 = _mm256_setr_epi32( b2 , 0, b3 , 0, b4 , 0, 0  , 0 );
	//   - multiply
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column12:
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a11, 0, a11, 0, a11, 0, a11, 0 );
	//   for i=11 to 14: r[i] += a[11] * b[i-11]
	b_8_11  = _mm256_setr_epi32( 0  , 0, 0  , 0, 0  , 0, b0 , 0 );
	b_12_15 = _mm256_setr_epi32( b1 , 0, b2 , 0, b3 , 0, 0  , 0 );
	//   - multiply
	tmp_8_11  = _mm256_mul_epi32( all_vals, b_8_11  );
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column13:
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a12, 0, a12, 0, a12, 0, a12, 0 );
	//   for i=12 to 14: r[i] += a[12] * b[i-12]
	b_12_15 = _mm256_setr_epi32( b0 , 0, b1 , 0, b2 , 0, 0  , 0 );
	//   - multiply
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column14:
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( a13, 0, a13, 0, a13, 0, a13, 0 );
	//   for i=13 to 14: r[i] += a[13] * b[i-13]
	b_12_15 = _mm256_setr_epi32( 0  , 0, b0 , 0, b1 , 0, 0  , 0 );
	//   - multiply
	tmp_12_15 = _mm256_mul_epi32( all_vals, b_12_15 );
	//   - add to result
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// store result
    _mm256_storeu_si256 ((__m256i *) &(r[0] ), r_0_3   );
    _mm256_storeu_si256 ((__m256i *) &(r[4] ), r_4_7   );
    _mm256_storeu_si256 ((__m256i *) &(r[8] ), r_8_11  );
    _mm256_storeu_si256 ((__m256i *) &(r[12]), r_12_15 );

    // column15:
    r[14] += (signed long) a14 * b0;
    
	// BOTTOM HALF
	// column1:
	//   for i=15 to 28: r[i] = a[i-14] * b[14]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b14, 0, b14, 0, b14, 0, b14, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a1 , 0, a2 , 0, a3 , 0, a4 , 0 );
	a_4_7   = _mm256_setr_epi32( a5 , 0, a6 , 0, a7 , 0, a8 , 0 );
	a_8_11  = _mm256_setr_epi32( a9 , 0, a10, 0, a11, 0, a12, 0 );
	a_12_15 = _mm256_setr_epi32( a13, 0, a14, 0, 0  , 0, 0  , 0 );
    //   - multiply
	r_0_3   = _mm256_mul_epi32( a_0_3  , all_vals );
	r_4_7   = _mm256_mul_epi32( a_4_7  , all_vals );
	r_8_11  = _mm256_mul_epi32( a_8_11 , all_vals );
	r_12_15 = _mm256_mul_epi32( a_12_15, all_vals );
	// column2:
	//   for i=15 to 28: r[i] = a[i-13] * b[13]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b13, 0, b13, 0, b13, 0, b13, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a2 , 0, a3 , 0, a4 , 0, a5 , 0 );
	a_4_7   = _mm256_setr_epi32( a6 , 0, a7 , 0, a8 , 0, a9 , 0 );
	a_8_11  = _mm256_setr_epi32( a10, 0, a11, 0, a12, 0, a13, 0 );
	a_12_15 = _mm256_setr_epi32( a14, 0, 0  , 0, 0  , 0, 0  , 0 );
    //   - multiply
	tmp_0_3   = _mm256_mul_epi32( a_0_3  , all_vals );
	tmp_4_7   = _mm256_mul_epi32( a_4_7  , all_vals );
	tmp_8_11  = _mm256_mul_epi32( a_8_11 , all_vals );
	tmp_12_15 = _mm256_mul_epi32( a_12_15, all_vals );
	//   - add to result
	r_0_3   = _mm256_add_epi64( r_0_3  , tmp_0_3   );
	r_4_7   = _mm256_add_epi64( r_4_7  , tmp_4_7   );
	r_8_11  = _mm256_add_epi64( r_8_11 , tmp_8_11  );
	r_12_15 = _mm256_add_epi64( r_12_15, tmp_12_15 );
	// column3:
	//   for i=15 to 28: r[i] = a[i-12] * b[12]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b12, 0, b12, 0, b12, 0, b12, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a3 , 0, a4 , 0, a5 , 0, a6 , 0 );
	a_4_7   = _mm256_setr_epi32( a7 , 0, a8 , 0, a9 , 0, a10, 0 );
	a_8_11  = _mm256_setr_epi32( a11, 0, a12, 0, a13, 0, a14, 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3  , all_vals );
	tmp_4_7  = _mm256_mul_epi32( a_4_7  , all_vals );
	tmp_8_11 = _mm256_mul_epi32( a_8_11 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3  , tmp_0_3  );
	r_4_7  = _mm256_add_epi64( r_4_7  , tmp_4_7  );
	r_8_11 = _mm256_add_epi64( r_8_11 , tmp_8_11 );
	// column4:
	//   for i=15 to 28: r[i] = a[i-11] * b[11]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b11, 0, b11, 0, b11, 0, b11, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a4 , 0, a5 , 0, a6 , 0, a7 , 0 );
	a_4_7   = _mm256_setr_epi32( a8 , 0, a9 , 0, a10, 0, a11, 0 );
	a_8_11  = _mm256_setr_epi32( a12, 0, a13, 0, a14, 0, 0  , 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3  , all_vals );
	tmp_4_7  = _mm256_mul_epi32( a_4_7  , all_vals );
	tmp_8_11 = _mm256_mul_epi32( a_8_11 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3  , tmp_0_3  );
	r_4_7  = _mm256_add_epi64( r_4_7  , tmp_4_7  );
	r_8_11 = _mm256_add_epi64( r_8_11 , tmp_8_11 );
	// column5:
	//   for i=15 to 28: r[i] = a[i-10] * b[10]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b10, 0, b10, 0, b10, 0, b10, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a5 , 0, a6 , 0, a7 , 0, a8 , 0 );
	a_4_7   = _mm256_setr_epi32( a9 , 0, a10, 0, a11, 0, a12, 0 );
	a_8_11  = _mm256_setr_epi32( a13, 0, a14, 0, 0  , 0, 0  , 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3  , all_vals );
	tmp_4_7  = _mm256_mul_epi32( a_4_7  , all_vals );
	tmp_8_11 = _mm256_mul_epi32( a_8_11 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3  , tmp_0_3  );
	r_4_7  = _mm256_add_epi64( r_4_7  , tmp_4_7  );
	r_8_11 = _mm256_add_epi64( r_8_11 , tmp_8_11 );
	// column6:
	//   for i=15 to 28: r[i] = a[i-9] * b[9]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b9, 0, b9, 0, b9, 0, b9, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a6 , 0, a7 , 0, a8 , 0, a9 , 0 );
	a_4_7   = _mm256_setr_epi32( a10, 0, a11, 0, a12, 0, a13, 0 );
	a_8_11  = _mm256_setr_epi32( a14, 0, 0  , 0, 0  , 0, 0  , 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3  , all_vals );
	tmp_4_7  = _mm256_mul_epi32( a_4_7  , all_vals );
	tmp_8_11 = _mm256_mul_epi32( a_8_11 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3  , tmp_0_3  );
	r_4_7  = _mm256_add_epi64( r_4_7  , tmp_4_7  );
	r_8_11 = _mm256_add_epi64( r_8_11 , tmp_8_11 );
	// column7:
	//   for i=15 to 28: r[i] = a[i-8] * b[8]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b8, 0, b8, 0, b8, 0, b8, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a7 , 0, a8 , 0, a9 , 0, a10, 0 );
	a_4_7   = _mm256_setr_epi32( a11, 0, a12, 0, a13, 0, a14, 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3 , all_vals );
	tmp_4_7  = _mm256_mul_epi32( a_4_7 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3 , tmp_0_3  );
	r_4_7  = _mm256_add_epi64( r_4_7 , tmp_4_7  );
	// column8:
	//   for i=15 to 28: r[i] = a[i-7] * b[7]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b7, 0, b7, 0, b7, 0, b7, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a8 , 0, a9 , 0, a10, 0, a11, 0 );
	a_4_7   = _mm256_setr_epi32( a12, 0, a13, 0, a14, 0, 0  , 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3 , all_vals );
	tmp_4_7  = _mm256_mul_epi32( a_4_7 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3 , tmp_0_3  );
	r_4_7  = _mm256_add_epi64( r_4_7 , tmp_4_7  );
	// column9:
	//   for i=15 to 28: r[i] = a[i-6] * b[6]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b6, 0, b6, 0, b6, 0, b6, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a9 , 0, a10, 0, a11, 0, a12, 0 );
	a_4_7   = _mm256_setr_epi32( a13, 0, a14, 0, 0  , 0, 0  , 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3 , all_vals );
	tmp_4_7  = _mm256_mul_epi32( a_4_7 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3 , tmp_0_3  );
	r_4_7  = _mm256_add_epi64( r_4_7 , tmp_4_7  );
	// column10:
	//   for i=15 to 28: r[i] = a[i-5] * b[5]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b5, 0, b5, 0, b5, 0, b5, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a10, 0, a11, 0, a12, 0, a13, 0 );
	a_4_7   = _mm256_setr_epi32( a14, 0, 0  , 0, 0  , 0, 0  , 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3 , all_vals );
	tmp_4_7  = _mm256_mul_epi32( a_4_7 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3 , tmp_0_3  );
	r_4_7  = _mm256_add_epi64( r_4_7 , tmp_4_7  );
	// column11:
	//   for i=15 to 28: r[i] = a[i-4] * b[4]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b4, 0, b4, 0, b4, 0, b4, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a11, 0, a12, 0, a13, 0, a14, 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3 , tmp_0_3  );
	// column12:
	//   for i=15 to 28: r[i] = a[i-3] * b[3]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b3, 0, b3, 0, b3, 0, b3, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a12, 0, a13, 0, a14, 0, 0, 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3 , tmp_0_3  );
	// column13:
	//   for i=15 to 28: r[i] = a[i-2] * b[2]
    //   - set all_vals
    all_vals  = _mm256_setr_epi32( b2, 0, b2, 0, b2, 0, b2, 0 );
	//   - set a
	a_0_3   = _mm256_setr_epi32( a13, 0, a14, 0, 0, 0, 0, 0 );
    //   - multiply
	tmp_0_3  = _mm256_mul_epi32( a_0_3 , all_vals );
	//   - add to result
	r_0_3  = _mm256_add_epi64( r_0_3 , tmp_0_3  );
	// store result
    _mm256_storeu_si256 ((__m256i *) &(r[15]), r_0_3   );
    _mm256_storeu_si256 ((__m256i *) &(r[19]), r_4_7   );
    _mm256_storeu_si256 ((__m256i *) &(r[23]), r_8_11  );
    _mm256_storeu_si256 ((__m256i *) &(r[27]), r_12_15 );
	
    // column14:
    r[15] += (signed long) a14 * b1;

    reduce_mul_512(r);

    reduce_mul_256(r);
    
    for(int i = 0; i<NO_LIMBS; i++){
        res->d[i] = r[i];
    }

    
    reduce_256(res);
}

int cmp_coord_256(radix17_t *x1, radix17_t *x2)
{
    //return 1 if x1 == x2

    int bool;
    bool = (x1->d[0] == x2->d[0]);

    for (size_t i = 1; i < NO_LIMBS; i++)
    {
        bool = bool && (x1->d[i] == x2->d[i]);
    }
    return bool;
}
int cmp_point_256(radix17_t *x1, radix17_t *y1, radix17_t *x2, radix17_t *y2)
{
    //return 1 if x1 == y1 and x2 == y2

    return cmp_coord_256(x1, y1) && cmp_coord_256(x2, y2);
}
int iszero_256(radix17_t *x)
{
    //return 1 if x == y == 0

    int bool;
    bool = (x->d[0] == 0);

    for (size_t i = 1; i < NO_LIMBS; i++)
    {
        bool = bool && (x->d[i] == 0);
    }
    return bool;
}
void set_ui_256(radix17_t *x, signed int ui)
{
    //set radix17_t to a unsigned int ui (lowest limb only, others 0)
    // maximum 17 bit numbers !!!

    x->d[0] = (signed int)ui;

    for (size_t i = 1; i < NO_LIMBS; i++)
    {
        x->d[i] = (signed int)0;
    }
}

void square_256(radix17_t *res, radix17_t *x)
{

    unsigned long r[30];

    signed int a0 = x->d[0];
    signed int a1 = x->d[1];
    signed int a2 = x->d[2];
    signed int a3 = x->d[3];
    signed int a4 = x->d[4];
    signed int a5 = x->d[5];
    signed int a6 = x->d[6];
    signed int a7 = x->d[7];
    signed int a8 = x->d[8];
    signed int a9 = x->d[9];
    signed int a10 = x->d[10];
    signed int a11 = x->d[11];
    signed int a12 = x->d[12];
    signed int a13 = x->d[13];
    signed int a14 = x->d[14];

    signed int _2a0 = a0 << 1;   // why does this work?? result could be > 32 bits (int) Update: works because of radix17
    signed int _2a1 = a1 << 1;
    signed int _2a2 = a2 << 1;
    signed int _2a3 = a3 << 1;
    signed int _2a4 = a4 << 1;
    signed int _2a5 = a5 << 1;
    signed int _2a6 = a6 << 1;
    signed int _2a7 = a7 << 1;
    signed int _2a8 = a8 << 1;
    signed int _2a9 = a9 << 1;
    signed int _2a10 = a10 << 1;
    signed int _2a11 = a11 << 1;
    signed int _2a12 = a12 << 1;
    signed int _2a13 = a13 << 1;

    // result variables
	__m256i r_0_3, r_4_7, r_8_11, r_12_15, r_16_19, r_20_23, r_24_27, r_28_31;
	
	__m256i col1_left, col1_right;
    __m256i tmp1, tmp2;

    //column 1
    col1_left  = _mm256_setr_epi32( _2a0, 0, _2a0, 0, _2a0, 0, _2a0, 0 );

    col1_right = _mm256_setr_epi32(0, 0, a1, 0, a2, 0, a3, 0);
    r_0_3 = _mm256_mul_epi32(col1_left, col1_right);
    //r_0_3 = _mm256_slli_epi64(r_0_3, 1);    // why does this not work? lshift(1) == *2 so it should not make a difference at what point we shift
    
    col1_right = _mm256_setr_epi32(a4, 0, a5, 0, a6, 0, a7, 0);
    r_4_7 = _mm256_mul_epi32(col1_left, col1_right);
    //r_4_7 = _mm256_slli_epi64(r_4_7, 1);

    col1_right = _mm256_setr_epi32(a8, 0, a9, 0, a10, 0, a11, 0);
    r_8_11 = _mm256_mul_epi32(col1_left, col1_right);
    //r_8_11 = _mm256_slli_epi64(r_8_11, 1);
    
    col1_left  = _mm256_setr_epi32( _2a0, 0, _2a0, 0, _2a0, 0, _2a1, 0 );
    col1_right = _mm256_setr_epi32(a12, 0, a13, 0, a14, 0, a14, 0);
    r_12_15 = _mm256_mul_epi32(col1_left, col1_right);
    //r_12_15 = _mm256_slli_epi64(r_12_15, 1);

    col1_left  = _mm256_setr_epi32( _2a2, 0, _2a3, 0, _2a4, 0, _2a5, 0 );
    col1_right = _mm256_setr_epi32(a14, 0, a14, 0, a14, 0, a14, 0);
    r_16_19 = _mm256_mul_epi32(col1_left, col1_right);
    //r_16_19 = _mm256_slli_epi64(r_16_19, 1);

    col1_left  = _mm256_setr_epi32( _2a6, 0, _2a7, 0, _2a8, 0, _2a9, 0 );
    col1_right = _mm256_setr_epi32(a14, 0, a14, 0, a14, 0, a14, 0);
    r_20_23 = _mm256_mul_epi32(col1_left, col1_right);
    //r_20_23 = _mm256_slli_epi64(r_20_23, 1);

    col1_left  = _mm256_setr_epi32( _2a10, 0, _2a11, 0, _2a12, 0, _2a13, 0 );
    col1_right = _mm256_setr_epi32(a14, 0, a14, 0, a14, 0, a14, 0);
    r_24_27 = _mm256_mul_epi32(col1_left, col1_right);
    //r_24_27 = _mm256_slli_epi64(r_24_27, 1);

    //column 2
    col1_left  = _mm256_setr_epi32( 0, 0, 0, 0, 0, 0, _2a1, 0 );
    col1_right = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, a2, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_0_3 = _mm256_add_epi64(r_0_3, tmp1);
    //r_0_3 = _mm256_slli_epi64(r_0_3, 1);

    col1_left  = _mm256_setr_epi32( _2a1, 0, _2a1, 0, _2a1, 0, _2a1, 0 );
    col1_right = _mm256_setr_epi32(a3, 0, a4, 0, a5, 0, a6, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_4_7 = _mm256_add_epi64(r_4_7, tmp1);
    //r_4_7 = _mm256_slli_epi64(r_4_7, 1);


    col1_left  = _mm256_setr_epi32( _2a1, 0, _2a1, 0, _2a1, 0, _2a1, 0 );
    col1_right = _mm256_setr_epi32(a7, 0, a8, 0, a9, 0, a10, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_8_11 = _mm256_add_epi64(r_8_11, tmp1);
    //r_8_11 = _mm256_slli_epi64(r_8_11, 1);

    col1_left  = _mm256_setr_epi32( _2a1, 0, _2a1, 0, _2a1, 0, _2a2, 0 );
    col1_right = _mm256_setr_epi32(a11, 0, a12, 0, a13, 0, a13, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_12_15 = _mm256_add_epi64(r_12_15, tmp1);
    //r_12_15 = _mm256_slli_epi64(r_12_15, 1);

    col1_left  = _mm256_setr_epi32( _2a3, 0, _2a4, 0, _2a5, 0, _2a6, 0 );
    col1_right = _mm256_setr_epi32(a13, 0, a13, 0, a13, 0, a13, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_16_19 = _mm256_add_epi64(r_16_19, tmp1);
    //r_16_19 = _mm256_slli_epi64(r_16_19, 1);

    col1_left  = _mm256_setr_epi32( _2a7, 0, _2a8, 0, _2a9, 0, _2a10, 0 );
    col1_right = _mm256_setr_epi32(a13, 0, a13, 0, a13, 0, a13, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_20_23 = _mm256_add_epi64(r_20_23, tmp1);
    //r_20_23 = _mm256_slli_epi64(r_20_23, 1);

    col1_left  = _mm256_setr_epi32( _2a11, 0, _2a12, 0, 0, 0, 0, 0 );
    col1_right = _mm256_setr_epi32(a13, 0, a13, 0, 0, 0, 0, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_24_27 = _mm256_add_epi64(r_24_27, tmp1);
    //r_24_27 = _mm256_slli_epi64(r_24_27, 1);


    //column 3
    col1_left  = _mm256_setr_epi32( 0, 0, _2a2, 0, _2a2, 0, _2a2, 0 );
    col1_right = _mm256_setr_epi32(0, 0, a3, 0, a4, 0, a5, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_4_7 = _mm256_add_epi64(r_4_7, tmp1);
    //r_4_7 = _mm256_slli_epi64(r_4_7, 1);

    col1_left  = _mm256_setr_epi32(_2a2, 0, _2a2, 0, _2a2, 0, _2a2, 0 );
    col1_right = _mm256_setr_epi32(a6, 0, a7, 0, a8, 0, a9, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_8_11 = _mm256_add_epi64(r_8_11, tmp1);
    //r_8_11 = _mm256_slli_epi64(r_8_11, 1);

    col1_left  = _mm256_setr_epi32(_2a2, 0, _2a2, 0, _2a2, 0, _2a3, 0 );
    col1_right = _mm256_setr_epi32(a10, 0, a11, 0, a12, 0, a12, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_12_15 = _mm256_add_epi64(r_12_15, tmp1);
    //r_12_15 = _mm256_slli_epi64(r_12_15, 1);

    col1_left  = _mm256_setr_epi32(_2a4, 0, _2a5, 0, _2a6, 0, _2a7, 0 );
    col1_right = _mm256_setr_epi32(a12, 0, a12, 0, a12, 0, a12, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_16_19 = _mm256_add_epi64(r_16_19, tmp1);
    //r_16_19 = _mm256_slli_epi64(r_16_19, 1);

    col1_left  = _mm256_setr_epi32(_2a8, 0, _2a9, 0, _2a10, 0, _2a11, 0 );
    col1_right = _mm256_setr_epi32(a12, 0, a12, 0, a12, 0, a12, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_20_23 = _mm256_add_epi64(r_20_23, tmp1);
    //r_20_23 = _mm256_slli_epi64(r_20_23, 1);

    //column 4 
    col1_left  = _mm256_setr_epi32( 0, 0, 0, 0, 0, 0, _2a3, 0 );
    col1_right = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, a4, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_4_7 = _mm256_add_epi64(r_4_7, tmp1);
    //r_4_7 = _mm256_slli_epi64(r_4_7, 1);

    col1_left  = _mm256_setr_epi32(_2a3, 0, _2a3, 0, _2a3, 0, _2a3, 0 );
    col1_right = _mm256_setr_epi32(a5, 0, a6, 0, a7, 0, a8, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_8_11 = _mm256_add_epi64(r_8_11, tmp1);
    //r_8_11 = _mm256_slli_epi64(r_8_11, 1);

    col1_left  = _mm256_setr_epi32(_2a3, 0, _2a3, 0, _2a3, 0, _2a4, 0 );
    col1_right = _mm256_setr_epi32(a9, 0, a10, 0, a11, 0, a11, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_12_15 = _mm256_add_epi64(r_12_15, tmp1);
    //r_12_15 = _mm256_slli_epi64(r_12_15, 1);

    col1_left  = _mm256_setr_epi32(_2a5, 0, _2a6, 0, _2a7, 0, _2a8, 0 );
    col1_right = _mm256_setr_epi32(a11, 0, a11, 0, a11, 0, a11, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_16_19 = _mm256_add_epi64(r_16_19, tmp1);
    //r_16_19 = _mm256_slli_epi64(r_16_19, 1);

    col1_left  = _mm256_setr_epi32(_2a9, 0, _2a10, 0, 0, 0, 0, 0 );
    col1_right = _mm256_setr_epi32(a11, 0, a11, 0, 0, 0, 0, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_20_23 = _mm256_add_epi64(r_20_23, tmp1);
    //r_20_23 = _mm256_slli_epi64(r_20_23, 1);

    //column 5
    col1_left  = _mm256_setr_epi32(0, 0, _2a4, 0, _2a4, 0, _2a4, 0 );
    col1_right = _mm256_setr_epi32(0, 0, a5, 0, a6, 0, a7, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_8_11 = _mm256_add_epi64(r_8_11, tmp1);
    //r_8_11 = _mm256_slli_epi64(r_8_11, 1);

    col1_left  = _mm256_setr_epi32(_2a4, 0, _2a4, 0, _2a4, 0, _2a5, 0 );
    col1_right = _mm256_setr_epi32(a8, 0, a9, 0, a10, 0, a10, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_12_15 = _mm256_add_epi64(r_12_15, tmp1);
    //r_12_15 = _mm256_slli_epi64(r_12_15, 1);

    col1_left  = _mm256_setr_epi32(_2a6, 0, _2a7, 0, _2a8, 0, _2a9, 0 );
    col1_right = _mm256_setr_epi32(a10, 0, a10, 0, a10, 0, a10, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_16_19 = _mm256_add_epi64(r_16_19, tmp1);
    //r_16_19 = _mm256_slli_epi64(r_16_19, 1);

    //column 6
    col1_left  = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, _2a5, 0 );
    col1_right = _mm256_setr_epi32(0, 0, 0, 0, 0, 0, a6, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_8_11 = _mm256_add_epi64(r_8_11, tmp1);
    //r_8_11 = _mm256_slli_epi64(r_8_11, 1);

    col1_left  = _mm256_setr_epi32(_2a5, 0, _2a5, 0, _2a5, 0, _2a6, 0 );
    col1_right = _mm256_setr_epi32(a7, 0, a8, 0, a9, 0, a9, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_12_15 = _mm256_add_epi64(r_12_15, tmp1);
    //r_12_15 = _mm256_slli_epi64(r_12_15, 1);

    col1_left  = _mm256_setr_epi32(_2a7, 0, _2a8, 0, 0, 0, 0, 0 );
    col1_right = _mm256_setr_epi32(a9, 0, a9, 0, 0, 0, 0, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_16_19 = _mm256_add_epi64(r_16_19, tmp1);
    //r_16_19 = _mm256_slli_epi64(r_16_19, 1);

    //column 7
    col1_left  = _mm256_setr_epi32(0, 0, _2a6, 0, _2a6, 0, _2a7, 0 );
    col1_right = _mm256_setr_epi32(0, 0, a7, 0, a8, 0, a8, 0);
    tmp1 = _mm256_mul_epi32(col1_left, col1_right);
    r_12_15 = _mm256_add_epi64(r_12_15, tmp1);
    //r_12_15 = _mm256_slli_epi64(r_12_15, 1);

    //column 8 (squares)
    tmp1 = _mm256_setr_epi32(a0, 0, 0, 0, a1, 0, 0, 0 );
    tmp2 = _mm256_mul_epi32(tmp1, tmp1);
    r_0_3 = _mm256_add_epi64(r_0_3, tmp2);

    tmp1 = _mm256_setr_epi32(a2, 0, 0, 0, a3, 0, 0, 0 );
    tmp2 = _mm256_mul_epi32(tmp1, tmp1);
    r_4_7 = _mm256_add_epi64(r_4_7, tmp2);

    tmp1 = _mm256_setr_epi32(a4, 0, 0, 0, a5, 0, 0, 0 );
    tmp2 = _mm256_mul_epi32(tmp1, tmp1);
    r_8_11 = _mm256_add_epi64(r_8_11, tmp2);

    tmp1 = _mm256_setr_epi32(a6, 0, 0, 0, a7, 0, 0, 0 );
    tmp2 = _mm256_mul_epi32(tmp1, tmp1);
    r_12_15 = _mm256_add_epi64(r_12_15, tmp2);

    tmp1 = _mm256_setr_epi32(a8, 0, 0, 0, a9, 0, 0, 0 );
    tmp2 = _mm256_mul_epi32(tmp1, tmp1);
    r_16_19 = _mm256_add_epi64(r_16_19, tmp2);

    tmp1 = _mm256_setr_epi32(a10, 0, 0, 0, a11, 0, 0, 0 );
    tmp2 = _mm256_mul_epi32(tmp1, tmp1);
    r_20_23 = _mm256_add_epi64(r_20_23, tmp2);

    tmp1 = _mm256_setr_epi32(a12, 0, 0, 0, a13, 0, 0, 0 );
    tmp2 = _mm256_mul_epi32(tmp1, tmp1);
    r_24_27 = _mm256_add_epi64(r_24_27, tmp2);

    r[28] = (signed long) a14*a14;

    _mm256_storeu_si256 ((__m256i *) r, r_0_3   );
    _mm256_storeu_si256 ((__m256i *) (r+4), r_4_7   );
    _mm256_storeu_si256 ((__m256i *) (r+8), r_8_11  );
    _mm256_storeu_si256 ((__m256i *) (r+12), r_12_15 );
    _mm256_storeu_si256 ((__m256i *) (r+16), r_16_19   );
    _mm256_storeu_si256 ((__m256i *) (r+20), r_20_23   );
    _mm256_storeu_si256 ((__m256i *) (r+24), r_24_27  );
    

/*
    r[0] =                                                                                     a0*a0;
    r[1] = _2a0*a1;
    r[2] = _2a0*a2                                                                            + a1*a1;
    r[3] = _2a0*a3    + _2a1*a2;
    
    r[4] = _2a0*a4    + _2a1*a3                                                               + a2*a2;
    r[5] = _2a0*a5    + _2a1*a4   + _2a2*a3;
    r[6] = _2a0*a6    + _2a1*a5   + _2a2*a4                                                   + a3*a3;
    r[7] = _2a0*a7    + _2a1*a6   + _2a2*a5     + _2a3*a4;
    
    r[8] = _2a0*a8    + _2a1*a7   + _2a2*a6     + _2a3*a5                                     + a4*a4;
    r[9] = _2a0*a9    + _2a1*a8   + _2a2*a7     + _2a3*a6   + _2a4*a5;
    r[10] = _2a0*a10  + _2a1*a9   + _2a2*a8     + _2a3*a7   + _2a4*a6                         + a5*a5;
    r[11] = _2a0*a11  + _2a1*a10  + _2a2*a9     + _2a3*a8   + _2a4*a7   + _2a5*a6;
    
    r[12] = _2a0*a12  + _2a1*a11  + _2a2*a10    + _2a3*a9   + _2a4*a8   + _2a5*a7             + a6*a6;
    r[13] = _2a0*a13  + _2a1*a12  + _2a2*a11    + _2a3*a10  + _2a4*a9   + _2a5*a8  + _2a6*a7;
    r[14] = _2a0*a14  + _2a1*a13  + _2a2*a12    + _2a3*a11  + _2a4*a10  + _2a5*a9  + _2a6*a8  + a7*a7;
    r[15] = _2a1*a14  + _2a2*a13  + _2a3*a12    + _2a4*a11  + _2a5*a10  + _2a6*a9  + _2a7*a8;
    
    r[16] = _2a2*a14  + _2a3*a13  + _2a4*a12    + _2a5*a11  + _2a6*a10  + _2a7*a9             + a8*a8;
    r[17] = _2a3*a14  + _2a4*a13  + _2a5*a12    + _2a6*a11  + _2a7*a10  + _2a8*a9;
    r[18] = _2a4*a14  + _2a5*a13  + _2a6*a12    + _2a7*a11  + _2a8*a10                       + a9*a9;
    r[19] = _2a5*a14  + _2a6*a13  + _2a7*a12    + _2a8*a11  + _2a9*a10;
    
    r[20] = _2a6*a14  + _2a7*a13  + _2a8*a12    + _2a9*a11                                   + a10*a10;
    r[21] = _2a7*a14  + _2a8*a13  + _2a9*a12    + _2a10*a11;
    r[22] = _2a8*a14  + _2a9*a13  + _2a10*a12                                                + a11*a11;
    r[23] = _2a9*a14  + _2a10*a13 + _2a11*a12;
    
    r[24] = _2a10*a14 + _2a11*a13                                                            + a12*a12;
    r[25] = _2a11*a14 + _2a12*a13;
    r[26] = _2a12*a14                                                                       +  a13*a13;
    r[27] = _2a13*a14;
    
    r[28] =                                                                                   a14*a14;
*/

    reduce_mul_512(r);

    reduce_mul_256(r);
    
     for(int i = 0; i<NO_LIMBS; i++){
        res->d[i] = r[i];
    }

    reduce_256(res);
}
void inverse_256(radix17_t *res, radix17_t *x)
{
    radix17_t z2, z9, z11, z2_5_0, z2_10_0, z2_20_0, z2_50_0, z2_100_0, t;
    int i;
    /* 2 */ square_256(&z2, x);
    /* 4 */ square_256(&t, &z2);
    /* 8 */ square_256(&t, &t);
    /* 9 */ mul_256(&z9, &t, x);
    /* 11 */ mul_256(&z11, &z9, &z2);
    /* 22 */ square_256(&t, &z11);
    /* 2^5 - 2^0 = 31 */ mul_256(&z2_5_0, &t, &z9);
    /* 2^6 - 2^1 */ square_256(&t, &z2_5_0);
    /* 2^10 - 2^5 */ for (i = 1; i < 5; i++)
    {
        square_256(&t, &t);
    }
    /* 2^10 - 2^0 */ mul_256(&z2_10_0, &t, &z2_5_0);
    /* 2^11 - 2^1 */ square_256(&t, &z2_10_0);
    /* 2^20 - 2^10 */ for (i = 1; i < 10; i++)
    {
        square_256(&t, &t);
    }
    /* 2^20 - 2^0 */ mul_256(&z2_20_0, &t, &z2_10_0);
    /* 2^21 - 2^1 */ square_256(&t, &z2_20_0);
    /* 2^40 - 2^20 */ for (i = 1; i < 20; i++)
    {
        square_256(&t, &t);
    }
    /* 2^40 - 2^0 */ mul_256(&t, &t, &z2_20_0);
    /* 2^41 - 2^1 */ square_256(&t, &t);
    /* 2^50 - 2^10 */ for (i = 1; i < 10; i++)
    {
        square_256(&t, &t);
    }
    /* 2^50 - 2^0 */ mul_256(&z2_50_0, &t, &z2_10_0);
    /* 2^51 - 2^1 */ square_256(&t, &z2_50_0);
    /* 2^100 - 2^50 */ for (i = 1; i < 50; i++)
    {
        square_256(&t, &t);
    }
    /* 2^100 - 2^0 */ mul_256(&z2_100_0, &t, &z2_50_0);
    /* 2^101 - 2^1 */ square_256(&t, &z2_100_0);
    /* 2^200 - 2^100 */ for (i = 1; i < 100; i++)
    {
        square_256(&t, &t);
    }
    /* 2^200 - 2^0 */ mul_256(&t, &t, &z2_100_0);
    /* 2^201 - 2^1 */ square_256(&t, &t);
    /* 2^250 - 2^50 */ for (i = 1; i < 50; i++)
    {
        square_256(&t, &t);
    }
    /* 2^250 - 2^0 */ mul_256(&t, &t, &z2_50_0);
    /* 2^251 - 2^1 */ square_256(&t, &t);
    /* 2^252 - 2^2 */ square_256(&t, &t);
    /* 2^253 - 2^3 */ square_256(&t, &t);
    /* 2^254 - 2^4 */ square_256(&t, &t);
    /* 2^255 - 2^5 */ square_256(&t, &t);
    /* 2^255 - 21 */ mul_256(res, &t, &z11);
    //print_radix17( res, "inv" );
}

void div_256(radix17_t *res, radix17_t *a, radix17_t *b)
{
    // inverse and mult
    radix17_t inv;
    inverse_256(&inv, b);
    mul_256(res, a, &inv);
}