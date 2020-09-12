#include <stdio.h>
#include "bignum.h"
#include <immintrin.h>
#include <emmintrin.h>
#include "../../_common_libs/radix17_vec.h"
#include "../../_common_libs/vec_utils.h"

const int NO_LIMBS = 15;
const int LIMB_SIZE = 17;

unsigned long times19(unsigned long a)
{
    /*
    multiply a by 19
    */
  return (a << 4) + (a << 1) + a;
}

void reduce_2_to_2(__m256i *xlo, __m256i *xhi)
{
    /*
    reduce xhi, xlo which might be overflowed by one limb
    */
    signed long carry;
    signed int d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15;

    int* i_v = (int*) xlo;
    int* i_v2 = (int*) xhi;
    d0 = i_v[0];
    d1 = i_v[1];
    d2 = i_v[2];
    d3 = i_v[3];
    d4 = i_v[4];
    d5 = i_v[5];
    d6 = i_v[6];
    d7 = i_v[7];
    d8 = i_v2[0];
    d9 = i_v2[1];
    d10 = i_v2[2];
    d11 = i_v2[3];
    d12 = i_v2[4];
    d13 = i_v2[5];
    d14 = i_v2[6];
    d15 = 0;
    
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
    d0 += times19(carry);
    carry <<= LIMB_SIZE;
    d14 -= carry;

    *xlo = _mm256_setr_epi32(d0, d1, d2, d3, d4, d5, d6, d7);
    *xhi = _mm256_setr_epi32(d8, d9, d10, d11, d12, d13, d14, d15); 

}

void reduce_4_to_2(__m256i *rv0, __m256i *rv1, __m256i *rv2, __m256i *rv3)
{
    /*
    reduce 4 vectors into 2
    */

    signed long carry;
    unsigned long r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14;
    long* i_v0 = (long*) rv0;
    long* i_v1 = (long*) rv1;
    long* i_v2 = (long*) rv2;
    long* i_v3 = (long*) rv3;
    
    r0 = i_v0[0];
    r1 = i_v0[1];
    r2 = i_v0[2];
    r3 = i_v0[3];
    r4 = i_v1[0];
    r5 = i_v1[1];
    r6 = i_v1[2];
    r7 = i_v1[3];
    r8 = i_v2[0];
    r9 = i_v2[1];
    r10 = i_v2[2];
    r11 = i_v2[3];
    r12 = i_v3[0];
    r13 = i_v3[1];
    r14 = i_v3[2];
    
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
    // write back the result
    *rv0 = _mm256_setr_epi32((int)r0,(int)r1,(int)r2,(int)r3,(int)r4,(int)r5,(int)r6,(int)r7);
    *rv1 = _mm256_setr_epi32((int)r8,(int)r9,(int)r10,(int)r11,(int)r12,(int)r13,(int)r14,0);
}

void reduce_8_to_4(__m256i *r0, __m256i *r1, __m256i *r2, __m256i *r3, __m256i *r4, __m256i *r5, __m256i *r6, __m256i *r7)
{
    /*
    reduce 8 vectors inot 4
    */
    __m256i r0_3, r4_7,r8_11, r12_15,
    r15_18, r19_22, r23_16, r27_30,
    r15_18x16, r19_22x16, r23_16x16, r27_30x16,
    r15_18x2, r19_22x2, r23_16x2, r27_30x2;
    r0_3 = *r0;
    r4_7 = *r1;
    r8_11 = *r2;
    r12_15 = *r3;
    r15_18 = *r4;
    r19_22 = *r5;
    r23_16 = *r6;
    r27_30 = *r7;
    
    // The following needs to be multiplyed with 19 --> 19x = 16x + 2x + x --> x<<4 + x<<1 +x
    r15_18 = _mm256_blend_epi32(*r4,*r3,192);
    r15_18 = _mm256_permute4x64_epi64(r15_18, 0b10010011);

    r19_22 = _mm256_blend_epi32(*r5,*r4,192);
    r19_22 = _mm256_permute4x64_epi64(r19_22, 0b10010011);

    r23_16 = _mm256_blend_epi32(*r6,*r5,192);
    r23_16 = _mm256_permute4x64_epi64(r23_16, 0b10010011);
    
    r27_30 = _mm256_blend_epi32(*r7,*r6,192);
    r27_30 = _mm256_permute4x64_epi64(r27_30, 0b10010011);
        
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
    *r0 = _mm256_add_epi64(r15_18, r0_3);
    *r1 = _mm256_add_epi64(r19_22, r4_7);
    *r2 = _mm256_add_epi64(r23_16, r8_11);
    *r3 = _mm256_add_epi64(r27_30, r12_15);
}


void add_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi)
{
    /*
    res = a + b
    */
    // a + b
    *reslo = _mm256_add_epi32(alo,blo);
    *reshi = _mm256_add_epi32(ahi,bhi);
    // reduce res
    //reduce_2_to_2(reslo, reshi); 
}

void sub_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi)
{
    /*
    res = a - b
    */
    
    // a - b
    *reslo = _mm256_sub_epi32(alo,blo);
    *reshi = _mm256_sub_epi32(ahi,bhi);
    int mask = _mm256_movemask_ps((__m256)*reshi);
    while(_mm256_movemask_ps((__m256)*reslo) != 0 || ( mask != 0 && mask != 0b10000000 )){ //reduce if there are negative numbers in the limbs only
        reduce_2_to_2(reslo, reshi); 
        mask = _mm256_movemask_ps((__m256)*reshi);
    }
    // reduce res
    // reduce_2_to_2(reslo, reshi);
    // reduce_2_to_2(reslo, reshi);
}

void mul_256(__m256i *reslo, __m256i *reshi, __m256i alo,__m256i ahi, __m256i blo,__m256i bhi)
{
    /*
    res = a * b
    */
	/*
	//       column1          column2          column3          column4          column5         column6       column7         column8        column9         column10       column11        column12       column13       column14      column15

    r[0]  =  a0  * b0;
    r[1]  =  a0 *  b1     +   a1 *  b0;
    r[2]  =  a0 *  b2     +   a1 *  b1     +   a2  * b0;
    r[3]  =  a0 *  b3     +   a1 *  b2     +   a2  * b1     +   a3 *  b0;
    r[4]  =  a0 *  b4     +   a1 *  b3     +   a2  * b2     +   a3 *  b1     +   a4  * b0;
    r[5]  =  a0 *  b5     +   a1 *  b4     +   a2  * b3     +   a3 *  b2     +   a4  * b1    +   a5 *  b0;
    r[6]  =  a0 *  b6     +   a1 *  b5     +   a2  * b4     +   a3 *  b3     +   a4  * b2    +   a5 *  b1   +   a6  * b0;
    r[7]  =  a0 *  b7     +   a1 *  b6     +   a2  * b5     +   a3 *  b4     +   a4  * b3    +   a5 *  b2   +   a6  * b1   +   a7 *  b0;
    r[8]  =  a0 *  b8     +   a1 *  b7     +   a2  * b6     +   a3 *  b5     +   a4  * b4    +   a5 *  b3   +   a6  * b2   +   a7 *  b1   +   a8  * b0;
    r[9]  =  a0 *  b9     +   a1 *  b8     +   a2  * b7     +   a3 *  b6     +   a4  * b5    +   a5 *  b4   +   a6  * b3   +   a7 *  b2   +   a8  * b1   +    a9 *  b0;
    r[10] =  a0 *  b10    +   a1 *  b9     +   a2  * b8     +   a3 *  b7     +   a4  * b6    +   a5 *  b5   +   a6  * b4   +   a7 *  b3   +   a8  * b2   +    a9 *  b1   +   a10 * b0;
    r[11] =  a0 *  b11    +   a1 *  b10    +   a2  * b9     +   a3 *  b8     +   a4  * b7    +   a5 *  b6   +   a6  * b5   +   a7 *  b4   +   a8  * b3   +    a9 *  b2   +   a10 * b1   +   a11 * b0;
    r[12] =  a0 *  b12    +   a1 *  b11    +   a2  * b10    +   a3 *  b9     +   a4  * b8    +   a5 *  b7   +   a6  * b6   +   a7 *  b5   +   a8  * b4   +    a9 *  b3   +   a10 * b2   +   a11 * b1   +   a12 * b0;
    r[13] =  a0 *  b13    +   a1 *  b12    +   a2  * b11    +   a3 *  b10    +   a4  * b9    +   a5 *  b8   +   a6  * b7   +   a7 *  b6   +   a8  * b5   +    a9 *  b4   +   a10 * b3   +   a11 * b2   +   a12 * b1   +   a13 * b0;
    r[14] =  a0 *  b14    +   a1 *  b13    +   a2  * b12    +   a3 *  b11    +   a4  * b10   +   a5 *  b9   +   a6  * b8   +   a7 *  b7   +   a8  * b6   +    a9 *  b5   +   a10 * b4   +   a11 * b3   +   a12 * b2   +   a13 * b1   +   a14 * b0;
    
    r[15] =               +   a1 * b14     +   a2  * b13    +   a3 *  b12    +   a4  * b11   +   a5 *  b10   +   a6 *  b9   +   a7 *  b8   +   a8 *  b7   +   a9 *  b6   +   a10 * b5   +   a11 * b4   +   a12 * b3   +   a13 * b2   +   a14 * b1;
    r[16] =                                +   a2  * b14    +   a3 *  b13    +   a4  * b12   +   a5 *  b11   +   a6 *  b10  +   a7 *  b9   +   a8 *  b8   +   a9 *  b7   +   a10 * b6   +   a11 * b5   +   a12 * b4   +   a13 * b3   +   a14 * b2;
    r[17] =                                                 +   a3 *  b14    +   a4  * b13   +   a5 *  b12   +   a6 *  b11  +   a7 *  b10  +   a8 *  b9   +   a9 *  b8   +   a10 * b7   +   a11 * b6   +   a12 * b5   +   a13 * b4   +   a14 * b3;
    r[18] =                                                                  +   a4  * b14   +   a5 *  b13   +   a6 *  b12  +   a7 *  b11  +   a8 *  b10  +   a9 *  b9   +   a10 * b8   +   a11 * b7   +   a12 * b6   +   a13 * b5   +   a14 * b4;
    r[19] =                                                                                  +   a5 *  b14   +   a6 *  b13  +   a7 *  b12  +   a8 *  b11  +   a9 *  b10  +   a10 * b9   +   a11 * b8   +   a12 * b7   +   a13 * b6   +   a14 * b5;
    r[20] =                                                                                                  +   a6 *  b14  +   a7 *  b13  +   a8 *  b12  +   a9 *  b11  +   a10 * b10  +   a11 * b9   +   a12 * b8   +   a13 * b7   +   a14 * b6;
    r[21] =                                                                                                                 +   a7 *  b14  +   a8 *  b13  +   a9 *  b12  +   a10 * b11  +   a11 * b10  +   a12 * b9   +   a13 * b8   +   a14 * b7;
    r[22] =                                                                                                                                +   a8 *  b14  +   a9 *  b13  +   a10 * b12  +   a11 * b11  +   a12 * b10  +   a13 * b9   +   a14 * b8;
    r[23] =                                                                                                                                               +   a9 *  b14  +   a10 * b13  +   a11 * b12  +   a12 * b11  +   a13 * b10  +   a14 * b9;
    r[24] =                                                                                                                                                              +   a10 * b14  +   a11 * b13  +   a12 * b12  +   a13 * b11  +   a14 * b10;
    r[25] =                                                                                                                                                                             +   a11 * b14  +   a12 * b13  +   a13 * b12  +   a14 * b11;
    r[26] =                                                                                                                                                                                            +   a12 * b14  +   a13 * b13  +   a14 * b12;
    r[27] =                                                                                                                                                                                                           +  a13 * b14    +   a14 * b13;
    r[28] =                                                                                                                                                                                                                            +  a14 * b14;                                                                                                

	*/
	
	__m256i b_0_3   = _mm256_permutevar8x32_epi32( blo, index1_vec );
	__m256i b_4_7   = _mm256_permutevar8x32_epi32( blo, index2_vec );
	__m256i b_8_11  = _mm256_permutevar8x32_epi32( bhi, index1_vec );
	__m256i b_12_15 = _mm256_permutevar8x32_epi32( bhi, index2_vec );
	
	__m256i a_0     = _mm256_broadcastq_epi64( _mm256_castsi256_si128( alo ) );
	__m256i a_1     = _mm256_permutevar8x32_epi32( a_0, index_alt_vec );
	__m256i a_2     = _mm256_broadcastq_epi64( _mm256_castsi256_si128( _mm256_unpackhi_epi64( alo, index1_vec ) ) );
	__m256i a_3     = _mm256_permutevar8x32_epi32( a_2, index_alt_vec );
	
	__m128i a_lo_hi = _mm256_extracti128_si256( alo, 1 );
	__m256i a_4     = _mm256_broadcastq_epi64( a_lo_hi );
	__m256i a_5     = _mm256_permutevar8x32_epi32( a_4, index_alt_vec );
	__m256i a_6     = _mm256_broadcastq_epi64( _mm256_castsi256_si128( _mm256_unpackhi_epi64( _mm256_castsi128_si256( a_lo_hi ), index1_vec ) ) );
	__m256i a_7     = _mm256_permutevar8x32_epi32( a_6 , index_alt_vec );

	__m256i a_8     = _mm256_broadcastq_epi64( _mm256_castsi256_si128( ahi ) );
	__m256i a_9     = _mm256_permutevar8x32_epi32( a_8, index_alt_vec );
	__m256i a_10    = _mm256_broadcastq_epi64( _mm256_castsi256_si128( _mm256_unpackhi_epi64( ahi, index1_vec ) ) );
	__m256i a_11    = _mm256_permutevar8x32_epi32( a_10, index_alt_vec );
	
	__m128i a_hi_hi = _mm256_extracti128_si256( ahi, 1 );
	__m256i a_12    = _mm256_broadcastq_epi64( a_hi_hi );
	__m256i a_13    = _mm256_permutevar8x32_epi32( a_12, index_alt_vec );
	__m256i a_14    = _mm256_broadcastq_epi64( _mm256_castsi256_si128( _mm256_unpackhi_epi64( _mm256_castsi128_si256( a_hi_hi ), index1_vec ) ) );
	__m256i a_15    = _mm256_permutevar8x32_epi32( a_14 , index_alt_vec );
	

    __m256i r_v0, r_v1, r_v2, r_v3, r_v4, r_v5, r_v6, r_v7;  // real result vectors

     __m256i c1_v0, c1_v1, c1_v2, c1_v3;  //column 1 results
     __m256i c2_v0, c2_v1, c2_v2, c2_v3;     
     __m256i c3_v0, c3_v1, c3_v2, c3_v3;
     __m256i c4_v0, c4_v1, c4_v2, c4_v3;
     __m256i c5_v0, c5_v1, c5_v2, c5_v3;
     __m256i c6_v0, c6_v1, c6_v2, c6_v3;
     __m256i c7_v0, c7_v1, c7_v2, c7_v3;
     __m256i c8_v0, c8_v1, c8_v2, c8_v3;
     __m256i c9_v0, c9_v1, c9_v2, c9_v3;
     __m256i c10_v0, c10_v1, c10_v2, c10_v3;
     __m256i c11_v0, c11_v1, c11_v2, c11_v3;
     __m256i c12_v0, c12_v1, c12_v2, c12_v3;
     __m256i c13_v0, c13_v1, c13_v2, c13_v3;
     __m256i c14_v0, c14_v1, c14_v2, c14_v3;
     __m256i c15_v0, c15_v1, c15_v2, c15_v3;

    c1_v0 = _mm256_mul_epi32(a_0, b_0_3);
    c1_v1 = _mm256_mul_epi32(a_0, b_4_7);
    c1_v2 = _mm256_mul_epi32(a_0, b_8_11);
    c1_v3 = _mm256_mul_epi32(a_0, b_12_15);
    c1_v3 = _mm256_blend_epi32(c1_v3, const_0_vec, 0b11000000);

    c2_v0 = _mm256_mul_epi32(a_1, b_0_3);
    c2_v1 = _mm256_mul_epi32(a_1, b_4_7);
    c2_v2 = _mm256_mul_epi32(a_1, b_8_11);
    c2_v3 = _mm256_mul_epi32(a_1, b_12_15);
    c2_v3 = _mm256_blend_epi32(c2_v3, const_0_vec, 0b11000000);

    c3_v0 = _mm256_mul_epi32(a_2, b_0_3);
    c3_v1 = _mm256_mul_epi32(a_2, b_4_7);
    c3_v2 = _mm256_mul_epi32(a_2, b_8_11);
    c3_v3 = _mm256_mul_epi32(a_2, b_12_15);
    c3_v3 = _mm256_blend_epi32(c3_v3, const_0_vec, 0b11000000);

    c4_v0 = _mm256_mul_epi32(a_3, b_0_3);
    c4_v1 = _mm256_mul_epi32(a_3, b_4_7);
    c4_v2 = _mm256_mul_epi32(a_3, b_8_11);
    c4_v3 = _mm256_mul_epi32(a_3, b_12_15);
    c4_v3 = _mm256_blend_epi32(c4_v3, const_0_vec, 0b11000000);

    c5_v0 = _mm256_mul_epi32(a_4, b_0_3);
    c5_v1 = _mm256_mul_epi32(a_4, b_4_7);
    c5_v2 = _mm256_mul_epi32(a_4, b_8_11);
    c5_v3 = _mm256_mul_epi32(a_4, b_12_15);
    c5_v3 = _mm256_blend_epi32(c5_v3, const_0_vec, 0b11000000);

    c6_v0 = _mm256_mul_epi32(a_5, b_0_3);
    c6_v1 = _mm256_mul_epi32(a_5, b_4_7);
    c6_v2 = _mm256_mul_epi32(a_5, b_8_11);
    c6_v3 = _mm256_mul_epi32(a_5, b_12_15);
    c6_v3 = _mm256_blend_epi32(c6_v3, const_0_vec, 0b11000000);

    c7_v0 = _mm256_mul_epi32(a_6, b_0_3);
    c7_v1 = _mm256_mul_epi32(a_6, b_4_7);
    c7_v2 = _mm256_mul_epi32(a_6, b_8_11);
    c7_v3 = _mm256_mul_epi32(a_6, b_12_15);
    c7_v3 = _mm256_blend_epi32(c7_v3, const_0_vec, 0b11000000);

    c8_v0 = _mm256_mul_epi32(a_7, b_0_3);
    c8_v1 = _mm256_mul_epi32(a_7, b_4_7);
    c8_v2 = _mm256_mul_epi32(a_7, b_8_11);
    c8_v3 = _mm256_mul_epi32(a_7, b_12_15);
    c8_v3 = _mm256_blend_epi32(c8_v3, const_0_vec, 0b11000000);

    c9_v0 = _mm256_mul_epi32(a_8, b_0_3);
    c9_v1 = _mm256_mul_epi32(a_8, b_4_7);
    c9_v2 = _mm256_mul_epi32(a_8, b_8_11);
    c9_v3 = _mm256_mul_epi32(a_8, b_12_15);
    c9_v3 = _mm256_blend_epi32(c9_v3, const_0_vec, 0b11000000);

    c10_v0 = _mm256_mul_epi32(a_9, b_0_3);
    c10_v1 = _mm256_mul_epi32(a_9, b_4_7);
    c10_v2 = _mm256_mul_epi32(a_9, b_8_11);
    c10_v3 = _mm256_mul_epi32(a_9, b_12_15);
    c10_v3 = _mm256_blend_epi32(c10_v3, const_0_vec, 0b11000000);

    c11_v0 = _mm256_mul_epi32(a_10, b_0_3);
    c11_v1 = _mm256_mul_epi32(a_10, b_4_7);
    c11_v2 = _mm256_mul_epi32(a_10, b_8_11);
    c11_v3 = _mm256_mul_epi32(a_10, b_12_15);
    c11_v3 = _mm256_blend_epi32(c11_v3, const_0_vec, 0b11000000);

    c12_v0 = _mm256_mul_epi32(a_11, b_0_3);
    c12_v1 = _mm256_mul_epi32(a_11, b_4_7);
    c12_v2 = _mm256_mul_epi32(a_11, b_8_11);
    c12_v3 = _mm256_mul_epi32(a_11, b_12_15);
    c12_v3 = _mm256_blend_epi32(c12_v3, const_0_vec, 0b11000000);

    c13_v0 = _mm256_mul_epi32(a_12, b_0_3);
    c13_v1 = _mm256_mul_epi32(a_12, b_4_7);
    c13_v2 = _mm256_mul_epi32(a_12, b_8_11);
    c13_v3 = _mm256_mul_epi32(a_12, b_12_15);
    c13_v3 = _mm256_blend_epi32(c13_v3, const_0_vec, 0b11000000);

    c14_v0 = _mm256_mul_epi32(a_13, b_0_3);
    c14_v1 = _mm256_mul_epi32(a_13, b_4_7);
    c14_v2 = _mm256_mul_epi32(a_13, b_8_11);
    c14_v3 = _mm256_mul_epi32(a_13, b_12_15);
    c14_v3 = _mm256_blend_epi32(c14_v3, const_0_vec, 0b11000000);

    c15_v0 = _mm256_mul_epi32(a_14, b_0_3);
    c15_v1 = _mm256_mul_epi32(a_14, b_4_7);
    c15_v2 = _mm256_mul_epi32(a_14, b_8_11);
    c15_v3 = _mm256_mul_epi32(a_14, b_12_15);
    c15_v3 = _mm256_blend_epi32(c15_v3, const_0_vec, 0b11000000);

    // offset 0 means:  0 1 2 3, 4 5 6 7, 8 9 10 11, ...
    __m256i offset0_v0, offset0_v1, offset0_v2, offset0_v3, offset0_v4, offset0_v5, offset0_v6, tmp_offset0;  // column 1 + 5 + 9 + 13

    offset0_v0 = c1_v0;
    offset0_v1 = _mm256_add_epi64(c1_v1, c5_v0);
    offset0_v2 = _mm256_add_epi64(c1_v2, c5_v1);
    offset0_v2 = _mm256_add_epi64(offset0_v2, c9_v0);
    offset0_v3 = _mm256_add_epi64(c1_v3, c5_v2);
    tmp_offset0 = _mm256_add_epi64(c9_v1, c13_v0);
    offset0_v3 = _mm256_add_epi64(offset0_v3, tmp_offset0);
    offset0_v4 = _mm256_add_epi64(c5_v3, c9_v2);
    offset0_v4 = _mm256_add_epi64(offset0_v4, c13_v1);
    offset0_v5 = _mm256_add_epi64(c9_v3, c13_v2);
    offset0_v6 = c13_v3;

    __m256i offset1_v0, offset1_v1, offset1_v2, offset1_v3, offset1_v4, offset1_v5, offset1_v6, tmp_offset1;  // column 2 + 6 + 10 + 14

    offset1_v0 = c2_v0;
    offset1_v1 = _mm256_add_epi64(c2_v1, c6_v0);
    offset1_v2 = _mm256_add_epi64(c2_v2, c6_v1);
    offset1_v2 = _mm256_add_epi64(offset1_v2, c10_v0);
    offset1_v3 = _mm256_add_epi64(c2_v3, c6_v2);
    tmp_offset1 = _mm256_add_epi64(c10_v1, c14_v0);
    offset1_v3 = _mm256_add_epi64(offset1_v3, tmp_offset1);
    offset1_v4 = _mm256_add_epi64(c6_v3, c10_v2);
    offset1_v4 = _mm256_add_epi64(offset1_v4, c14_v1);
    offset1_v5 = _mm256_add_epi64(c10_v3, c14_v2);
    offset1_v6 = c14_v3;

    __m256i offset2_v0, offset2_v1, offset2_v2, offset2_v3, offset2_v4, offset2_v5, offset2_v6, tmp_offset2;  // column 3 + 7 + 11 + 15

    offset2_v0 = c3_v0;
    offset2_v1 = _mm256_add_epi64(c3_v1, c7_v0);
    offset2_v2 = _mm256_add_epi64(c3_v2, c7_v1);
    offset2_v2 = _mm256_add_epi64(offset2_v2, c11_v0);
    offset2_v3 = _mm256_add_epi64(c3_v3, c7_v2);
    tmp_offset2 = _mm256_add_epi64(c11_v1, c15_v0);
    offset2_v3 = _mm256_add_epi64(offset2_v3, tmp_offset2);
    offset2_v4 = _mm256_add_epi64(c7_v3, c11_v2);
    offset2_v4 = _mm256_add_epi64(offset2_v4, c15_v1);
    offset2_v5 = _mm256_add_epi64(c11_v3, c15_v2);
    offset2_v6 = c15_v3;


    __m256i offset3_v0, offset3_v1, offset3_v2, offset3_v3, offset3_v4, offset3_v5;  // column 3 + 7 + 11 + 15

    offset3_v0 = c4_v0;
    offset3_v1 = _mm256_add_epi64(c4_v1, c8_v0);
    offset3_v2 = _mm256_add_epi64(c4_v2, c8_v1);
    offset3_v2 = _mm256_add_epi64(offset3_v2, c12_v0);
    offset3_v3 = _mm256_add_epi64(c4_v3, c8_v2);
    offset3_v3 = _mm256_add_epi64(offset3_v3, c12_v1);
    offset3_v4 = _mm256_add_epi64(c8_v3, c12_v2);
    offset3_v5 = c12_v3;

    //  results :      v0      v1        v2         v3            v4            v5          v6            v7
    //  results :   0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 19, 20 21 22 23, 24 25 26 27,  28 29 30 31
    
    // offset 0 :      v0      v1        v2         v3            v4            v5          v6        
    // offset 0 :   0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 19, 20 21 22 23, 24 25 26 -27-

    // offset 1 :      v0      v1        v2          v3            v4            v5          v6            
    // offset 1 :   1 2 3 4, 5 6 7 8, 9 10 11 12, 13 14 15 16, 17 18 19 20, 21 22 23 24, 25 26 27 -28-
    // permuted :   4 1 2 3, 8 5 6 7, 12 9 10 11, 16 13 14 15, 20 17 18 19, 24 21 22 23, -28- 25 26 27

    // offset 2 :      v0      v1        v2          v3            v4            v5          v6            
    // offset 2 :   2 3 4 5, 6 7 8 9, 10 11 12 13, 14 15 16 17, 18 19 20 21, 22 23 24 25, 26 27 28 -29-
    // permuted :   4 5 2 3, 8 9 6 7, 12 13 10 11, 16 17 14 15, 20 21 18 19, 24 25 22 23, 28 -29- 26 27

    // offset 3 :      v0      v1        v2          v3            v4            v5                      
    // offset 3 :   3 4 5 6, 7 8 9 10, 11 12 13 14, 15 16 17 18, 19 20 21 22, 23 24 25 -26-
    // permuted :   4 5 6 3, 8 9 10 7, 12 13 14 11, 16 17 18 15, 20 21 22 19, 24 25 -26- 23 

	// permute all offset vectors
	
	// offset 1
    __m256i offset1_v0_perm, offset1_v1_perm, offset1_v2_perm, offset1_v3_perm, offset1_v4_perm, offset1_v5_perm, offset1_v6_perm;
	offset1_v0_perm = _mm256_permute4x64_epi64( offset1_v0, 0b10010011 );
	offset1_v1_perm = _mm256_permute4x64_epi64( offset1_v1, 0b10010011 );
	offset1_v2_perm = _mm256_permute4x64_epi64( offset1_v2, 0b10010011 );
	offset1_v3_perm = _mm256_permute4x64_epi64( offset1_v3, 0b10010011 );
	offset1_v4_perm = _mm256_permute4x64_epi64( offset1_v4, 0b10010011 );
	offset1_v5_perm = _mm256_permute4x64_epi64( offset1_v5, 0b10010011 );
	offset1_v6_perm = _mm256_permute4x64_epi64( offset1_v6, 0b10010011 );

	// offset 2
    __m256i offset2_v0_perm, offset2_v1_perm, offset2_v2_perm, offset2_v3_perm, offset2_v4_perm, offset2_v5_perm, offset2_v6_perm;
	offset2_v0_perm = _mm256_permute4x64_epi64( offset2_v0, 0b01001110 );
	offset2_v1_perm = _mm256_permute4x64_epi64( offset2_v1, 0b01001110 );
	offset2_v2_perm = _mm256_permute4x64_epi64( offset2_v2, 0b01001110 );
	offset2_v3_perm = _mm256_permute4x64_epi64( offset2_v3, 0b01001110 );
	offset2_v4_perm = _mm256_permute4x64_epi64( offset2_v4, 0b01001110 );
	offset2_v5_perm = _mm256_permute4x64_epi64( offset2_v5, 0b01001110 );
	offset2_v6_perm = _mm256_permute4x64_epi64( offset2_v6, 0b01001110 );

	// offset 3
    __m256i offset3_v0_perm, offset3_v1_perm, offset3_v2_perm, offset3_v3_perm, offset3_v4_perm, offset3_v5_perm;
	offset3_v0_perm = _mm256_permute4x64_epi64( offset3_v0, 0b00111001 );
	offset3_v1_perm = _mm256_permute4x64_epi64( offset3_v1, 0b00111001 );
	offset3_v2_perm = _mm256_permute4x64_epi64( offset3_v2, 0b00111001 );
	offset3_v3_perm = _mm256_permute4x64_epi64( offset3_v3, 0b00111001 );
	offset3_v4_perm = _mm256_permute4x64_epi64( offset3_v4, 0b00111001 );
	offset3_v5_perm = _mm256_permute4x64_epi64( offset3_v5, 0b00111001 );
	
	// add offset 0 and offset 1
	r_v0 = _mm256_add_epi64( offset0_v0, _mm256_blend_epi32( const_0_vec    , offset1_v0_perm, 0b11111100 ) );
	r_v1 = _mm256_add_epi64( offset0_v1, _mm256_blend_epi32( offset1_v0_perm, offset1_v1_perm, 0b11111100 ) );
	r_v2 = _mm256_add_epi64( offset0_v2, _mm256_blend_epi32( offset1_v1_perm, offset1_v2_perm, 0b11111100 ) );
	r_v3 = _mm256_add_epi64( offset0_v3, _mm256_blend_epi32( offset1_v2_perm, offset1_v3_perm, 0b11111100 ) );
	r_v4 = _mm256_add_epi64( offset0_v4, _mm256_blend_epi32( offset1_v3_perm, offset1_v4_perm, 0b11111100 ) );
	r_v5 = _mm256_add_epi64( offset0_v5, _mm256_blend_epi32( offset1_v4_perm, offset1_v5_perm, 0b11111100 ) );
	r_v6 = _mm256_add_epi64( offset0_v6, _mm256_blend_epi32( offset1_v5_perm, offset1_v6_perm, 0b11111100 ) );
	
	// add offset 2 to result
	r_v0 = _mm256_add_epi64( r_v0, _mm256_blend_epi32( const_0_vec    , offset2_v0_perm, 0b11110000 ) );
	r_v1 = _mm256_add_epi64( r_v1, _mm256_blend_epi32( offset2_v0_perm, offset2_v1_perm, 0b11110000 ) );
	r_v2 = _mm256_add_epi64( r_v2, _mm256_blend_epi32( offset2_v1_perm, offset2_v2_perm, 0b11110000 ) );
	r_v3 = _mm256_add_epi64( r_v3, _mm256_blend_epi32( offset2_v2_perm, offset2_v3_perm, 0b11110000 ) );
	r_v4 = _mm256_add_epi64( r_v4, _mm256_blend_epi32( offset2_v3_perm, offset2_v4_perm, 0b11110000 ) );
	r_v5 = _mm256_add_epi64( r_v5, _mm256_blend_epi32( offset2_v4_perm, offset2_v5_perm, 0b11110000 ) );
	r_v6 = _mm256_add_epi64( r_v6, _mm256_blend_epi32( offset2_v5_perm, offset2_v6_perm, 0b11110000 ) );
	r_v7 = _mm256_blend_epi32( offset2_v6_perm, const_0_vec, 0b11111100 );

	// add offset 3 to result
	r_v0 = _mm256_add_epi64( r_v0, _mm256_blend_epi32( const_0_vec    , offset3_v0_perm, 0b11000000 ) );
	r_v1 = _mm256_add_epi64( r_v1, _mm256_blend_epi32( offset3_v0_perm, offset3_v1_perm, 0b11000000 ) );
	r_v2 = _mm256_add_epi64( r_v2, _mm256_blend_epi32( offset3_v1_perm, offset3_v2_perm, 0b11000000 ) );
	r_v3 = _mm256_add_epi64( r_v3, _mm256_blend_epi32( offset3_v2_perm, offset3_v3_perm, 0b11000000 ) );
	r_v4 = _mm256_add_epi64( r_v4, _mm256_blend_epi32( offset3_v3_perm, offset3_v4_perm, 0b11000000 ) );
	r_v5 = _mm256_add_epi64( r_v5, _mm256_blend_epi32( offset3_v4_perm, offset3_v5_perm, 0b11000000 ) );
	r_v6 = _mm256_add_epi64( r_v6, _mm256_blend_epi32( offset3_v5_perm, const_0_vec    , 0b11110000 ) );
	    
    // reduce
	reduce_8_to_4(&r_v0, &r_v1, &r_v2, &r_v3, &r_v4, &r_v5, &r_v6, &r_v7);
    reduce_4_to_2(&r_v0, &r_v1, &r_v2, &r_v3);
    reduce_2_to_2(&r_v0, &r_v1);
    // copy back the result vector
    copy_256(reslo, reshi, r_v0, r_v1);
}

void square_256_old(__m256i *reslo, __m256i *reshi, __m256i xlo,__m256i xhi)
{
    /*
    res = x^2
    */

	/*        cl1       2             3            4            5           6           7           8          9           10           11        12         13          14                  15 
    r[0] =   a0*a0;
    r[1] = _2a0*a1;
    r[2] = _2a0*a2    + a1*a1;
    r[3] = _2a0*a3    + _2a1*a2;
    r[4] = _2a0*a4    + _2a1*a3     + a2*a2;
    r[5] = _2a0*a5    + _2a1*a4   + _2a2*a3;
    r[6] = _2a0*a6    + _2a1*a5   + _2a2*a4      + a3*a3;
    r[7] = _2a0*a7    + _2a1*a6   + _2a2*a5     + _2a3*a4;
    r[8] = _2a0*a8    + _2a1*a7   + _2a2*a6     + _2a3*a5                                                                                               + a4*a4;
    r[9] = _2a0*a9    + _2a1*a8   + _2a2*a7     + _2a3*a6   + _2a4*a5;
    r[10] = _2a0*a10  + _2a1*a9   + _2a2*a8     + _2a3*a7   + _2a4*a6                                                                                                   + a5*a5;
    r[11] = _2a0*a11  + _2a1*a10  + _2a2*a9     + _2a3*a8   + _2a4*a7   + _2a5*a6;
    r[12] = _2a0*a12  + _2a1*a11  + _2a2*a10    + _2a3*a9   + _2a4*a8   + _2a5*a7                                                                                           + a6*a6;
    r[13] = _2a0*a13  + _2a1*a12  + _2a2*a11    + _2a3*a10  + _2a4*a9   + _2a5*a8  + _2a6*a7;
    r[14] = _2a0*a14  + _2a1*a13  + _2a2*a12    + _2a3*a11  + _2a4*a10  + _2a5*a9  + _2a6*a8                                                                                                +  a7*a7;
    r[15] =           + _2a1*a14  + _2a2*a13  + _2a3*a12    + _2a4*a11  + _2a5*a10  + _2a6*a9  + _2a7*a8;
    r[16] =                       + _2a2*a14  + _2a3*a13  + _2a4*a12    + _2a5*a11  + _2a6*a10  + _2a7*a9                                                                                                           + a8*a8;
    r[17] =                                   + _2a3*a14  + _2a4*a13  + _2a5*a12    + _2a6*a11  + _2a7*a10  + _2a8*a9;
    r[18] =                                               + _2a4*a14  + _2a5*a13  + _2a6*a12    + _2a7*a11  + _2a8*a10                                                                                                      + a9*a9;
    r[19] =                                                           +  _2a5*a14  + _2a6*a13  + _2a7*a12    + _2a8*a11  + _2a9*a10;
    r[20] =                                                                        + _2a6*a14  + _2a7*a13  + _2a8*a12    + _2a9*a11                                                 + a10*a10;
    r[21] =                                                                                    + _2a7*a14  + _2a8*a13  + _2a9*a12    + _2a10*a11;
    r[22] =                                                                                                + _2a8*a14  + _2a9*a13  + _2a10*a12  + a11*a11;
    r[23] =                                                                                                            + _2a9*a14  + _2a10*a13 + _2a11*a12;
    r[24] =                                                                                                                         +_2a10*a14 + _2a11*a13    + a12*a12;
    r[25] =                                                                                                                                     + _2a11*a14   +  _2a12*a13
    r[26] =                                                                                                                                                   + _2a12*a14   +  a13*a13;
    r[27] =                                                                                                                                                                _2a13*a14;
    
    r[28] =                                                                                                                                                                                   a14*a14;
*/  
    
    // xlo = [0,1,2,3,4,5,6,7] xhi = [8,9,10,11,12,13,14,15]
	__m256i a_11_14   = _mm256_permutevar8x32_epi32( xhi, index3_vec );
    a_11_14 = _mm256_blend_epi32(const_0_vec, a_11_14, 0b01010101);
	__m256i a_7_10   = _mm256_permutevar8x32_epi32( xhi, index4_vec );
    a_7_10 = _mm256_blend_epi32(const_0_vec, a_7_10, 0b01010101);
    __m256i tmp = _mm256_permutevar8x32_epi32(xlo, index5_vec);
    a_7_10 = _mm256_blend_epi32(tmp, a_7_10, 0b11111110);
	__m256i a_3_6  = _mm256_permutevar8x32_epi32( xlo, index3_vec );
    a_3_6 = _mm256_blend_epi32(const_0_vec, a_3_6, 0b01010101);
	__m256i a_m1_2 = _mm256_permutevar8x32_epi32( xlo, index4_vec );
    a_m1_2 = _mm256_blend_epi32(const_0_vec, a_m1_2, 0b01010101);
    tmp = _mm256_permutevar8x32_epi32(xlo, index5_vec);
    a_m1_2 = _mm256_blend_epi32(tmp, a_m1_2, 0b11111110);
    a_m1_2 = _mm256_blend_epi32(const_0_vec, a_m1_2, 0b11111100);
    
    __m256i a2_0_7   = _mm256_slli_epi64( xlo, 1 );
	__m256i a2_8_14   = _mm256_slli_epi64( xhi, 1 );
	
	__m256i a2_0     = _mm256_broadcastq_epi64( _mm256_castsi256_si128( a2_0_7 ) );
	__m256i a2_1     = _mm256_permutevar8x32_epi32( a2_0, index_alt_vec );
	__m256i a2_2     = _mm256_broadcastq_epi64( _mm256_castsi256_si128( _mm256_unpackhi_epi64( a2_0_7, index1_vec ) ) );
	__m256i a2_3     = _mm256_permutevar8x32_epi32( a2_2, index_alt_vec );
	
	__m128i a_lo_hi = _mm256_extracti128_si256( a2_0_7, 1 );
	__m256i a2_4     = _mm256_broadcastq_epi64( a_lo_hi );
	__m256i a2_5     = _mm256_permutevar8x32_epi32( a2_4, index_alt_vec );
	__m256i a2_6     = _mm256_broadcastq_epi64( _mm256_castsi256_si128( _mm256_unpackhi_epi64( _mm256_castsi128_si256( a_lo_hi ), index1_vec ) ) );
	__m256i a2_7     = _mm256_permutevar8x32_epi32( a2_6 , index_alt_vec );

	__m256i a2_8     = _mm256_broadcastq_epi64( _mm256_castsi256_si128( a2_8_14 ) );
	__m256i a2_9     = _mm256_permutevar8x32_epi32( a2_8, index_alt_vec );
	__m256i a2_10    = _mm256_broadcastq_epi64( _mm256_castsi256_si128( _mm256_unpackhi_epi64( a2_8_14, index1_vec ) ) );
	__m256i a2_11    = _mm256_permutevar8x32_epi32( a2_10, index_alt_vec );
	
	__m128i a_hi_hi = _mm256_extracti128_si256( a2_8_14, 1 );
	__m256i a2_12    = _mm256_broadcastq_epi64( a_hi_hi );
	__m256i a2_13    = _mm256_permutevar8x32_epi32( a2_12, index_alt_vec );

    // column vectors
     __m256i c1_v0, c1_v1, c1_v2, c1_v3;  //column 1 results
     __m256i c2_v0, c2_v1, c2_v2, c2_v3;     
     __m256i c3_v0, c3_v1, c3_v2;
     __m256i c4_v0, c4_v1, c4_v2;
     __m256i c5_v0, c5_v1, c5_v2;
     __m256i c6_v0, c6_v1, c6_v2;
     __m256i c7_v0, c7_v1;
     __m256i c8_v0, c8_v1;
     __m256i c9_v0, c9_v1;
     __m256i c10_v0, c10_v1;
     __m256i c11_v0;
     __m256i c12_v0;
     __m256i c13_v0;
     __m256i c14_v0;
     __m256i c15_v0;
     // col 1
     c1_v0 = _mm256_mul_epi32(a2_0, a_11_14);
     c1_v1 = _mm256_mul_epi32(a2_0, a_7_10);
     c1_v2 = _mm256_mul_epi32(a2_0, a_3_6);
     c1_v3 = _mm256_mul_epi32(a2_0, a_m1_2);
     c1_v3 = _mm256_blend_epi32(const_0_vec, c1_v3, 0b11110000);
     // col 2
     c2_v0 = _mm256_mul_epi32(a2_1, a_11_14);
     c2_v1 = _mm256_mul_epi32(a2_1, a_7_10);
     c2_v2 = _mm256_mul_epi32(a2_1, a_3_6);
     c2_v3 = _mm256_mul_epi32(a2_1, a_m1_2);
     c2_v3 = _mm256_blend_epi32(const_0_vec, c2_v3, 0b11000000);
     // col 3
     c3_v0 = _mm256_mul_epi32(a2_2, a_11_14);
     c3_v1 = _mm256_mul_epi32(a2_2, a_7_10);
     c3_v2 = _mm256_mul_epi32(a2_2, a_3_6);
     // col 4
     c4_v0 = _mm256_mul_epi32(a2_3, a_11_14);
     c4_v1 = _mm256_mul_epi32(a2_3, a_7_10);
     c4_v2 = _mm256_mul_epi32(a2_3, a_3_6);
     c4_v2 = _mm256_blend_epi32(const_0_vec, c4_v2, 0b11111100);
     // col 5
     c5_v0 = _mm256_mul_epi32(a2_4, a_11_14);
     c5_v1 = _mm256_mul_epi32(a2_4, a_7_10);
     c5_v2 = _mm256_mul_epi32(a2_4, a_3_6);
     c5_v2 = _mm256_blend_epi32(const_0_vec, c5_v2, 0b11110000);
     // col 6
     c6_v0 = _mm256_mul_epi32(a2_5, a_11_14);
     c6_v1 = _mm256_mul_epi32(a2_5, a_7_10);
     c6_v2 = _mm256_mul_epi32(a2_5, a_3_6);
     c6_v2 = _mm256_blend_epi32(const_0_vec, c6_v2, 0b11000000);
     // col 7
     c7_v0 = _mm256_mul_epi32(a2_6, a_11_14);
     c7_v1 = _mm256_mul_epi32(a2_6, a_7_10);
     // col 8
     c8_v0 = _mm256_mul_epi32(a2_7, a_11_14);
     c8_v1 = _mm256_mul_epi32(a2_7, a_7_10);
     c8_v1 = _mm256_blend_epi32(const_0_vec, c8_v1, 0b11111100);
     // col 9
     c9_v0 = _mm256_mul_epi32(a2_8, a_11_14);
     c9_v1 = _mm256_mul_epi32(a2_8, a_7_10);
     c9_v1 = _mm256_blend_epi32(const_0_vec, c9_v1, 0b11110000);
     // col 10
     c10_v0 = _mm256_mul_epi32(a2_9, a_11_14);
     c10_v1 = _mm256_mul_epi32(a2_9, a_7_10);
     c10_v1 = _mm256_blend_epi32(const_0_vec, c10_v1, 0b11000000);
     // col 11
     c11_v0 = _mm256_mul_epi32(a2_10, a_11_14);
     // col 12
     c12_v0 = _mm256_mul_epi32(a2_11, a_11_14);
     c12_v0 = _mm256_blend_epi32(const_0_vec, c12_v0, 0b11111100);
     // col 13
     c13_v0 = _mm256_mul_epi32(a2_12, a_11_14);
     c13_v0 = _mm256_blend_epi32(const_0_vec, c13_v0, 0b11110000);
     // col 14
     c14_v0 = _mm256_mul_epi32(a2_13, a_11_14);
     c14_v0 = _mm256_blend_epi32(const_0_vec, c14_v0, 0b11000000);

    // build the shifted result vectors
     __m256i    r0_3,   r4_7,   r8_11,  r12_15, r16_19, r20_23, r24_27, r28,
                rm1_2,  r3_6,   r7_10,  r11_14, r15_18, r19_22, r23_26,
                                r6_9,   r10_13, r14_17, r18_21, r22_25,
                                r5_8,   r9_12,  r13_16, r17_20, r21_24;
    //  r0_3+=c2v3      r4_7+=c2v2      r8_11+=c2v1+c6v2        r12_15+=c2v0+c6v1       r16_19+=c6v0+c10v1  r20_23+=c10v0       r24_27+=c14v0       
    //  rm1_2+=c1v3     r3_6+=c1v2      r7_10+=c1v1+c5v2        r11_14+=c1v0+c5v1       r15_18+=c5v0+c9v1   r19_22+=c9v0        r23_26+=c13v0       
    //                                  r6_9+=c4v2              r10_13+=c4v1            r14_17+=c4v0+c8v1   r18_21+=c8v0        r22_25+=c12v0       
    //                                  r5_8+=c3v2              r9_12+=c3v1             r13_16+=c3v0+c7v1   r17_20+=c7v0        r21_24+=c11v0       
    r28 = const_0_vec;
    // Offset 0
    r0_3 = c2_v3;
    r4_7 = c2_v2;
    r8_11 = _mm256_add_epi64(c2_v1,c6_v2);
    r12_15 = _mm256_add_epi64(c2_v0,c6_v1);
    r16_19 = _mm256_add_epi64(c6_v0, c10_v1);
    r20_23 = c10_v0;
    r24_27 = c14_v0;
    // Offset -1
    rm1_2 = c1_v3;
    r3_6 = c1_v2;
    r7_10 = _mm256_add_epi64(c1_v1, c5_v2);
    r11_14 = _mm256_add_epi64(c1_v0, c5_v1);
    r15_18 = _mm256_add_epi64(c5_v0, c9_v1);
    r19_22 = c9_v0;
    r23_26 = c13_v0;
    // Offset -2
    r6_9 = c4_v2;
    r10_13 = c4_v1;
    r14_17 = _mm256_add_epi64(c4_v0, c8_v1);
    r18_21 = c8_v0;
    r22_25 = c12_v0;
    // Offset -3
    r5_8 = c3_v2;
    r9_12 = c3_v1;
    r13_16 = _mm256_add_epi64(c3_v0, c7_v1);
    r17_20 = c7_v0;
    r21_24 = c11_v0;

    // add up the different shifted result vectors

    // add offset -1 to offset 0
    // 0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 19, 20 21 22 23, 24 25 26 27
    //-1 0 1 2, 3 4 5 6, 7 8 9 10, 11 12 13 14, 15 16 17 18, 19 20 21 22, 23 24 25 26
    // permute offset -1
    // 0 1 2 -1,4 5 6 3, 8 9 10 7, 12 13 14 11, 16 17 18 15, 20 21 22 19, 24 25 26 23
    rm1_2 = _mm256_permute4x64_epi64(rm1_2, 0b00111001);
    r3_6 = _mm256_permute4x64_epi64(r3_6, 0b00111001);
    r7_10 = _mm256_permute4x64_epi64(r7_10, 0b00111001);
    r11_14 = _mm256_permute4x64_epi64(r11_14, 0b00111001);
    r15_18 = _mm256_permute4x64_epi64(r15_18, 0b00111001);
    r19_22 = _mm256_permute4x64_epi64(r19_22, 0b00111001);
    r23_26 = _mm256_permute4x64_epi64(r23_26, 0b00111001);
    // blend and add
    r0_3 = _mm256_add_epi64(r0_3, _mm256_blend_epi32(r3_6, rm1_2, 0b00111111));
    r4_7 = _mm256_add_epi64(r4_7, _mm256_blend_epi32(r7_10, r3_6, 0b00111111));
    r8_11 = _mm256_add_epi64(r8_11, _mm256_blend_epi32(r11_14, r7_10, 0b00111111));
    r12_15 = _mm256_add_epi64(r12_15, _mm256_blend_epi32(r15_18, r11_14, 0b00111111));
    r16_19 = _mm256_add_epi64(r16_19, _mm256_blend_epi32(r19_22, r15_18, 0b00111111));
    r20_23 = _mm256_add_epi64(r20_23, _mm256_blend_epi32(r23_26, r19_22, 0b00111111));
    r24_27 = _mm256_add_epi64(r24_27, _mm256_blend_epi32(const_0_vec, r23_26, 0b00111111));

    // add offset -2 to offset 0
    // 0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 19, 20 21 22 23, 24 25 26 27
    //          6 7 8 9, 10 11 12 13, 14 15 16 17, 18 19 20 21, 22 23 24 25
    // permute offset -2
    //          8 9 6 7, 12 13 10 11, 16 17 14 15, 20 21 18 19, 24 25 22 23
    r6_9 = _mm256_permute4x64_epi64(r6_9, 0b01001110);
    r10_13 = _mm256_permute4x64_epi64(r10_13, 0b01001110);
    r14_17 = _mm256_permute4x64_epi64(r14_17, 0b01001110);
    r18_21 = _mm256_permute4x64_epi64(r18_21, 0b01001110);
    r22_25 = _mm256_permute4x64_epi64(r22_25, 0b01001110);
    // blend and add
    __m256i test = _mm256_blend_epi32(const_0_vec, r4_7, 0b11110000);
    r4_7 = _mm256_add_epi64(r4_7, _mm256_blend_epi32(const_0_vec, r6_9, 0b11110000));
    r8_11 = _mm256_add_epi64(r8_11, _mm256_blend_epi32(r6_9, r10_13, 0b11110000));
    r12_15 = _mm256_add_epi64(r12_15, _mm256_blend_epi32(r10_13, r14_17, 0b11110000));
    r16_19 = _mm256_add_epi64(r16_19, _mm256_blend_epi32(r14_17, r18_21, 0b11110000));
    r20_23 = _mm256_add_epi64(r20_23, _mm256_blend_epi32(r18_21, r22_25, 0b11110000));
    r24_27 = _mm256_add_epi64(r24_27,_mm256_blend_epi32(r22_25, const_0_vec, 0b11110000));

    // add offset -3 to offset 0
    // 0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 19, 20 21 22 23, 24 25 26 27
    //          5 6 7 8, 9 10 11 12, 13 14 15 16, 17 18 19 20, 21 22 23 24
    // permute offset -3
    //          8 5 6 7, 12 9 10 11, 16 13 14 15, 20 17 18 19, 24 21 22 23 
    r5_8 = _mm256_permute4x64_epi64(r5_8, 0b10010011);
    r9_12 = _mm256_permute4x64_epi64(r9_12, 0b10010011);
    r13_16 = _mm256_permute4x64_epi64(r13_16, 0b10010011);
    r17_20 = _mm256_permute4x64_epi64(r17_20, 0b10010011);
    r21_24 = _mm256_permute4x64_epi64(r21_24, 0b10010011); 
    // blend and add
    r4_7 = _mm256_add_epi64(r4_7, _mm256_blend_epi32(const_0_vec, r5_8, 0b11111100));
    r8_11 = _mm256_add_epi64(r8_11, _mm256_blend_epi32(r5_8, r9_12, 0b11111100));
    r12_15 = _mm256_add_epi64(r12_15, _mm256_blend_epi32(r9_12, r13_16, 0b11111100));
    r16_19 = _mm256_add_epi64(r16_19, _mm256_blend_epi32(r13_16, r17_20, 0b11111100));
    r20_23 = _mm256_add_epi64(r20_23, _mm256_blend_epi32(r17_20, r21_24, 0b11111100));
    r24_27 = _mm256_add_epi64(r24_27,_mm256_blend_epi32(r21_24, const_0_vec, 0b11111100));
    
    // create the square vectors
    // 0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 19, 20 21 22 23, 24 25 26 27
   
    __m256i sqr_a_11_14, sqr_a_7_10, sqr_a_3_6, sqr_a_m1_2; 
    sqr_a_m1_2 = _mm256_mul_epi32(a_m1_2,a_m1_2); // [a-1^2, a0^2, a1^2, a2^2]
    sqr_a_3_6 = _mm256_mul_epi32(a_3_6,a_3_6);
    sqr_a_7_10 = _mm256_mul_epi32(a_7_10,a_7_10);
    sqr_a_11_14 = _mm256_mul_epi32(a_11_14,a_11_14);
    
    //[a0^2, 0, a1^2, 0]  
    __m256i sqrv0, sqrv1, sqrv2, sqrv3, sqrv4, sqrv5, sqrv6, sqrv7;
    sqrv0 = _mm256_permute4x64_epi64(sqr_a_m1_2, 0b11100001);       //a0, a1
    sqrv0 = _mm256_blend_epi32(const_0_vec, sqrv0, 0b00110011);

    sqrv1 = _mm256_blend_epi32(sqr_a_m1_2,sqr_a_3_6,0b00000011);    //a2, a3
    sqrv1 = _mm256_permute4x64_epi64(sqrv1, 0b01001011);
    sqrv1 = _mm256_blend_epi32(const_0_vec, sqrv1, 0b00110011);

    sqrv2 = _mm256_permute4x64_epi64(sqr_a_3_6, 0b11100001);        //a4, a5
    sqrv2 = _mm256_blend_epi32(const_0_vec, sqrv2, 0b00110011);

    sqrv3 = _mm256_blend_epi32(sqr_a_3_6,sqr_a_7_10,0b00000011);    //a6, a7
    sqrv3 = _mm256_permute4x64_epi64(sqrv3, 0b01001011);
    sqrv3 = _mm256_blend_epi32(const_0_vec, sqrv3, 0b00110011);

    sqrv4 = _mm256_permute4x64_epi64(sqr_a_7_10, 0b11100001);        //a8, a9
    sqrv4 = _mm256_blend_epi32(const_0_vec, sqrv4, 0b00110011);

    sqrv5 = _mm256_blend_epi32(sqr_a_7_10,sqr_a_11_14,0b00000011);    //a10, a11
    sqrv5 = _mm256_permute4x64_epi64(sqrv5, 0b01001011);
    sqrv5 = _mm256_blend_epi32(const_0_vec, sqrv5, 0b00110011);

    sqrv6 = _mm256_permute4x64_epi64(sqr_a_11_14, 0b11100001);       //a12, a13
    sqrv6 = _mm256_blend_epi32(const_0_vec, sqrv6, 0b00110011);

    sqrv7 = _mm256_permute4x64_epi64(sqr_a_11_14, 0b11111111);       //a14
    sqrv7 = _mm256_blend_epi32(const_0_vec, sqrv7, 0b00000011);
    
    // add squares to result
    r0_3 = _mm256_add_epi64(r0_3, sqrv0);
    r4_7 = _mm256_add_epi64(r4_7, sqrv1);
    r8_11 = _mm256_add_epi64(r8_11, sqrv2);
    r12_15 = _mm256_add_epi64(r12_15, sqrv3);
    r16_19 = _mm256_add_epi64(r16_19, sqrv4);
    r20_23 = _mm256_add_epi64(r20_23, sqrv5);
    r24_27 = _mm256_add_epi64(r24_27, sqrv6);
    r28 = _mm256_add_epi64(r28, sqrv7);

    // reduce the result
    reduce_8_to_4(&r0_3, &r4_7, &r8_11, &r12_15, &r16_19, &r20_23, &r24_27, &r28);
    reduce_4_to_2(&r0_3, &r4_7, &r8_11, &r12_15);
    reduce_2_to_2(&r0_3, &r4_7);
    copy_256(reslo, reshi, r0_3,r4_7);
}

void square_256(__m256i *reslo, __m256i *reshi, __m256i xlo,__m256i xhi)
{
    /*
    res = x^2
    */
    
    // xlo = [0,1,2,3,4,5,6,7] xhi = [8,9,10,11,12,13,14,15]            //new
	__m256i a_0_3   = _mm256_permutevar8x32_epi32( xlo, index1_vec );
	__m256i a_4_7   = _mm256_permutevar8x32_epi32( xlo, index2_vec );
	__m256i a_8_11  = _mm256_permutevar8x32_epi32( xhi, index1_vec );
	__m256i a_12_15 = _mm256_permutevar8x32_epi32( xhi, index2_vec );
    
    __m256i a2_0_7  = _mm256_slli_epi32( xlo, 1 );      //changed this from 64
	__m256i a2_8_14 = _mm256_slli_epi32( xhi, 1 );
	
	__m256i a2_0    = _mm256_broadcastq_epi64( _mm256_castsi256_si128( a2_0_7 ) );
	__m256i a2_1    = _mm256_permutevar8x32_epi32( a2_0, index_alt_vec );
	__m256i a2_2    = _mm256_broadcastq_epi64( _mm256_castsi256_si128( _mm256_unpackhi_epi64( a2_0_7, index1_vec ) ) );
	__m256i a2_3    = _mm256_permutevar8x32_epi32( a2_2, index_alt_vec );
	
	__m128i a_lo_hi = _mm256_extracti128_si256( a2_0_7, 1 );
	__m256i a2_4    = _mm256_broadcastq_epi64( a_lo_hi );
	__m256i a2_5    = _mm256_permutevar8x32_epi32( a2_4, index_alt_vec );
	__m256i a2_6    = _mm256_broadcastq_epi64( _mm256_castsi256_si128( _mm256_unpackhi_epi64( _mm256_castsi128_si256( a_lo_hi ), index1_vec ) ) );
	__m256i a2_7    = _mm256_permutevar8x32_epi32( a2_6 , index_alt_vec );

	__m256i a2_8    = _mm256_broadcastq_epi64( _mm256_castsi256_si128( a2_8_14 ) );
	__m256i a2_9    = _mm256_permutevar8x32_epi32( a2_8, index_alt_vec );
	__m256i a2_10   = _mm256_broadcastq_epi64( _mm256_castsi256_si128( _mm256_unpackhi_epi64( a2_8_14, index1_vec ) ) );
	__m256i a2_11   = _mm256_permutevar8x32_epi32( a2_10, index_alt_vec );
	
	__m128i a_hi_hi = _mm256_extracti128_si256( a2_8_14, 1 );
	__m256i a2_12   = _mm256_broadcastq_epi64( a_hi_hi );
	__m256i a2_13   = _mm256_permutevar8x32_epi32( a2_12, index_alt_vec );

    // column vectors
     __m256i c1_v0 , c1_v1 , c1_v2 , c1_v3;  //column 1 results
     __m256i c2_v0 , c2_v1 , c2_v2 , c2_v3;     
     __m256i c3_v0 , c3_v1 , c3_v2 , c3_v3;
     __m256i c4_v0 , c4_v1 , c4_v2 , c4_v3;
     __m256i /*c5_v0 ,*/ c5_v1 , c5_v2 , c5_v3;
     __m256i /*c6_v0 ,*/ c6_v1 , c6_v2 , c6_v3;
     __m256i /*c7_v0 ,*/ c7_v1 , c7_v2 , c7_v3;
     __m256i /*c8_v0 ,*/ c8_v1 , c8_v2 , c8_v3;
     __m256i /*c9_v0 , c9_v1 ,*/ c9_v2 , c9_v3;
     __m256i /*c10_v0, c10_v1,*/ c10_v2, c10_v3;
     __m256i /*c11_v0, c11_v1,*/ c11_v2, c11_v3;
     __m256i /*c12_v0, c12_v1,*/ c12_v2, c12_v3;
     __m256i /*c13_v0, c13_v1, c13_v2,*/ c13_v3;
     __m256i /*c14_v0, c14_v1, c14_v2,*/ c14_v3;
     __m256i /*c15_v0, c15_v1, c15_v2,*/ c15_v3;

	/*        cl1       2             3            4            5           6           7           8          9           10           11        12         13          14                  15 
    r[0]  =   a0*a0;     
    r[1]  = _2a0*a1;     0             
    r[2]  = _2a0*a2    + a1*a1;         0           
    r[3]  = _2a0*a3    + _2a1*a2;       1            0        
    r[4]  = _2a0*a4    + _2a1*a3     + a2*a2;        1           0
    r[5]  = _2a0*a5    + _2a1*a4   + _2a2*a3;        2           1           0
    r[6]  = _2a0*a6    + _2a1*a5   + _2a2*a4      + a3*a3;       2           1           0
    r[7]  = _2a0*a7    + _2a1*a6   + _2a2*a5     + _2a3*a4;      3           2           1           0
    r[8]  = _2a0*a8    + _2a1*a7   + _2a2*a6     + _2a3*a5   + a4*a4;        3           2           1          0
    r[9]  = _2a0*a9    + _2a1*a8   + _2a2*a7     + _2a3*a6   + _2a4*a5;      4           3           2          1           0
    r[10] = _2a0*a10  + _2a1*a9   + _2a2*a8     + _2a3*a7   + _2a4*a6     + a5*a5;      4           3           2           1           0
    r[11] = _2a0*a11  + _2a1*a10  + _2a2*a9     + _2a3*a8   + _2a4*a7   + _2a5*a6;      5           4           3           2           1           0
    r[12] = _2a0*a12  + _2a1*a11  + _2a2*a10    + _2a3*a9   + _2a4*a8   + _2a5*a7    + a6*a6;       5           4           3           2           1           0
    r[13] = _2a0*a13  + _2a1*a12  + _2a2*a11    + _2a3*a10  + _2a4*a9   + _2a5*a8  + _2a6*a7;       6           5           4           3           2           1           0          
    r[14] = _2a0*a14  + _2a1*a13  + _2a2*a12    + _2a3*a11  + _2a4*a10  + _2a5*a9  + _2a6*a8       a7*a7;       6           5           4           3           2           1           0
    r[15] =           + _2a1*a14  + _2a2*a13  + _2a3*a12    + _2a4*a11  + _2a5*a10  + _2a6*a9  + _2a7*a8;       7           6           5           4           3           2           1           0
    r[16] =                       + _2a2*a14  + _2a3*a13  + _2a4*a12    + _2a5*a11  + _2a6*a10  + _2a7*a9     + a8*a8;      7           6           5           4           3           2           1
    r[17] =                                   + _2a3*a14  + _2a4*a13  + _2a5*a12    + _2a6*a11  + _2a7*a10  + _2a8*a9;      8           7           6           5           4           3           2           
    r[18] =                                               + _2a4*a14  + _2a5*a13  + _2a6*a12    + _2a7*a11  + _2a8*a10       + a9*a9;   8           7           6           5           4           3
    r[19] =                                                           +  _2a5*a14  + _2a6*a13  + _2a7*a12    + _2a8*a11  + _2a9*a10;    9           8           7           6           5           4   
    r[20] =                                                                        + _2a6*a14  + _2a7*a13  + _2a8*a12    + _2a9*a11    + a10*a10;   9           8           7           6           5
    r[21] =                                                                                    + _2a7*a14  + _2a8*a13  + _2a9*a12    + _2a10*a11;   10          9           8           7           6
    r[22] =                                                                                                + _2a8*a14  + _2a9*a13  + _2a10*a12  + a11*a11;      10          9           8           7
    r[23] =                                                                                                            + _2a9*a14  + _2a10*a13 + _2a11*a12;     11          10          9           8
    r[24] =                                                                                                                         +_2a10*a14 + _2a11*a13    + a12*a12;     11         10              
    r[25] =                                                                                                                                     + _2a11*a14   +  _2a12*a13    12           11
    r[26] =                                                                                                                                                   + _2a12*a14   +  a13*a13;     12
    r[27] =                                                                                                                                                                _2a13*a14;       13
    r[28] =                                                                                                                                                                            a14*a14;
*/  

    c1_v0 = _mm256_mul_epi32(_mm256_blend_epi32(a2_0, a_0_3, 0b00000011), a_0_3);
    c1_v1 = _mm256_mul_epi32(a2_0, a_4_7);
    c1_v2 = _mm256_mul_epi32(a2_0, a_8_11);
    c1_v3 = _mm256_mul_epi32(a2_0, a_12_15);
    c1_v3 = _mm256_blend_epi32(c1_v3, const_0_vec, 0b11000000);
 
    c2_v0 = _mm256_mul_epi32(_mm256_blend_epi32(a2_1, a_0_3, 0b00001100), _mm256_blend_epi32(a_0_3, const_0_vec, 0b00000011));
    c2_v1 = _mm256_mul_epi32(a2_1, a_4_7);
    c2_v2 = _mm256_mul_epi32(a2_1, a_8_11);
    c2_v3 = _mm256_mul_epi32(a2_1, a_12_15);
    c2_v3 = _mm256_blend_epi32(c2_v3, const_0_vec, 0b11000000);

    c3_v0 = _mm256_mul_epi32(_mm256_blend_epi32(a2_2, a_0_3, 0b00110000), _mm256_blend_epi32(a_0_3, const_0_vec, 0b00001111));
    c3_v1 = _mm256_mul_epi32(a2_2, a_4_7);
    c3_v2 = _mm256_mul_epi32(a2_2, a_8_11);
    c3_v3 = _mm256_mul_epi32(a2_2, a_12_15);
    c3_v3 = _mm256_blend_epi32(c3_v3, const_0_vec, 0b11000000); 

    c4_v0 = _mm256_mul_epi32(_mm256_blend_epi32(a2_3, a_0_3, 0b11000000), _mm256_blend_epi32(a_0_3, const_0_vec, 0b00111111));
    c4_v1 = _mm256_mul_epi32(a2_3, a_4_7);
    c4_v2 = _mm256_mul_epi32(a2_3, a_8_11);
    c4_v3 = _mm256_mul_epi32(a2_3, a_12_15);
    c4_v3 = _mm256_blend_epi32(c4_v3, const_0_vec, 0b11000000); 

    c5_v1 = _mm256_mul_epi32(_mm256_blend_epi32(a2_4, a_4_7, 0b00000011), a_4_7);
    c5_v2 = _mm256_mul_epi32(a2_4, a_8_11);
    c5_v3 = _mm256_mul_epi32(a2_4, a_12_15);
    c5_v3 = _mm256_blend_epi32(c5_v3, const_0_vec, 0b11000000); 

    c6_v1 = _mm256_mul_epi32(_mm256_blend_epi32(a2_5, a_4_7, 0b00001100), _mm256_blend_epi32(a_4_7, const_0_vec, 0b00000011));
    c6_v2 = _mm256_mul_epi32(a2_5, a_8_11);
    c6_v3 = _mm256_mul_epi32(a2_5, a_12_15);
    c6_v3 = _mm256_blend_epi32(c6_v3, const_0_vec, 0b11000000); 

    c7_v1 = _mm256_mul_epi32(_mm256_blend_epi32(a2_6, a_4_7, 0b00110000), _mm256_blend_epi32(a_4_7, const_0_vec, 0b00001111));
    c7_v2 = _mm256_mul_epi32(a2_6, a_8_11);
    c7_v3 = _mm256_mul_epi32(a2_6, a_12_15);
    c7_v3 = _mm256_blend_epi32(c7_v3, const_0_vec, 0b11000000); 

    c8_v1 = _mm256_mul_epi32(_mm256_blend_epi32(a2_7, a_4_7, 0b11000000), _mm256_blend_epi32(a_4_7, const_0_vec, 0b00111111));
    c8_v2 = _mm256_mul_epi32(a2_7, a_8_11);
    c8_v3 = _mm256_mul_epi32(a2_7, a_12_15);
    c8_v3 = _mm256_blend_epi32(c8_v3, const_0_vec, 0b11000000); 

    c9_v2 = _mm256_mul_epi32(_mm256_blend_epi32(a2_8, a_8_11, 0b00000011), a_8_11);
    c9_v3 = _mm256_mul_epi32(a2_8, a_12_15);
    c9_v3 = _mm256_blend_epi32(c9_v3, const_0_vec, 0b11000000); 

    c10_v2 = _mm256_mul_epi32(_mm256_blend_epi32(a2_9, a_8_11, 0b00001100), _mm256_blend_epi32(a_8_11, const_0_vec, 0b00000011));
    c10_v3 = _mm256_mul_epi32(a2_9, a_12_15);
    c10_v3 = _mm256_blend_epi32(c10_v3, const_0_vec, 0b11000000); 

    c11_v2 = _mm256_mul_epi32(_mm256_blend_epi32(a2_10, a_8_11, 0b00110000), _mm256_blend_epi32(a_8_11, const_0_vec, 0b00001111));
    c11_v3 = _mm256_mul_epi32(a2_10, a_12_15);
    c11_v3 = _mm256_blend_epi32(c11_v3, const_0_vec, 0b11000000);
     
    c12_v2 = _mm256_mul_epi32(_mm256_blend_epi32(a2_11, a_8_11, 0b11000000), _mm256_blend_epi32(a_8_11, const_0_vec, 0b00111111));
    c12_v3 = _mm256_mul_epi32(a2_11, a_12_15);
    c12_v3 = _mm256_blend_epi32(c12_v3, const_0_vec, 0b11000000);
    
    c13_v3 = _mm256_mul_epi32(_mm256_blend_epi32(a2_12, a_12_15, 0b00000011), a_12_15);
    c13_v3 = _mm256_blend_epi32(c13_v3, const_0_vec, 0b11000000);
    
    c14_v3 = _mm256_mul_epi32(_mm256_blend_epi32(a2_13, a_12_15, 0b00001100), _mm256_blend_epi32(a_12_15, const_0_vec, 0b00000011));
    c14_v3 = _mm256_blend_epi32(c14_v3, const_0_vec, 0b11000000);
    
    c15_v3 = _mm256_mul_epi32(a_12_15, _mm256_blend_epi32(a_12_15, const_0_vec, 0b11001111));

    __m256i r_v0, r_v1, r_v2, r_v3, r_v4, r_v5, r_v6, r_v7;  // real result vectors

    // offset 0 means:  0 1 2 3, 4 5 6 7, 8 9 10 11, ...
    __m256i offset0_v0, offset0_v1, offset0_v2, offset0_v3, offset0_v4, offset0_v5, offset0_v6;  // column 1 + 5 + 9 + 13

    offset0_v0 = c1_v0;
    offset0_v1 = c1_v1;
    offset0_v2 = _mm256_add_epi64(c1_v2, c5_v1);
    offset0_v3 = _mm256_add_epi64(c1_v3, c5_v2);
    offset0_v4 = _mm256_add_epi64(c5_v3, c9_v2);
    offset0_v5 = c9_v3;
    offset0_v6 = c13_v3;

    __m256i offset1_v0, offset1_v1, offset1_v2, offset1_v3, offset1_v4, offset1_v5, offset1_v6, tmp_offset1;  // column 2 + 6 + 10 + 14

    offset1_v0 = c2_v0;
    offset1_v1 = c2_v1;
    offset1_v2 = _mm256_add_epi64(c2_v2, c6_v1);
    offset1_v3 = _mm256_add_epi64(c2_v3, c6_v2);
    offset1_v4 = _mm256_add_epi64(c6_v3, c10_v2);
    offset1_v5 = c10_v3;
    offset1_v6 = c14_v3;

    __m256i offset2_v0, offset2_v1, offset2_v2, offset2_v3, offset2_v4, offset2_v5, offset2_v6, tmp_offset2;  // column 3 + 7 + 11 + 15

    offset2_v0 = c3_v0;
    offset2_v1 = c3_v1;
    offset2_v2 = _mm256_add_epi64(c3_v2, c7_v1);
    offset2_v3 = _mm256_add_epi64(c3_v3, c7_v2);
    offset2_v4 = _mm256_add_epi64(c7_v3, c11_v2);
    offset2_v5 = c11_v3;
    offset2_v6 = c15_v3;


    __m256i offset3_v0, offset3_v1, offset3_v2, offset3_v3, offset3_v4, offset3_v5;  // column 4 8 12 16

    offset3_v0 = c4_v0;
    offset3_v1 = c4_v1;
    offset3_v2 = _mm256_add_epi64(c4_v2, c8_v1);
    offset3_v3 = _mm256_add_epi64(c4_v3, c8_v2);
    offset3_v4 = _mm256_add_epi64(c8_v3, c12_v2);
    offset3_v5 = c12_v3;

    //  results :      v0      v1        v2         v3            v4            v5          v6            v7
    //  results :   0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 19, 20 21 22 23, 24 25 26 27,  28 29 30 31
    
    // offset 0 :      v0      v1        v2         v3            v4            v5          v6           v7           
    // offset 0 :   0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 19, 20 21 22 23, 24 25 26 27  28 29 30 31

    // offset 1 :      v0      v1        v2          v3            v4            v5          v6          
    // offset 1 :   1 2 3 4, 5 6 7 8, 9 10 11 12, 13 14 15 16, 17 18 19 20, 21 22 23 24, 25 26 27 28  
    // permuted :   4 1 2 3, 8 5 6 7, 12 9 10 11, 16 13 14 15, 20 17 18 19, 24 21 22 23, 28 25 26 27

    // offset 2 :      v0      v1        v2          v3            v4            v5          v6            
    // offset 2 :   2 3 4 5, 6 7 8 9, 10 11 12 13, 14 15 16 17, 18 19 20 21, 22 23 24 25, 26 27 28 29
    // permuted :   4 5 2 3, 8 9 6 7, 12 13 10 11, 16 17 14 15, 20 21 18 19, 24 25 22 23, 28 29 26 27

    // offset 3 :      v0      v1        v2          v3            v4            v5                      
    // offset 3 :   3 4 5 6, 7 8 9 10, 11 12 13 14, 15 16 17 18, 19 20 21 22, 23 24 25 26
    // permuted :   4 5 6 3, 8 9 10 7, 12 13 14 11, 16 17 18 15, 20 21 22 19, 24 25 26 23 

	// permute all offset vectors
	
	// offset 1
    __m256i offset1_v0_perm, offset1_v1_perm, offset1_v2_perm, offset1_v3_perm, offset1_v4_perm, offset1_v5_perm, offset1_v6_perm;
	offset1_v0_perm = _mm256_permute4x64_epi64( offset1_v0, 0b10010011 );
	offset1_v1_perm = _mm256_permute4x64_epi64( offset1_v1, 0b10010011 );
	offset1_v2_perm = _mm256_permute4x64_epi64( offset1_v2, 0b10010011 );
	offset1_v3_perm = _mm256_permute4x64_epi64( offset1_v3, 0b10010011 );
	offset1_v4_perm = _mm256_permute4x64_epi64( offset1_v4, 0b10010011 );
	offset1_v5_perm = _mm256_permute4x64_epi64( offset1_v5, 0b10010011 );
	offset1_v6_perm = _mm256_permute4x64_epi64( offset1_v6, 0b10010011 );

	// offset 2
    __m256i offset2_v0_perm, offset2_v1_perm, offset2_v2_perm, offset2_v3_perm, offset2_v4_perm, offset2_v5_perm, offset2_v6_perm;
	offset2_v0_perm = _mm256_permute4x64_epi64( offset2_v0, 0b01001110 );
	offset2_v1_perm = _mm256_permute4x64_epi64( offset2_v1, 0b01001110 );
	offset2_v2_perm = _mm256_permute4x64_epi64( offset2_v2, 0b01001110 );
	offset2_v3_perm = _mm256_permute4x64_epi64( offset2_v3, 0b01001110 );
	offset2_v4_perm = _mm256_permute4x64_epi64( offset2_v4, 0b01001110 );
	offset2_v5_perm = _mm256_permute4x64_epi64( offset2_v5, 0b01001110 );
	offset2_v6_perm = _mm256_permute4x64_epi64( offset2_v6, 0b01001110 );

	// offset 3
    __m256i offset3_v0_perm, offset3_v1_perm, offset3_v2_perm, offset3_v3_perm, offset3_v4_perm, offset3_v5_perm;
	offset3_v0_perm = _mm256_permute4x64_epi64( offset3_v0, 0b00111001 );
	offset3_v1_perm = _mm256_permute4x64_epi64( offset3_v1, 0b00111001 );
	offset3_v2_perm = _mm256_permute4x64_epi64( offset3_v2, 0b00111001 );
	offset3_v3_perm = _mm256_permute4x64_epi64( offset3_v3, 0b00111001 );
	offset3_v4_perm = _mm256_permute4x64_epi64( offset3_v4, 0b00111001 );
	offset3_v5_perm = _mm256_permute4x64_epi64( offset3_v5, 0b00111001 );
	
	// add offset 0 and offset 1
	r_v0 = _mm256_add_epi64( offset0_v0, _mm256_blend_epi32( const_0_vec    , offset1_v0_perm, 0b11111100 ) );
	r_v1 = _mm256_add_epi64( offset0_v1, _mm256_blend_epi32( offset1_v0_perm, offset1_v1_perm, 0b11111100 ) );
	r_v2 = _mm256_add_epi64( offset0_v2, _mm256_blend_epi32( offset1_v1_perm, offset1_v2_perm, 0b11111100 ) );
	r_v3 = _mm256_add_epi64( offset0_v3, _mm256_blend_epi32( offset1_v2_perm, offset1_v3_perm, 0b11111100 ) );
	r_v4 = _mm256_add_epi64( offset0_v4, _mm256_blend_epi32( offset1_v3_perm, offset1_v4_perm, 0b11111100 ) );
	r_v5 = _mm256_add_epi64( offset0_v5, _mm256_blend_epi32( offset1_v4_perm, offset1_v5_perm, 0b11111100 ) );
	r_v6 = _mm256_add_epi64( offset0_v6, _mm256_blend_epi32( offset1_v5_perm, offset1_v6_perm, 0b11111100 ) );
	
	// add offset 2 to result
	r_v0 = _mm256_add_epi64( r_v0, _mm256_blend_epi32( const_0_vec    , offset2_v0_perm, 0b11110000 ) );
	r_v1 = _mm256_add_epi64( r_v1, _mm256_blend_epi32( offset2_v0_perm, offset2_v1_perm, 0b11110000 ) );
	r_v2 = _mm256_add_epi64( r_v2, _mm256_blend_epi32( offset2_v1_perm, offset2_v2_perm, 0b11110000 ) );
	r_v3 = _mm256_add_epi64( r_v3, _mm256_blend_epi32( offset2_v2_perm, offset2_v3_perm, 0b11110000 ) );
	r_v4 = _mm256_add_epi64( r_v4, _mm256_blend_epi32( offset2_v3_perm, offset2_v4_perm, 0b11110000 ) );
	r_v5 = _mm256_add_epi64( r_v5, _mm256_blend_epi32( offset2_v4_perm, offset2_v5_perm, 0b11110000 ) );
	r_v6 = _mm256_add_epi64( r_v6, _mm256_blend_epi32( offset2_v5_perm, offset2_v6_perm, 0b11110000 ) );
	r_v7 = _mm256_blend_epi32( offset2_v6_perm, const_0_vec, 0b11111100 );

	// add offset 3 to result
	r_v0 = _mm256_add_epi64( r_v0, _mm256_blend_epi32( const_0_vec    , offset3_v0_perm, 0b11000000 ) );
	r_v1 = _mm256_add_epi64( r_v1, _mm256_blend_epi32( offset3_v0_perm, offset3_v1_perm, 0b11000000 ) );
	r_v2 = _mm256_add_epi64( r_v2, _mm256_blend_epi32( offset3_v1_perm, offset3_v2_perm, 0b11000000 ) );
	r_v3 = _mm256_add_epi64( r_v3, _mm256_blend_epi32( offset3_v2_perm, offset3_v3_perm, 0b11000000 ) );
	r_v4 = _mm256_add_epi64( r_v4, _mm256_blend_epi32( offset3_v3_perm, offset3_v4_perm, 0b11000000 ) );
	r_v5 = _mm256_add_epi64( r_v5, _mm256_blend_epi32( offset3_v4_perm, offset3_v5_perm, 0b11000000 ) );
	r_v6 = _mm256_add_epi64( r_v6, _mm256_blend_epi32( offset3_v5_perm, const_0_vec    , 0b11110000 ) );
	    
    // reduce
	reduce_8_to_4(&r_v0, &r_v1, &r_v2, &r_v3, &r_v4, &r_v5, &r_v6, &r_v7);
    reduce_4_to_2(&r_v0, &r_v1, &r_v2, &r_v3);
    reduce_2_to_2(&r_v0, &r_v1);
    // copy back the result vector
    copy_256(reslo, reshi, r_v0, r_v1);
	
}


void inverse_256(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi)
{
    /*
    res = x^-1 (the multiplicative inverse modulo p)
    */
    __m256i z2lo, z2hi, z9lo, z9hi, z11lo, z11hi, z2_5_0lo, z2_5_0hi, z2_10_0lo, z2_10_0hi, z2_20_0lo, z2_20_0hi, z2_50_0lo, z2_50_0hi, z2_100_0lo, z2_100_0hi, tlo, thi;
    int i;
    /* 2 */ square_256(&z2lo, &z2hi, xlo, xhi);
    /* 4 */ square_256(&tlo, &thi, z2lo, z2hi);
    /* 8 */ square_256(&tlo, &thi, tlo, thi);
    /* 9 */ mul_256(&z9lo, &z9hi,tlo, thi, xlo, xhi);
    /* 11 */ mul_256(&z11lo, &z11hi, z9lo, z9hi, z2lo, z2hi);
    /* 22 */ square_256(&tlo, &thi, z11lo, z11hi);
    /* 2^5 - 2^0 = 31 */ mul_256(&z2_5_0lo, &z2_5_0hi, tlo, thi, z9lo, z9hi);
    /* 2^6 - 2^1 */ square_256(&tlo, &thi, z2_5_0lo, z2_5_0hi);
    /* 2^10 - 2^5 */ for (i = 1; i < 5; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
    }
    /* 2^10 - 2^0 */ mul_256(&z2_10_0lo, &z2_10_0hi, tlo, thi, z2_5_0lo, z2_5_0hi);
    /* 2^11 - 2^1 */ square_256(&tlo, &thi, z2_10_0lo, z2_10_0hi);
    /* 2^20 - 2^10 */ for (i = 1; i < 10; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
    }
    /* 2^20 - 2^0 */ mul_256(&z2_20_0lo, &z2_20_0hi, tlo, thi, z2_10_0lo, z2_10_0hi);
    /* 2^21 - 2^1 */ square_256(&tlo, &thi, z2_20_0lo,z2_20_0hi);
    /* 2^40 - 2^20 */ for (i = 1; i < 20; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
    }
    /* 2^40 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_20_0lo, z2_20_0hi);
    /* 2^41 - 2^1 */ square_256(&tlo, &thi, tlo, thi);
    /* 2^50 - 2^10 */ for (i = 1; i < 10; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
    }
    /* 2^50 - 2^0 */ mul_256(&z2_50_0lo, &z2_50_0hi, tlo, thi, z2_10_0lo, z2_10_0hi);
    /* 2^51 - 2^1 */ square_256(&tlo, &thi, z2_50_0lo, z2_50_0hi);
    /* 2^100 - 2^50 */ for (i = 1; i < 50; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
    }
    /* 2^100 - 2^0 */ mul_256(&z2_100_0lo, &z2_100_0hi, tlo, thi, z2_50_0lo,z2_50_0hi);
    /* 2^101 - 2^1 */ square_256(&tlo, &thi, z2_100_0lo, z2_100_0hi);
    /* 2^200 - 2^100 */ for (i = 1; i < 100; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
    }
    /* 2^200 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_100_0lo, z2_100_0hi);
    /* 2^201 - 2^1 */ square_256(&tlo, &thi, tlo, thi);
    /* 2^250 - 2^50 */ for (i = 1; i < 50; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
    }
    /* 2^250 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_50_0lo, z2_50_0hi);
    /* 2^251 - 2^1 */ square_256(&tlo, &thi, tlo, thi);
    /* 2^252 - 2^2 */ square_256(&tlo, &thi, tlo, thi);
    /* 2^253 - 2^3 */ square_256(&tlo, &thi, tlo, thi);
    /* 2^254 - 2^4 */ square_256(&tlo, &thi, tlo, thi);
    /* 2^255 - 2^5 */ square_256(&tlo, &thi, tlo, thi);
    /* 2^255 - 21 */ mul_256(reslo, reshi, tlo, thi, z11lo, z11hi);  
}

void div_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi)
{
    /*
    res = a/b
    */
    // inverse and mult
    __m256i invlo, invhi;
    inverse_256(&invlo, &invhi, blo, bhi);
    mul_256(reslo, reshi, alo, ahi, invlo, invhi);
}


