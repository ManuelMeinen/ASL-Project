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


void square_256_non_vec(signed int reslo[8], signed int reshi[8], signed int xlo[8], signed int xhi[8]){
    /*
    res = x^2
    */
    unsigned long r[29]; 

    signed long a0 = xlo[0];
    signed long a1 = xlo[1];
    signed long a2 = xlo[2];
    signed long a3 = xlo[3];
    signed long a4 = xlo[4];
    signed long a5 = xlo[5];
    signed long a6 = xlo[6];
    signed long a7 = xlo[7];
    signed long a8 = xhi[0];
    signed long a9 = xhi[1];
    signed long a10 = xhi[2];
    signed long a11 = xhi[3];
    signed long a12 = xhi[4];
    signed long a13 = xhi[5];
    signed long a14 = xhi[6];

    signed long _2a0 = a0 << 1;
    signed long _2a1 = a1 << 1;
    signed long _2a2 = a2 << 1;
    signed long _2a3 = a3 << 1;
    signed long _2a4 = a4 << 1;
    signed long _2a5 = a5 << 1;
    signed long _2a6 = a6 << 1;
    signed long _2a7 = a7 << 1;
    signed long _2a8 = a8 << 1;
    signed long _2a9 = a9 << 1;
    signed long _2a10 = a10 << 1;
    signed long _2a11 = a11 << 1;
    signed long _2a12 = a12 << 1;
    signed long _2a13 = a13 << 1;


    r[0] = a0*a0;
    r[1] = _2a0*a1;
    r[2] = _2a0*a2 + a1*a1;
    r[3] = _2a0*a3 + _2a1*a2;
    r[4] = _2a0*a4 + _2a1*a3 + a2*a2;
    r[5] = _2a0*a5 + _2a1*a4 + _2a2*a3;
    r[6] = _2a0*a6 + _2a1*a5 + _2a2*a4 + a3*a3;
    r[7] = _2a0*a7 + _2a1*a6 + _2a2*a5 + _2a3*a4;
    r[8] = _2a0*a8 + _2a1*a7 + _2a2*a6 + _2a3*a5 + a4*a4;
    r[9] = _2a0*a9 + _2a1*a8 + _2a2*a7 + _2a3*a6 + _2a4*a5;
    r[10] = _2a0*a10 + _2a1*a9 + _2a2*a8 + _2a3*a7 + _2a4*a6 + a5*a5;
    r[11] = _2a0*a11 + _2a1*a10 + _2a2*a9 + _2a3*a8 + _2a4*a7 + _2a5*a6;
    r[12] = _2a0*a12 + _2a1*a11 + _2a2*a10 + _2a3*a9 + _2a4*a8 + _2a5*a7 + a6*a6;
    r[13] = _2a0*a13 + _2a1*a12 + _2a2*a11 + _2a3*a10 + _2a4*a9 + _2a5*a8 + _2a6*a7;
    r[14] = _2a0*a14 + _2a1*a13 + _2a2*a12 + _2a3*a11 + _2a4*a10 + _2a5*a9 + _2a6*a8 + a7*a7;
    r[15] = _2a1*a14 + _2a2*a13 + _2a3*a12 + _2a4*a11 + _2a5*a10 + _2a6*a9 + _2a7*a8;
    r[16] = _2a2*a14 + _2a3*a13 + _2a4*a12 + _2a5*a11 + _2a6*a10 + _2a7*a9 + a8*a8;
    r[17] = _2a3*a14 + _2a4*a13 + _2a5*a12 + _2a6*a11 + _2a7*a10 + _2a8*a9;
    r[18] = _2a4*a14 + _2a5*a13 + _2a6*a12 + _2a7*a11 + _2a8*a10 + a9*a9;
    r[19] = _2a5*a14 + _2a6*a13 + _2a7*a12 + _2a8*a11 + _2a9*a10;
    r[20] = _2a6*a14 + _2a7*a13 + _2a8*a12 + _2a9*a11 + a10*a10;
    r[21] = _2a7*a14 + _2a8*a13 + _2a9*a12 + _2a10*a11;
    r[22] = _2a8*a14 + _2a9*a13 + _2a10*a12 + a11*a11;
    r[23] = _2a9*a14 + _2a10*a13 + _2a11*a12;
    r[24] = _2a10*a14 + _2a11*a13 + a12*a12;
    r[25] = _2a11*a14 + _2a12*a13;
    r[26] = _2a12*a14 + a13*a13;
    r[27] = _2a13*a14;
    r[28] = a14*a14;

    reduce_8_to_4_non_vec(r);

	reduce_4_to_2_non_vec(r);
	
    for(int i = 0; i<8; i++){
        reslo[i] = r[i];
        reshi[i] = r[i+8];
    }
    reshi[7] = 0;
    
    reduce_2_to_2_non_vec(reslo, reshi);
}
void reduce_8_to_4_non_vec(unsigned long r[29]){
    /*
    reduce r
    */
    for(int i = 0; i<NO_LIMBS-1; i++){
        //r[i] += r[i+NO_LIMBS] = times19(r[i+NO_LIMBS]);
        r[i] += times19(r[i+NO_LIMBS]);
    }
}

void reduce_4_to_2_non_vec(unsigned long r[29]){
    /*
    reduce r
    */
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
    r0 += times19(carry);
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

void reduce_2_to_2_non_vec(signed int xlo[8], signed int xhi[8]){
    /*
    reduce x
    */
    signed long carry;
    signed int d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14;

    d0 = xlo[0];
    d1 = xlo[1];
    d2 = xlo[2];
    d3 = xlo[3];
    d4 = xlo[4];
    d5 = xlo[5];
    d6 = xlo[6];
    d7 = xlo[7];
    d8 = xhi[0];
    d9 = xhi[1];
    d10 = xhi[2];
    d11 = xhi[3];
    d12 = xhi[4];
    d13 = xhi[5];
    d14 = xhi[6];
    
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

    xlo[0] = d0;
    xlo[1] = d1;
    xlo[2] = d2;
    xlo[3] = d3;
    xlo[4] = d4;
    xlo[5] = d5;
    xlo[6] = d6;
    xlo[7] = d7;
    xhi[0] = d8;
    xhi[1] = d9;
    xhi[2] = d10;
    xhi[3] = d11;
    xhi[4] = d12;
    xhi[5] = d13;
    xhi[6] = d14; 
}

void square_256(__m256i *reslo, __m256i *reshi, __m256i xlo,__m256i xhi) {
    
    signed int reslo_arr[8], reshi_arr[8];

    // read input vectors into arrays
    signed int* xlo_arr = (int*) &xlo;
    signed int* xhi_arr = (int*) &xhi;
    
    // square arrays
    square_256_non_vec( reslo_arr, reshi_arr, xlo_arr, xhi_arr );
    
    // store results back into vectors
    *reslo = _mm256_loadu_si256( (const __m256i*) reslo_arr); 
    *reshi = _mm256_loadu_si256( (const __m256i*) reshi_arr); 
    
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
    // to array
    signed int* alo = (int*) &tlo;
    signed int* ahi = (int*) &thi;
    /* 2^10 - 2^5 */ for (i = 1; i < 5; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256( (const __m256i*) alo); 
    thi = _mm256_loadu_si256( (const __m256i*) ahi); 

    /* 2^10 - 2^0 */ mul_256(&z2_10_0lo, &z2_10_0hi, tlo, thi, z2_5_0lo, z2_5_0hi);
    /* 2^11 - 2^1 */ square_256(&tlo, &thi, z2_10_0lo, z2_10_0hi);
    // to array
    alo = (int*) &tlo;
    ahi = (int*) &thi;
    /* 2^20 - 2^10 */ for (i = 1; i < 10; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256( (const __m256i*) alo); 
    thi = _mm256_loadu_si256( (const __m256i*) ahi); 
    /* 2^20 - 2^0 */ mul_256(&z2_20_0lo, &z2_20_0hi, tlo, thi, z2_10_0lo, z2_10_0hi);
    /* 2^21 - 2^1 */ square_256(&tlo, &thi, z2_20_0lo,z2_20_0hi);
    // to array
    alo = (int*) &tlo;
    ahi = (int*) &thi;
    /* 2^40 - 2^20 */ for (i = 1; i < 20; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256( (const __m256i*) alo); 
    thi = _mm256_loadu_si256( (const __m256i*) ahi); 
    /* 2^40 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_20_0lo, z2_20_0hi);
    /* 2^41 - 2^1 */ //square_256(&tlo, &thi, tlo, thi);
    // to array
    alo = (int*) &tlo;
    ahi = (int*) &thi;
    /* 2^50 - 2^10 */ for (i = 1; i < 11; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    /* 2^50 - 2^0 */ mul_256(&z2_50_0lo, &z2_50_0hi, tlo, thi, z2_10_0lo, z2_10_0hi);
    /* 2^51 - 2^1 */ square_256(&tlo, &thi, z2_50_0lo, z2_50_0hi);
    // to array
    alo = (int*) &tlo;
    ahi = (int*) &thi;
    /* 2^100 - 2^50 */ for (i = 1; i < 50; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256( (const __m256i*) alo); 
    thi = _mm256_loadu_si256( (const __m256i*) ahi); 
    /* 2^100 - 2^0 */ mul_256(&z2_100_0lo, &z2_100_0hi, tlo, thi, z2_50_0lo,z2_50_0hi);
    /* 2^101 - 2^1 */ square_256(&tlo, &thi, z2_100_0lo, z2_100_0hi);
    // to array
    alo = (int*) &tlo;
    ahi = (int*) &thi;
    /* 2^200 - 2^100 */ for (i = 1; i < 100; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256( (const __m256i*) alo); 
    thi = _mm256_loadu_si256( (const __m256i*) ahi); 
    /* 2^200 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_100_0lo, z2_100_0hi);
    /* 2^201 - 2^1 */ //square_256(&tlo, &thi, tlo, thi);
    // to array
    alo = (int*) &tlo;
    ahi = (int*) &thi;
    /* 2^250 - 2^50 */ for (i = 1; i < 51; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256( (const __m256i*) alo); 
    thi = _mm256_loadu_si256( (const __m256i*) ahi); 
    /* 2^250 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_50_0lo, z2_50_0hi);
    // to array
    alo = (int*) &tlo;
    ahi = (int*) &thi;
    /* 2^10 - 2^5 */ for (i = 0; i < 5; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256( (const __m256i*) alo); 
    thi = _mm256_loadu_si256( (const __m256i*) ahi); 
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


