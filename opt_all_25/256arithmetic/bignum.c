#include <stdio.h>
#include "bignum.h"
#include <immintrin.h>
#include <emmintrin.h>
#include "../../_common_libs/radix25_vec.h"
#include "../../_common_libs/vec_utils.h"

const int NO_LIMBS = 10;
const int LIMB_SIZE = 25;

unsigned long times19(unsigned long a)
{
    /*
    multiply a by 19
    */
    return (a << 4) + (a << 1) + a;
}

void reduce_2_to_2(__m256i *xlo, __m256i *xhi)
{

    signed long carry;
    unsigned long d0, d1, d2, d3, d4, d5, d6, d7, d8, d9;

    int *d = (int *)xlo;
    int *i_v2 = (int *)xhi;
    // d0 = i_v[0];
    // d1 = i_v[1];
    // d2 = i_v[2];
    // d3 = i_v[3];
    // d4 = i_v[4];
    // d5 = i_v[5];
    // d6 = i_v[6];
    // d7 = i_v[7];
    d[8] = i_v2[0];
    d[9] = i_v2[1];

    for (int i = 0; i < NO_LIMBS - 2; i += 2)
    {
        carry = d[i] >> (LIMB_SIZE + 1); // 26
        d[i + 1] += carry;
        carry <<= (LIMB_SIZE + 1);
        d[i] -= carry;

        carry = d[i + 1] >> LIMB_SIZE; // 25
        d[i + 2] += carry;
        carry <<= LIMB_SIZE;
        d[i + 1] -= carry;
    }

    carry = d[8] >> (LIMB_SIZE + 1); // 26
    d[9] += carry;
    carry <<= (LIMB_SIZE + 1);
    d[8] -= carry;

    carry = d[NO_LIMBS - 1] >> LIMB_SIZE; //25
    d[0] += 19 * carry;
    carry <<= LIMB_SIZE;
    d[NO_LIMBS - 1] -= carry;

    *xlo = _mm256_setr_epi32(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]);
    *xhi = _mm256_setr_epi32(d[8], d[9], 0, 0, 0, 0, 0, 0);
}

void reduce_3_to_2(__m256i *rv0, __m256i *rv1, __m256i *rv2)
{
    /*
    reduce 3 vectors into 2
    */

    unsigned long carry;
    unsigned long r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14;
    unsigned long *i_v0 = (unsigned long *)rv0;
    unsigned long *i_v1 = (unsigned long *)rv1;
    unsigned long *i_v2 = (unsigned long *)rv2;
    unsigned long r[20];
    r[0] = i_v0[0];
    r[1] = i_v0[1];
    r[2] = i_v0[2];
    r[3] = i_v0[3];
    r[4] = i_v1[0];
    r[5] = i_v1[1];
    r[6] = i_v1[2];
    r[7] = i_v1[3];
    r[8] = i_v2[0];
    r[9] = i_v2[1];
    // for(int i =0;i<10;i++){
    //     printf("%ld, ", r[i]);
    // }
    // printf("\n");
     int n = 2;
     while(n>0){
        for (int i = 0; i < NO_LIMBS - 2; i += 2)
        {
            carry = r[i] >> (LIMB_SIZE + 1); //26
            r[i + 1] += carry;
            carry <<= (LIMB_SIZE + 1);
            r[i] -= carry;

            carry = r[i + 1] >> LIMB_SIZE; //25
            r[i + 2] += carry;
            carry <<= LIMB_SIZE;
            r[i + 1] -= carry;

        }

        carry = r[8] >> (LIMB_SIZE + 1); //26
        r[9] += carry;
        carry <<= (LIMB_SIZE + 1);
        r[8] -= carry;

        carry = r[NO_LIMBS - 1] >> LIMB_SIZE; //25
        r[0] += 19 * carry;
        carry <<= LIMB_SIZE;
        r[NO_LIMBS - 1] -= carry;
        n--;
        
     }

// printf("sdasdadadads");
//     for(int i =0;i<10;i++){
//         printf("%ld, ", r[i]);
//     }
    // printf("\n");
   // printf("*********************%ud \n",(unsigned int) r[0]);
    *rv0 = _mm256_setr_epi32((unsigned  int)r[0], (unsigned  int)r[1], (unsigned int)r[2], (unsigned int)r[3], (unsigned int)r[4], (unsigned int)r[5], (unsigned int)r[6], (unsigned int)r[7]);
     *rv1 = _mm256_setr_epi32((unsigned int)r[8], (unsigned int)r[9], 0, 0, 0, 0, 0, 0);
//    print_vec_int(_mm256_setr_epi32((unsigned  int)r[0], (unsigned int)r[1], (unsigned int)r[2], (unsigned int)r[3], (unsigned int)r[4], (unsigned int)r[5], (unsigned int)r[6], (unsigned int)r[7]), "rv0");
//    print_vec_int(_mm256_setr_epi32((unsigned int)r[8], (unsigned int)r[9], 0, 0, 0, 0, 0, 0),"rv1");
}

void add_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi)
{
    /*
    res = a + b
    */
    // a + b
    *reslo = _mm256_add_epi32(alo, blo);
    *reshi = _mm256_add_epi32(ahi, bhi);
    // reduce res
    reduce_2_to_2(reslo, reshi);
}

void sub_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi)
{
    /*
    res = a - b
    */

    // a - b
    *reslo = _mm256_sub_epi32(alo, blo);
    *reshi = _mm256_sub_epi32(ahi, bhi);
    int mask = _mm256_movemask_ps((__m256)*reshi);
    while(_mm256_movemask_ps((__m256)*reslo) != 0 || ( mask != 0 && mask != 0b10000000 )){ //reduce if there are negative numbers in the limbs only
        reduce_2_to_2(reslo, reshi); 
        mask = _mm256_movemask_ps((__m256)*reshi);
    }
    // reduce res
    // reduce_2_to_2(reslo, reshi);
    // reduce_2_to_2(reslo, reshi);
}

void mul_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi)
{

    /*r[0]  =     a[0] * b[0]                                                                                                                                                                                        
r[1]  =       a[0] * b[1]    +       a[1] * b[0]                                                                                                                                                                          
r[2]  =       a[0] * b[2]    +  2 *  a[1] * b[1]    +       a[2] * b[0]                                                                                                                                                      
r[3]  =       a[0] * b[3]    +       a[1] * b[2]    +       a[2] * b[1]   +       a[3] * b[0]                                                                                                                                
r[4]  =       a[0] * b[4]    +  2 *  a[1] * b[3]    +       a[2] * b[2]   +  2 *  a[3] * b[1]   +       a[4] * b[0]                                                                                                          
r[5]  =       a[0] * b[5]    +       a[1] * b[4]    +       a[2] * b[3]   +       a[3] * b[2]   +       a[4] * b[1]    +       a[5] * b[0]                                                                                   
r[6]  =       a[0] * b[6]    +  2 *  a[1] * b[5]    +       a[2] * b[4]   +  2 *  a[3] * b[3]   +       a[4] * b[2]    +  2  * a[5] * b[1]    +       a[6] * b[0]                                                              
r[7]  =       a[0] * b[7]    +       a[1] * b[6]    +       a[2] * b[5]   +       a[3] * b[4]   +       a[4] * b[3]    +       a[5] * b[2]    +       a[6] * b[1]    +       a[7] * b[0]                                         
r[8]  =       a[0] * b[8]    +  2 *  a[1] * b[7]    +       a[2] * b[6]   +  2 *  a[3] * b[5]   +       a[4] * b[4]    +  2  * a[5] * b[3]    +       a[6] * b[2]    +  2  * a[7] * b[1]     +       a[8] * b[0]                    
r[9]  =       a[0] * b[9]    +       a[1] * b[8]    +       a[2] * b[7]   +       a[3] * b[6]   +       a[4] * b[5]    +       a[5] * b[4]    +       a[6] * b[3]    +       a[7] * b[2]     +       a[8] * b[1]     +       a[9] * b[0];
r[10] =                         38 * a[1] * b[9]    +  19 * a[2] * b[8]   +  38 * a[3] * b[7]   +  19 * a[4] * b[6]    +  38 * a[5] * b[5]    +  19 * a[6] * b[4]    +  38 * a[7] * b[3]     +  19 * a[8] * b[2]     +  38 * a[9] * b[1];
r[11] =                                                19 * a[2] * b[9]   +  19 * a[3] * b[8]   +  19 * a[4] * b[7]    +  19 * a[5] * b[6]    +  19 * a[6] * b[5]    +  19 * a[7] * b[4]     +  19 * a[8] * b[3]     +  19 * a[9] * b[2];
r[12] =                                                                      38 * a[3] * b[9]   +  19 * a[4] * b[8]    +  38 * a[5] * b[7]    +  19 * a[6] * b[6]    +  38 * a[7] * b[5]     +  19 * a[8] * b[4]     +  38 * a[9] * b[3];
r[13] =                                                                                            19 * a[4] * b[9]    +  19 * a[5] * b[8]    +  19 * a[6] * b[7]    +  19 * a[7] * b[6]     +  19 * a[8] * b[5]     +  19 * a[9] * b[4];
r[14] =                                                                                                                   38 * a[5] * b[9]    +  19 * a[6] * b[8]    +  38 * a[7] * b[7]     +  19 * a[8] * b[6]     +  38 * a[9] * b[5];
r[15] =                                                                                                                                          19 * a[6] * b[9]    +  19 * a[7] * b[8]     +  19 * a[8] * b[7]     +  19 * a[9] * b[6];
r[16] =                                                                                                                                                                 38 * a[7] * b[9]     +  19 * a[8] * b[8]     +  38 * a[9] * b[7];
r[17] =                                                                                                                                                                                         19 * a[8] * b[9]     +  19 * a[9] * b[8];
r[18] =                                                                                                                                                                                                                 38 * a[9] * b[9];
*/

    __m256i b_0_3 = _mm256_permutevar8x32_epi32(blo, index1_vec);
    __m256i b_4_7 = _mm256_permutevar8x32_epi32(blo, index2_vec);
    __m256i b_8_11 = _mm256_permutevar8x32_epi32(bhi, index1_vec);

    __m256i a_0 = _mm256_broadcastq_epi64(_mm256_castsi256_si128(alo));
    __m256i a_1 = _mm256_permutevar8x32_epi32(a_0, index_alt_vec);
    __m256i a_2 = _mm256_broadcastq_epi64(_mm256_castsi256_si128(_mm256_unpackhi_epi64(alo, index1_vec)));
    __m256i a_3 = _mm256_permutevar8x32_epi32(a_2, index_alt_vec);

    __m128i a_lo_hi = _mm256_extracti128_si256(alo, 1);
    __m256i a_4 = _mm256_broadcastq_epi64(a_lo_hi);
    __m256i a_5 = _mm256_permutevar8x32_epi32(a_4, index_alt_vec);
    __m256i a_6 = _mm256_broadcastq_epi64(_mm256_castsi256_si128(_mm256_unpackhi_epi64(_mm256_castsi128_si256(a_lo_hi), index1_vec)));
    __m256i a_7 = _mm256_permutevar8x32_epi32(a_6, index_alt_vec);

    __m256i a_8 = _mm256_broadcastq_epi64(_mm256_castsi256_si128(ahi));
    __m256i a_9 = _mm256_permutevar8x32_epi32(a_8, index_alt_vec);

    __m256i r_v0, r_v1, r_v2, r_v3, r_v4, r_v5, r_v6, r_v7; // real result vectors

    // Naming conventions: cx_vy -> column x result vector y
    __m256i c1_v0, c1_v1, c1_v2;
    __m256i c2_v0, c2_v1, c2_v2;
    __m256i c3_v0, c3_v1, c3_v2;
    __m256i c4_v0, c4_v1, c4_v2;
    __m256i c5_v0, c5_v1, c5_v2;
    __m256i c6_v0, c6_v1, c6_v2;
    __m256i c7_v0, c7_v1, c7_v2;
    __m256i c8_v0, c8_v1, c8_v2;
    __m256i c9_v0, c9_v1, c9_v2;
    __m256i c10_v0, c10_v1, c10_v2;

    // precomputations
    __m256i shifted_b_0_3 = _mm256_sllv_epi64(b_0_3, c_alt_0_1__64_vec);
    __m256i mul_b_8_11__2x19_2x0 = _mm256_mul_epi32(b_8_11, c_2x19_2x0__64_vec);
    __m256i mul_b_8_11__19_38_2x0 = _mm256_mul_epi32(b_8_11, c_19_38_2x0__64_vec);
    __m256i mul_b_4_7__4x19 = _mm256_mul_epi32(b_4_7, c_all_19__64_vec);
    __m256i mul_b_4_7__alt_19_38 = _mm256_mul_epi32(b_4_7, c_alt_19_38__64_vec);

    /* DEBUGGING HELP
	// FOR DEBUGGING
    signed int* alo_arr = (int*) &alo;
    signed int* ahi_arr = (int*) &ahi;
    signed int* blo_arr = (int*) &blo;
    signed int* bhi_arr = (int*) &bhi;

	unsigned long a0 = (unsigned long) alo_arr[0];
	unsigned long a1 = (unsigned long) alo_arr[1];
	unsigned long a2 = (unsigned long) alo_arr[2];
	unsigned long a3 = (unsigned long) alo_arr[3];
	unsigned long a4 = (unsigned long) alo_arr[4];
	unsigned long a5 = (unsigned long) alo_arr[5];
	unsigned long a6 = (unsigned long) alo_arr[6];
	unsigned long a7 = (unsigned long) alo_arr[7];
	unsigned long a8 = (unsigned long) ahi_arr[0];
	unsigned long a9 = (unsigned long) ahi_arr[1];
	unsigned long b0 = (unsigned long) blo_arr[0];
	unsigned long b1 = (unsigned long) blo_arr[1];
	unsigned long b2 = (unsigned long) blo_arr[2];
	unsigned long b3 = (unsigned long) blo_arr[3];
	unsigned long b4 = (unsigned long) blo_arr[4];
	unsigned long b5 = (unsigned long) blo_arr[5];
	unsigned long b6 = (unsigned long) blo_arr[6];
	unsigned long b7 = (unsigned long) blo_arr[7];
	unsigned long b8 = (unsigned long) bhi_arr[0];
	unsigned long b9 = (unsigned long) bhi_arr[1];
	// END DEBUGGING VARIABLES ASSIGNMENTS (TODO: DELETE ABOVE)
	
	unsigned long r0  = a0*b0;
	unsigned long r1  = a0*b1 +   a1*b0;
	unsigned long r2  = a0*b2 + 2*a1*b1 +    a2*b0;
	unsigned long r3  = a0*b3 +   a1*b2 +    a2*b1 +   a3*b0;
	unsigned long r4  = a0*b4 + 2*a1+b3 +    a2*b2 + 2*a3*b1 +    a4*b0;
	unsigned long r5  = a0*b5 +   a1*b4 +    a2*b3 +   a3*b2 +    a4*b1 + a5*b0;
	unsigned long r6  = a0*b6 + 2*a1*b5 +    a2*b4 + 2*a3*b3 +    a4*b2 + 2*a5*b1 +    a6*b0;
	unsigned long r7  = a0*b7 +   a1*b6 +    a2*b5 +   a3*b4 +    a4*b3 +   a5*b2 +    a6*b1 + a7*b0;
	unsigned long r8  = a0*b8 + 2*a1*b7 +    a2*b6 + 2*a3*b5 +    a4*b4 + 2*a5*b3 +    a6*b2 + 2*a7*b1 +    a8*b0;
	unsigned long r9  = a0*b9 +   a1*b8 +    a2*b7 +   a3*b6 +    a4*b5 +   a5*b4 +    a6*b3 +   a7*b2 +    a8*b1 +    a9*b0;
	unsigned long r10 =        38*a1*b9 + 19*a2*b8 +38*a3*b7 + 19*a4*b6 +38*a5*b5 + 19*a6*b4 +38*a7*b3 + 19*a8*b2 + 38*a9*b1;
	unsigned long r11 =                   19*a2*b9 +19*a3*b8 + 19*a4*b7 +19*a5*b6 + 19*a6*b5 +19*a7*b4 + 19*a8*b3 + 19*a9*b2;
	unsigned long r12 =                             38*a3*b9 + 19*a4*b8 +38*a5*b7 + 19*a6*b6 +38*a7*b5 + 19*a8*b4 + 38*a9*b3;
	unsigned long r13 =                                        19*a4*b9 +19*a5*b8 + 19*a6*b7 +19*a7*b6 + 19*a8*b5 + 19*a9*b4;
	unsigned long r14 =                                                  38*a5*b9 + 19*a6*b8 +38*a7*b7 + 19*a8*b6 + 38*a9*b5;
	unsigned long r15 =                                                             19*a6*b9 +19*a7*b8 + 19*a8*b7 + 19*a9*b6;
	unsigned long r16 =                                                                       38*a7*b9 + 19*a8*b8 + 38*a9*b7;
	unsigned long r17 =                                                                                  19*a8*b9 + 19*a9*b8;
	unsigned long r18 =                                                                                             38*a9*b9;
	printf( "r_v0 = %ld, %ld, %ld, %ld\n", r0, r1, r2, r3 );
	printf( "r_v1 = %ld, %ld, %ld, %ld\n", r4, r5, r6, r7 );
	printf( "r_v2 = %ld, %ld, %ld, %ld\n", r8, r9, r10, r11 );
	printf( "r_v3 = %ld, %ld, %ld, %ld\n", r12, r13, r14, r15 );
	printf( "r_v4 = %ld, %ld, %ld, %ld\n", r16, r17, r18, (long) 0 );
	*/

    c1_v0 = _mm256_mul_epi32(a_0, b_0_3);
    c1_v1 = _mm256_mul_epi32(a_0, b_4_7);
    c1_v2 = _mm256_mul_epi32(a_0, b_8_11);
    c1_v2 = _mm256_blend_epi32(c1_v2, const_0_vec, 0b11110000);
    /* --> OK
    check_vec_long(c1_v0, (unsigned long)alo_arr[0]*blo_arr[0],(unsigned long)alo_arr[0]*blo_arr[1],(unsigned long)alo_arr[0]*blo_arr[2],(unsigned long)alo_arr[0]*blo_arr[3]);
    check_vec_long(c1_v1, (unsigned long)alo_arr[0]*blo_arr[4],(unsigned long)alo_arr[0]*blo_arr[5],(unsigned long)alo_arr[0]*blo_arr[6],(unsigned long)alo_arr[0]*blo_arr[7]);
    check_vec_long(c1_v2, (unsigned long)alo_arr[0]*bhi_arr[0],(unsigned long)alo_arr[0]*bhi_arr[1],0,0);
    */

    c2_v1 = _mm256_sllv_epi64(b_4_7, c_alt_0_1__64_vec);
    c2_v2 = _mm256_mul_epi32(b_8_11, c_1_38_2x0__64_vec);
    c2_v0 = _mm256_mul_epi32(a_1, shifted_b_0_3);
    c2_v1 = _mm256_mul_epi32(a_1, c2_v1);
    c2_v2 = _mm256_mul_epi32(a_1, c2_v2);
    /* --> OK
    check_vec_long(c2_v0, (unsigned long)alo_arr[1]*blo_arr[0],(unsigned long)2*alo_arr[1]*blo_arr[1],(unsigned long)alo_arr[1]*blo_arr[2],(unsigned long)2*alo_arr[1]*blo_arr[3]);
    check_vec_long(c2_v1, (unsigned long)alo_arr[1]*blo_arr[4],(unsigned long)2*alo_arr[1]*blo_arr[5],(unsigned long)alo_arr[1]*blo_arr[6],(unsigned long)2*alo_arr[1]*blo_arr[7]);
    check_vec_long(c2_v2, (unsigned long)alo_arr[1]*bhi_arr[0],(unsigned long)38*alo_arr[1]*bhi_arr[1],0,0);
    */

    c3_v0 = _mm256_mul_epi32(a_2, b_0_3);
    c3_v1 = _mm256_mul_epi32(a_2, b_4_7);
    c3_v2 = _mm256_mul_epi32(a_2, mul_b_8_11__2x19_2x0);
    /* --> OK
    check_vec_long(c3_v0, (unsigned long)alo_arr[2]*blo_arr[0],(unsigned long)alo_arr[2]*blo_arr[1],(unsigned long)alo_arr[2]*blo_arr[2],(unsigned long)alo_arr[2]*blo_arr[3]);
    check_vec_long(c3_v1, (unsigned long)alo_arr[2]*blo_arr[4],(unsigned long)alo_arr[2]*blo_arr[5],(unsigned long)alo_arr[2]*blo_arr[6],(unsigned long)alo_arr[2]*blo_arr[7]);
    check_vec_long(c3_v2, (unsigned long)19*alo_arr[2]*bhi_arr[0],(unsigned long)19*alo_arr[2]*bhi_arr[1],0,0);
    */

    c4_v0 = _mm256_mul_epi32(a_3, shifted_b_0_3);
    c4_v1 = _mm256_mul_epi32(b_4_7, c_1_2_1_38__64_vec);
    c4_v1 = _mm256_mul_epi32(a_3, c4_v1);
    c4_v2 = _mm256_mul_epi32(a_3, mul_b_8_11__19_38_2x0);
    /* --> OK
    check_vec_long(c4_v0, (unsigned long)alo_arr[3]*blo_arr[0],(unsigned long)2*alo_arr[3]*blo_arr[1],(unsigned long)alo_arr[3]*blo_arr[2],(unsigned long)2*alo_arr[3]*blo_arr[3]);
    check_vec_long(c4_v1, (unsigned long)alo_arr[3]*blo_arr[4],(unsigned long)2*alo_arr[3]*blo_arr[5],(unsigned long)alo_arr[3]*blo_arr[6],(unsigned long)38*alo_arr[3]*blo_arr[7]);
    check_vec_long(c4_v2, (unsigned long)19*alo_arr[3]*bhi_arr[0],(unsigned long)38*alo_arr[3]*bhi_arr[1],0,0);
	*/

    c5_v0 = _mm256_mul_epi32(a_4, b_0_3);
    c5_v1 = _mm256_mul_epi32(b_4_7, c_2x1_2x19__64_vec);
    c5_v1 = _mm256_mul_epi32(a_4, c5_v1);
    c5_v2 = _mm256_mul_epi32(a_4, mul_b_8_11__2x19_2x0);
    /* --> OK
    check_vec_long(c5_v0, (unsigned long)alo_arr[4]*blo_arr[0],(unsigned long)alo_arr[4]*blo_arr[1],(unsigned long)alo_arr[4]*blo_arr[2],(unsigned long)alo_arr[4]*blo_arr[3]);
    check_vec_long(c5_v1, (unsigned long)alo_arr[4]*blo_arr[4],(unsigned long)alo_arr[4]*blo_arr[5],(unsigned long)19*alo_arr[4]*blo_arr[6],(unsigned long)19*alo_arr[4]*blo_arr[7]);
    check_vec_long(c5_v2, (unsigned long)19*alo_arr[4]*bhi_arr[0],(unsigned long)19*alo_arr[4]*bhi_arr[1],0,0);
	*/

    c6_v1 = _mm256_mul_epi32(b_4_7, c_1_38_19_38__64_vec);
    c6_v0 = _mm256_mul_epi32(a_5, shifted_b_0_3);
    c6_v1 = _mm256_mul_epi32(a_5, c6_v1);
    c6_v2 = _mm256_mul_epi32(a_5, mul_b_8_11__19_38_2x0);
    /* --> OK
    check_vec_long(c6_v0, (unsigned long)alo_arr[5]*blo_arr[0],(unsigned long)2*alo_arr[5]*blo_arr[1],(unsigned long)alo_arr[5]*blo_arr[2],(unsigned long)2*alo_arr[5]*blo_arr[3]);
    check_vec_long(c6_v1, (unsigned long)alo_arr[5]*blo_arr[4],(unsigned long)38*alo_arr[5]*blo_arr[5],(unsigned long)19*alo_arr[5]*blo_arr[6],(unsigned long)38*alo_arr[5]*blo_arr[7]);
    check_vec_long(c6_v2, (unsigned long)19*alo_arr[5]*bhi_arr[0],(unsigned long)38*alo_arr[5]*bhi_arr[1],0,0);
	*/

    c7_v0 = _mm256_mul_epi32(a_6, b_0_3);
    c7_v1 = _mm256_mul_epi32(a_6, mul_b_4_7__4x19);
    c7_v2 = _mm256_mul_epi32(a_6, mul_b_8_11__2x19_2x0);
    /* --> OK
    check_vec_long(c7_v0, (unsigned long)alo_arr[6]*blo_arr[0],(unsigned long)alo_arr[6]*blo_arr[1],(unsigned long)alo_arr[6]*blo_arr[2],(unsigned long)alo_arr[6]*blo_arr[3]);
    check_vec_long(c7_v1, (unsigned long)19*alo_arr[6]*blo_arr[4],(unsigned long)19*alo_arr[6]*blo_arr[5],(unsigned long)19*alo_arr[6]*blo_arr[6],(unsigned long)19*alo_arr[6]*blo_arr[7]);
    check_vec_long(c7_v2, (unsigned long)19*alo_arr[6]*bhi_arr[0],(unsigned long)19*alo_arr[6]*bhi_arr[1],0,0);
	*/

    c8_v0 = _mm256_mul_epi32(b_0_3, c_1_2_1_38__64_vec);
    c8_v0 = _mm256_mul_epi32(a_7, c8_v0);
    c8_v1 = _mm256_mul_epi32(a_7, mul_b_4_7__alt_19_38);
    c8_v2 = _mm256_mul_epi32(a_7, mul_b_8_11__19_38_2x0);
    /* --> OK
    check_vec_long(c8_v0, (unsigned long)alo_arr[7]*blo_arr[0],(unsigned long)2*alo_arr[7]*blo_arr[1],(unsigned long)alo_arr[7]*blo_arr[2],(unsigned long)38*alo_arr[7]*blo_arr[3]);
    check_vec_long(c8_v1, (unsigned long)19*alo_arr[7]*blo_arr[4],(unsigned long)38*alo_arr[7]*blo_arr[5],(unsigned long)19*alo_arr[7]*blo_arr[6],(unsigned long)38*alo_arr[7]*blo_arr[7]);
    check_vec_long(c8_v2, (unsigned long)19*alo_arr[7]*bhi_arr[0],(unsigned long)38*alo_arr[7]*bhi_arr[1],0,0);
	*/

    c9_v0 = _mm256_mul_epi32(b_0_3, c_2x1_2x19__64_vec);
    c9_v0 = _mm256_mul_epi32(a_8, c9_v0);
    c9_v1 = _mm256_mul_epi32(a_8, mul_b_4_7__4x19);
    c9_v2 = _mm256_mul_epi32(a_8, mul_b_8_11__2x19_2x0);
    /* --> OK
    check_vec_long(c9_v0, (unsigned long)ahi_arr[0]*blo_arr[0],(unsigned long)ahi_arr[0]*blo_arr[1],(unsigned long)19*ahi_arr[0]*blo_arr[2],(unsigned long)19*ahi_arr[0]*blo_arr[3]);
    check_vec_long(c9_v1, (unsigned long)19*ahi_arr[0]*blo_arr[4],(unsigned long)19*ahi_arr[0]*blo_arr[5],(unsigned long)19*ahi_arr[0]*blo_arr[6],(unsigned long)19*ahi_arr[0]*blo_arr[7]);
    check_vec_long(c9_v2, (unsigned long)19*ahi_arr[0]*bhi_arr[0],(unsigned long)19*ahi_arr[0]*bhi_arr[1],0,0);
    */

    c10_v0 = _mm256_mul_epi32(b_0_3, c_1_38_19_38__64_vec);
    c10_v0 = _mm256_mul_epi32(a_9, c10_v0);
    c10_v1 = _mm256_mul_epi32(a_9, mul_b_4_7__alt_19_38);
    c10_v2 = _mm256_mul_epi32(a_9, mul_b_8_11__19_38_2x0);
    /* --> OK
    check_vec_long(c10_v0, (unsigned long)ahi_arr[1]*blo_arr[0],(unsigned long)38*ahi_arr[1]*blo_arr[1],(unsigned long)19*ahi_arr[1]*blo_arr[2],(unsigned long)38*ahi_arr[1]*blo_arr[3]);
    check_vec_long(c10_v1, (unsigned long)19*ahi_arr[1]*blo_arr[4],38*(unsigned long)ahi_arr[1]*blo_arr[5],19*(unsigned long)ahi_arr[1]*blo_arr[6],(unsigned long)38*ahi_arr[1]*blo_arr[7]);
    check_vec_long(c10_v2, (unsigned long)19*ahi_arr[1]*bhi_arr[0],38*(unsigned long)ahi_arr[1]*bhi_arr[1],0,0);
    */

    __m256i offset0_v0, offset0_v1, offset0_v2, offset0_v3, offset0_v4; // column 1 + 5 + 9 + 13

    offset0_v0 = c1_v0;
    offset0_v1 = _mm256_add_epi64(c1_v1, c5_v0);
    offset0_v2 = _mm256_add_epi64(c1_v2, c5_v1);
    offset0_v2 = _mm256_add_epi64(offset0_v2, c9_v0);
    offset0_v3 = _mm256_add_epi64(c9_v1, c5_v2);
    offset0_v4 = c9_v2;
    /* --> OK
 	check_vec_long( offset0_v0, a0*b0, a0*b1, a0*b2, a0*b3 );
 	check_vec_long( offset0_v1, a0*b4 + a4*b0, a0*b5 + a4*b1, a0*b6 + a4*b2, a0*b7 + a4*b3 );
 	check_vec_long( offset0_v2, a0*b8 + a4*b4 + a8*b0, a0*b9 + a4*b5 + a8*b1, 19*a4*b6 + 19*a8*b2, 19*a4*b7 + 19*a8*b3 );
 	check_vec_long( offset0_v3, 19*a4*b8 + 19*a8*b4, 19*a4*b9 + 19*a8*b5, 19*a8*b6, 19*a8*b7 );
 	check_vec_long( offset0_v4, 19*a8*b8, 19*a8*b9, 0, 0 );
 	*/

    __m256i offset1_v0, offset1_v1, offset1_v2, offset1_v3, offset1_v4; // column 2 + 6 + 10 + 14

    offset1_v0 = c2_v0;
    offset1_v1 = _mm256_add_epi64(c2_v1, c6_v0);
    offset1_v2 = _mm256_add_epi64(c2_v2, c6_v1);
    offset1_v2 = _mm256_add_epi64(offset1_v2, c10_v0);
    offset1_v3 = _mm256_add_epi64(c10_v1, c6_v2);
    offset1_v4 = c10_v2;
    /* --> OK
 	check_vec_long( offset1_v0, a1*b0, 2*a1*b1, a1*b2, 2*a1*b3 );
 	check_vec_long( offset1_v1, a1*b4 + a5*b0, 2*a1*b5 + 2*a5*b1, a1*b6 + a5*b2, 2*a1*b7 + 2*a5*b3 );
 	check_vec_long( offset1_v2, a1*b8 + a5*b4 + a9*b0, 38*a1*b9 + 38*a5*b5 + 38*a9*b1, 19*a5*b6 + 19*a9*b2, 38*a5*b7 + 38*a9*b3 );
 	check_vec_long( offset1_v3, 19*a5*b8 + 19*a9*b4, 38*a5*b9 + 38*a9*b5, 19*a9*b6, 38*a9*b7 );
 	check_vec_long( offset1_v4, 19*a9*b8, 38*a9*b9, 0, 0 );
 	*/

    __m256i offset2_v0, offset2_v1, offset2_v2, offset2_v3; // column 3 + 7 + 11 + 15

    offset2_v0 = c3_v0;
    offset2_v1 = _mm256_add_epi64(c3_v1, c7_v0);
    offset2_v2 = _mm256_add_epi64(c3_v2, c7_v1);
    offset2_v3 = c7_v2;
    /* --> OK
 	check_vec_long( offset2_v0, a2*b0, a2*b1, a2*b2, a2*b3 );
 	check_vec_long( offset2_v1, a2*b4 + a6*b0, a2*b5 + a6*b1, a2*b6 + a6*b2, a2*b7 + a6*b3 );
 	check_vec_long( offset2_v2, 19*a2*b8 + 19*a6*b4, 19*a2*b9 + 19*a6*b5, 19*a6*b6, 19*a6*b7 );
 	check_vec_long( offset2_v3, 19*a6*b8, 19*a6*b9, 0, 0 );
 	*/

    __m256i offset3_v0, offset3_v1, offset3_v2, offset3_v3; // column 3 + 7 + 11 + 15

    offset3_v0 = c4_v0;
    offset3_v1 = _mm256_add_epi64(c4_v1, c8_v0);
    offset3_v2 = _mm256_add_epi64(c4_v2, c8_v1);
    offset3_v3 = c8_v2;
    /* --> OK
 	check_vec_long( offset3_v0, a3*b0, 2*a3*b1, a3*b2, 2*a3*b3 );
 	check_vec_long( offset3_v1, a3*b4 + a7*b0, 2*a3*b5 + 2*a7*b1, a3*b6 + a7*b2, 38*a3*b7 + 38*a7*b3 );
 	check_vec_long( offset3_v2, 19*a3*b8 + 19*a7*b4, 38*a3*b9 + 38*a7*b5, 19*a7*b6, 38*a7*b7 );
 	check_vec_long( offset3_v3, 19*a7*b8, 38*a7*b9, 0, 0 );
 	*/

    //  results :      v0      v1        v2         v3            v4
    //  results :   0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 19,

    // offset 0 :      v0      v1        v2         v3            v4
    // offset 0 :   0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 19,

    // offset 1 :      v0      v1        v2          v3            v4
    // offset 1 :   1 2 3 4, 5 6 7 8, 9 10 11 12, 13 14 15 16, 17 18 -19- -20-,
    // permuted :   4 1 2 3, 8 5 6 7, 12 9 10 11, 16 13 14 15, -20- 17 18 -19-,

    // offset 2 :      v0      v1        v2          v3
    // offset 2 :   2 3 4 5, 6 7 8 9, 10 11 12 13, 14 15 -16- -17-,
    // permuted :   4 5 2 3, 8 9 6 7, 12 13 10 11, -16- -17- 14 15,

    // offset 3 :      v0      v1        v2          v3
    // offset 3 :   3 4 5 6, 7 8 9 10, 11 12 13 14, 15 16 -17- -18-,
    // permuted :   4 5 6 3, 8 9 10 7, 12 13 14 11, 16 -17- -18- 15,

    // permute all offset vectors

    // offset 1
    __m256i offset1_v0_perm, offset1_v1_perm, offset1_v2_perm, offset1_v3_perm, offset1_v4_perm;
    offset1_v0_perm = _mm256_permute4x64_epi64(offset1_v0, 0b10010011);
    offset1_v1_perm = _mm256_permute4x64_epi64(offset1_v1, 0b10010011);
    offset1_v2_perm = _mm256_permute4x64_epi64(offset1_v2, 0b10010011);
    offset1_v3_perm = _mm256_permute4x64_epi64(offset1_v3, 0b10010011);
    offset1_v4_perm = _mm256_permute4x64_epi64(offset1_v4, 0b10010011);

    // offset 2
    __m256i offset2_v0_perm, offset2_v1_perm, offset2_v2_perm, offset2_v3_perm;
    offset2_v0_perm = _mm256_permute4x64_epi64(offset2_v0, 0b01001110);
    offset2_v1_perm = _mm256_permute4x64_epi64(offset2_v1, 0b01001110);
    offset2_v2_perm = _mm256_permute4x64_epi64(offset2_v2, 0b01001110);
    offset2_v3_perm = _mm256_permute4x64_epi64(offset2_v3, 0b01001110);

    // offset 3
    __m256i offset3_v0_perm, offset3_v1_perm, offset3_v2_perm, offset3_v3_perm;
    offset3_v0_perm = _mm256_permute4x64_epi64(offset3_v0, 0b00111001);
    offset3_v1_perm = _mm256_permute4x64_epi64(offset3_v1, 0b00111001);
    offset3_v2_perm = _mm256_permute4x64_epi64(offset3_v2, 0b00111001);
    offset3_v3_perm = _mm256_permute4x64_epi64(offset3_v3, 0b00111001);

    // add offset 0 and offset 1
    r_v0 = _mm256_add_epi64(offset0_v0, _mm256_blend_epi32(const_0_vec, offset1_v0_perm, 0b11111100));
    r_v1 = _mm256_add_epi64(offset0_v1, _mm256_blend_epi32(offset1_v0_perm, offset1_v1_perm, 0b11111100));
    r_v2 = _mm256_add_epi64(offset0_v2, _mm256_blend_epi32(offset1_v1_perm, offset1_v2_perm, 0b11111100));
    r_v3 = _mm256_add_epi64(offset0_v3, _mm256_blend_epi32(offset1_v2_perm, offset1_v3_perm, 0b11111100));
    r_v4 = _mm256_add_epi64(offset0_v4, _mm256_blend_epi32(offset1_v3_perm, offset1_v4_perm, 0b11111100));
    /* DEBUGGING HELP
	unsigned long tmp0  = a0*b0;
	unsigned long tmp1  = a0*b1 + a1*b0;
	unsigned long tmp2  = a0*b2 + 2*a1*b1;
	unsigned long tmp3  = a0*b3 + a1*b2;
	unsigned long tmp4  = a0*b4 + 2*a1*b3 + a4*b0;
	unsigned long tmp5  = a0*b5 + a1*b4 + a4*b1 + a5*b0;
	unsigned long tmp6  = a0*b6 + 2*a1*b5 + a4*b2 + 2*a5*b1;
	unsigned long tmp7  = a0*b7 + a1*b6 + a4*b3 + a5*b2;
	unsigned long tmp8  = a0*b8 + 2*a1*b7 + a4*b4 + 2*a5*b3 + a8*b0;
	unsigned long tmp9  = a0*b9 + a1*b8 + a4*b5 + a5*b4 + a8*b1 + a9*b0;
	unsigned long tmp10 = 38*a1*b9 + 19*a4*b6 + 38*a5*b5 + 19*a8*b2 + 38*a9*b1;
	unsigned long tmp11 = 19*a4*b7 + 19*a5*b6 + 19*a8*b3 + 19*a9*b2;
	unsigned long tmp12 = 19*a4*b8 + 38*a5*b7 + 19*a8*b4 + 38*a9*b3;
	unsigned long tmp13 = 19*a4*b9 + 19*a5*b8 + 19*a8*b5 + 19*a9*b4;
	unsigned long tmp14 = 38*a5*b9 + 19*a8*b6 + 38*a9*b5;
	unsigned long tmp15 = 19*a8*b7 + 19*a9*b6;
	unsigned long tmp16 = 19*a8*b8 + 38*a9*b7;
	unsigned long tmp17 = 19*a8*b9 + 19*a9*b8;
	unsigned long tmp18 = 38*a9*b9;
	check_vec_long( r_v0, tmp0, tmp1, tmp2, tmp3 );
	check_vec_long( r_v1, tmp4, tmp5, tmp6, tmp7 );
	check_vec_long( r_v2, tmp8, tmp9, tmp10, tmp11 );
	check_vec_long( r_v3, tmp12, tmp13, tmp14, tmp15 );
	check_vec_long( r_v4, tmp16, tmp17, tmp18, 0 );
	*/

    // add offset 2 to result
    r_v0 = _mm256_add_epi64(r_v0, _mm256_blend_epi32(const_0_vec, offset2_v0_perm, 0b11110000));
    r_v1 = _mm256_add_epi64(r_v1, _mm256_blend_epi32(offset2_v0_perm, offset2_v1_perm, 0b11110000));
    r_v2 = _mm256_add_epi64(r_v2, _mm256_blend_epi32(offset2_v1_perm, offset2_v2_perm, 0b11110000));
    r_v3 = _mm256_add_epi64(r_v3, _mm256_blend_epi32(offset2_v2_perm, offset2_v3_perm, 0b11110000));
    /* DEBUGGING HELP
	tmp2  += a2*b0;
	tmp3  += a2*b1;
	tmp4  += a2*b2;
	tmp5  += a2*b3;
	tmp6  += a2*b4 + a6*b0;
	tmp7  += a2*b5 + a6*b1;
	tmp8  += a2*b6 + a6*b2;
	tmp9  += a2*b7 + a6*b3;
	tmp10 += 19*a2*b8 + 19*a6*b4;
	tmp11 += 19*a2*b9 + 19*a6*b5;
	tmp12 += 19*a6*b6;
	tmp13 += 19*a6*b7;
	tmp14 += 19*a6*b8;
	tmp15 += 19*a6*b9;
	--> OK
	check_vec_long( r_v0, tmp0, tmp1, tmp2, tmp3 );
	check_vec_long( r_v1, tmp4, tmp5, tmp6, tmp7 );
	check_vec_long( r_v2, tmp8, tmp9, tmp10, tmp11 );
	check_vec_long( r_v3, tmp12, tmp13, tmp14, tmp15 );
	check_vec_long( r_v4, tmp16, tmp17, tmp18, 0 );
	*/

    // add offset 3 to result
    r_v0 = _mm256_add_epi64(r_v0, _mm256_blend_epi32(const_0_vec, offset3_v0_perm, 0b11000000));
    r_v1 = _mm256_add_epi64(r_v1, _mm256_blend_epi32(offset3_v0_perm, offset3_v1_perm, 0b11000000));
    r_v2 = _mm256_add_epi64(r_v2, _mm256_blend_epi32(offset3_v1_perm, offset3_v2_perm, 0b11000000));
    r_v3 = _mm256_add_epi64(r_v3, _mm256_blend_epi32(offset3_v2_perm, offset3_v3_perm, 0b11000000));
    r_v4 = _mm256_add_epi64(r_v4, _mm256_blend_epi32(offset3_v3_perm, const_0_vec, 0b11111100));
    /* DEBUGGING HELP
	tmp3  += a3*b0;
	tmp4  += 2*a3*b1;
	tmp5  += a3*b2;
	tmp6  += 2*a3*b3;
	tmp7  += a3*b4 + a7*b0;
	tmp8  += 2*a3*b5 + 2*a7*b1;
	tmp9  += a3*b6 + a7*b2;
	tmp10 += 38*a3*b7 + 38*a7*b3;
	tmp11 += 19*a3*b8 + 19*a7*b4;
	tmp12 += 38*a3*b9 + 38*a7*b5;
	tmp13 += 19*a7*b6;
	tmp14 += 38*a7*b7;
	tmp15 += 19*a7*b8;
	tmp16 += 38*a7*b9;
	check_vec_long( r_v0, tmp0, tmp1, tmp2, tmp3 );
	check_vec_long( r_v1, tmp4, tmp5, tmp6, tmp7 );
	check_vec_long( r_v2, tmp8, tmp9, tmp10, tmp11 );
	check_vec_long( r_v3, tmp12, tmp13, tmp14, tmp15 );
	check_vec_long( r_v4, tmp16, tmp17, tmp18, 0 );

	print_vec_long( r_v0, "r_v0" );
	print_vec_long( r_v1, "r_v1" );
	print_vec_long( r_v2, "r_v2" );
	print_vec_long( r_v3, "r_v3" );
	print_vec_long( r_v4, "r_v4" );
	*/

    // v0 : 0,1,2,3 / v1: 4,5,6,7 / v2: 8,9,10,11 / v3: 12, 13, 14, 15 / v4: 16,17,18,19

    // v0 : 0,1,2,3 / v1: 4,5,6,7 / v2: 8, 9 - - /
    // ++++++++++++++++++++++++++++++++++++++++++
    //   10,11,12,13 /  14,15,16,17 / 18
    // r_v0 = _mm256_setr_epi64x(0,1,2,3);
    // r_v1 = _mm256_setr_epi64x(4,5,6,7);
    // r_v2 = _mm256_setr_epi64x(8,9,10,11);
    // r_v3 = _mm256_setr_epi64x(12,13,14,15);
    // r_v4 = _mm256_setr_epi64x(16,17,18,19);

    __m256i tmp_10_13, tmp_14_17, tmp_18_;

    tmp_10_13 = _mm256_blend_epi32(r_v2, r_v3, 0b00001111); // 12 13 10 11
    tmp_10_13 = _mm256_permute4x64_epi64(tmp_10_13, 0b01001110);

    tmp_14_17 = _mm256_blend_epi32(r_v3, r_v4, 0b00001111); // 16 17 14 15
    tmp_14_17 = _mm256_permute4x64_epi64(tmp_14_17, 0b01001110);

    tmp_18_ = _mm256_permute4x64_epi64(r_v4, 0b10101010);
    tmp_18_ = _mm256_blend_epi32(tmp_18_, const_0_vec, 0b11111100);

    r_v0 = _mm256_add_epi64(r_v0, tmp_10_13);
    r_v1 = _mm256_add_epi64(r_v1, tmp_14_17);
    r_v2 = _mm256_add_epi64(r_v2, tmp_18_);
    r_v2 = _mm256_blend_epi32(r_v2, const_0_vec, 0b11110000);

    /*
	print_vec_long( r_v0, "r_v0" );
	print_vec_long( r_v1, "r_v1" );
	print_vec_long( r_v2, "r_v2" );
	*/

    reduce_3_to_2(&r_v0, &r_v1, &r_v2);
    reduce_2_to_2(&r_v0, &r_v1);

    copy_256(reslo, reshi, r_v0, r_v1);
}

void square_256_vec_to_non_vec(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi)
{

    signed int reslo_arr[8], reshi_arr[8];

    // read input vectors into arrays
    signed int *xlo_arr = (int *)&xlo;
    signed int *xhi_arr = (int *)&xhi;

    // square arrays
    square_256_non_vec(reslo_arr, reshi_arr, xlo_arr, xhi_arr);

    // store results back into vectors
    *reslo = _mm256_loadu_si256((const __m256i *)reslo_arr);
    *reshi = _mm256_loadu_si256((const __m256i *)reshi_arr);
}

void square_256(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi)
{
    
    // FOR DEBUGGING
    signed int* alo_arr = (int*) &xlo;
    signed int* ahi_arr = (int*) &xhi;

	unsigned long a0 = (unsigned long) alo_arr[0];
	unsigned long a1 = (unsigned long) alo_arr[1];
	unsigned long a2 = (unsigned long) alo_arr[2];
	unsigned long a3 = (unsigned long) alo_arr[3];
	unsigned long a4 = (unsigned long) alo_arr[4];
	unsigned long a5 = (unsigned long) alo_arr[5];
	unsigned long a6 = (unsigned long) alo_arr[6];
	unsigned long a7 = (unsigned long) alo_arr[7];
	unsigned long a8 = (unsigned long) ahi_arr[0];
	unsigned long a9 = (unsigned long) ahi_arr[1];
	// END DEBUGGING VARIABLES ASSIGNMENTS (TODO: DELETE ABOVE)
    

    // xlo = [0,1,2,3,4,5,6,7] xhi = [8,9,10,11,12,13,14,15]
    __m256i a_0_3 = _mm256_permutevar8x32_epi32(xlo, index1_vec);
    //check_vec_int(a_0_3, a0, a0, a1, a1, a2, a2, a3, a3); // OK
    __m256i a_4_7 = _mm256_permutevar8x32_epi32(xlo, index2_vec);
    //check_vec_int(a_4_7, a4, a4, a5, a5, a6, a6, a7, a7); // OK
    __m256i a_8_11 = _mm256_permutevar8x32_epi32(xhi, index1_vec);
    //check_vec_int(a_8_11, a8, a8, a9, a9, 0, 0, 0, 0); // OK

    __m256i a2_0_7 = _mm256_slli_epi32(xlo, 1);
    //check_vec_int(a2_0_7, 2*a0, 2*a1, 2*a2, 2*a3, 2*a4, 2*a5, 2*a6, 2*a7); // OK
    __m256i a2_8_14 = _mm256_slli_epi32(xhi, 1);
    //check_vec_int(a2_8_14, 2*a8, 2*a9, 0, 0, 0, 0, 0, 0); // OK

    __m256i a2_0 = _mm256_broadcastq_epi64(_mm256_castsi256_si128(a2_0_7));
    //check_vec_int(a2_0, 2*a0, 2*a1, 2*a0, 2*a1, 2*a0, 2*a1, 2*a0, 2*a1); // OK
    __m256i a2_1 = _mm256_permutevar8x32_epi32(a2_0, index_alt_vec);
    //check_vec_int(a2_1, 2*a1, 2*a0, 2*a1, 2*a0, 2*a1, 2*a0, 2*a1, 2*a0); // OK
    __m256i a2_2 = _mm256_broadcastq_epi64(_mm256_castsi256_si128(_mm256_unpackhi_epi64(a2_0_7, index1_vec)));
    //check_vec_int(a2_2, 2*a2, 2*a3, 2*a2, 2*a3, 2*a2, 2*a3, 2*a2, 2*a3); // OK
    __m256i a2_3 = _mm256_permutevar8x32_epi32(a2_2, index_alt_vec);
    //check_vec_int(a2_3, 2*a3, 2*a2, 2*a3, 2*a2, 2*a3, 2*a2, 2*a3, 2*a2); // OK

    __m128i a_lo_hi = _mm256_extracti128_si256(a2_0_7, 1);
    __m256i a2_4 = _mm256_broadcastq_epi64(a_lo_hi);
    //check_vec_int(a2_4, 2*a4, 2*a5, 2*a4, 2*a5, 2*a4, 2*a5, 2*a4, 2*a5); // OK
    __m256i a2_5 = _mm256_permutevar8x32_epi32(a2_4, index_alt_vec);
    //check_vec_int(a2_5, 2*a5, 2*a4, 2*a5, 2*a4, 2*a5, 2*a4, 2*a5, 2*a4); // OK
    __m256i a2_6 = _mm256_broadcastq_epi64(_mm256_castsi256_si128(_mm256_unpackhi_epi64(_mm256_castsi128_si256(a_lo_hi), index1_vec)));
    //check_vec_int(a2_6, 2*a6, 2*a7, 2*a6, 2*a7, 2*a6, 2*a7, 2*a6, 2*a7); // OK
    __m256i a2_7 = _mm256_permutevar8x32_epi32(a2_6, index_alt_vec);
    //check_vec_int(a2_7, 2*a7, 2*a6, 2*a7, 2*a6, 2*a7, 2*a6, 2*a7, 2*a6); // OK
    __m256i a2_8 = _mm256_broadcastq_epi64(_mm256_castsi256_si128(a2_8_14));
    //check_vec_int(a2_8, 2*a8, 2*a9, 2*a8, 2*a9, 2*a8, 2*a9, 2*a8, 2*a9); // OK

    __m256i r_v0, r_v1, r_v2, r_v3, r_v4, r_v5, r_v6, r_v7; // real result vectors

    __m256i c1_v0, c1_v1, c1_v2; //column 1 results
    __m256i c2_v0, c2_v1, c2_v2;
    __m256i c3_v0, c3_v1, c3_v2;
    __m256i c4_v0, c4_v1, c4_v2;
    __m256i c5_v1, c5_v2;
    __m256i c6_v1, c6_v2;
    __m256i c7_v1, c7_v2;
    __m256i c8_v1, c8_v2;
    __m256i c9_v2;
    __m256i c10_v2;

    /*
                c1 of 0     c2 of1        c3 of 2         c4 of 3       c5 of 0          c6             c7              c8                  c9           c10
    r[0]  =     a0*a0   
    r[1]  =   _2a0*a1          0           
    r[2]  =   _2a0*a2 +     _2a1*a1             0        
    r[3]  =   _2a0*a3 +     _2a1*a2             1              0
    r[4]  =   _2a0*a4 + 2 * _2a1*a3 +         a2*a2            1               0
    r[5]  =   _2a0*a5 +     _2a1*a4 +       _2a2*a3            2               1              0
    r[6]  =   _2a0*a6 + 2 * _2a1*a5 +       _2a2*a4   +     _2a3*a3            2              1               0
    r[7]  =   _2a0*a7 +     _2a1*a6 +       _2a2*a5   +     _2a3*a4            3              2               1            0
    r[8]  =   _2a0*a8 + 2 * _2a1*a7 +       _2a2*a6   + 2 * _2a3*a5  +       a4*a4            3               2            1                 0
    r[9]  =   _2a0*a9 +     _2a1*a8 +       _2a2*a7   +     _2a3*a6  +     _2a4*a5            4               3            2                 1         0
    r[10] =             38* _2a1*a9 +  19 * _2a2*a8   + 38* _2a3*a7  + 19* _2a4*a6  + 38*   a5*a5             4            3                 2         1
    r[11] =                            19 * _2a2*a9   + 19* _2a3*a8  + 19* _2a4*a7  + 19* _2a5*a6             5            4                 3         2
    r[12] =                                             38* _2a3*a9  + 19* _2a4*a8  + 38* _2a5*a7  +  19*   a6*a6          5                 4         3
    r[13] =                                                            19* _2a4*a9  + 19* _2a5*a8  +  19* _2a6*a7          6                 5         4
    r[14] =                                                                           38* _2a5*a9  +  19* _2a6*a8  + 38*   a7*a7             6         5
    r[15] =                                                                                           19* _2a6*a9  + 19* _2a7*a8             7         6
    r[16] =                                                                                                          38* _2a7*a9  + 19*   a8*a8        7
    r[17] =                                                                                                                         19* _2a8*a9        8
    r[18] =                                                                                                                                       38*a9*a9;

	*/

    c1_v0 = _mm256_mul_epi32(_mm256_blend_epi32(a2_0, a_0_3, 0b00000011), a_0_3);
    c1_v1 = _mm256_mul_epi32(a2_0, a_4_7);
    c1_v2 = _mm256_mul_epi32(a2_0, a_8_11);
    c1_v2 = _mm256_blend_epi32(c1_v2, const_0_vec, 0b11110000);
    //check_vec_long(c1_v0, a0*a0, 2*a0*a1, 2*a0*a2, 2*a0*a3); // OK
    //check_vec_long(c1_v1, 2*a0*a4, 2*a0*a5, 2*a0*a6, 2*a0*a7); // OK
    //check_vec_long(c1_v2, 2*a0*a8, 2*a0*a9, 0, 0); // OK

    c2_v0 = _mm256_mul_epi32(a_0_3, c_0_1_1_2__64_vec);
    c2_v0 = _mm256_mul_epi32(a2_1, c2_v0);
    c2_v1 = _mm256_sllv_epi64(a_4_7, c_alt_0_1__64_vec);
    c2_v1 = _mm256_mul_epi32(a2_1, c2_v1);
    c2_v2 = _mm256_mul_epi32(a_8_11, c_1_38_2x0__64_vec);
    c2_v2 = _mm256_mul_epi32(a2_1, c2_v2);
    c2_v2 = _mm256_blend_epi32(c2_v2, const_0_vec, 0b11110000);
    //check_vec_long(c2_v0, 0, 2*a1*a1, 2*a1*a2, 2*2*a1*a3); // OK
    //check_vec_long(c2_v1, 2*a1*a4, 2*2*a1*a5, 2*a1*a6, 2*2*a1*a7); // OK
    //check_vec_long(c2_v2, 2*a1*a8, 38*2*a1*a9, 0, 0); // OK

    c3_v0 = _mm256_mul_epi32(_mm256_blend_epi32(a_0_3, const_0_vec, 0b00001111), _mm256_blend_epi32(a_0_3, a2_2, 0b11000000));
    c3_v1 = _mm256_mul_epi32(a2_2, a_4_7);
    c3_v2 = _mm256_mul_epi32(a_8_11, c_2x19_2x0__64_vec);
    c3_v2 = _mm256_mul_epi32(a2_2, c3_v2);
    //check_vec_long(c3_v0, 0, 0, a2*a2, 2*a2*a3); // OK
    //check_vec_long(c3_v1, 2*a2*a4, 2*a2*a5, 2*a2*a6, 2*a2*a7); // OK
    //check_vec_long(c3_v2, 19*2*a2*a8, 19*2*a2*a9, 0, 0); // OK

    c4_v0 = _mm256_mul_epi32(a2_3, _mm256_blend_epi32(a_0_3, const_0_vec, 0b00111111));
    c4_v1 = _mm256_mul_epi32(a_4_7, c_1_2_1_38__64_vec);
    c4_v1 = _mm256_mul_epi32(a2_3, c4_v1);
    c4_v2 = _mm256_mul_epi32(a_8_11, c_19_38_2x0__64_vec);
    c4_v2 = _mm256_mul_epi32(a2_3, c4_v2);
    //check_vec_long(c4_v0, 0, 0, 0, 2*a3*a3); // OK
    //check_vec_long(c4_v1, 2*a3*a4, 2*2*a3*a5, 2*a3*a6, 38*2*a3*a7); // OK
    //check_vec_long(c4_v2, 19*2*a3*a8, 38*2*a3*a9, 0, 0); // OK

    c5_v1 = _mm256_mul_epi32(a_4_7, c_2x1_2x19__64_vec);
    c5_v1 = _mm256_mul_epi32(_mm256_blend_epi32(a2_4, a_4_7, 0b00000011), c5_v1);
    c5_v2 = _mm256_mul_epi32(a_8_11, c_2x19_2x0__64_vec);
    c5_v2 = _mm256_mul_epi32(a2_4, c5_v2);
    //check_vec_long(c5_v1, a4*a4, 2*a4*a5, 19*2*a4*a6, 19*2*a4*a7); // OK
    //check_vec_long(c5_v2, 19*2*a4*a8, 19*2*a4*a9, 0, 0); // OK

    c6_v1 = _mm256_mul_epi32(a_4_7, c_0_38_19_38__64_vec);
    c6_v1 = _mm256_mul_epi32(_mm256_blend_epi32(a2_5, a_4_7, 0b00001111), c6_v1);
    c6_v2 = _mm256_mul_epi32(a_8_11, c_19_38_2x0__64_vec);
    c6_v2 = _mm256_mul_epi32(a2_5, c6_v2);
    //check_vec_long(c6_v1, 0, 38*a5*a5, 19*2*a5*a6, 38*2*a5*a7); // OK
    //check_vec_long(c6_v2, 19*2*a5*a8, 38*2*a5*a9, 0, 0); // OK

    c7_v1 = _mm256_mul_epi32(a_4_7, c_2x0_2x19__64_vec);
    c7_v1 = _mm256_mul_epi32(_mm256_blend_epi32(a2_6, a_4_7, 0b00111111), c7_v1);
    c7_v2 = _mm256_mul_epi32(a_8_11, c_2x19_2x0__64_vec);
    c7_v2 = _mm256_mul_epi32(a2_6, c7_v2);
    //check_vec_long(c7_v1, 0, 0, 19*a6*a6, 19*2*a6*a7); // OK
    //check_vec_long(c7_v2, 19*2*a6*a8, 19*2*a6*a9, 0, 0); // OK

    c8_v1 = _mm256_mul_epi32(a_4_7, c_3x0_38__64_vec);
    c8_v1 = _mm256_mul_epi32(a_4_7, c8_v1);
    c8_v2 = _mm256_mul_epi32(a_8_11, c_19_38_2x0__64_vec);
    c8_v2 = _mm256_mul_epi32(a2_7, c8_v2);
    //check_vec_long(c8_v1, 0, 0, 0, 38*a7*a7); // OK
    //check_vec_long(c8_v2, 19*2*a7*a8, 38*2*a7*a9, 0, 0); // OK

    c9_v2 = _mm256_mul_epi32(a_8_11, c_2x19_2x0__64_vec);
    c9_v2 = _mm256_mul_epi32(_mm256_blend_epi32(a2_8, a_8_11, 0b00000011), c9_v2);
    //check_vec_long(c9_v2, 19*a8*a8, 19*2*a8*a9, 0, 0); // OK

    c10_v2 = _mm256_mul_epi32(a_8_11, c_0_38_2x0__64_vec);
    c10_v2 = _mm256_mul_epi32(a_8_11, c10_v2);
    //check_vec_long(c10_v2, 0, 38*a9*a9, 0, 0); // OK

    // offset 0 means:  0 1 2 3, 4 5 6 7, 8 9 10 11, ...
    __m256i offset0_v0, offset0_v1, offset0_v2, offset0_v3, offset0_v4; // column 1 + 5 + 9
    offset0_v0 = c1_v0;
    offset0_v1 = c1_v1;
    offset0_v2 = _mm256_add_epi64(c1_v2, c5_v1);
    offset0_v3 = c5_v2;
    offset0_v4 = c9_v2;

    __m256i offset1_v0, offset1_v1, offset1_v2, offset1_v3, offset1_v4; // column 2 + 6 + 10

    offset1_v0 = c2_v0;
    offset1_v1 = c2_v1;
    offset1_v2 = _mm256_add_epi64(c2_v2, c6_v1);
    offset1_v3 = c6_v2;
    offset1_v4 = c10_v2;

    __m256i offset2_v0, offset2_v1, offset2_v2, offset2_v3; // column 3 + 7

    offset2_v0 = c3_v0;
    offset2_v1 = c3_v1;
    offset2_v2 = _mm256_add_epi64(c3_v2, c7_v1);
    offset2_v3 = c7_v2;

    __m256i offset3_v0, offset3_v1, offset3_v2, offset3_v3; // column 4 + 8
    offset3_v0 = c4_v0;
    offset3_v1 = c4_v1;
    offset3_v2 = _mm256_add_epi64(c4_v2, c8_v1);
    offset3_v3 = c8_v2;

    //  results :      v0      v1        v2         v3            v4
    //  results :   0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 -19-,

    // offset 0 :      v0      v1        v2         v3            v4
    // offset 0 :   0 1 2 3, 4 5 6 7, 8 9 10 11, 12 13 14 15, 16 17 18 -19-,

    // offset 1 :      v0      v1        v2          v3            v4
    // offset 1 :   1 2 3 4, 5 6 7 8, 9 10 11 12, 13 14 15 16, 17 18 -19- -20-
    // permuted :   4 1 2 3, 8 5 6 7, 12 9 10 11, 16 13 14 15, -20- 17 18 -19-

    // offset 2 :      v0      v1        v2          v3
    // offset 2 :   2 3 4 5, 6 7 8 9, 10 11 12 13, 14 15 16 17,
    // permuted :   4 5 2 3, 8 9 6 7, 12 13 10 11, 16 17 14 15,

    // offset 3 :      v0      v1        v2          v3
    // offset 3 :   3 4 5 6, 7 8 9 10, 11 12 13 14, 15 16 17 18
    // permuted :   4 5 6 3, 8 9 10 7, 12 13 14 11, 16 17 18 15

    // permute all offset vectors

    // offset 1
    __m256i offset1_v0_perm, offset1_v1_perm, offset1_v2_perm, offset1_v3_perm, offset1_v4_perm;
    offset1_v0_perm = _mm256_permute4x64_epi64(offset1_v0, 0b10010011);
    offset1_v1_perm = _mm256_permute4x64_epi64(offset1_v1, 0b10010011);
    offset1_v2_perm = _mm256_permute4x64_epi64(offset1_v2, 0b10010011);
    offset1_v3_perm = _mm256_permute4x64_epi64(offset1_v3, 0b10010011);
    offset1_v4_perm = _mm256_permute4x64_epi64(offset1_v4, 0b10010011);

    // offset 2
    __m256i offset2_v0_perm, offset2_v1_perm, offset2_v2_perm, offset2_v3_perm;
    offset2_v0_perm = _mm256_permute4x64_epi64(offset2_v0, 0b01001110);
    offset2_v1_perm = _mm256_permute4x64_epi64(offset2_v1, 0b01001110);
    offset2_v2_perm = _mm256_permute4x64_epi64(offset2_v2, 0b01001110);
    offset2_v3_perm = _mm256_permute4x64_epi64(offset2_v3, 0b01001110);

    // offset 3
    __m256i offset3_v0_perm, offset3_v1_perm, offset3_v2_perm, offset3_v3_perm;
    offset3_v0_perm = _mm256_permute4x64_epi64(offset3_v0, 0b00111001);
    offset3_v1_perm = _mm256_permute4x64_epi64(offset3_v1, 0b00111001);
    offset3_v2_perm = _mm256_permute4x64_epi64(offset3_v2, 0b00111001);
    offset3_v3_perm = _mm256_permute4x64_epi64(offset3_v3, 0b00111001);

    //check_vec_long(offset0_v2, 2*a0*a8+a4*a4, 2*a0*a9+2*a4*a5, 19*2*a4*a6, 19*2*a4*a7);

    // add offset 0 and offset 1
    r_v0 = _mm256_add_epi64(offset0_v0, _mm256_blend_epi32(const_0_vec, offset1_v0_perm, 0b11111100));
    r_v1 = _mm256_add_epi64(offset0_v1, _mm256_blend_epi32(offset1_v0_perm, offset1_v1_perm, 0b11111100));
    r_v2 = _mm256_add_epi64(offset0_v2, _mm256_blend_epi32(offset1_v1_perm, offset1_v2_perm, 0b11111100));
    r_v3 = _mm256_add_epi64(offset0_v3, _mm256_blend_epi32(offset1_v2_perm, offset1_v3_perm, 0b11111100));
    r_v4 = _mm256_add_epi64(offset0_v4, _mm256_blend_epi32(offset1_v3_perm, offset1_v4_perm, 0b11111100));

    ////check_vec_long(r_v2, 2*a0*a8+2*2*a1*a7+a4*a4, 2*a0*a9+2*a1*a8+2*a4*a5, 38*2*a1*a9+19*2*a4*a6+38*a5*a5, 19*2*a4*a7+19*2*a5*a6);

    // add offset 2 to result
    r_v0 = _mm256_add_epi64(r_v0, _mm256_blend_epi32(const_0_vec, offset2_v0_perm, 0b11110000));
    r_v1 = _mm256_add_epi64(r_v1, _mm256_blend_epi32(offset2_v0_perm, offset2_v1_perm, 0b11110000));
    r_v2 = _mm256_add_epi64(r_v2, _mm256_blend_epi32(offset2_v1_perm, offset2_v2_perm, 0b11110000));
    r_v3 = _mm256_add_epi64(r_v3, _mm256_blend_epi32(offset2_v2_perm, offset2_v3_perm, 0b11110000));
    r_v4 = _mm256_add_epi64(r_v4, _mm256_blend_epi32(offset2_v3_perm, const_0_vec, 0b11110000));

    //check_vec_long(r_v2, 2*a0*a8+2*2*a1*a7+2*a2*a6+a4*a4, 2*a0*a9+2*a1*a8+2*a2*a7+2*a4*a5, 38*2*a1*a9+19*2*a2*a8+19*2*a4*a6+38*a5*a5, 19*2*a2*a9+19*2*a4*a7+19*2*a5*a6);

    // add offset 3 to result
    r_v0 = _mm256_add_epi64(r_v0, _mm256_blend_epi32(const_0_vec, offset3_v0_perm, 0b11000000));
    r_v1 = _mm256_add_epi64(r_v1, _mm256_blend_epi32(offset3_v0_perm, offset3_v1_perm, 0b11000000));
    r_v2 = _mm256_add_epi64(r_v2, _mm256_blend_epi32(offset3_v1_perm, offset3_v2_perm, 0b11000000));
    r_v3 = _mm256_add_epi64(r_v3, _mm256_blend_epi32(offset3_v2_perm, offset3_v3_perm, 0b11000000));
    r_v4 = _mm256_add_epi64(r_v4, _mm256_blend_epi32(offset3_v3_perm, const_0_vec, 0b11000000));

    //check_vec_long(r_v0, a0*a0, 2*a0*a1, 2*a0*a2+2*a1*a1, 2*a0*a3+2*a1*a2);
    // Ok
    //check_vec_long(r_v1, 2*a0*a4+2*2*a1*a3+a2*a2, 2*a0*a5+2*a1*a4+2*a2*a3, 2*a0*a6+2*2*a1*a5+2*a2*a4+2*a3*a3, 2*a0*a7+2*a1*a6+2*a2*a5+2*a3*a4);
    // Ok
    //check_vec_long(r_v2, 2*a0*a8+2*2*a1*a7+2*a2*a6+2*2*a3*a5+a4*a4, 2*a0*a9+2*a1*a8+2*a2*a7+2*a3*a6+2*a4*a5, 38*2*a1*a9+19*2*a2*a8+38*2*a3*a7+19*2*a4*a6+38*a5*a5, 19*2*a2*a9+19*2*a3*a8+19*2*a4*a7+19*2*a5*a6);
    // Ok
    //check_vec_long(r_v3, 38*2*a3*a9+19*2*a4*a8+38*2*a5*a7+19*a6*a6, 19*2*a4*a9+19*2*a5*a8+19*2*a6*a7, 38*2*a5*a9+19*2*a6*a8+38*a7*a7, 19*2*a6*a9+19*2*a7*a8);
    // Ok
    //check_vec_long(r_v4, 38*2*a7*a9+19*a8*a8, 19*2*a8*a9, 38*a9*a9,0);
    // Ok

    __m256i tmp_10_13, tmp_14_17, tmp_18_;

    tmp_10_13 = _mm256_blend_epi32(r_v2, r_v3, 0b00001111); // 12 13 10 11
    tmp_10_13 = _mm256_permute4x64_epi64(tmp_10_13, 0b01001110);

    tmp_14_17 = _mm256_blend_epi32(r_v3, r_v4, 0b00001111); // 16 17 14 15
    tmp_14_17 = _mm256_permute4x64_epi64(tmp_14_17, 0b01001110);

    tmp_18_ = _mm256_permute4x64_epi64(r_v4, 0b10101010);
    tmp_18_ = _mm256_blend_epi32(tmp_18_, const_0_vec, 0b11111100);

    r_v0 = _mm256_add_epi64(r_v0, tmp_10_13);
    r_v1 = _mm256_add_epi64(r_v1, tmp_14_17);
    r_v2 = _mm256_add_epi64(r_v2, tmp_18_);
    r_v2 = _mm256_blend_epi32(r_v2, const_0_vec, 0b11110000);

    reduce_3_to_2(&r_v0, &r_v1, &r_v2);
    reduce_2_to_2(&r_v0, &r_v1);

    copy_256(reslo, reshi, r_v0, r_v1);
}

void square_256_non_vec(signed int reslo[8], signed int reshi[8], signed int xlo[8], signed int xhi[8])
{

    unsigned long r[10];

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
    signed long _2a0 = a0 << 1;
    signed long _2a1 = a1 << 1;
    signed long _2a2 = a2 << 1;
    signed long _2a3 = a3 << 1;
    signed long _2a4 = a4 << 1;
    signed long _2a5 = a5 << 1;
    signed long _2a6 = a6 << 1;
    signed long _2a7 = a7 << 1;
    signed long _2a8 = a8 << 1;

    r[0] = a0 * a0 + 38 * _2a1 * a9 + 19 * _2a2 * a8 + 38 * _2a3 * a7 + 19 * _2a4 * a6 + 38 * a5 * a5;
    r[1] = _2a0 * a1 + 19 * _2a2 * a9 + 19 * _2a3 * a8 + 19 * _2a4 * a7 + 19 * _2a5 * a6;
    r[2] = _2a0 * a2 + _2a1 * a1 + 38 * _2a3 * a9 + 19 * _2a4 * a8 + 38 * _2a5 * a7 + 19 * a6 * a6;
    r[3] = _2a0 * a3 + _2a1 * a2 + 19 * _2a4 * a9 + 19 * _2a5 * a8 + 19 * _2a6 * a7;
    r[4] = _2a0 * a4 + 2 * _2a1 * a3 + a2 * a2 + 38 * _2a5 * a9 + 19 * _2a6 * a8 + 38 * a7 * a7;
    r[5] = _2a0 * a5 + _2a1 * a4 + _2a2 * a3 + 19 * _2a6 * a9 + 19 * _2a7 * a8;
    r[6] = _2a0 * a6 + 2 * _2a1 * a5 + _2a2 * a4 + _2a3 * a3 + 38 * _2a7 * a9 + 19 * a8 * a8;
    r[7] = _2a0 * a7 + _2a1 * a6 + _2a2 * a5 + _2a3 * a4 + 19 * _2a8 * a9;
    r[8] = _2a0 * a8 + 2 * _2a1 * a7 + _2a2 * a6 + 2 * _2a3 * a5 + a4 * a4 + 38 * a9 * a9;
    r[9] = _2a0 * a9 + _2a1 * a8 + _2a2 * a7 + _2a3 * a6 + _2a4 * a5;

    reduce_3_to_2_non_vec(r);
    reduce_3_to_2_non_vec(r);

    for (int i = 0; i < 8; i++)
    {
        reslo[i] = r[i];
        reshi[i] = r[i + 8];
    }
    reshi[2] = 0;
    reshi[3] = 0;
    reshi[4] = 0;
    reshi[5] = 0;
    reshi[6] = 0;
    reshi[7] = 0;

    reduce_2_to_2_non_vec(reslo, reshi);
    //print_radix17( res, "square" );
}
void reduce_2_to_2_non_vec(signed int xlo[8], signed int xhi[8])
{

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

    carry = d0 >> (LIMB_SIZE + 1);
    d1 += carry;
    carry <<= (LIMB_SIZE + 1);
    d0 -= carry;

    carry = d1 >> LIMB_SIZE;
    d2 += carry;
    carry <<= LIMB_SIZE;
    d1 -= carry;

    carry = d2 >> (LIMB_SIZE + 1);
    d3 += carry;
    carry <<= (LIMB_SIZE + 1);
    d2 -= carry;

    carry = d3 >> LIMB_SIZE;
    d4 += carry;
    carry <<= LIMB_SIZE;
    d3 -= carry;

    carry = d4 >> (LIMB_SIZE + 1);
    d5 += carry;
    carry <<= (LIMB_SIZE + 1);
    d4 -= carry;

    carry = d5 >> LIMB_SIZE;
    d6 += carry;
    carry <<= LIMB_SIZE;
    d5 -= carry;

    carry = d6 >> (LIMB_SIZE + 1);
    d7 += carry;
    carry <<= (LIMB_SIZE + 1);
    d6 -= carry;

    carry = d7 >> LIMB_SIZE;
    d8 += carry;
    carry <<= LIMB_SIZE;
    d7 -= carry;

    carry = d8 >> (LIMB_SIZE + 1);
    d9 += carry;
    carry <<= (LIMB_SIZE + 1);
    d8 -= carry;

    carry = d9 >> LIMB_SIZE;
    d0 += 19 * carry;
    carry <<= LIMB_SIZE;
    d9 -= carry;

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
}
void reduce_3_to_2_non_vec(unsigned long r[10])
{

    signed long carry;
    unsigned long r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;

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

    carry = r0 >> (LIMB_SIZE + 1);
    r1 += carry;
    carry <<= (LIMB_SIZE + 1);
    r0 -= carry;

    carry = r1 >> LIMB_SIZE;
    r2 += carry;
    carry <<= LIMB_SIZE;
    r1 -= carry;

    carry = r2 >> (LIMB_SIZE + 1);
    r3 += carry;
    carry <<= (LIMB_SIZE + 1);
    r2 -= carry;

    carry = r3 >> LIMB_SIZE;
    r4 += carry;
    carry <<= LIMB_SIZE;
    r3 -= carry;

    carry = r4 >> (LIMB_SIZE + 1);
    r5 += carry;
    carry <<= (LIMB_SIZE + 1);
    r4 -= carry;

    carry = r5 >> LIMB_SIZE;
    r6 += carry;
    carry <<= LIMB_SIZE;
    r5 -= carry;

    carry = r6 >> (LIMB_SIZE + 1);
    r7 += carry;
    carry <<= (LIMB_SIZE + 1);
    r6 -= carry;

    carry = r7 >> LIMB_SIZE;
    r8 += carry;
    carry <<= LIMB_SIZE;
    r7 -= carry;

    carry = r8 >> (LIMB_SIZE + 1);
    r9 += carry;
    carry <<= (LIMB_SIZE + 1);
    r8 -= carry;

    carry = r9 >> LIMB_SIZE; //25
    r0 += times19(carry);
    carry <<= LIMB_SIZE;
    r9 -= carry;

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
}

//For debugging
void print(__m256i lo, __m256i hi, char* str){
    /*
    radix 25 from vec
    */
    int* i_lo = (int*) &lo;
    int* i_hi = (int*) &hi;
    printf("%s:\n%d, %d, %d, %d, %d, %d, %d, %d, %d, %d,\n",str, i_lo[0], i_lo[1], i_lo[2], i_lo[3], i_lo[4], i_lo[5], i_lo[6], i_lo[7], i_hi[0], i_hi[1]);
}


void inverse_256(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi)
{
     /*
    res = x^-1 (the multiplicative inverse modulo p)
    */
    __m256i z2lo, z2hi, z9lo, z9hi, z11lo, z11hi, z2_5_0lo, z2_5_0hi, z2_10_0lo, z2_10_0hi, z2_20_0lo, z2_20_0hi, z2_50_0lo, z2_50_0hi, z2_100_0lo, z2_100_0hi, tlo, thi;
    int i;
    //print(xlo, xhi, "x");
    /* 2 */ square_256(&z2lo, &z2hi, xlo, xhi);
    //print(z2lo, z2hi, "z2");
    /* 4 */ square_256(&tlo, &thi, z2lo, z2hi);
    //print(tlo, thi, "t");
    /* 8 */ square_256(&tlo, &thi, tlo, thi);
    //print(tlo, thi, "t");
    /* 9 */ mul_256(&z9lo, &z9hi,tlo, thi, xlo, xhi);
    //print(z9lo, z9hi, "z9");
    /* 11 */ mul_256(&z11lo, &z11hi, z9lo, z9hi, z2lo, z2hi);
    //print(z11lo, z11hi, "z11");
    /* 22 */ square_256(&tlo, &thi, z11lo, z11hi);
    //print(tlo, thi, "t");
    /* 2^5 - 2^0 = 31 */ mul_256(&z2_5_0lo, &z2_5_0hi, tlo, thi, z9lo, z9hi);
    //print(z2_5_0lo, z2_5_0hi, "z2_5_0");
    /* 2^6 - 2^1 */ square_256(&tlo, &thi, z2_5_0lo, z2_5_0hi);
    //print(tlo, thi, "t");
    /* 2^10 - 2^5 */ for (i = 1; i < 5; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
        //print(tlo, thi, "t");
    }
    /* 2^10 - 2^0 */ mul_256(&z2_10_0lo, &z2_10_0hi, tlo, thi, z2_5_0lo, z2_5_0hi);
    //print(z2_10_0lo, z2_10_0hi, "z2_10_0");
    /* 2^11 - 2^1 */ square_256(&tlo, &thi, z2_10_0lo, z2_10_0hi);
    //print(tlo, thi, "t");
    /* 2^20 - 2^10 */ for (i = 1; i < 10; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
        //print(tlo, thi, "t");
    }
    /* 2^20 - 2^0 */ mul_256(&z2_20_0lo, &z2_20_0hi, tlo, thi, z2_10_0lo, z2_10_0hi);
    //print(z2_20_0lo, z2_20_0hi, "z2_20_0");
    /* 2^21 - 2^1 */ square_256(&tlo, &thi, z2_20_0lo,z2_20_0hi);
    //print(tlo, thi, "t");
    /* 2^40 - 2^20 */ for (i = 1; i < 20; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
        //print(tlo, thi, "t");
    }
    /* 2^40 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_20_0lo, z2_20_0hi);
    //print(tlo, thi, "t");
    /* 2^41 - 2^1 */ square_256(&tlo, &thi, tlo, thi);
    //print(tlo, thi, "t");
    /* 2^50 - 2^10 */ for (i = 1; i < 10; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
        //print(tlo, thi, "t");
    }
    /* 2^50 - 2^0 */ mul_256(&z2_50_0lo, &z2_50_0hi, tlo, thi, z2_10_0lo, z2_10_0hi);
    //print(z2_50_0lo, z2_50_0hi, "z2_50_0");
    /* 2^51 - 2^1 */ square_256(&tlo, &thi, z2_50_0lo, z2_50_0hi);
    //print(tlo, thi, "t");
    /* 2^100 - 2^50 */ for (i = 1; i < 50; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
        //print(tlo, thi, "t");
    }
    /* 2^100 - 2^0 */ mul_256(&z2_100_0lo, &z2_100_0hi, tlo, thi, z2_50_0lo,z2_50_0hi);
    //print(z2_100_0lo, z2_100_0hi, "z2_100_0");
    /* 2^101 - 2^1 */ square_256(&tlo, &thi, z2_100_0lo, z2_100_0hi);
    //print(tlo, thi, "t");
    /* 2^200 - 2^100 */ for (i = 1; i < 100; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
        //print(tlo, thi, "t");
    }
    /* 2^200 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_100_0lo, z2_100_0hi);
    //print(tlo, thi, "t");
    /* 2^201 - 2^1 */ square_256(&tlo, &thi, tlo, thi);
    //print(tlo, thi, "t");
    /* 2^250 - 2^50 */ for (i = 1; i < 50; i++)
    {
        square_256(&tlo, &thi, tlo, thi);
        //print(tlo, thi, "t");
    }
    /* 2^250 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_50_0lo, z2_50_0hi);
    //print(tlo, thi, "t");
    /* 2^251 - 2^1 */ square_256(&tlo, &thi, tlo, thi);
    //print(tlo, thi, "t");
    /* 2^252 - 2^2 */ square_256(&tlo, &thi, tlo, thi);
    //print(tlo, thi, "t");
    /* 2^253 - 2^3 */ square_256(&tlo, &thi, tlo, thi);
    //print(tlo, thi, "t");
    /* 2^254 - 2^4 */ square_256(&tlo, &thi, tlo, thi);
    //print(tlo, thi, "t");
    /* 2^255 - 2^5 */ square_256(&tlo, &thi, tlo, thi);
    //print(tlo, thi, "t");
    /* 2^255 - 21 */ mul_256(reslo, reshi, tlo, thi, z11lo, z11hi);
    //print(*reslo, *reshi, "res");
    //__m256i z2lo, z2hi, z9lo, z9hi, z11lo, z11hi, z2_5_0lo, z2_5_0hi, z2_10_0lo, z2_10_0hi, z2_20_0lo, z2_20_0hi, z2_50_0lo, z2_50_0hi, z2_100_0lo, z2_100_0hi, tlo, thi;
    //int i;
    ///* 2 */ square_256(&z2lo, &z2hi, xlo, xhi);
    ///* 4 */ square_256(&tlo, &thi, z2lo, z2hi);
    ///* 8 */ square_256(&tlo, &thi, tlo, thi);
    ///* 9 */ mul_256(&z9lo, &z9hi, tlo, thi, xlo, xhi);
    ///* 11 */ mul_256(&z11lo, &z11hi, z9lo, z9hi, z2lo, z2hi);
    ///* 22 */ square_256(&tlo, &thi, z11lo, z11hi);
    ///* 2^5 - 2^0 = 31 */ mul_256(&z2_5_0lo, &z2_5_0hi, tlo, thi, z9lo, z9hi);
    ///* 2^6 - 2^1 */ square_256(&tlo, &thi, z2_5_0lo, z2_5_0hi);
    ///* 2^10 - 2^5 */ for (i = 1; i < 5; i++)
    //{
    //    square_256(&tlo, &thi, tlo, thi);
    //}
    ///* 2^10 - 2^0 */ mul_256(&z2_10_0lo, &z2_10_0hi, tlo, thi, z2_5_0lo, z2_5_0hi);
    ///* 2^11 - 2^1 */ square_256(&tlo, &thi, z2_10_0lo, z2_10_0hi);
    ///* 2^20 - 2^10 */ for (i = 1; i < 10; i++)
    //{
    //    square_256(&tlo, &thi, tlo, thi);
    //}
    ///* 2^20 - 2^0 */ mul_256(&z2_20_0lo, &z2_20_0hi, tlo, thi, z2_10_0lo, z2_10_0hi);
    ///* 2^21 - 2^1 */ square_256(&tlo, &thi, z2_20_0lo, z2_20_0hi);
    ///* 2^40 - 2^20 */ for (i = 1; i < 20; i++)
    //{
    //    square_256(&tlo, &thi, tlo, thi);
    //}
    ///* 2^40 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_20_0lo, z2_20_0hi);
    ///* 2^41 - 2^1 */ square_256(&tlo, &thi, tlo, thi);
    ///* 2^50 - 2^10 */ for (i = 1; i < 10; i++)
    //{
    //    square_256(&tlo, &thi, tlo, thi);
    //}
    ///* 2^50 - 2^0 */ mul_256(&z2_50_0lo, &z2_50_0hi, tlo, thi, z2_10_0lo, z2_10_0hi);
    ///* 2^51 - 2^1 */ square_256(&tlo, &thi, z2_50_0lo, z2_50_0hi);
    ///* 2^100 - 2^50 */ for (i = 1; i < 50; i++)
    //{
    //    square_256(&tlo, &thi, tlo, thi);
    //}
    ///* 2^100 - 2^0 */ mul_256(&z2_100_0lo, &z2_100_0hi, tlo, thi, z2_50_0lo, z2_50_0hi);
    ///* 2^101 - 2^1 */ square_256(&tlo, &thi, z2_100_0lo, z2_100_0hi);
    ///* 2^200 - 2^100 */ for (i = 1; i < 100; i++)
    //{
    //    square_256(&tlo, &thi, tlo, thi);
    //}
    ///* 2^200 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_100_0lo, z2_100_0hi);
    ///* 2^201 - 2^1 */ square_256(&tlo, &thi, tlo, thi);
    ///* 2^250 - 2^50 */ for (i = 1; i < 50; i++)
    //{
    //    square_256(&tlo, &thi, tlo, thi);
    //}
    ///* 2^250 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_50_0lo, z2_50_0hi);
    ///* 2^251 - 2^1 */ square_256(&tlo, &thi, tlo, thi);
    ///* 2^252 - 2^2 */ square_256(&tlo, &thi, tlo, thi);
    ///* 2^253 - 2^3 */ square_256(&tlo, &thi, tlo, thi);
    ///* 2^254 - 2^4 */ square_256(&tlo, &thi, tlo, thi);
    ///* 2^255 - 2^5 */ square_256(&tlo, &thi, tlo, thi);
    ///* 2^255 - 21 */ mul_256(reslo, reshi, tlo, thi, z11lo, z11hi);
}

void inverse_256_non_vec(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi)
{
    /*
    res = x^-1 (the multiplicative inverse modulo p)
    */
    __m256i z2lo, z2hi, z9lo, z9hi, z11lo, z11hi, z2_5_0lo, z2_5_0hi, z2_10_0lo, z2_10_0hi, z2_20_0lo, z2_20_0hi, z2_50_0lo, z2_50_0hi, z2_100_0lo, z2_100_0hi, tlo, thi;
    int i;
    /* 2 */ square_256(&z2lo, &z2hi, xlo, xhi);
    /* 4 */ square_256(&tlo, &thi, z2lo, z2hi);
    /* 8 */ square_256(&tlo, &thi, tlo, thi);
    /* 9 */ mul_256(&z9lo, &z9hi, tlo, thi, xlo, xhi);
    /* 11 */ mul_256(&z11lo, &z11hi, z9lo, z9hi, z2lo, z2hi);
    /* 22 */ square_256(&tlo, &thi, z11lo, z11hi);
    /* 2^5 - 2^0 = 31 */ mul_256(&z2_5_0lo, &z2_5_0hi, tlo, thi, z9lo, z9hi);
    /* 2^6 - 2^1 */ square_256(&tlo, &thi, z2_5_0lo, z2_5_0hi);
    // to array
    signed int *alo = (int *)&tlo;
    signed int *ahi = (int *)&thi;
    /* 2^10 - 2^5 */ for (i = 1; i < 5; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256((const __m256i *)alo);
    thi = _mm256_loadu_si256((const __m256i *)ahi);

    /* 2^10 - 2^0 */ mul_256(&z2_10_0lo, &z2_10_0hi, tlo, thi, z2_5_0lo, z2_5_0hi);
    /* 2^11 - 2^1 */ square_256(&tlo, &thi, z2_10_0lo, z2_10_0hi);
    // to array
    alo = (int *)&tlo;
    ahi = (int *)&thi;
    /* 2^20 - 2^10 */ for (i = 1; i < 10; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256((const __m256i *)alo);
    thi = _mm256_loadu_si256((const __m256i *)ahi);
    /* 2^20 - 2^0 */ mul_256(&z2_20_0lo, &z2_20_0hi, tlo, thi, z2_10_0lo, z2_10_0hi);
    /* 2^21 - 2^1 */ square_256(&tlo, &thi, z2_20_0lo, z2_20_0hi);
    // to array
    alo = (int *)&tlo;
    ahi = (int *)&thi;
    /* 2^40 - 2^20 */ for (i = 1; i < 20; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256((const __m256i *)alo);
    thi = _mm256_loadu_si256((const __m256i *)ahi);
    /* 2^40 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_20_0lo, z2_20_0hi);
    /* 2^41 - 2^1 */ //square_256(&tlo, &thi, tlo, thi);
    // to array
    alo = (int *)&tlo;
    ahi = (int *)&thi;
    /* 2^50 - 2^10 */ for (i = 1; i < 11; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    /* 2^50 - 2^0 */ mul_256(&z2_50_0lo, &z2_50_0hi, tlo, thi, z2_10_0lo, z2_10_0hi);
    /* 2^51 - 2^1 */ square_256(&tlo, &thi, z2_50_0lo, z2_50_0hi);
    // to array
    alo = (int *)&tlo;
    ahi = (int *)&thi;
    /* 2^100 - 2^50 */ for (i = 1; i < 50; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256((const __m256i *)alo);
    thi = _mm256_loadu_si256((const __m256i *)ahi);
    /* 2^100 - 2^0 */ mul_256(&z2_100_0lo, &z2_100_0hi, tlo, thi, z2_50_0lo, z2_50_0hi);
    /* 2^101 - 2^1 */ square_256(&tlo, &thi, z2_100_0lo, z2_100_0hi);
    // to array
    alo = (int *)&tlo;
    ahi = (int *)&thi;
    /* 2^200 - 2^100 */ for (i = 1; i < 100; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256((const __m256i *)alo);
    thi = _mm256_loadu_si256((const __m256i *)ahi);
    /* 2^200 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_100_0lo, z2_100_0hi);
    /* 2^201 - 2^1 */ //square_256(&tlo, &thi, tlo, thi);
    // to array
    alo = (int *)&tlo;
    ahi = (int *)&thi;
    /* 2^250 - 2^50 */ for (i = 1; i < 51; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256((const __m256i *)alo);
    thi = _mm256_loadu_si256((const __m256i *)ahi);
    /* 2^250 - 2^0 */ mul_256(&tlo, &thi, tlo, thi, z2_50_0lo, z2_50_0hi);
    // to array
    alo = (int *)&tlo;
    ahi = (int *)&thi;
    /* 2^10 - 2^5 */ for (i = 0; i < 5; i++)
    {
        square_256_non_vec(alo, ahi, alo, ahi);
    }
    //to vector
    tlo = _mm256_loadu_si256((const __m256i *)alo);
    thi = _mm256_loadu_si256((const __m256i *)ahi);
    /* 2^255 - 21 */ mul_256(reslo, reshi, tlo, thi, z11lo, z11hi);
}

void div_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi)
{
    // inverse and mult
    __m256i invlo, invhi;
    inverse_256(&invlo, &invhi, blo, bhi);
    mul_256(reslo, reshi, alo, ahi, invlo, invhi);
}