#include "vec_utils.h"
#include <immintrin.h>
#include <emmintrin.h>
#include <stdio.h>

void init_vec_constants() {
	index1_vec = _mm256_setr_epi32( 0, 0, 1, 1, 2, 2, 3, 3 );
	index2_vec = _mm256_setr_epi32( 4, 4, 5, 5, 6, 6, 7, 7 );
	index3_vec = _mm256_setr_epi32( 3, 3, 4, 4, 5, 5, 6, 6 );
    index4_vec = _mm256_setr_epi32( 0, 0, 0, 0, 1, 1, 2, 2 );
    index5_vec = _mm256_setr_epi32( 7, 7, 0, 0, 0, 0, 0, 0 );
	index_alt_vec = _mm256_setr_epi32( 1, 0, 1, 0, 1, 0, 1, 0 );
	const_0_vec = _mm256_set1_epi32( 0 );
	full_vec = _mm256_set1_epi32(-1);
	ones1_vec = _mm256_set1_epi32(-1);
	ones2_vec = _mm256_setr_epi32(-1,-1,-1,-1,-1,-1,-1,0);
	
	c_alt_0_1__64_vec    = _mm256_setr_epi64x(0,1,0,1);
	c_all_19__64_vec     = _mm256_set1_epi64x(19);
	c_2x19_2x0__64_vec   = _mm256_setr_epi64x(19,19,0,0);
	c_1_38_19_38__64_vec = _mm256_setr_epi64x(1,38,19,38);
	c_1_2_1_38__64_vec   = _mm256_setr_epi64x(1,2,1,38);
	c_2x1_2x19__64_vec   = _mm256_setr_epi64x(1,1,19,19);
	c_alt_19_38__64_vec  = _mm256_setr_epi64x(19,38,19,38);
	c_19_38_2x0__64_vec  = _mm256_setr_epi64x(19,38,0,0);
	c_1_38_2x0__64_vec   = _mm256_setr_epi64x(1,38,0,0);
    c_0_1_1_2__64_vec = _mm256_setr_epi64x(0,1,1,2);
    c_0_38_19_38__64_vec = _mm256_setr_epi64x(0,38,19,38);
    c_2x0_2x19__64_vec = _mm256_setr_epi64x(0,0,19,19);
    c_3x0_38__64_vec = _mm256_setr_epi64x(0,0,0,38);
    c_0_38_2x0__64_vec = _mm256_setr_epi64x(0,38,0,0);
}

// print functions
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

// check functions (for debugging)
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
        printf("right values: %ld %ld %ld %ld\n", v0, v1, v2, v3);
    }
}

void check_vec_int(__m256i vec, int v0, int v1, int v2, int v3, int v4, int v5, int v6, int v7){
    /*
    check if vec (8xint) is equal to (v0, v1, v2, v3, v4, v5, v6, v7)
    */
    int* i_v = (int*) &vec;
    if(i_v[0] == v0 && i_v[1] == v1 && i_v[2] == v2 && i_v[3] == v3 && i_v[4] == v4 && i_v[5] == v5 && i_v[6] == v6 && i_v[7] == v7){
        printf("OK\n");
    }else{
        printf("Not OK!\n");
        print_vec_int(vec, "wrong vector");
        printf("right values: %d %d %d %d %d %d %d %d\n", v0, v1, v2, v3, v4, v5, v6, v7);
    }
}

// point vector copy
void copy_256(__m256i *dst_vec_lo, __m256i *dst_vec_hi, __m256i src_vec_lo, __m256i src_vec_hi)
{
    /*
    set dest to src
    */
    __m256i dummy;
    *dst_vec_hi = _mm256_blend_epi32(src_vec_hi,dummy, 0);
    *dst_vec_lo = _mm256_blend_epi32(src_vec_lo,dummy, 0);
}

// point vector compare functions
int cmp_coord_256(__m256i x1_lo, __m256i x1_hi, __m256i x2_lo,__m256i x2_hi)
{
    /*
    return 1 if x1 == x2
    */

    __m256i res0 = _mm256_cmpeq_epi32(x1_lo, x2_lo);
    __m256i res1 = _mm256_cmpeq_epi32(x1_hi, x2_hi);
    
    res0 = _mm256_xor_si256(res0, full_vec);
    res1 = _mm256_xor_si256(res1, full_vec);

    int r0 = _mm256_movemask_epi8(res0);
    int r1 = _mm256_movemask_epi8(res1);
   
    return !r0&&!r1;
}

int cmp_point_256(__m256i x1_lo, __m256i x1_hi,__m256i y1_lo, __m256i y1_hi,  __m256i x2_lo,__m256i x2_hi, __m256i y2_lo, __m256i y2_hi)
{
    /*
    return 1 if x1 == y1 and x2 == y2
    */
    return cmp_coord_256(x1_lo, x1_hi, y1_lo, y1_hi) && cmp_coord_256(x2_lo, x2_hi, y2_lo, y2_hi);
}

int iszero_256(__m256i xlo, __m256i xhi)
{
    /*
    return 1 if x == y == 0
    */
    int cmp1, cmp2;
    cmp1 = _mm256_testz_si256(ones1_vec, xlo);
    cmp2 = _mm256_testz_si256(ones2_vec, xhi);    
    return cmp1&&cmp2;
}

// point vector set functions
void set_ui_256(__m256i *xlo, __m256i *xhi, signed int ui)
{
    /*
    set radix17_t to a unsigned int ui (lowest limb only, others 0)
    maximum 17 bit numbers !!!
    */
    *xhi = _mm256_set1_epi32(0);
    *xlo = _mm256_setr_epi32(ui, 0,0,0,0,0,0,0);
}

void set_zero(__m256i *xlo, __m256i *xhi)
{
    /*
    set radix17_t to a unsigned int ui (lowest limb only, others 0)
    maximum 17 bit numbers !!!
    */
    *xhi = _mm256_set1_epi32(0);
    *xlo = _mm256_set1_epi32(0);
}

