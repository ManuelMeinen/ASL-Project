#ifndef __VEC_UTILS_H__
#include <immintrin.h>
#define __VEC_UTILS_H__

// vector constants
// - mask vectors
__m256i index1_vec, index2_vec, index3_vec, index4_vec, index5_vec, index_alt_vec, full_vec, ones1_vec, ones2_vec;
__m256i c_alt_0_1__64_vec;
// - constant vectors
__m256i const_0_vec;
__m256i c_all_19__64_vec;
__m256i c_2x19_2x0__64_vec;
__m256i c_1_38_19_38__64_vec;
__m256i c_1_2_1_38__64_vec;
__m256i c_2x1_2x19__64_vec;
__m256i c_alt_19_38__64_vec;
__m256i c_19_38_2x0__64_vec;
__m256i c_1_38_2x0__64_vec;
__m256i c_0_1_1_2__64_vec;
__m256i c_0_38_19_38__64_vec; 
__m256i c_2x0_2x19__64_vec;
__m256i c_3x0_38__64_vec;
__m256i c_0_38_2x0__64_vec;

void init_vec_constants();

// vector print functions
void print_vec_int(__m256i vec, char* str);
void print_vec_long(__m256i vec, char* str);

// check functions (for debugging)
void check_vec_long(__m256i vec, long v0, long v1, long v2, long v3);
void check_vec_int(__m256i vec, int v0, int v1, int v2, int v3, int v4, int v5, int v6, int v7);

// point vector copy
void copy_256(__m256i *dst_vec_lo, __m256i *dst_vec_hi, __m256i src_vec_lo, __m256i src_vec_hi);

// point vector compare functions
int cmp_coord_256(__m256i x1_lo, __m256i x1_hi, __m256i x2_lo,__m256i x2_hi);
int cmp_point_256(__m256i x1_lo, __m256i x1_hi,__m256i y1_lo, __m256i y1_hi,  __m256i x2_lo,__m256i x2_hi, __m256i y2_lo, __m256i y2_hi);
int iszero_256(__m256i xlo, __m256i xhi);

// point vector set functions
void set_ui_256(__m256i *xlo, __m256i *xhi, signed int ui);
void set_zero(__m256i *xlo, __m256i *xhi);

#endif /* __VEC_UTILS_H__ */
