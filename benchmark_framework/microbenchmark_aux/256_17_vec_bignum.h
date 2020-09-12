#ifndef __BIGNUM_RAD17_VEC_H__
#include <stdio.h>
#include <string.h>
#include "../../_common_libs/radix17_vec.h"
#define __BIGNUM_RAD17_VEC_H__

void add_256_17_vec(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void reduce_2_to_2_vec(__m256i *xlo, __m256i *xhi);
void reduce_mul_512_17_vec(__m256i *r0, __m256i *r1, __m256i *r2, __m256i *r3, __m256i *r4, __m256i *r5, __m256i *r6, __m256i *r7);
void reduce_4_to_2_vec(__m256i *rv0, __m256i *rv1, __m256i *rv2, __m256i *rv3);
void sub_256_17_vec(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void mul_256_17_vec(__m256i *reslo, __m256i *reshi, __m256i alo,__m256i ahi, __m256i blo,__m256i bhi);
void div_256_17_vec(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
int cmp_point_256_17_vec(__m256i x1_lo, __m256i x1_hi,__m256i y1_lo, __m256i y1_hi,  __m256i x2_lo,__m256i x2_hi, __m256i y2_lo, __m256i y2_hi);
int cmp_coord_256_17_vec(__m256i x1_lo, __m256i x1_hi, __m256i x2_lo,__m256i x2_hi);
int iszero_256_17_vec(__m256i xlo, __m256i xhi);
void set_ui_256_17_vec(__m256i *xlo, __m256i *xhi, signed int ui);
void set_zero(__m256i *xlo, __m256i *xhi);
void copy_256_17_vec(__m256i *dst_vec_lo, __m256i *dst_vec_hi, __m256i src_vec_lo, __m256i src_vec_hi);
void square_256_17_vec(__m256i *reslo, __m256i *reshi, __m256i xlo,__m256i xhi);
void square_256_17_vec_alt(__m256i *reslo, __m256i *reshi, __m256i xlo,__m256i xhi);
void inverse_256_17_vec(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi);
void print_vec_int(__m256i vec, char* str);
#endif /* __BIGNUM_RAD17_VEC_H__ */