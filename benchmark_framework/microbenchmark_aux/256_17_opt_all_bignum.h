#ifndef __BIGNUM_RAD17_OPT_ALL_H__
#include <stdio.h>
#include <string.h>
#include "../../_common_libs/radix17_vec.h"
#define __BIGNUM_RAD17_OPT_ALL_H__

void add_256_17_all(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void reduce_2_to_2_17_all(__m256i *xlo, __m256i *xhi);
void reduce_4_to_2_17_all(__m256i *rv0, __m256i *rv1, __m256i *rv2, __m256i *rv3);
void sub_256_17_all(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void mul_256_17_all(__m256i *reslo, __m256i *reshi, __m256i alo,__m256i ahi, __m256i blo,__m256i bhi);
void div_256_17_all(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void div_256_non_vec(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void div_256_non_vec_alt(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
int cmp_point_256_17_all(__m256i x1_lo, __m256i x1_hi,__m256i y1_lo, __m256i y1_hi,  __m256i x2_lo,__m256i x2_hi, __m256i y2_lo, __m256i y2_hi);
int cmp_coord_256_17_all(__m256i x1_lo, __m256i x1_hi, __m256i x2_lo,__m256i x2_hi);
int iszero_256_17_all(__m256i xlo, __m256i xhi);
void set_ui_256_17_all(__m256i *xlo, __m256i *xhi, signed int ui);
void set_zero_17_all(__m256i *xlo, __m256i *xhi);
void copy_256_17_all(__m256i *dst_vec_lo, __m256i *dst_vec_hi, __m256i src_vec_lo, __m256i src_vec_hi);
void square_256_17_all(__m256i *reslo, __m256i *reshi, __m256i xlo,__m256i xhi);
void square_256_17_all_new(__m256i *reslo, __m256i *reshi, __m256i xlo,__m256i xhi);
void inverse_256_17_all(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi);
void reduce_2_to_2_non_vec(signed int xlo[8], signed int xhi[8]);
void reduce_4_to_2_non_vec(unsigned long r[29]);
void reduce_8_to_4_non_vec(unsigned long r[29]);
void square_256_17_all_non_vec(signed int reslo[8], signed int reshi[8], signed int xlo[8], signed int xhi[8]);
void square_256_opt_vec(__m256i *reslo, __m256i *reshi, __m256i xlo,__m256i xhi);
void inverse_256_non_vec(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi);
void inverse_256_non_vec_alt(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi);

#endif /* __BIGNUM_RAD17_OPT_ALL_H__ */