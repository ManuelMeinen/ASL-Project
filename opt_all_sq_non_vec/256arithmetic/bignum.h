#ifndef __BIGNUM_RAD17_H__
#include <stdio.h>
#include <string.h>
#include "../../_common_libs/radix17_vec.h"
#define __BIGNUM_RAD17_H__

void add_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void sub_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void mul_256(__m256i *reslo, __m256i *reshi, __m256i alo,__m256i ahi, __m256i blo,__m256i bhi);
void div_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void square_256(__m256i *reslo, __m256i *reshi, __m256i xlo,__m256i xhi);
void inverse_256(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi);

void reduce_2_to_2(__m256i *xlo, __m256i *xhi);
void reduce_mul_512(__m256i *r0, __m256i *r1, __m256i *r2, __m256i *r3, __m256i *r4, __m256i *r5, __m256i *r6, __m256i *r7);
void reduce_4_to_2(__m256i *rv0, __m256i *rv1, __m256i *rv2, __m256i *rv3);

void square_256_non_vec(signed int reslo[8], signed int reshi[8], signed int xlo[8], signed int xhi[8]);

void reduce_2_to_2_non_vec(signed int xlo[8], signed int xhi[8]);
void reduce_8_to_4_non_vec(unsigned long r[29]);
void reduce_4_to_2_non_vec(unsigned long r[29]);


#endif /* __BIGNUM_RAD17_H__ */