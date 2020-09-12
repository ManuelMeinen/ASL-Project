#ifndef __BIGNUM_RAD25_H__
#include <stdio.h>
#include <string.h>
#include "../../_common_libs/radix25_vec.h"
#define __BIGNUM_RAD25_H__

void reduce_2_to_2(__m256i *reslo, __m256i *reshi);
void reduce_256(unsigned long* r);
void reduce_mul_512(unsigned long r[29]);
void reduce_mul_256(unsigned long r[29]);

void add_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void sub_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void mul_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void div_256(__m256i *reslo, __m256i *reshi, __m256i alo, __m256i ahi, __m256i blo, __m256i bhi);
void square_256(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi);
void inverse_256(__m256i *reslo, __m256i *reshi, __m256i xlo, __m256i xhi);

void square_256_non_vec(signed int reslo[8], signed int reshi[8], signed int xlo[8], signed int xhi[8]);
void reduce_3_to_2_non_vec(unsigned long r[10]);
void reduce_2_to_2_non_vec(signed int xlo[8], signed int xhi[8]);

#endif /* __BIGNUM_RAD25_H__ */