#ifndef __BIGNUM_RAD25_H__
#include <stdio.h>
#include <string.h>
#include "../../_common_libs/radix25.h"
#define __BIGNUM_RAD25_H__

void add_256(radix25_t *res, radix25_t *a, radix25_t *b);
void reduce_256(radix25_t *x);
void reduce_mul_512(unsigned long r[29]);
void reduce_mul_256(unsigned long r[29]);
void sub_256(radix25_t *res, radix25_t *a, radix25_t *b);
void mul_256(radix25_t *res, radix25_t *a, radix25_t *b);
void div_256(radix25_t *res, radix25_t *a, radix25_t *b);
int cmp_point_256(radix25_t *x1, radix25_t *y1, radix25_t *x2, radix25_t *y2);
int cmp_coord_256(radix25_t *x1, radix25_t *x2);
int iszero_256(radix25_t *x);
void set_ui_256(radix25_t *x, signed int ui);
void copy_256(radix25_t *dest, radix25_t src);
void square_256(radix25_t *res, radix25_t *x);
void inverse_256(radix25_t *res, radix25_t *x);


#endif /* __BIGNUM_RAD25_H__ */