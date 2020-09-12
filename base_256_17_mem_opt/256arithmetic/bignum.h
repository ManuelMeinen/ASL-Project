#ifndef __BIGNUM_RAD17_H__
#include <stdio.h>
#include <string.h>
#include "../../_common_libs/radix17.h"
#define __BIGNUM_RAD17_H__

void add_256(radix17_t *res, radix17_t *a, radix17_t *b);
void reduce_2_to_2(radix17_t *x);
void reduce_mul_512(unsigned long r[29]);
void reduce_4_to_2(unsigned long r[29]);
void sub_256(radix17_t *res, radix17_t *a, radix17_t *b);
void mul_256(radix17_t *res, radix17_t *a, radix17_t *b);
void div_256(radix17_t *res, radix17_t *a, radix17_t *b);
int cmp_point_256(radix17_t *x1, radix17_t *y1, radix17_t *x2, radix17_t *y2);
int cmp_coord_256(radix17_t *x1, radix17_t *x2);
int iszero_256(radix17_t *x);
void set_ui_256(radix17_t *x, signed int ui);
void copy_256(radix17_t *dest, radix17_t src);
void square_256(radix17_t *res, radix17_t *x);
void inverse_256(radix17_t *res, radix17_t *x);


#endif /* __BIGNUM_RAD17_H__ */