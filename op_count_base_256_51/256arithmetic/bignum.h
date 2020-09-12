#ifndef __BIGNUM_RAD51_H__
#include <stdio.h>
#include <string.h>
#include "../../_common_libs/radix51.h"
#define __BIGNUM_RAD51_H__

void add_256(radix51_t *res, radix51_t *a, radix51_t *b);
void reduce_2_to_2(radix51_t *x);
void reduce_mul_512(__int128 r[9]);
void reduce_4_to_2(__int128 r[9]);
void sub_256(radix51_t *res, radix51_t *a, radix51_t *b);
void mul_256(radix51_t *res, radix51_t *a, radix51_t *b);
void div_256(radix51_t *res, radix51_t *a, radix51_t *b);
int cmp_point_256(radix51_t *x1, radix51_t *y1, radix51_t *x2, radix51_t *y2);
int cmp_coord_256(radix51_t *x1, radix51_t *x2);
int iszero_256(radix51_t *x);
void set_ui_256(radix51_t *x, unsigned int ui);
void copy_256(radix51_t *dest, radix51_t src);
void square_256(radix51_t *res, radix51_t *x);
void inverse_256(radix51_t *res, radix51_t *x);
int get_total_ops();

#endif /* __BIGNUM_RAD51_H__ */