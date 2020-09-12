#ifndef __BIGNUM_RAD17_H__
#include <stdio.h>
#include <string.h>
#include "../../_common_libs/radix17.h"
#define __BIGNUM_RAD17_H__

void add_256_17(radix17_t *res, radix17_t *a, radix17_t *b);
void add_256_17_no_reduce(radix17_t *res, radix17_t *a, radix17_t *b);
void sub_256_17(radix17_t *res, radix17_t *a, radix17_t *b);
void sub_256_17_no_reduce(radix17_t *res, radix17_t *a, radix17_t *b);
void mul_256_17(radix17_t *res, radix17_t *a, radix17_t *b);
void div_256_17(radix17_t *res, radix17_t *a, radix17_t *b);
void square_256_17(radix17_t *res, radix17_t *x);
void inverse_256_17(radix17_t *res, radix17_t *x);


#endif /* __BIGNUM_RAD17_H__ */