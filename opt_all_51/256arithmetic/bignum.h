#ifndef __BIGNUM_RAD51_H__
#include <stdio.h>
#include <string.h>
#define __BIGNUM_RAD51_H__

void set_ui_256(long x[5], signed int ui);
void copy_256(long dest[5], long src[5]);

void add_256(long res[5], long a[5], long b[5]);
void add_1_256(long a[5]);
void sub_256(long res[5], long a[5], long b[5]);
void mul_256(long res[5], long x[5], long y[5]);
void div_256(long res[5], long a[5], long b[5]);
void square_256 (long res[5], long a[5]);
void inverse_256(long res[5], long x[5]);

int cmp_point_256(long x1[5], long y1[5], long x2[5], long y2[5]);
int cmp_coord_256(long x1[5], long x2[5]);
int iszero_256(long x[5]);

void reduce_2_to_2(long x[5]);


#endif /* __BIGNUM_RAD51_H__ */