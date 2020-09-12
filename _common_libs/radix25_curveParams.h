#ifndef __radix25_curveParams_H__
#include <stdio.h>
#include <string.h>
#include "radix25.h"
#define __radix25_curveParams_H__


// points
radix25_t A_rad25, B_rad25, Gx_rad25, Gy_rad25;
// constants
radix25_t const_1_rad25, const_2_rad25, const_3_rad25, const_4_rad25;

void print_radix25(radix25_t *x, char *string);
void radix25_toString( char* dest, radix25_t x );
void init_curve25519_params_rad25();

#endif /* __radix25_curveParams_H__ */