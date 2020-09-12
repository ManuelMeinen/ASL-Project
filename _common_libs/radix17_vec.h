#ifndef __RADIX17_VEC_H__
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <immintrin.h>
#include <emmintrin.h>
#define __RADIX17_VEC_H__

void print_radix17_vec(__m256i xlo, __m256i xhi, char *string);
void radix17_vec_toString( char* dest, __m256i xlo, __m256i xhi ) ;


#endif /* __RADIX17_H__ */