#ifndef __radix25_VEC_H__
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <immintrin.h>
#include <emmintrin.h>
#define __radix25_VEC_H__

void print_radix25_vec(__m256i xlo, __m256i xhi, char *string);
void radix25_vec_toString( char* dest, __m256i xlo, __m256i xhi ) ;


#endif /* __radix25_H__ */