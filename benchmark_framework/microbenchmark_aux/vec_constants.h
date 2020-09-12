#ifndef __VEC_CONSTS_H__
#include <immintrin.h>
#define __VEC_CONSTS_H__

// vector constants
// - index vectors
__m256i index1_vec, index2_vec, index3_vec, index4_vec, index5_vec, index_alt_vec, full_vec, ones1_vec, ones2_vec;
// - constant vectors
__m256i const_0_vec;

void init_vec_constants();

#endif /* __VEC_CONSTS_H__ */
