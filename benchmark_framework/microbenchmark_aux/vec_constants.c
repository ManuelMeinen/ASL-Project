#include "vec_constants.h"
#include <immintrin.h>
#include <emmintrin.h>
#include <stdio.h>

void init_vec_constants() {
	index1_vec = _mm256_setr_epi32( 0, 0, 1, 1, 2, 2, 3, 3 );
	index2_vec = _mm256_setr_epi32( 4, 4, 5, 5, 6, 6, 7, 7 );
	index3_vec = _mm256_setr_epi32( 3, 3, 4, 4, 5, 5, 6, 6 );
    index4_vec = _mm256_setr_epi32( 0, 0, 0, 0, 1, 1, 2, 2 );
    index5_vec = _mm256_setr_epi32( 7, 7, 0, 0, 0, 0, 0, 0 );
	index_alt_vec = _mm256_setr_epi32( 1, 0, 1, 0, 1, 0, 1, 0 );
	const_0_vec = _mm256_set1_epi32( 0 );
	full_vec = _mm256_set1_epi32(-1);
	ones1_vec = _mm256_set1_epi32(-1);
	ones2_vec = _mm256_setr_epi32(-1,-1,-1,-1,-1,-1,-1,0);
}