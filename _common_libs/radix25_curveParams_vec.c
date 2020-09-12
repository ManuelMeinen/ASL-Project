#include "radix25_curveParams_vec.h"
#include "vec_utils.h"


__m256i A_rad25lo, A_rad25hi, B_rad25lo, B_rad25hi, Gx_rad25lo, Gx_rad25hi, Gy_rad25lo, Gy_rad25hi;

void init_curve25519_params_rad25() {
    
    /*
    init intrinsic vectors
    	A_rad25lo, A_rad25hi, B_rad25lo, B_rad25hi, Gx_rad25lo, Gx_rad25hi, Gy_rad25lo, Gy_rad25hi
    */
    A_rad25lo = _mm256_setr_epi32( 486662, 0, 0, 0, 0, 0, 0, 0 );
    A_rad25hi = _mm256_set1_epi32( 0 );
    
    B_rad25lo = _mm256_setr_epi32( 1, 0, 0, 0, 0, 0, 0, 0 );
    B_rad25hi = _mm256_set1_epi32( 0 );
    
    Gx_rad25lo = _mm256_setr_epi32( 9, 0, 0, 0, 0, 0, 0, 0 );
    Gx_rad25hi = _mm256_set1_epi32( 0 );

    Gy_rad25lo = _mm256_setr_epi32( 47109081, 7432351, 3556669, 32730083, 38335797, 21549388, 57579067, 14064643);
    Gy_rad25hi = _mm256_setr_epi32( 35359240, 8566886,0,0,0,0,0,0);
    
	/*
    init intrinsic vector constants
    	const_1_rad25lo, const_1_rad25hi, const_2_rad25lo, const_2_rad25hi, const_3_rad25lo, const_3_rad25hi, const_4_rad25lo, const_4_rad25hi
    */
    const_1_rad25lo = _mm256_setr_epi32( 1, 0, 0, 0, 0, 0, 0, 0 );
    const_1_rad25hi = _mm256_set1_epi32( 0 );

    const_2_rad25lo = _mm256_setr_epi32( 2, 0, 0, 0, 0, 0, 0, 0 );
    const_2_rad25hi = _mm256_set1_epi32( 0 );
    
    const_3_rad25lo = _mm256_setr_epi32( 3, 0, 0, 0, 0, 0, 0, 0 );
    const_3_rad25hi = _mm256_set1_epi32( 0 );

    const_4_rad25lo = _mm256_setr_epi32( 4, 0, 0, 0, 0, 0, 0, 0 );
    const_4_rad25hi = _mm256_set1_epi32( 0 );
    
    /*
    init constant index vectors
    */
	init_vec_constants();
}	
