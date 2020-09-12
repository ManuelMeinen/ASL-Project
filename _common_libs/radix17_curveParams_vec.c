#include "radix17_curveParams_vec.h"
#include "vec_utils.h"


__m256i A_rad17lo, A_rad17hi, B_rad17lo, B_rad17hi, Gx_rad17lo, Gx_rad17hi, Gy_rad17lo, Gy_rad17hi;

void init_curve25519_params_rad17() {
    
    /*
    init intrinsic vectors
    	A_rad17lo, A_rad17hi, B_rad17lo, B_rad17hi, Gx_rad17lo, Gx_rad17hi, Gy_rad17lo, Gy_rad17hi
    */
    A_rad17lo = _mm256_setr_epi32( 93446, 3, 0, 0, 0, 0, 0, 0 );
    A_rad17hi = _mm256_set1_epi32( 0 );
    
    B_rad17lo = _mm256_setr_epi32( 1, 0, 0, 0, 0, 0, 0, 0 );
    B_rad17hi = _mm256_set1_epi32( 0 );
    
    Gx_rad17lo = _mm256_setr_epi32( 9, 0, 0, 0, 0, 0, 0, 0 );
    Gx_rad17hi = _mm256_set1_epi32( 0 );
    
    Gy_rad17lo = _mm256_setr_epi32( 54233, 81767, 29032, 17725, 116251, 127851, 62773, 39204 );
    Gy_rad17hi = _mm256_setr_epi32( 84177, 38459, 1975 , 54940, 100872, 52493 , 33464, 0 );
    
	/*
    init intrinsic vector constants
    	const_1_rad17lo, const_1_rad17hi, const_2_rad17lo, const_2_rad17hi, const_3_rad17lo, const_3_rad17hi, const_4_rad17lo, const_4_rad17hi
    */
    const_1_rad17lo = _mm256_setr_epi32( 1, 0, 0, 0, 0, 0, 0, 0 );
    const_1_rad17hi = _mm256_set1_epi32( 0 );

    const_2_rad17lo = _mm256_setr_epi32( 2, 0, 0, 0, 0, 0, 0, 0 );
    const_2_rad17hi = _mm256_set1_epi32( 0 );
    
    const_3_rad17lo = _mm256_setr_epi32( 3, 0, 0, 0, 0, 0, 0, 0 );
    const_3_rad17hi = _mm256_set1_epi32( 0 );

    const_4_rad17lo = _mm256_setr_epi32( 4, 0, 0, 0, 0, 0, 0, 0 );
    const_4_rad17hi = _mm256_set1_epi32( 0 );
    
    /*
    init constant index vectors
    */
	init_vec_constants();
}	
