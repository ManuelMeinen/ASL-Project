#ifndef __radix25_curveParams_vec_H__
#include "radix25_vec.h"
#define __radix25_curveParams_vec_H__

// points
__m256i A_rad25lo, A_rad25hi, B_rad25lo, B_rad25hi, Gx_rad25lo, Gx_rad25hi, Gy_rad25lo, Gy_rad25hi;
// constants
__m256i const_1_rad25lo, const_1_rad25hi, const_2_rad25lo, const_2_rad25hi, const_3_rad25lo, const_3_rad25hi, const_4_rad25lo, const_4_rad25hi;


void init_curve25519_params_rad25();

#endif  //__radix25_curveParams_H__


