#ifndef __radix17_curveParams_vec_H__
#include "radix17_vec.h"
#define __radix17_curveParams_vec_H__

// points
__m256i A_rad17lo, A_rad17hi, B_rad17lo, B_rad17hi, Gx_rad17lo, Gx_rad17hi, Gy_rad17lo, Gy_rad17hi;
// constants
__m256i const_1_rad17lo, const_1_rad17hi, const_2_rad17lo, const_2_rad17hi, const_3_rad17lo, const_3_rad17hi, const_4_rad17lo, const_4_rad17hi;


void init_curve25519_params_rad17();

#endif  //__radix17_curveParams_H__


