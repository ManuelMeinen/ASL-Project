#ifndef __radix51_curveParams_H__
#include "radix51.h"
#define __radix51_curveParams_H__

// points
radix51_t A_rad51, B_rad51, Gx_rad51, Gy_rad51;
// constants
radix51_t const_1_rad51, const_2_rad51, const_3_rad51, const_4_rad51;

void init_curve25519_params_rad51();

#endif  /* __radix51_curveParams_H__ */
