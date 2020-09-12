#ifndef __gmp_curveParams_H__
#include <gmp.h>
#define __gmp_curveParams_H__

// points
mpz_t A_gmp, B_gmp, Gx_gmp, Gy_gmp, p_gmp;

void init_curve25519_params_gmp();

#endif  /* __gmp_curveParams_H__ */
