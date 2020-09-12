#include <gmp.h>
#include "gmp_curveParams.h"

mpz_t A_gmp, B_gmp, Gx_gmp, Gy_gmp, p_gmp;


void init_curve25519_params_gmp()
{
    // curve parameters for curve25519
    mpz_init(A_gmp);
    mpz_init(B_gmp);
    mpz_init(p_gmp);

    mpz_set_ui(A_gmp, 486662); // 486662
    mpz_set_ui(B_gmp, 1);      // 1
    mpz_set_str(p_gmp, "7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed", 16);

    mpz_init(Gx_gmp);
    mpz_set_ui(Gx_gmp, 0x9);

    mpz_init(Gy_gmp);
    mpz_set_str(Gy_gmp, "20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9", 16);
}
