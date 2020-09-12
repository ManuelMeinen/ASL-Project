#ifndef _curve25519_base_256_25_h_
#include "../../_common_libs/radix25.h"
#define _curve25519_base_256_25_h_

/*
curve25519_t A, B, Gx, Gy, p;
//void add_mod_p(curve25519_t res, curve25519_t a, curve25519_t b);
//void sub_mod_p(curve25519_t res, const curve25519_t a, const curve25519_t b);
//void mul_mod_p(curve25519_t res, const curve25519_t a, const curve25519_t b);
//void div_mod_p(curve25519_t res, const curve25519_t a, const curve25519_t b);
void print_mpz(char name[], curve25519_t var);
int iszero(curve25519_t x, curve25519_t y);
int isequal(curve25519_t x1, curve25519_t x2, curve25519_t y1, curve25519_t y2);
void point_double(curve25519_t *resX, curve25519_t *resY, curve25519_t *x, curve25519_t *y);
void point_add(curve25519_t *resX, curve25519_t *resY, curve25519_t *op1x, curve25519_t *op1y, curve25519_t *op2x, curve25519_t *op2y);
*/

/* general functions (to be called by benchmarking) */
void scalar_mult_base_256_25(radix25_t *base_point_x, radix25_t *base_point_y, char *sk, radix25_t *Pkx, radix25_t *Pky);
void public_key_gen_base_256_25(char *sk, radix25_t *x, radix25_t *y);
int validate_base_256_25(char **sk_A, char **sk_B, int num_test_cases);
int checkPk_base_256_25(char **sk, char **pk, int num_test_cases);

#endif // _curve25519_base_256_25_h_