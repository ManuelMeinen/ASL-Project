#ifndef _curve25519_base_256_17_h_
#include "../../_common_libs/radix17.h"
#define _curve25519_base_256_17_h_

/* general functions (to be called by benchmarking) */
void scalar_mult_base_256_17(radix17_t *base_point_x, radix17_t *base_point_y, char *sk, radix17_t *Pkx, radix17_t *Pky);
void public_key_gen_base_256_17(char *sk, radix17_t *x, radix17_t *y);
int validate_base_256_17(char **sk_A, char **sk_B, int num_test_cases);
int checkPk_base_256_17(char **sk, char **pk, int num_test_cases);
void point_double_base_256_17(radix17_t *resX, radix17_t *resY, radix17_t *x, radix17_t *y);
void point_add_base_256_17(radix17_t *resX, radix17_t *resY, radix17_t *op1x, radix17_t *op1y, radix17_t *op2x, radix17_t *op2y);
#endif // _curve25519_base_256_17_h_