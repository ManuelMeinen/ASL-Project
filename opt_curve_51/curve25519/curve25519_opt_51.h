#ifndef _curve25519_base_256_51_h_
#include "../../_common_libs/radix51.h"
#define _curve25519_base_256_51_h_

/* general functions (to be called by benchmarking) */
void scalar_mult_base_256_51(radix51_t *base_point_x, radix51_t *base_point_y, char *sk, radix51_t *Pkx, radix51_t *Pky);
void public_key_gen_base_256_51(char *sk, radix51_t *x, radix51_t *y);
int validate_base_256_51(char **sk_A, char **sk_B, int num_test_cases);
int checkPk_base_256_51(char **sk, char **pk, int num_test_cases);

#endif // _curve25519_base_256_51_h_