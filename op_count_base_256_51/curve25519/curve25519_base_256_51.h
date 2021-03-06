#ifndef _curve25519_base_256_h_
#include "../../_common_libs/radix51.h"
#define _curve25519_base_256_h_

/* general functions (to be called by benchmarking) */
void pub_key_gen(char *sk, radix51_t *x, radix51_t *y);
int validate(char **sk_A, char **sk_B, int num_test_cases);
int checkPk(char **sk_A, char **pk, int num_test_cases);

#endif// _curve25519_base_256_h_
