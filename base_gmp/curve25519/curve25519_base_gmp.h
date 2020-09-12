#ifndef __CURVE25519_BASE_GMP_H__
#include <gmp.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#define __CURVE25519_BASE_GMP_H__

void scalar_mult_base_gmp(mpz_t base_point_x, mpz_t base_point_y, char *sk, mpz_t Pkx, mpz_t Pky);
void public_key_gen_base_gmp(char *sk);
int validate_base_gmp( char** sk_A, char** sk_B, int num_test_cases );
int checkPk_base_gmp(char **sk, char **pk, int num_test_cases);
int test_validate();

#endif /* __CURVE25519_BASE_GMP_H__ */