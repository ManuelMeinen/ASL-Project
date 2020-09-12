#ifndef __CURVE25519_H__
#include <gmp.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#define __CURVE25519_H__

void gmp_public_key_gen(char *sk);
void init_pubkeygen_gmp();
int validate_pubkeygen_gmp( char* skA, char* skB );

#endif /* __CURVE25519_H__ */