/* Copyright (c) 2007 Michele Bini
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Derived from public domain code (original notice follows)
 *
 * test-curve25519 version 20050915
 * D. J. Bernstein
 * Public domain.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "curve25519.h"

void doit(curve25519key_t *ek,curve25519key_t *e,curve25519key_t *k)
{
  int i;
#if 0
  char s[64];
  base32_encode(s, e); printf("%s ", s);
  base32_encode(s, k); printf("%s ", s);
  curve25519(ek,e,k);
  base32_encode(s, ek); printf("%s\n", s);
#else
  for (i = 0;i < 32;++i) printf("%02x",(unsigned int) curve25519key_getbyte(e,i)); printf(" ");
  for (i = 0;i < 32;++i) printf("%02x",(unsigned int) curve25519key_getbyte(k, i)); printf(" ");
  curve25519(ek,e,k);
  for (i = 0;i < 32;++i) printf("%02x",(unsigned int) curve25519key_getbyte(ek, i)); printf("\n");
#endif
}

curve25519key_t e1k, e2k, e1e2k, e2e1k;
curve25519key_t e1 = {3};
curve25519key_t e2 = {5};
curve25519key_t k = {9};

int
main() {
  int loop;
  int i;

  for (loop = 0;loop < 1000000000;++loop) {
    doit(&e1k,&e1,&k);
    doit(&e2e1k,&e2,&e1k);
    doit(&e2k,&e2,&k);
    doit(&e1e2k,&e1,&e2k);
    for (i = 0;i < C25519N; ++i) {
      if (e1e2k[i] != e2e1k[i]) {
	printf("fail\n");
      }
    }
    for (i = 0;i < C25519N; ++i) e1[i] ^= e2k[i];
    for (i = 0;i < C25519N; ++i) e2[i] ^= e1k[i];
    for (i = 0;i < C25519N; ++i) k[i] ^= e1e2k[i];
  }
  exit(EXIT_SUCCESS);
}
