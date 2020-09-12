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
 */


#include <gmp.h>
#include "curve25519.h"
#include "base32.h"

static char base32chars[32] = {
  'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
  'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
  'q', 'r', 's', 't', 'u', 'v', 'w', 'x',
  'y', 'z', '2', '3', '4', '5', '6', '7'
};

extern void
base32_encode(char*s, curve25519key_t *x) {
  int c; int i = 0; int v = 0; int a = 0;
  for (c = 0; c < C25519BITS; c++) {
    v |= curve25519key_getbit(x, c) << i;
    /* fprintf(stderr, "v: %x\n", v); */
    i++;
    if (i > 4) {
      s[a] = base32chars[v];
      v = 0;
      i = 0;
      a++;
    }
  }
  if (v) {
    s[a] = base32chars[v];
    a++;
  }
  s[a] = 0;
  a--;
  i=0;
  while (i < a) {
    v = s[a];
    s[a] = s[i];
    s[i] = v;
    i++; a--;
  }
}

extern void
base32_decode(const char*s, curve25519key_t *x) {
  int c; int i = 0; int v = 0; int a = 0;
  while (s[a]) {
    a++;
  }
  a--;
  while (a >= 0) {
    c = s[a];
    switch (c) {
    case '2': case '3': case '4': case '5': case '6': case '7':
      v = (c - '2')+26;
      break;
    default:
      v = c - 'a';
    }
    curve25519key_setbit(x, i, v&1); v>>=1; i++; if (i >= C25519BITS) { break; };
    curve25519key_setbit(x, i, v&1); v>>=1; i++; if (i >= C25519BITS) { break; };
    curve25519key_setbit(x, i, v&1); v>>=1; i++; if (i >= C25519BITS) { break; };
    curve25519key_setbit(x, i, v&1); v>>=1; i++; if (i >= C25519BITS) { break; };
    curve25519key_setbit(x, i, v&1); v>>=1; i++; if (i >= C25519BITS) { break; };
    a--;
  }
  while (i < C25519BITS) {
    curve25519key_setbit(x, i, 0);
    i++;
  }
}
