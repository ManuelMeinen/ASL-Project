/* Copyright (c) 2007, 2013 Michele Bini
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

#if GMP_LIMB_BITS == 32
static curve25519key_t p25519 = { 0xffffffed, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0x7fffffff };
static curve25519key_t zerocmp = { 0, 0, 0, 0, 0, 0, 0, 0 };
static curve25519key_t onecmp = { 1, 0, 0, 0, 0, 0, 0, 0 };
static curve25519key_t unsafe[12] =
  {{ 0 },
   { 1 },
   { 0x7C7AEBE0, 0xAEB8413B, 0xFAE35616, 0x6AC49FF1, 0xEB8D09DA, 0xFDB1329C, 0x16056286, 0xB8495F },
   { 0xBC959C5F, 0x248C50A3, 0x55B1D0B1, 0x5BEF839C, 0xC45C4404, 0x868E1C58, 0xDD4E22D8, 0x57119FD0 },
   { 0xFFFFFFEC, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x7FFFFFFF },
   { 0xFFFFFFED, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x7FFFFFFF },
   { 0xFFFFFFEE, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x7FFFFFFF },
   { 0x7C7AEBCD, 0xAEB8413B, 0xFAE35616, 0x6AC49FF1, 0xEB8D09DA, 0xFDB1329C, 0x16056286, 0x80B8495F },
   { 0xBC959C4C, 0x248C50A3, 0x55B1D0B1, 0x5BEF839C, 0xC45C4404, 0x868E1C58, 0xDD4E22D8, 0xD7119FD0 },
   { 0xFFFFFFD9, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF },
   { 0xFFFFFFDA, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF },
   { 0xFFFFFFDB, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF },
  };
#elif GMP_LIMB_BITS == 64
static curve25519key_t p25519 = { 0xffffffffffffffed, 0xffffffffffffffff, 0xffffffffffffffff, 0x7fffffffffffffff };
static curve25519key_t zerocmp = { 0, 0, 0, 0 };
static curve25519key_t onecmp = { 1, 0, 0, 0 };
static curve25519key_t unsafe[12] =
  {{ 0, 0, 0, 0 },
   { 1, 0, 0, 0 },
   { 0xAEB8413B7C7AEBE0, 0x6AC49FF1FAE35616, 0xFDB1329CEB8D09DA, 0xB8495F16056286  },
   { 0x248C50A3BC959C5F, 0x5BEF839C55B1D0B1, 0x868E1C58C45C4404, 0x57119FD0DD4E22D8 },
   { 0xFFFFFFFFFFFFFFEC, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF },
   { 0xFFFFFFFFFFFFFFED, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF },
   { 0xFFFFFFFFFFFFFFEE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF },
   { 0xAEB8413B7C7AEBCD, 0x6AC49FF1FAE35616, 0xFDB1329CEB8D09DA, 0x80B8495F16056286 },
   { 0x248C50A3BC959C4C, 0x5BEF839C55B1D0B1, 0x868E1C58C45C4404, 0xD7119FD0DD4E22D8 },
   { 0xFFFFFFFFFFFFFFD9, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF },
   { 0xFFFFFFFFFFFFFFDA, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF },
   { 0xFFFFFFFFFFFFFFDB, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF },
  };
#else
#error "GMP_LIMBS_BITS not supported for this architecture"
#endif

#define CMP(a, b) mpn_cmp((mp_limb_t*)(a), (mp_limb_t*)(b), C25519N)

extern int
curve25519key_getbit(curve25519key_t *x, unsigned int n) {
  unsigned int d = GMP_LIMB_BITS;
  return (x[0][n / d] >> (n % d)) & 1;
}

extern void
curve25519key_setbit(curve25519key_t *x, unsigned int n, int v) {
  unsigned int d = GMP_LIMB_BITS;
  unsigned int i = n / d;
  mp_limb_t l = x[0][i];
  if (v) {
    l |= ((mp_limb_t)1)<<(n % d);
  } else {
    l &= ~(((mp_limb_t)(1))<<(n % d));
  }
  x[0][i] = l;
}

extern unsigned int
curve25519key_getbyte(curve25519key_t *x, unsigned int n) {
  unsigned int d = GMP_LIMB_BITS;
  n *= 8;
  return (x[0][n / d] >> (n % d)) & 0xff;
}

extern void
curve25519key_setbyte(curve25519key_t *x, unsigned int n, unsigned int v) {
  n *= 8;
  {
    unsigned int d = GMP_LIMB_BITS;
    unsigned int i = n / d;
    mp_limb_t l = x[0][i];
    l = ~l;
    l |= 0xff << (n % d);
    l = ~l;
    l |= v << (n % d);
    x[0][i] = l;
  }
}

#if GMP_LIMB_BITS == 32

extern unsigned int
curve25519key_getuint32(curve25519key_t *x, unsigned int n) {
  return x[0][n];
}

extern void
curve25519key_setuint32(curve25519key_t *x, unsigned int n, unsigned int v) {
  x[0][n] = v;
}

#elif GMP_LIMB_BITS == 64

extern unsigned int
curve25519key_getuint32(curve25519key_t *x, unsigned int n) {
  return (x[0][n>>1] >> ((n&1)*32))&0xffffffff;
}

extern void
curve25519key_setuint32(curve25519key_t *x, unsigned int n, unsigned int v) {
  x[0][n>>1] = (n&1)
  ? (x[0][n>>1] & 0xffffffff) | (((mp_limb_t)v)<<32)
  : (x[0][n>>1] & 0xffffffff00000000) | v;
}

#endif

extern int
curve25519key_validate(curve25519key_t *x) {
  int r;
  if ((r = CMP(x, unsafe + 5)) > 0) {
    if ((r = CMP(x, unsafe + 8)) > 0) {
      if ((r = CMP(x, unsafe + 10)) > 0) {
	if (CMP(x, unsafe + 11) == 0) {
	  return 0;
	}
      } else if (r < 0) {
	if (CMP(x, unsafe + 9) == 0) {
	  return 0;
	}
      } else {
	return 0;
      }
    } else if (r < 0) {
      if ((r = CMP(x, unsafe + 6)) > 0) {
	if (CMP(x, unsafe + 7) == 0) {
	  return 0;
	}
      } else if (r == 0) {
	return 0;
      }
    } else {
      return 0;
    }
  } else if (r < 0) {
    if ((r = CMP(x, unsafe + 2)) > 0) {
      if ((r = CMP(x, unsafe + 3)) > 0) {
	if (CMP(x, unsafe + 4) == 0) {
	  return 0;
	}
      } else if (r == 0) {
	return 0;
      }
    } else if (r < 0) {
      if ((r = CMP(x, unsafe + 1)) < 0) {
	if (CMP(x, unsafe) == 0) {
	  return 0;
	}
      } else if (r == 0) {
	return 0;
      }
    } else {
      return 0;
    }
  } else {
    return 0;
  }
  return 1;
}

#if 1
#include <stdio.h>
#include "base32.h"
static void
tracev(char*m, curve25519key_t *x) {
  char s[(C25519BITS/4)+2];
  base32_encode(s, x);
  fprintf(stderr, "%s: %s\n", m, s);
}
#endif

static void
copykey(curve25519key_t *n, curve25519key_t *x) {
  int c;
  for (c = 0; c < C25519N; c++) {
    n[0][c] = x[0][c];
  }
}

static
int zeromodp(curve25519key_t *x) {
  return (mpn_cmp((mp_limb_t*)x, (mp_limb_t*)&zerocmp, C25519N) == 0);
}

static
void addmodp(curve25519key_t *a, curve25519key_t *b) {
  mpn_add_n((mp_limb_t*)a, (mp_limb_t*)a, (mp_limb_t*)b, C25519N);
  if (mpn_cmp((mp_limb_t*)a, (mp_limb_t*)&p25519, C25519N) >= 0) {
    mpn_sub_n((mp_limb_t*)a, (mp_limb_t*)a, (mp_limb_t*)&p25519, C25519N);
  }
}

static
void submodp(curve25519key_t *a, curve25519key_t *b) {
  if (mpn_cmp((mp_limb_t*)b, (mp_limb_t*)a, C25519N) > 0) {
    mpn_add_n((mp_limb_t*)a, (mp_limb_t*)&p25519, (mp_limb_t*)a, C25519N);
  }
  mpn_sub_n((mp_limb_t*)a, (mp_limb_t*)a, (mp_limb_t*)b, C25519N);
}

static void
mulmodp(curve25519key_t *a, curve25519key_t *b) {
  mp_limb_t d[C25519N*2]; 
  mpn_mul_n(d, (mp_limb_t*)a, (mp_limb_t*)b, C25519N);
  if (0) {
    // unoptimized, this makes the curve25519 function ~150% slower
    mp_limb_t r[C25519N+1];
    mpn_tdiv_qr(r, (mp_limb_t*)a, 0, d, C25519N*2, (mp_limb_t*)&p25519, C25519N);
  } else {
    mp_limb_t r = mpn_addmul_1(d, d+C25519N, C25519N, 19*2);
    r = mpn_add_1((mp_limb_t*)a, d, C25519N, r * (19*2));
    r <<= 1;
#if GMP_LIMB_BITS == 32
    if (((mp_limb_t*)a)[C25519N - 1] & 0x80000000) {
      r |= 1;
      ((mp_limb_t*)a)[C25519N - 1] &= 0x7fffffff;
    }
#elif GMP_LIMB_BITS == 64
    if (((mp_limb_t*)a)[C25519N - 1] & 0x8000000000000000) {
      r |= 1;
      ((mp_limb_t*)a)[C25519N - 1] &= 0x7fffffffffffffff;
    }
#endif
    mpn_add_1((mp_limb_t*)a, (mp_limb_t*)a, C25519N, r * 19);
    if (mpn_cmp((mp_limb_t*)a, (mp_limb_t*)&p25519, C25519N) >= 0) {
      mpn_sub_n((mp_limb_t*)a, (mp_limb_t*)a, (mp_limb_t*)&p25519, C25519N);
    }
  }
}

static void
sqrmodp(curve25519key_t *a) {
  mulmodp(a, a);
}

static void
invmodp(curve25519key_t *a) {
  /* a = a ** (p-2)
     0111 + (1111) x 7 + (1111) x (8*6) + (1111) x 6 + 1110 + 1011
     0 . 1 x (3 + 4*7 + 4 * 8*6 + 4*6 + 3) . 0 . 1011
     0 . 1 x (250) . 0 . 1011
  */
  curve25519key_t c; copykey(&c, a);
  int i = 250;
  while (--i) {
    sqrmodp(a);
    //if (i > 240) { tracev("invmodp a", a); }
    mulmodp(a, &c);
    //if (i > 240) { tracev("invmodp a 2", a); }
  }
  sqrmodp(a);
  sqrmodp(a); mulmodp(a, &c);
  sqrmodp(a);
  sqrmodp(a); mulmodp(a, &c);
  sqrmodp(a); mulmodp(a, &c);
}

static mp_limb_t asmall = 121665; /* (486662 - 2) / 4; */

static
void mulasmall(curve25519key_t *a) {
  if (0) {
    // unoptimized: this makes the function ~5 % slower
    mp_limb_t d[C25519N+1]; mp_limb_t r[2];
    //tracev("mulasmall a", a);
    d[C25519N] = mpn_mul_1(d, (mp_limb_t*)a, C25519N, asmall);
    //tracev("mulasmall d", d);
    //tracev("mulasmall a", a);
    mpn_tdiv_qr(r, (mp_limb_t*)a, 0, d, C25519N+1, (mp_limb_t*)&p25519, C25519N);
    //tracev("mulasmall a", a);
  } else {
    mp_limb_t r = mpn_mul_1((mp_limb_t*)a, (mp_limb_t*)a, C25519N, asmall);
    // Limb size must be at least 32-bits for this to work
    // r = mpn_mul_1((mp_limb_t*)a, (mp_limb_t*)a, C25519N, r*19*2);
    r = mpn_add_1((mp_limb_t*)a, (mp_limb_t*)a, C25519N, r * (19*2));
    r <<= 1;
#if GMP_LIMB_BITS == 32
    if (((mp_limb_t*)a)[C25519N - 1] & 0x80000000) {
      r |= 1;
      ((mp_limb_t*)a)[C25519N - 1] &= 0x7fffffff;
    }
#elif GMP_LIMB_BITS == 64
    if (((mp_limb_t*)a)[C25519N - 1] & 0x8000000000000000) {
      r |= 1;
      ((mp_limb_t*)a)[C25519N - 1] &= 0x7fffffffffffffff;
    }
#endif
    mpn_add_1((mp_limb_t*)a, (mp_limb_t*)a, C25519N, r * 19);
    if (mpn_cmp((mp_limb_t*)a, (mp_limb_t*)&p25519, C25519N) >= 0) {
      mpn_sub_n((mp_limb_t*)a, (mp_limb_t*)a, (mp_limb_t*)&p25519, C25519N);
    }
  }
}

static
void dbl(curve25519key_t *x_2, curve25519key_t *z_2, curve25519key_t *x, curve25519key_t *z) {
  curve25519key_t m, n, o;
  //tracev("dbl x", x);
  //tracev("dbl z", z);
  copykey(&m, x); addmodp(&m, z); sqrmodp(&m);
  //tracev("dbl m", &m);
  copykey(&n, x); submodp(&n, z); sqrmodp(&n);
  //tracev("dbl n", &n);
  copykey(&o, &m); submodp(&o, &n);
  //tracev("dbl o", &o);
  copykey(x_2, &n); mulmodp(x_2, &m);
  //tracev("dbl x_2", x_2);
  copykey(z_2, &o); mulasmall(z_2); addmodp(z_2, &m); mulmodp(z_2, &o);
  //tracev("dbl z_2", z_2);
}

static
void sum(curve25519key_t *x_3, curve25519key_t *z_3, curve25519key_t *x, curve25519key_t *z, curve25519key_t *x_p, curve25519key_t *z_p, curve25519key_t *x_1) {
  curve25519key_t k, l, p, q;
  //tracev("sum x", x);
  //tracev("sum z", z);
  copykey(&p, x); submodp(&p, z); copykey(&k, x_p); addmodp(&k, z_p); mulmodp(&p, &k);
  copykey(&q, x); addmodp(&q, z); copykey(&l, x_p); submodp(&l, z_p); mulmodp(&q, &l);
  //tracev("sum p", &p);
  //tracev("sum q", &q);
  copykey(x_3, &p); addmodp(x_3, &q); sqrmodp(x_3);
  copykey(z_3, &p); submodp(z_3, &q); sqrmodp(z_3); mulmodp(z_3, x_1);
}

extern void
curve25519(curve25519key_t *r, curve25519key_t *f, curve25519key_t *c) {
  curve25519key_t x_1, x_a, z_a, x, z;

  //tracev("f", f);
  if (zeromodp(f)) {
    copykey(r, &zerocmp);
    return;
  }
  copykey(&x_1, c);
  //tracev("c", c);
  //tracev("x_1", x_1);
  dbl(&x_a, &z_a, &x_1, &onecmp);
  //tracev("x_a", &x_a);
  //tracev("z_a", &z_a);
  copykey(&x, &x_1);
  copykey(&z, &onecmp);

  int n = C25519BITS-1;

  while (curve25519key_getbit(f, n) == 0) {
    n--;
  }
  n--;

  while (n >= 0) {
    curve25519key_t nx, nz, nx_a, nz_a;
    int b = curve25519key_getbit(f, n);
    //fprintf(stderr, "b: %d\n", b);
    if (b == 0) {
      dbl(&nx, &nz, &x, &z);
      sum(&nx_a, &nz_a, &x_a, &z_a, &x, &z, &x_1);
    } else {
      sum(&nx, &nz, &x_a, &z_a, &x, &z, &x_1);
      dbl(&nx_a, &nz_a, &x_a, &z_a);
    }
    copykey(&x, &nx); copykey(&z, &nz); copykey(&x_a, &nx_a); copykey(&z_a, &nz_a);
    //tracev("xn", &x);
    //tracev("zn", &z);
    //tracev("x_a", &x_a);
    //tracev("z_a", &z_a);
    n--;
  }

  //tracev("x", &x);
  //tracev("z", &z);
  invmodp(&z);
  //tracev("1/z", &z);
  mulmodp(&x, &z);
  copykey(r, &x);
}
