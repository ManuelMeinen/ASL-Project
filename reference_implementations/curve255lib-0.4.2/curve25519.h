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

#ifndef __CURVE25519_H__
#define __CURVE25519_H__

#include <gmp.h>

#define C25519BITS 256
#define C25519USEDBITS (C25519BITS - 1)
#define C25519N (C25519BITS/GMP_LIMB_BITS)
typedef mp_limb_t curve25519key_t[C25519N];

extern void curve25519(curve25519key_t *r, curve25519key_t *f, curve25519key_t *c);
extern int curve25519key_validate(curve25519key_t *x);
extern int curve25519key_getbit(curve25519key_t *x, unsigned int n);
extern void curve25519key_setbit(curve25519key_t *x, unsigned int n, int v);
extern unsigned int curve25519key_getbyte(curve25519key_t *x, unsigned int n);
extern void curve25519key_setbyte(curve25519key_t *x, unsigned int n, unsigned int v);
extern unsigned int curve25519key_getuint32(curve25519key_t *x, unsigned int n);
extern void curve25519key_setuint32(curve25519key_t *x, unsigned int n, unsigned int v);

#endif /* __CURVE25519_H__ */
