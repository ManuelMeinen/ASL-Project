#ifndef __CURVE25519LIB_BASE32_H__
#define __CURVE25519LIB_BASE32_H__

extern void base32_encode(char*s, curve25519key_t *x);
extern void base32_decode(const char*s, curve25519key_t *x);

#endif /* __CURVE25519LIB_BASE32_H__ */
