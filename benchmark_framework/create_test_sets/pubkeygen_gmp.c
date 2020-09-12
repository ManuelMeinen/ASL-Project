#include <gmp.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

mpz_t A, B, Gx, Gy, p;

void add_mod_p(mpz_t res, const mpz_t a, const mpz_t b)
{
    mpz_add(res, a, b);
    mpz_mod(res, res, p);
}
void sub_mod_p(mpz_t res, const mpz_t a, const mpz_t b)
{
    mpz_sub(res, a, b);
    mpz_mod(res, res, p);
}
void mul_mod_p(mpz_t res, const mpz_t a, const mpz_t b)
{
    mpz_mul(res, a, b);
    mpz_mod(res, res, p);
}
void div_mod_p(mpz_t res, const mpz_t a, const mpz_t b)
{
    // a / b -> a * inv(b)
    mpz_t inv;
    mpz_init(inv);
    mpz_invert(inv, b, p);
    mul_mod_p(res, a, inv);
    mpz_clear(inv);
}

int iszero(mpz_t x, mpz_t y)
{
    return (mpz_sgn(x) == 0 && mpz_sgn(y) == 0);
}

int isequal(mpz_t x1, mpz_t x2, mpz_t y1, mpz_t y2)
{
    return (!mpz_cmp(x1, x2) && !mpz_cmp(y1, y2));
}

void point_double_gmp(mpz_t resX, mpz_t resY, mpz_t x, mpz_t y)
{
    // x <- pointdouble(x), y <- pointdouble(y)

    if (iszero(x, y))
    {
        mpz_set(resX, x);
        mpz_set(resY, y);
        return;
    }

    mpz_t upperfrac, lowerfrac, tmp1, tmp2, l;
    mpz_init(upperfrac);
    mpz_init(lowerfrac);
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(l);

    //x - coordinate

    mul_mod_p(tmp1, A, x);
    mpz_mul_ui(tmp1, tmp1, 2);
    mpz_mod(tmp1, tmp1, p); // tmp1 = 2Ax
    mpz_add_ui(tmp1, tmp1, 1);
    mpz_mod(tmp1, tmp1, p); // tmp1 = 2Ax + 1

    mpz_powm_ui(tmp2, x, 2, p); // tmp2 = x^2
    mpz_mul_ui(tmp2, tmp2, 3);
    mpz_mod(tmp2, tmp2, p); // tmp2 = 3x^2

    add_mod_p(upperfrac, tmp2, tmp1);

    mul_mod_p(lowerfrac, B, y); // By
    mpz_mul_ui(lowerfrac, lowerfrac, 2);
    mpz_mod(lowerfrac, lowerfrac, p); // 2by

    div_mod_p(l, upperfrac, lowerfrac);

    mpz_powm_ui(l, l, 2, p);  // l = l^2
    mul_mod_p(resX, B, l);    //  Bl^2
    sub_mod_p(resX, resX, A); // -A
    sub_mod_p(resX, resX, x); // -x
    sub_mod_p(resX, resX, x); // -x

    //y - coordinate

    mpz_init(upperfrac);
    mpz_init(lowerfrac);
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(l);

    // l = 3x^2 + 2 Ax + 1 / 2By
    mul_mod_p(tmp1, A, x);
    mpz_mul_ui(tmp1, tmp1, 2);
    mpz_mod(tmp1, tmp1, p); // tmp1 = 2Ax
    mpz_add_ui(tmp1, tmp1, 1);
    mpz_mod(tmp1, tmp1, p); // tmp1 = 2Ax + 1

    mpz_powm_ui(tmp2, x, 2, p); // tmp2 = x^2
    mpz_mul_ui(tmp2, tmp2, 3);
    mpz_mod(tmp2, tmp2, p); // tmp2 = 3x^2

    add_mod_p(upperfrac, tmp2, tmp1); // l = 3x^2 + 2Ax + 1

    mul_mod_p(lowerfrac, B, y); // By
    mpz_mul_ui(lowerfrac, lowerfrac, 2);
    mpz_mod(lowerfrac, lowerfrac, p); // 2by

    div_mod_p(l, upperfrac, lowerfrac);

    mpz_mul_ui(tmp1, x, 3);
    mpz_mod(tmp1, tmp1, p);   // 3x
    add_mod_p(tmp1, tmp1, A); // 3x + A

    mul_mod_p(tmp1, tmp1, l); // (3x+A) * l

    mpz_powm_ui(l, l, 3, p); // l^3

    mul_mod_p(l, l, B); // l^3 B

    sub_mod_p(resY, tmp1, l); // (3x+A)* l - Bl^3

    sub_mod_p(resY, resY, y); // -y

    //clear
    mpz_clear(upperfrac);
    mpz_clear(lowerfrac);
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(l);
}

void point_add_gmp(mpz_t resX, mpz_t resY, mpz_t op1x, mpz_t op1y, mpz_t op2x, mpz_t op2y)
{
    // resX <- op1x + op2 x, resY <- op1y + op2 y
    // https://en.wikipedia.org/wiki/Montgomery_curve#Addition

    if (!mpz_cmp(op1x, op2x) && mpz_cmp(op1y, op2y) != 0)
    {
        mpz_set_ui(resX, 0);
        mpz_set_ui(resY, 0);
        return;
    }

    if (iszero(op1x, op1y))
    { // probably not side channel safe anymore ?
        mpz_set(resX, op2x);
        mpz_set(resY, op2y);
        return;
    }
    if (iszero(op2x, op2y))
    { // point_add(0,P) = P = point_add(P,0)
        mpz_set(resX, op1x);
        mpz_set(resY, op1y);
        return;
    }
    if (isequal(op1x, op2x, op1y, op2y))
    { // point_add(p1, p1) = double(p1), otherwise exception div by 0
        point_double_gmp(resX, resY, op1x, op1y);
        return;
    }

    mpz_t diffx, diffy, tmp1, tmp2;
    mpz_init(diffx);
    mpz_init(diffy);
    mpz_init(tmp1);
    mpz_init(tmp2);

    // x - coordinate

    sub_mod_p(diffy, op2y, op1y); // y2 - y1
    sub_mod_p(diffx, op2x, op1x); // x2 - x1

    mpz_powm_ui(diffy, diffy, 2, p); // diffy = (y2 - y1) ^2

    mpz_powm_ui(diffx, diffx, 2, p); // diffx = (x2 - x1) ^2
    mul_mod_p(diffy, diffy, B);      //  diffx = B * (y2 - y1) ^2

    div_mod_p(resX, diffy, diffx); // B * (y2 - y1) ^2 / (x2 - x1) ^2

    sub_mod_p(resX, resX, A);    // - A
    sub_mod_p(resX, resX, op1x); // - x1
    sub_mod_p(resX, resX, op2x); // - x2

    // y - coordinate

    mpz_init(diffx);
    mpz_init(diffy);

    add_mod_p(tmp1, op1x, op2x); // x1 + x2
    add_mod_p(tmp1, tmp1, op1x); // + x1
    add_mod_p(tmp1, tmp1, A);    // + A  // tmp1 =  2x1 + x2 + A

    sub_mod_p(diffy, op2y, op1y); // y2 - y1

    mul_mod_p(tmp1, tmp1, diffy); // (2x1 + x2 + A) (y2-y1)

    sub_mod_p(diffx, op2x, op1x); // x2-x1

    div_mod_p(tmp1, tmp1, diffx); // tmp1 = (2x1 + x2 + A) (y2-y1) / (x2-x1)

    mpz_powm_ui(diffx, diffx, 3, p);
    mpz_powm_ui(diffy, diffy, 3, p);
    mul_mod_p(diffy, diffy, B);

    div_mod_p(tmp2, diffy, diffx);
    sub_mod_p(resY, tmp1, tmp2);
    sub_mod_p(resY, resY, op1y);

    //clear
    mpz_clear(diffx);
    mpz_clear(diffy);
    mpz_clear(tmp1);
    mpz_clear(tmp2);
}

void scalar_mult_gmp(mpz_t base_point_x, mpz_t base_point_y, char *sk, mpz_t Pkx, mpz_t Pky)
{
    int n = strlen(sk);

    mpz_t R0x, R0y, R1x, R1y, resX, resY;
    mpz_init(R0x);
    mpz_init(R0y);
    mpz_init(R1x);
    mpz_init(R1y);
    mpz_init(resX);
    mpz_init(resY);

    mpz_set(R1x, base_point_x); // R1 <- G
    mpz_set(R1y, base_point_y);

    //montgomery ladder
    // https://en.wikipedia.org/wiki/Elliptic_curve_point_multiplication#Montgomery_ladder

    for (int i = 0; i < n; i++)
    {
        mpz_init(resX);
        mpz_init(resY);

        if (sk[i] == '0')
        {
            //R1 <- add (R0, R1)

            point_add_gmp(resX, resY, R0x, R0y, R1x, R1y);
            mpz_set(R1x, resX);
            mpz_set(R1y, resY);

            mpz_init(resX);
            mpz_init(resY);

            // R0 <- double(R0)

            point_double_gmp(resX, resY, R0x, R0y);
            mpz_set(R0x, resX);
            mpz_set(R0y, resY);
        }
        else
        {
            point_add_gmp(resX, resY, R0x, R0y, R1x, R1y);

            mpz_set(R0x, resX);
            mpz_set(R0y, resY);

            mpz_init(resX);
            mpz_init(resY);

            point_double_gmp(resX, resY, R1x, R1y);
            mpz_set(R1x, resX);
            mpz_set(R1y, resY);
        }
    }

    // result: R0 (public key)
    mpz_set(Pkx, R0x);
    mpz_set(Pky, R0y);
}

void gmp_public_key_gen(char *sk)
{
	/*
    //base point (x,y)
    mpz_init(Gx);
    mpz_set_ui(Gx, 0x9);

    mpz_init(Gy);
    mpz_set_str(Gy, "20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9", 16);
    */

    //public key
    mpz_t Pkx, Pky;
    mpz_init(Pkx);
    mpz_init(Pky);

    scalar_mult_gmp(Gx, Gy, sk, Pkx, Pky);
    gmp_printf("%Zd\n", Pkx);
}

int validate_pubkeygen_gmp( char* skA, char* skB ) {

	// 1. generate Alice's public key
	mpz_t Pkx, Pky;
	mpz_init(Pkx);
	mpz_init(Pky);
	
	scalar_mult_gmp(Gx, Gy, skA, Pkx, Pky);
	
	// 2. generate Bob's public key
	// Bob's public key
	mpz_t Pk2x, Pk2y;
	mpz_init(Pk2x);
	mpz_init(Pk2y);
	
	scalar_mult_gmp(Gx, Gy, skB, Pk2x, Pk2y);
	
	// two shared secret instances
	mpz_t shared_secret_alice, shared_secret_bob, shared_secret_alice_y, shared_secret_bob_y;
	mpz_init(shared_secret_alice);
	mpz_init(shared_secret_bob);
	mpz_init(shared_secret_alice_y);
	mpz_init(shared_secret_bob_y);
	
	// check correctness
	scalar_mult_gmp(Pk2x, Pk2y, skA, shared_secret_alice, shared_secret_alice_y);
	scalar_mult_gmp(Pkx , Pky , skB, shared_secret_bob  , shared_secret_bob_y);
	
	if ( mpz_cmp( shared_secret_alice, shared_secret_bob ) ) {
		return 0;
	}
	else {
		return 1;
	}
}

void init_pubkeygen_gmp()
{
    // curve parameters for curve25519
    mpz_init(A);
    mpz_init(B);
    mpz_init(p);

    mpz_set_ui(A, 486662); // 486662
    mpz_set_ui(B, 1);      // 1
    mpz_set_str(p, "7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed", 16);    

	mpz_init(Gx);
	mpz_set_ui(Gx, 0x9);
	mpz_init(Gy);
	mpz_set_str(Gy, "20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9", 16);
}