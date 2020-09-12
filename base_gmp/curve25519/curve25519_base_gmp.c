#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include "../../_common_libs/gmp_curveParams.h"
#include "../../_common_libs/common.h"

void add_mod_p(mpz_t res, const mpz_t a, const mpz_t b)
{
    mpz_add(res, a, b);
    mpz_mod(res, res, p_gmp);

}
void sub_mod_p(mpz_t res, const mpz_t a, const mpz_t b)
{
    mpz_sub(res, a, b);
    mpz_mod(res, res, p_gmp);
}
void mul_mod_p(mpz_t res, const mpz_t a, const mpz_t b)
{
    mpz_mul(res, a, b);
    mpz_mod(res, res, p_gmp);
}
void div_mod_p(mpz_t res, const mpz_t a, const mpz_t b)
{
    // a / b -> a * inv(b)
    mpz_t inv;
    mpz_init(inv);
    mpz_invert(inv, b, p_gmp);
    mul_mod_p(res, a, inv);
    mpz_clear(inv);
}

void print_mpz(char name[], mpz_t var)
{
    gmp_printf("%s: %Zd\n", name, var);
}

int iszero(mpz_t x, mpz_t y)
{
    return (mpz_sgn(x) == 0 && mpz_sgn(y) == 0);
}

int isequal(mpz_t x1, mpz_t x2, mpz_t y1, mpz_t y2)
{
    return (!mpz_cmp(x1, x2) && !mpz_cmp(y1, y2));
}

void point_double_base_gmp(mpz_t resX, mpz_t resY, mpz_t x, mpz_t y)
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

    mul_mod_p(tmp1, A_gmp, x);
    mpz_mul_ui(tmp1, tmp1, 2);
    mpz_mod(tmp1, tmp1, p_gmp); // tmp1 = 2Ax
    mpz_add_ui(tmp1, tmp1, 1);
    mpz_mod(tmp1, tmp1, p_gmp); // tmp1 = 2Ax + 1

    mpz_powm_ui(tmp2, x, 2, p_gmp); // tmp2 = x^2
    mpz_mul_ui(tmp2, tmp2, 3);
    mpz_mod(tmp2, tmp2, p_gmp); // tmp2 = 3x^2

    add_mod_p(upperfrac, tmp2, tmp1);

    mul_mod_p(lowerfrac, B_gmp, y); // By
    mpz_mul_ui(lowerfrac, lowerfrac, 2);
    mpz_mod(lowerfrac, lowerfrac, p_gmp); // 2by

    div_mod_p(l, upperfrac, lowerfrac);

    mpz_powm_ui(l, l, 2, p_gmp);  // l = l^2
    mul_mod_p(resX, B_gmp, l);    //  Bl^2
    sub_mod_p(resX, resX, A_gmp); // -A
    sub_mod_p(resX, resX, x); // -x
    sub_mod_p(resX, resX, x); // -x

    //y - coordinate

    mpz_init(upperfrac);
    mpz_init(lowerfrac);
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(l);

    // l = 3x^2 + 2 Ax + 1 / 2By
    mul_mod_p(tmp1, A_gmp, x);
    mpz_mul_ui(tmp1, tmp1, 2);
    mpz_mod(tmp1, tmp1, p_gmp); // tmp1 = 2Ax
    mpz_add_ui(tmp1, tmp1, 1);
    mpz_mod(tmp1, tmp1, p_gmp); // tmp1 = 2Ax + 1

    mpz_powm_ui(tmp2, x, 2, p_gmp); // tmp2 = x^2
    mpz_mul_ui(tmp2, tmp2, 3);
    mpz_mod(tmp2, tmp2, p_gmp); // tmp2 = 3x^2

    add_mod_p(upperfrac, tmp2, tmp1); // l = 3x^2 + 2Ax + 1

    mul_mod_p(lowerfrac, B_gmp, y); // By
    mpz_mul_ui(lowerfrac, lowerfrac, 2);
    mpz_mod(lowerfrac, lowerfrac, p_gmp); // 2by

    div_mod_p(l, upperfrac, lowerfrac);

    mpz_mul_ui(tmp1, x, 3);
    mpz_mod(tmp1, tmp1, p_gmp);   // 3x
    add_mod_p(tmp1, tmp1, A_gmp); // 3x + A

    mul_mod_p(tmp1, tmp1, l); // (3x+A) * l

    mpz_powm_ui(l, l, 3, p_gmp); // l^3

    mul_mod_p(l, l, B_gmp); // l^3 B

    sub_mod_p(resY, tmp1, l); // (3x+A)* l - Bl^3

    sub_mod_p(resY, resY, y); // -y

    //clear
    mpz_clear(upperfrac);
    mpz_clear(lowerfrac);
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(l);
}

void point_add_base_gmp(mpz_t resX, mpz_t resY, mpz_t op1x, mpz_t op1y, mpz_t op2x, mpz_t op2y)
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
    { // point_add_base_gmp(0,P) = P = point_add_base_gmp(P,0)
        mpz_set(resX, op1x);
        mpz_set(resY, op1y);
        return;
    }
    if (isequal(op1x, op2x, op1y, op2y))
    { // point_add_base_gmp(p1, p1) = double(p1), otherwise exception div by 0
        point_double_base_gmp(resX, resY, op1x, op1y);
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

    mpz_powm_ui(diffy, diffy, 2, p_gmp); // diffy = (y2 - y1) ^2

    mpz_powm_ui(diffx, diffx, 2, p_gmp); // diffx = (x2 - x1) ^2
    mul_mod_p(diffy, diffy, B_gmp);      //  diffx = B * (y2 - y1) ^2

    div_mod_p(resX, diffy, diffx); // B * (y2 - y1) ^2 / (x2 - x1) ^2

    sub_mod_p(resX, resX, A_gmp);    // - A
    sub_mod_p(resX, resX, op1x); // - x1
    sub_mod_p(resX, resX, op2x); // - x2

    // y - coordinate

    mpz_init(diffx);
    mpz_init(diffy);

    add_mod_p(tmp1, op1x, op2x); // x1 + x2
    add_mod_p(tmp1, tmp1, op1x); // + x1
    add_mod_p(tmp1, tmp1, A_gmp);    // + A  // tmp1 =  2x1 + x2 + A

    sub_mod_p(diffy, op2y, op1y); // y2 - y1

    mul_mod_p(tmp1, tmp1, diffy); // (2x1 + x2 + A) (y2-y1)

    sub_mod_p(diffx, op2x, op1x); // x2-x1

    div_mod_p(tmp1, tmp1, diffx); // tmp1 = (2x1 + x2 + A) (y2-y1) / (x2-x1)

    mpz_powm_ui(diffx, diffx, 3, p_gmp);
    mpz_powm_ui(diffy, diffy, 3, p_gmp);
    mul_mod_p(diffy, diffy, B_gmp);

    div_mod_p(tmp2, diffy, diffx);
    sub_mod_p(resY, tmp1, tmp2);
    sub_mod_p(resY, resY, op1y);

    //clear
    mpz_clear(diffx);
    mpz_clear(diffy);
    mpz_clear(tmp1);
    mpz_clear(tmp2);
}

void scalar_mult_base_gmp(mpz_t base_point_x, mpz_t base_point_y, char *sk, mpz_t Pkx, mpz_t Pky)
{
    /*
    Pk = sk*base_point
    */
    int n = sk_length;

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

            point_add_base_gmp(resX, resY, R0x, R0y, R1x, R1y);
            mpz_set(R1x, resX);
            mpz_set(R1y, resY);

            mpz_init(resX);
            mpz_init(resY);

            // R0 <- double(R0)

            point_double_base_gmp(resX, resY, R0x, R0y);
            mpz_set(R0x, resX);
            mpz_set(R0y, resY);
        }
        else
        {
            point_add_base_gmp(resX, resY, R0x, R0y, R1x, R1y);

            mpz_set(R0x, resX);
            mpz_set(R0y, resY);

            mpz_init(resX);
            mpz_init(resY);

            point_double_base_gmp(resX, resY, R1x, R1y);
            mpz_set(R1x, resX);
            mpz_set(R1y, resY);
        }
    }

    // result: R0 (public key)
    mpz_set(Pkx, R0x);
    mpz_set(Pky, R0y);
}

void public_key_gen_base_gmp(char *sk)
{
    /*
    (x,y) = sk*G
    */
    //public key
    mpz_t Pkx, Pky;
    mpz_init(Pkx);
    mpz_init(Pky);

    scalar_mult_base_gmp(Gx_gmp, Gy_gmp, sk, Pkx, Pky);
}

int validate_base_gmp( char** sk_A, char** sk_B, int num_test_cases )
{
    /*
    validate by running an instance of DH
    */
	mpz_t Pkx, Pky;
	mpz_init(Pkx);
	mpz_init(Pky);

	mpz_t Pk2x, Pk2y;
	mpz_init(Pk2x);
	mpz_init(Pk2y);
	
	mpz_t shared_secret_alice, shared_secret_bob, shared_secret_alice_y, shared_secret_bob_y;
	mpz_init(shared_secret_alice);
	mpz_init(shared_secret_bob);
	mpz_init(shared_secret_alice_y);
	mpz_init(shared_secret_bob_y);
	
	for ( int i = 0; i < num_test_cases; ++i ) {
    	// 1. generate Alice's public key
	    scalar_mult_base_gmp(Gx_gmp, Gy_gmp, sk_A[i], Pkx, Pky);

    	// 2. generate Bob's public key
	    scalar_mult_base_gmp(Gx_gmp, Gy_gmp, sk_B[i], Pk2x, Pk2y);
		
	    // check correctness
	    scalar_mult_base_gmp(Pk2x, Pk2y, sk_A[i], shared_secret_alice, shared_secret_alice_y);
	    scalar_mult_base_gmp(Pkx , Pky , sk_B[i], shared_secret_bob  , shared_secret_bob_y);
	
	    if ( mpz_cmp( shared_secret_alice, shared_secret_bob ) ) {
	    	return 0;
	    }
	}
	
	return 1;
}

int checkPk_base_gmp(char **sk, char **pk, int num_test_cases) {
	// return dummy result as this computation yields the ground truth
	return 1;
}

int test_validate()
{
    printf("VALIDATION USING ECDH\n");

	// base point
	mpz_init(Gx_gmp);
	mpz_set_ui(Gx_gmp, 0x9);
	
	mpz_init(Gy_gmp);
	mpz_set_str(Gy_gmp, "20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9", 16);
	    
	// Alices's secret key
	char skA[257] = "0110100000111100010111001110100101110000011001011011010011111000100110000111011100010001001011110110110001100100010110010000101000101101001011000000110010101011000110100111101010001110000010100001010000011100100110001011000110010100110100101100111100010110"; //change this
	
	//assert(strlen(sk)==256);
	printf("\nAlice's secret key: \n%s\n\n", skA);
		
	// Alice's public key
	mpz_t Pkx, Pky;
	mpz_init(Pkx);
	mpz_init(Pky);
	
	scalar_mult_base_gmp(Gx_gmp, Gy_gmp, skA, Pkx, Pky);
	
	printf("Alice's public key: \n");
	print_mpz("Pk(x)", Pkx);
	print_mpz("Pk(y)", Pky);
	
	// 2. generate Bob's public key
	// base point
	mpz_init(Gx_gmp);
	mpz_set_ui(Gx_gmp, 0x9);
	
	mpz_init(Gy_gmp);
	mpz_set_str(Gy_gmp, "20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9", 16);
	
	// Bob's secret key
	char skB[257] = "0000000001100010110001011101101101111001001010000100011010111010110001100000011100001111100111011001000100010111000000110000100001010011010101000011011111100010000110011001011011111100011110010000010011010010001001101000111111110011110001100011110100111011";
	printf("\nBob's secret key: \n%s\n\n", skB);

	// Bob's public key
	mpz_t Pk2x, Pk2y;
	mpz_init(Pk2x);
	mpz_init(Pk2y);
	
	scalar_mult_base_gmp(Gx_gmp, Gy_gmp, skB, Pk2x, Pk2y);
	
	printf("Bob's public key: \n");
	print_mpz("Pk2(x)", Pk2x);
	print_mpz("Pk2(y)", Pk2y);
	
	// two shared secret instances
	mpz_t shared_secret_alice, shared_secret_bob, shared_secret_alice_y, shared_secret_bob_y;
	mpz_init(shared_secret_alice);
	mpz_init(shared_secret_bob);
	mpz_init(shared_secret_alice_y);
	mpz_init(shared_secret_bob_y);
	
	// check correctness
	scalar_mult_base_gmp(Pk2x, Pk2y, skA, shared_secret_alice, shared_secret_alice_y);
	scalar_mult_base_gmp(Pkx , Pky , skB, shared_secret_bob  , shared_secret_bob_y);
	
	printf("\nShared Secrets: \n");
	print_mpz("Alice (skA * Pk2) [x-coordinate]", shared_secret_alice);
	print_mpz("Bob   (skB * Pk ) [x-coordinate]", shared_secret_bob);
	printf("\n");

    if (!mpz_cmp(shared_secret_alice, shared_secret_bob))
    {
        printf("CORRECT!\n");
    }
    else
    {
        printf("FALSE!\n");
    }
    return 1;
}
