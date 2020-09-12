#include "curve25519_base_gmp.h"
#include "../../_common_libs/gmp_curveParams.h"

int main(int argc, char const *argv[])
{
    init_curve25519_params_gmp();
    char sk[257] = "0110100000111100010111001110100101110000011001011011010011111000100110000111011100010001001011110110110001100100010110010000101000101101001011000000110010101011000110100111101010001110000010100001010000011100100110001011000110010100110100101100111100010110"; //change this

    public_key_gen_base_gmp(sk);
    
    test_validate();

    return 0;
}