#include "radix51_curveParams.h"
#include "radix51.h"

radix51_t A_rad51, B_rad51, Gx_rad51, Gy_rad51;
radix51_t const_1_rad51, const_2_rad51, const_3_rad51, const_4_rad51;

void init_curve25519_params_rad51() {
    /*
    init A_rad51, B_rad51, Gx_rad51, Gy_rad51
    */
    A_rad51 .d[0] = 486662;
    B_rad51 .d[0] = 1;
    Gx_rad51.d[0] = 9;
    for (int i = 1; i < 5; i++){
    	A_rad51 .d[i] = 0;
    	B_rad51 .d[i] = 0;
    	Gx_rad51.d[i] = 0;
    }

    Gy_rad51.d[0] = (signed long long)498776679568345;
    Gy_rad51.d[1] = (signed long long)2196478692312381;
    Gy_rad51.d[2] = (signed long long)1446154986911029;
    Gy_rad51.d[3] = (signed long long)943862271874619;
    Gy_rad51.d[4] = (signed long long)574914022836744;
    
    /*
    init radix51_t constant values
    */
    const_1_rad51.d[0] = 1;
    const_2_rad51.d[0] = 2;
    const_3_rad51.d[0] = 3;
    const_4_rad51.d[0] = 4;
    for(int i = 1; i < 5; i++) {
    	const_1_rad51.d[i] = 0;
    	const_2_rad51.d[i] = 0;
    	const_3_rad51.d[i] = 0;
    	const_4_rad51.d[i] = 0;
    }
}