#include "radix25_curveParams.h"
#include "radix25.h"

radix25_t A_rad25, B_rad25, Gx_rad25, Gy_rad25;
radix25_t const_1_rad25, const_2_rad25, const_3_rad25, const_4_rad25;

void init_curve25519_params_rad25() {
     /*
    init A_rad25, B_rad25, Gx_rad25, Gy_rad25
    */
    A_rad25.d[0] = 486662;
    for(int i = 1; i<10; i++){
        A_rad25.d[i] = 0;
    }

	B_rad25.d[0] = 1;
	for(int i = 1; i < 10; i++){
		B_rad25.d[i] = 0;
	}

	Gx_rad25.d[0] = 9;
	for(int i = 1; i < 10; i++){
		Gx_rad25.d[i] = 0;
	}
    
    signed int arr[] = {47109081, 7432351, 3556669, 32730083, 38335797, 21549388, 57579067, 14064643, 35359240, 8566886};
    for (int i =0; i<10; i++){
        Gy_rad25.d[i] = arr[i];
    }

    /*
    init radix17_t constant values
    */
    const_1_rad25.d[0] = 1;
    const_2_rad25.d[0] = 2;
    const_3_rad25.d[0] = 3;
    const_4_rad25.d[0] = 4;
    for(int i = 1; i < 10; i++) {
    	const_1_rad25.d[i] = 0;
    	const_2_rad25.d[i] = 0;
    	const_3_rad25.d[i] = 0;
    	const_4_rad25.d[i] = 0;
    }
}