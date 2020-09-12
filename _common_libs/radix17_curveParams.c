#include "radix17_curveParams.h"
#include "radix17.h"

radix17_t A_rad17, B_rad17, Gx_rad17, Gy_rad17;
radix17_t const_1_rad17, const_2_rad17, const_3_rad17, const_4_rad17;


void init_curve25519_params_rad17() {
    
    /*
    init A_rad17, B_rad17, Gx_rad17, Gy_rad17
    */
    A_rad17.d[0] = 93446;
    A_rad17.d[1] = 3;
    for(int i = 2; i<15; i++){
        A_rad17.d[i] = 0;
    }

	B_rad17.d[0] = 1;
	for(int i = 1; i < 15; i++){
		B_rad17.d[i] = 0;
	}

	Gx_rad17.d[0] = 9;
	for(int i = 1; i < 15; i++){
		Gx_rad17.d[i] = 0;
	}
    
    signed int arr[] = {54233, 81767, 29032, 17725, 116251, 127851, 62773, 39204, 84177, 38459, 1975, 54940, 100872, 52493, 33464};
    for (int i =0; i<15; i++){
        Gy_rad17.d[i] = arr[i];
    }

    /*
    init radix17_t constant values
    */
    const_1_rad17.d[0] = 1;
    const_2_rad17.d[0] = 2;
    const_3_rad17.d[0] = 3;
    const_4_rad17.d[0] = 4;
    for(int i = 1; i < 15; i++) {
    	const_1_rad17.d[i] = 0;
    	const_2_rad17.d[i] = 0;
    	const_3_rad17.d[i] = 0;
    	const_4_rad17.d[i] = 0;
    }
    

}