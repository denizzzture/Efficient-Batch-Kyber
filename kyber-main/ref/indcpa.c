#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "symmetric.h"
#include "randombytes.h"
#include <stdio.h>

void poly_sum(poly *result, poly *poly_array[], unsigned int n) {
    unsigned int i, j;
    
    for (i = 0; i < KYBER_N; ++i) {
        result->coeffs[i] = 0; 
        for (j = 0; j < n; ++j) {
            result->coeffs[i] += poly_array[j]->coeffs[i]; 
        }
    }
	 
}

int mult_p4_str(polyvec C[BATCH_SIZE], polyvec mat[KYBER_K], polyvec s1hats[BATCH_SIZE]){

  poly addition, addition2, addition3, addition4, subtraction, subtraction2;
  poly p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,p33,p34,p35,p36,p37,p38,p39,p40,p41,p42,p43,p44,p45,p46,p47,p48,p49;


	// Computation of p_i's
	
	poly *poly_array[] = {&mat[0].vec[0], &mat[2].vec[2], &mat[1].vec[1], &mat[3].vec[3]};
	poly_sum(&addition, poly_array, 4);
	poly *poly_array2[] = {&s1hats[0].vec[0], &s1hats[2].vec[2], &s1hats[1].vec[1], &s1hats[3].vec[3]};
	poly_sum(&addition2, poly_array2, 4);
	poly_basemul_montgomery(&p1,  &addition, &addition2);
	
	poly *poly_array3[] = {&mat[1].vec[0], &mat[3].vec[2], &mat[1].vec[1], &mat[3].vec[3]};
	poly_sum(&addition, poly_array3, 4);
	poly_add(&addition2, &s1hats[0].vec[0], &s1hats[2].vec[2]);
	poly_basemul_montgomery(&p2,  &addition, &addition2);
	
	poly_add(&addition, &mat[0].vec[0], &mat[2].vec[2]);
	poly_add(&addition2, &s1hats[1].vec[0], &s1hats[3].vec[2]);
	poly_add(&addition3, &s1hats[1].vec[1], &s1hats[3].vec[3]);
	poly_sub(&subtraction,  &addition2, &addition3);
	poly_basemul_montgomery(&p3, &addition, &subtraction);
	
	poly_add(&addition, &mat[1].vec[1], &mat[3].vec[3]);
	poly_add(&addition2, &s1hats[0].vec[1], &s1hats[2].vec[3]);
	poly_add(&addition3, &s1hats[0].vec[0], &s1hats[2].vec[2]);
	poly_sub(&subtraction,  &addition2, &addition3);
	poly_basemul_montgomery(&p4, &addition, &subtraction);
	
	poly *poly_array4[] = {&mat[0].vec[0], &mat[2].vec[2], &mat[0].vec[1], &mat[2].vec[3]};
	poly_sum(&addition, poly_array4, 4);
	poly_add(&addition2, &s1hats[1].vec[1], &s1hats[3].vec[3]);
	poly_basemul_montgomery(&p5,  &addition, &addition2);
	
	poly_add(&addition, &mat[1].vec[0], &mat[3].vec[2]);
	poly_add(&addition2, &mat[0].vec[0], &mat[2].vec[2]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array5[] = {&s1hats[0].vec[0], &s1hats[2].vec[2], &s1hats[1].vec[0], &s1hats[3].vec[2]};
	poly_sum(&addition3, poly_array5, 4);
	poly_basemul_montgomery(&p6,  &addition3, &subtraction);
	
	poly_add(&addition, &mat[0].vec[1], &mat[2].vec[3]);
	poly_add(&addition2, &mat[1].vec[1], &mat[3].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array6[] = {&s1hats[0].vec[1], &s1hats[2].vec[3], &s1hats[1].vec[1], &s1hats[3].vec[3]};
	poly_sum(&addition3, poly_array6, 4);
	poly_basemul_montgomery(&p7,  &addition3, &subtraction);
	
	poly *poly_array7[] = {&mat[2].vec[0], &mat[2].vec[2], &mat[3].vec[1], &mat[3].vec[3]};
	poly_sum(&addition, poly_array7, 4);
	poly_add(&addition2, &s1hats[0].vec[0], &s1hats[1].vec[1]);
	poly_basemul_montgomery(&p8,  &addition, &addition2);
	
	poly *poly_array8[] = {&mat[3].vec[0], &mat[3].vec[2], &mat[3].vec[1], &mat[3].vec[3]};
	poly_sum(&addition, poly_array8, 4);
	poly_basemul_montgomery(&p9,  &addition, &s1hats[0].vec[0]);
	
	poly_add(&addition, &mat[2].vec[0], &mat[2].vec[2]);
	poly_sub(&subtraction,  &s1hats[1].vec[0], &s1hats[1].vec[1]);
	poly_basemul_montgomery(&p10,  &addition, &subtraction);
	
	poly_add(&addition, &mat[3].vec[1], &mat[3].vec[3]);
	poly_sub(&subtraction,  &s1hats[0].vec[1], &s1hats[0].vec[0]);
	poly_basemul_montgomery(&p11,  &addition, &subtraction);
	
	poly *poly_array9[] = {&mat[2].vec[0], &mat[2].vec[2], &mat[2].vec[1], &mat[2].vec[3]};
	poly_sum(&addition, poly_array9, 4);
	poly_basemul_montgomery(&p12,  &addition, &s1hats[1].vec[1]);
	
	poly_add(&addition, &mat[3].vec[0], &mat[3].vec[2]);
	poly_add(&addition2, &mat[2].vec[0], &mat[2].vec[2]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[0].vec[0], &s1hats[1].vec[0]);
	poly_basemul_montgomery(&p13, &subtraction, &addition3);
	
	poly_add(&addition, &mat[2].vec[1], &mat[2].vec[3]);
	poly_add(&addition2, &mat[3].vec[1], &mat[3].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[0].vec[1], &s1hats[1].vec[1]);
	poly_basemul_montgomery(&p14, &subtraction, &addition3);
	
	poly_add(&addition, &mat[0].vec[0], &mat[1].vec[1]);
	poly_add(&addition2, &s1hats[2].vec[0], &s1hats[3].vec[1]);
	poly_add(&addition3, &s1hats[2].vec[2], &s1hats[3].vec[3]);
	poly_sub(&subtraction,  &addition2, &addition3);
	poly_basemul_montgomery(&p15, &addition, &subtraction);
	
	poly_add(&addition, &mat[1].vec[0], &mat[1].vec[1]);
	poly_sub(&subtraction,  &s1hats[2].vec[0], &s1hats[2].vec[2]);
	poly_basemul_montgomery(&p16,  &addition, &subtraction);
	
	poly_add(&addition, &s1hats[3].vec[0], &s1hats[3].vec[3]);
	poly_add(&addition2, &s1hats[3].vec[2], &s1hats[3].vec[1]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_basemul_montgomery(&p17,  &mat[0].vec[0], &subtraction);
	
	poly_add(&addition, &s1hats[2].vec[1], &s1hats[2].vec[2]);
	poly_add(&addition2, &s1hats[2].vec[3], &s1hats[2].vec[0]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_basemul_montgomery(&p18,  &mat[1].vec[1], &subtraction);
	
	poly_add(&addition, &mat[0].vec[0], &mat[0].vec[1]);
	poly_sub(&subtraction,  &s1hats[3].vec[1], &s1hats[3].vec[3]);
	poly_basemul_montgomery(&p19,  &addition, &subtraction);
	
	poly_sub(&subtraction,  &mat[1].vec[0], &mat[0].vec[0]);
	poly_add(&addition, &s1hats[2].vec[0], &s1hats[3].vec[0]);
	poly_add(&addition2, &s1hats[2].vec[2], &s1hats[3].vec[2]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_basemul_montgomery(&p20,  &subtraction, &subtraction2);
	
	poly_sub(&subtraction,  &mat[0].vec[1], &mat[1].vec[1]);
	poly_add(&addition, &s1hats[2].vec[1], &s1hats[3].vec[1]);
	poly_add(&addition2, &s1hats[2].vec[3], &s1hats[3].vec[3]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_basemul_montgomery(&p21,  &subtraction, &subtraction2);
	
	poly_add(&addition,  &mat[2].vec[2], &mat[3].vec[3]);
	poly_add(&addition2, &s1hats[0].vec[2], &s1hats[1].vec[3]);
	poly_add(&addition3, &s1hats[0].vec[0], &s1hats[1].vec[1]);
	poly_sub(&subtraction,  &addition2, &addition3);
	poly_basemul_montgomery(&p22,  &addition, &subtraction);
	
	poly_add(&addition, &mat[3].vec[2], &mat[3].vec[3]);
	poly_sub(&subtraction,  &s1hats[0].vec[2], &s1hats[0].vec[0]);
	poly_basemul_montgomery(&p23,  &addition, &subtraction);
	
	poly_add(&addition, &s1hats[1].vec[2], &s1hats[1].vec[1]);
	poly_add(&addition2, &s1hats[1].vec[0], &s1hats[1].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_basemul_montgomery(&p24,  &mat[2].vec[2], &subtraction);
	
	poly_add(&addition, &s1hats[0].vec[3], &s1hats[0].vec[0]);
	poly_add(&addition2, &s1hats[0].vec[1], &s1hats[0].vec[2]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_basemul_montgomery(&p25,  &mat[3].vec[3], &subtraction);
	
	poly_add(&addition, &mat[2].vec[2], &mat[2].vec[3]);
	poly_sub(&subtraction,  &s1hats[1].vec[3], &s1hats[1].vec[1]);
	poly_basemul_montgomery(&p26,  &addition, &subtraction);
	
	poly_sub(&subtraction,  &mat[3].vec[2], &mat[2].vec[2]);
	poly_add(&addition, &s1hats[0].vec[2], &s1hats[1].vec[2]);
	poly_add(&addition2, &s1hats[0].vec[0], &s1hats[1].vec[0]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_basemul_montgomery(&p27,  &subtraction, &subtraction2);
	
	poly_sub(&subtraction,  &mat[2].vec[3], &mat[3].vec[3]);
	poly_add(&addition, &s1hats[0].vec[3], &s1hats[1].vec[3]);
	poly_add(&addition2, &s1hats[0].vec[1], &s1hats[1].vec[1]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_basemul_montgomery(&p28,  &subtraction, &subtraction2);
	
	poly *poly_array10[] = {&mat[0].vec[0], &mat[0].vec[2], &mat[1].vec[1], &mat[1].vec[3]};
	poly_sum(&addition, poly_array10, 4);
	poly_add(&addition2, &s1hats[2].vec[2], &s1hats[3].vec[3]);
	poly_basemul_montgomery(&p29,  &addition, &addition2);
	
	poly *poly_array11[] = {&mat[1].vec[0], &mat[1].vec[2], &mat[1].vec[1], &mat[1].vec[3]};
	poly_sum(&addition, poly_array11, 4);
	poly_basemul_montgomery(&p30,  &addition, &s1hats[2].vec[2]);
	
	poly_add(&addition, &mat[0].vec[0], &mat[0].vec[2]);
	poly_sub(&subtraction,  &s1hats[3].vec[2], &s1hats[3].vec[3]);
	poly_basemul_montgomery(&p31,  &addition, &subtraction);
	
	poly_add(&addition, &mat[1].vec[1], &mat[1].vec[3]);
	poly_sub(&subtraction,  &s1hats[2].vec[3], &s1hats[2].vec[2]);
	poly_basemul_montgomery(&p32,  &addition, &subtraction);
	
	poly *poly_array12[] = {&mat[0].vec[0], &mat[0].vec[2], &mat[0].vec[1], &mat[0].vec[3]};
	poly_sum(&addition, poly_array12, 4);
	poly_basemul_montgomery(&p33,  &addition, &s1hats[3].vec[3]);
	
	poly_add(&addition, &mat[1].vec[0], &mat[1].vec[2]);
	poly_add(&addition2, &mat[0].vec[0], &mat[0].vec[2]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[2].vec[2], &s1hats[3].vec[2]);
	poly_basemul_montgomery(&p34, &subtraction, &addition3);
	
	poly_add(&addition, &mat[0].vec[1], &mat[0].vec[3]);
	poly_add(&addition2, &mat[1].vec[1], &mat[1].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[2].vec[3], &s1hats[3].vec[3]);
	poly_basemul_montgomery(&p35, &subtraction, &addition3);
	
	poly_add(&addition, &mat[2].vec[0], &mat[3].vec[1]);
	poly_add(&addition2, &mat[0].vec[0], &mat[1].vec[1]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array122[] = {&s1hats[0].vec[0], &s1hats[2].vec[0], &s1hats[1].vec[1], &s1hats[3].vec[1]};
	poly_sum(&addition3, poly_array122, 4);
	poly_basemul_montgomery(&p36,  &addition3, &subtraction);
	
	poly_add(&addition, &mat[3].vec[0], &mat[3].vec[1]);
	poly_add(&addition2, &mat[1].vec[0], &mat[1].vec[1]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[0].vec[0], &s1hats[2].vec[0]);
	poly_basemul_montgomery(&p37, &subtraction, &addition3);
	
	poly_sub(&subtraction,  &mat[2].vec[0], &mat[0].vec[0]);
	poly_add(&addition, &s1hats[1].vec[0], &s1hats[3].vec[0]);
	poly_add(&addition2, &s1hats[1].vec[1], &s1hats[3].vec[1]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_basemul_montgomery(&p38,  &subtraction, &subtraction2);
	
	poly_sub(&subtraction,  &mat[3].vec[1], &mat[1].vec[1]);
	poly_add(&addition, &s1hats[0].vec[1], &s1hats[2].vec[1]);
	poly_add(&addition2, &s1hats[0].vec[0], &s1hats[2].vec[0]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_basemul_montgomery(&p39,  &subtraction, &subtraction2);
	
	poly_add(&addition, &mat[2].vec[0], &mat[2].vec[1]);
	poly_add(&addition2, &mat[0].vec[0], &mat[0].vec[1]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[1].vec[1], &s1hats[3].vec[1]);
	poly_basemul_montgomery(&p40, &subtraction, &addition3);
	
	poly_add(&addition, &mat[3].vec[0], &mat[0].vec[0]);
	poly_add(&addition2, &mat[1].vec[0], &mat[2].vec[0]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array13[] = {&s1hats[0].vec[0], &s1hats[2].vec[0], &s1hats[1].vec[0], &s1hats[3].vec[0]};
	poly_sum(&addition3, poly_array13, 4);
	poly_basemul_montgomery(&p41,  &addition3, &subtraction);
	
	poly_add(&addition, &mat[2].vec[1], &mat[1].vec[1]);
	poly_add(&addition2, &mat[0].vec[1], &mat[3].vec[1]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array132[] = {&s1hats[0].vec[1], &s1hats[2].vec[1], &s1hats[1].vec[1], &s1hats[3].vec[1]};
	poly_sum(&addition3, poly_array132, 4);
	poly_basemul_montgomery(&p42,  &addition3, &subtraction);

	poly_add(&addition, &mat[0].vec[2], &mat[1].vec[3]);
	poly_add(&addition2, &mat[2].vec[2], &mat[3].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array14[] = {&s1hats[0].vec[2], &s1hats[2].vec[2], &s1hats[1].vec[3], &s1hats[3].vec[3]};
	poly_sum(&addition3, poly_array14, 4);
	poly_basemul_montgomery(&p43,  &addition3, &subtraction);

	poly_add(&addition, &mat[1].vec[2], &mat[1].vec[3]);
	poly_add(&addition2, &mat[3].vec[2], &mat[3].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[0].vec[2], &s1hats[2].vec[2]);
	poly_basemul_montgomery(&p44, &subtraction, &addition3);

	poly_sub(&subtraction,  &mat[0].vec[2], &mat[2].vec[2]);
	poly_add(&addition, &s1hats[1].vec[2], &s1hats[3].vec[2]);
	poly_add(&addition2, &s1hats[1].vec[3], &s1hats[3].vec[3]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_basemul_montgomery(&p45,  &subtraction, &subtraction2);
	
	poly_sub(&subtraction,  &mat[1].vec[3], &mat[3].vec[3]);
	poly_add(&addition, &s1hats[0].vec[3], &s1hats[2].vec[3]);
	poly_add(&addition2, &s1hats[0].vec[2], &s1hats[2].vec[2]);
	poly_sub(&subtraction2,  &addition, &addition2);
	poly_basemul_montgomery(&p46,  &subtraction, &subtraction2);

	poly_add(&addition, &mat[0].vec[2], &mat[0].vec[3]);
	poly_add(&addition2, &mat[2].vec[2], &mat[2].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly_add(&addition3, &s1hats[1].vec[3], &s1hats[3].vec[3]);
	poly_basemul_montgomery(&p47, &subtraction, &addition3);

	poly_add(&addition, &mat[1].vec[2], &mat[2].vec[2]);
	poly_add(&addition2, &mat[3].vec[2], &mat[0].vec[2]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array15[] = {&s1hats[0].vec[2], &s1hats[2].vec[2], &s1hats[1].vec[2], &s1hats[3].vec[2]};
	poly_sum(&addition3, poly_array15, 4);
	poly_basemul_montgomery(&p48,  &addition3, &subtraction);
	
	poly_add(&addition, &mat[0].vec[3], &mat[3].vec[3]);
	poly_add(&addition2, &mat[2].vec[3], &mat[1].vec[3]);
	poly_sub(&subtraction,  &addition, &addition2);
	poly *poly_array16[] = {&s1hats[0].vec[3], &s1hats[2].vec[3], &s1hats[1].vec[3], &s1hats[3].vec[3]};
	poly_sum(&addition3, poly_array16, 4);
	poly_basemul_montgomery(&p49,  &addition3, &subtraction);



	//Computation of the transpose of C=(A).(S1_hat)
	poly *poly_array17[] = {&p1, &p4, &p7, &p22, &p25, &p28, &p33, &p43, &p46, &p49};
	poly_sum(&addition, poly_array17, 10);
	poly *poly_array18[] = {&p5, &p26, &p29, &p32, &p35, &p47};
	poly_sum(&addition2, poly_array18, 6);
	poly_sub(&C[0].vec[0], &addition, &addition2);
		
	poly *poly_array19[] = {&p3, &p5, &p24, &p26, &p45, &p47};
	poly_sum(&addition, poly_array19, 6);
	poly_add(&addition2, &p31, &p33);
	poly_sub(&C[1].vec[0], &addition, &addition2);
		
	poly *poly_array21[] = {&p15, &p18, &p21, &p29, &p32, &p35};
	poly_sum(&addition, poly_array21, 6);
	poly_add(&addition2, &p19, &p33);
	poly_sub(&C[2].vec[0], &addition, &addition2);
		
	poly *poly_array23[] = {&p17, &p19, &p31, &p33};
	poly_sum(&C[3].vec[0], poly_array23, 4);



	poly *poly_array25[] = {&p2, &p4, &p23, &p25, &p44, &p46};
	poly_sum(&addition, poly_array25, 6);
	poly_add(&addition2, &p30, &p32);
	poly_sub(&C[0].vec[1], &addition, &addition2);
		
	poly *poly_array27[] = {&p1, &p3, &p6, &p22, &p24, &p27, &p30, &p43, &p45, &p48};
	poly_sum(&addition, poly_array27, 10);
	poly *poly_array28[] = {&p2, &p23, &p29, &p31, &p34, &p44};
	poly_sum(&addition2, poly_array28, 6); 
	poly_sub(&C[1].vec[1], &addition, &addition2);
		
	poly *poly_array29[] = {&p16, &p18, &p30, &p32};
	poly_sum(&C[2].vec[1], poly_array29, 4);
	
	poly *poly_array31[] = {&p15, &p17, &p20, &p29, &p31, &p34};
	poly_sum(&addition, poly_array31, 6);
	poly_add(&addition2, &p16, &p30);
	poly_sub(&C[3].vec[1], &addition, &addition2);
	

	
	poly *poly_array33[] = {&p8, &p11, &p14, &p22, &p25, &p28};
	poly_sum(&addition, poly_array33, 6);
	poly_add(&addition2, &p12, &p26);
	poly_sub(&C[0].vec[2], &addition, &addition2);
		
	poly *poly_array35[] = {&p10, &p12, &p24, &p26};
	poly_sum(&C[1].vec[2], poly_array35, 4);
			
	poly *poly_array37[] = {&p1, &p4, &p7, &p12, &p15, &p18, &p21, &p36, &p39, &p42};
	poly_sum(&addition, poly_array37, 10);
	poly *poly_array38[] = {&p5, &p8, &p11, &p14, &p19, &p40};
	poly_sum(&addition2, poly_array38, 6);
	poly_sub(&C[2].vec[2], &addition, &addition2);
	
	poly *poly_array39[] = {&p3, &p5, &p17, &p19, &p38, &p40};
	poly_sum(&addition, poly_array39, 6);
	poly_add(&addition2, &p10, &p12);
	poly_sub(&C[3].vec[2], &addition, &addition2);
	
	
	
	poly *poly_array41[] = {&p9, &p11, &p23, &p25};
	poly_sum(&C[0].vec[3], poly_array41, 4);

		
	poly *poly_array43[] = {&p8, &p10, &p13, &p22, &p24, &p27};
	poly_sum(&addition, poly_array43, 6);
	poly_add(&addition2, &p9, &p23);
	poly_sub(&C[1].vec[3], &addition, &addition2);
		
	poly *poly_array45[] = {&p2, &p4, &p16, &p18, &p37, &p39};
	poly_sum(&addition, poly_array45, 6);
	poly_add(&addition2, &p9, &p11);
	poly_sub(&C[2].vec[3], &addition, &addition2);
	
	poly *poly_array47[] = {&p1, &p3, &p6, &p9, &p15, &p17, &p20, &p36, &p38, &p41};
	poly_sum(&addition, poly_array47, 10);
	poly *poly_array48[] = {&p2, &p8, &p10, &p13, &p16, &p37};
	poly_sum(&addition2, poly_array48, 6);
	poly_sub(&C[3].vec[3], &addition, &addition2);

  return 0;
}
int mult_p4(polyvec C[BATCH_SIZE], polyvec mat[KYBER_K], polyvec s1hats[BATCH_SIZE], poly p[]){
    // C transpose is taken here
  poly addition, addition2, addition3, addition4, subtraction;
  for (int i = 0; i < 4; i++) {
      poly_add(&addition, &s1hats[0].vec[0], &mat[i].vec[1]);
      poly_basemul_montgomery(&p[2 * i], &mat[i].vec[0], &addition);

      poly_add(&addition, &s1hats[0].vec[2], &mat[i].vec[3]);
      poly_basemul_montgomery(&p[2 * i + 1], &mat[i].vec[2], &addition);
  }


  int index = 8;
  for (int i = 1; i <= 3; i++) {
      poly_add(&addition, &s1hats[0].vec[0], &s1hats[i].vec[0]);
      poly_basemul_montgomery(&p[index], &s1hats[i].vec[1], &addition);
      index++;

      poly_add(&addition, &s1hats[0].vec[2], &s1hats[i].vec[2]);
      poly_basemul_montgomery(&p[index], &s1hats[i].vec[3], &addition);
      index++;
  }
  index = 14;
  for (int i = 0; i < 4; i++) {
      poly_sub(&subtraction, &s1hats[0].vec[1], &mat[i].vec[0]);
      poly_basemul_montgomery(&p[index], &mat[i].vec[1], &subtraction);
      index++;

      poly_sub(&subtraction, &s1hats[0].vec[3], &mat[i].vec[2]);
      poly_basemul_montgomery(&p[index], &mat[i].vec[3], &subtraction);
      index++;
  }
  index = 22;

  for (int j = 0;j<4;j++){
    for (int i=0;i<3;i++){
      poly_add(&addition, &mat[j].vec[0],&s1hats[i+1].vec[1]);
      poly *poly_array[] = {&mat[j].vec[1], &s1hats[0].vec[0], &s1hats[i+1].vec[0]};
      poly_sum(&addition2, poly_array, 3);
      poly_basemul_montgomery(&p[index], &addition, &addition2);
      index++;
      poly_add(&addition, &mat[j].vec[2],&s1hats[i+1].vec[3]);
      poly *poly_array1[] = {&mat[j].vec[3], &s1hats[0].vec[2], &s1hats[i+1].vec[2]};
      poly_sum(&addition2, poly_array1, 3);
      poly_basemul_montgomery(&p[index], &addition, &addition2);
      index++;
    }
  }

	poly *poly_array24[] = {&p[0], &p[1], &p[14], &p[15]};
	poly_sum(&C[0].vec[0], poly_array24, 4);
	
	
	poly_add(&addition, &p[22],&p[23]);
	poly *poly_array25[] = { &p[0], &p[1], &p[8], &p[9]};
	poly_sum(&addition4, poly_array25, 4);

	poly_sub(&C[1].vec[0], &addition, &addition4);
	
	poly_add(&addition, &p[24],&p[25]);
	poly *poly_array26[] = { &p[0], &p[1], &p[10], &p[11]};
	poly_sum(&addition4, poly_array26, 4);

	poly_sub(&C[2].vec[0], &addition, &addition4);
	
	poly_add(&addition, &p[26],&p[27]);
	poly *poly_array27[] = { &p[0], &p[1], &p[12], &p[13]};
	poly_sum(&addition4, poly_array27, 4);

	poly_sub(&C[3].vec[0], &addition, &addition4);
 
 
	poly *poly_array28[] = {&p[2], &p[3], &p[16], &p[17]};
	poly_sum(&C[0].vec[1], poly_array28, 4);
	
	
	poly_add(&addition, &p[28],&p[29]);
	poly *poly_array29[] = { &p[2], &p[3], &p[8], &p[9]};
	poly_sum(&addition4, poly_array29, 4);

	poly_sub(&C[1].vec[1], &addition, &addition4);

	poly_add(&addition, &p[30],&p[31]);
	poly *poly_array30[] = { &p[2], &p[3], &p[10], &p[11]};
	poly_sum(&addition4, poly_array30, 4);

	poly_sub(&C[2].vec[1], &addition, &addition4);
    
	poly_add(&addition, &p[32],&p[33]);
	poly *poly_array31[] = { &p[2], &p[3], &p[12], &p[13]};
	poly_sum(&addition4, poly_array31, 4);

	poly_sub(&C[3].vec[1], &addition, &addition4);


	poly *poly_array32[] = {&p[4], &p[5], &p[18], &p[19]};
	poly_sum(&C[0].vec[2], poly_array32, 4);
	
	poly_add(&addition, &p[34],&p[35]);
	poly *poly_array33[] = { &p[4], &p[5], &p[8], &p[9]};
	poly_sum(&addition4, poly_array33, 4);

	poly_sub(&C[1].vec[2], &addition, &addition4);
	
	poly_add(&addition, &p[36],&p[37]);
	poly *poly_array34[] = { &p[4], &p[5], &p[10], &p[11]};
	poly_sum(&addition4, poly_array34, 4);

	poly_sub(&C[2].vec[2], &addition, &addition4);
	
	poly_add(&addition, &p[38],&p[39]);
	poly *poly_array35[] = { &p[4], &p[5], &p[12], &p[13]};
	poly_sum(&addition4, poly_array35, 4);

	poly_sub(&C[3].vec[2], &addition, &addition4);

	poly *poly_array36[] = {&p[6], &p[7], &p[20], &p[21]};
	poly_sum(&C[0].vec[3], poly_array36, 4);

	poly_add(&addition, &p[40],&p[41]);
	poly *poly_array37[] = { &p[6], &p[7], &p[8], &p[9]};
	poly_sum(&addition4, poly_array37, 4);

	poly_sub(&C[1].vec[3], &addition, &addition4);
	
	poly_add(&addition, &p[42],&p[43]);
	poly *poly_array38[] = { &p[6], &p[7], &p[10], &p[11]};
	poly_sum(&addition4, poly_array38, 4);

	poly_sub(&C[2].vec[3], &addition, &addition4);
    
	poly_add(&addition, &p[44],&p[45]);
	poly *poly_array39[] = { &p[6], &p[7], &p[12], &p[13]};
	poly_sum(&addition4, poly_array39, 4);

	poly_sub(&C[3].vec[3], &addition, &addition4);

  return 0;
}
int mult_p2(polyvec C[BATCH_SIZE], polyvec mat[KYBER_K], polyvec s1hats[BATCH_SIZE]){

  poly addition, addition2, subtraction, subtraction2;
  poly p1,p2,p3,p4,p5,p6,p7;

	// Computation of p_i's
	
	poly_add(&addition, &mat[0].vec[0], &mat[1].vec[1]);
	poly_add(&addition2, &s1hats[0].vec[0], &s1hats[1].vec[1]);
	poly_basemul_montgomery(&p1,  &addition, &addition2);
	
	poly_add(&addition, &mat[1].vec[0], &mat[1].vec[1]);
	poly_basemul_montgomery(&p2,  &addition, &s1hats[0].vec[0]);
	
	poly_sub(&subtraction,  &s1hats[1].vec[0], &s1hats[1].vec[1]);
	poly_basemul_montgomery(&p3,  &mat[0].vec[0], &subtraction);
	
	poly_sub(&subtraction,  &s1hats[0].vec[1], &s1hats[0].vec[0]);
	poly_basemul_montgomery(&p4,  &mat[1].vec[1], &subtraction);
	
	poly_add(&addition, &mat[0].vec[0], &mat[0].vec[1]);
	poly_basemul_montgomery(&p5,  &addition, &s1hats[1].vec[1]);
	
	poly_sub(&subtraction, &mat[1].vec[0], &mat[0].vec[0]);
	poly_add(&addition, &s1hats[0].vec[0], &s1hats[1].vec[0]);
	poly_basemul_montgomery(&p6,  &subtraction, &addition);
	
	poly_sub(&subtraction, &mat[0].vec[1], &mat[1].vec[1]);
	poly_add(&addition, &s1hats[0].vec[1], &s1hats[1].vec[1]);
	poly_basemul_montgomery(&p7,  &subtraction, &addition);

	//Computation of the transpose of C=(A).(S1_hat)
	
	poly *poly_array[] = {&p1, &p4, &p7};
	poly_sum(&addition, poly_array, 3);
	poly_sub(&C[0].vec[0], &addition, &p5);
	
	poly_add(&C[1].vec[0], &p3, &p5);
	
	poly_add(&C[0].vec[1], &p2, &p4);
	
	poly *poly_array2[] = {&p1, &p3, &p6};
	poly_sum(&addition, poly_array2, 3);
	poly_sub(&C[1].vec[1], &addition, &p2);

  return 0;
}


int mult_p3(polyvec C[BATCH_SIZE], polyvec mat[KYBER_K], polyvec s1hats[BATCH_SIZE], poly p[]){

  poly addition, addition2, subtraction, subtraction2;
  int j;
  
	// Computation of p_i's

	int index=0;

	for (int j=1; j<=2; j++){
		for(int i=0; i<=2; i++){
			poly_add(&addition, &mat[i].vec[0], &s1hats[0].vec[j]);
			poly_add(&addition2, &mat[i].vec[j], &s1hats[j].vec[0]);
			poly_basemul_montgomery(&p[index],  &addition, &addition2);
			index++;
		}

	}

	poly_basemul_montgomery(&p[6],  &s1hats[1].vec[0], &s1hats[0].vec[1]);
	
	poly_basemul_montgomery(&p[7],  &s1hats[2].vec[0], &s1hats[0].vec[2]);
	
	poly_basemul_montgomery(&p[8],  &s1hats[2].vec[1], &s1hats[1].vec[2]);


	index=9;

	for (int i=0; i<=2; i++){
		poly_add(&addition, &mat[i].vec[1], &s1hats[1].vec[2]);
		poly_add(&addition2, &mat[i].vec[2], &s1hats[2].vec[1]);
		poly_basemul_montgomery(&p[index],  &addition, &addition2);
		index++;
	}

	poly_sub(&subtraction , &s1hats[0].vec[1], &s1hats[3].vec[1]);
	poly_sub(&subtraction2, &s1hats[3].vec[0], &s1hats[1].vec[0]);
	poly_basemul_montgomery(&p[12],  &subtraction, &subtraction2);

		index=13;


	for(j=0; j<=2; j++){
			poly *poly_array[] = {&s1hats[1].vec[0], &s1hats[2].vec[0], &mat[j].vec[1], &mat[j].vec[2]};
			poly_sum(&addition, poly_array, 4);
			poly_sub(&subtraction, &s1hats[0].vec[0], &addition);
			poly_basemul_montgomery(&p[index],  &mat[j].vec[0], &subtraction);
			index++;
	}

	index=16;

	for(int i=0; i<=2; i++){
		poly *poly_array4[] = {&s1hats[0].vec[1], &s1hats[2].vec[1], &mat[i].vec[0], &mat[i].vec[2]};
		poly_sum(&addition, poly_array4, 4);
		poly_sub(&subtraction, &s1hats[1].vec[1], &addition);
		poly_basemul_montgomery(&p[index],  &mat[i].vec[1], &subtraction);
		index++;
	}

	index=19;

	for(int i=0; i<=2; i++){
		poly *poly_array7[] = {&s1hats[0].vec[2], &s1hats[1].vec[2], &mat[i].vec[0], &mat[i].vec[1]};
		poly_sum(&addition, poly_array7, 4);
		poly_sub(&subtraction, &s1hats[2].vec[2], &addition);
		poly_basemul_montgomery(&p[index],  &mat[i].vec[2], &subtraction);
		index++;
	}

	index=22;

	for(int i=0; i<=2; i++){
		poly_add(&addition, &mat[i].vec[0], &s1hats[0].vec[1]);
		poly_sub(&subtraction,  &addition, &s1hats[3].vec[1]);
		poly_add(&addition2, &mat[i].vec[1], &s1hats[1].vec[0]);
		poly_sub(&subtraction2, &s1hats[3].vec[0], &addition2);
		poly_basemul_montgomery(&p[index],  &subtraction, &subtraction2);
		index++;
	}

	index=25;

	for(int i=0; i<=2; i++){
		poly_basemul_montgomery(&p[index],  &mat[i].vec[2], &s1hats[3].vec[2]);
		index++;
	}
	
	//Computation of the transpose of C=(A).(S1_hat)
	
	poly *poly_array10[] = {&p[0], &p[3], &p[13]};
	poly_sum(&addition, poly_array10, 3);
	poly_add(&addition2, &p[6], &p[7]);
	poly_sub(&C[0].vec[0],  &addition, &addition2);
	
	poly *poly_array11[] = {&p[0], &p[9], &p[16]};
	poly_sum(&addition, poly_array11, 3);
	poly_add(&addition2, &p[6], &p[8]);
	poly_sub(&C[1].vec[0],  &addition, &addition2);
	
	poly *poly_array12[] = {&p[3], &p[9], &p[19]};
	poly_sum(&addition, poly_array12, 3);
	poly_add(&addition2, &p[7], &p[8]);
	poly_sub(&C[2].vec[0],  &addition, &addition2);
	
	poly *poly_array13[] = {&p[0], &p[22], &p[25]};
	poly_sum(&addition, poly_array13, 3);
	poly_add(&addition2, &p[6], &p[12]);
	poly_sub(&C[3].vec[0],  &addition, &addition2);
	
	
	poly *poly_array14[] = {&p[1], &p[4], &p[14]};
	poly_sum(&addition, poly_array14, 3);
	poly_add(&addition2, &p[6], &p[7]);
	poly_sub(&C[0].vec[1],  &addition, &addition2);
	
	poly *poly_array15[] = {&p[1], &p[10], &p[17]};
	poly_sum(&addition, poly_array15, 3);
	poly_add(&addition2, &p[6], &p[8]);
	poly_sub(&C[1].vec[1],  &addition, &addition2);
	
	poly *poly_array16[] = {&p[4], &p[10], &p[20]};
	poly_sum(&addition, poly_array16, 3);
	poly_add(&addition2, &p[7], &p[8]);
	poly_sub(&C[2].vec[1],  &addition, &addition2);
	
	poly *poly_array17[] = {&p[1], &p[23], &p[26]};
	poly_sum(&addition, poly_array17, 3);
	poly_add(&addition2, &p[6], &p[12]);
	poly_sub(&C[3].vec[1],  &addition, &addition2);
	
	
	
	poly *poly_array18[] = {&p[2], &p[5], &p[15]};
	poly_sum(&addition, poly_array18, 3);
	poly_add(&addition2, &p[6], &p[7]);
	poly_sub(&C[0].vec[2],  &addition, &addition2);
	
	poly *poly_array19[] = {&p[2], &p[11], &p[18]};
	poly_sum(&addition, poly_array19, 3);
	poly_add(&addition2, &p[6], &p[8]);
	poly_sub(&C[1].vec[2],  &addition, &addition2);
	
	poly *poly_array20[] = {&p[5], &p[11], &p[21]};
	poly_sum(&addition, poly_array20, 3);
	poly_add(&addition2, &p[7], &p[8]);
	poly_sub(&C[2].vec[2],  &addition, &addition2);
	
	poly *poly_array21[] = {&p[2], &p[24], &p[27]};
	poly_sum(&addition, poly_array21, 3);
	poly_add(&addition2, &p[6], &p[12]);
	poly_sub(&C[3].vec[2],  &addition, &addition2);
	
  return 0;
}
/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk
*              and the public seed used to generate the matrix A.
*
* Arguments:   uint8_t *r: pointer to the output serialized public key
*              polyvec *pk: pointer to the input public-key polyvec
*              const uint8_t *seed: pointer to the input public seed
**************************************************/
static void pack_pk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES])
{
  size_t i;
  polyvec_tobytes(r, pk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    r[i+KYBER_POLYVECBYTES] = seed[i];
}

/*************************************************
* Name:        unpack_pk
*
* Description: De-serialize public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk: pointer to output public-key polynomial vector
*              - uint8_t *seed: pointer to output seed to generate matrix A
*              - const uint8_t *packedpk: pointer to input serialized public key
**************************************************/
static void unpack_pk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES])
{
  size_t i;
  polyvec_frombytes(pk, packedpk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    seed[i] = packedpk[i+KYBER_POLYVECBYTES];
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key
*
* Arguments:   - uint8_t *r: pointer to output serialized secret key
*              - polyvec *sk: pointer to input vector of polynomials (secret key)
**************************************************/
static void pack_sk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  polyvec_tobytes(r, sk);
}

/*************************************************
* Name:        unpack_sk
*
* Description: De-serialize the secret key; inverse of pack_sk
*
* Arguments:   - polyvec *sk: pointer to output vector of polynomials (secret key)
*              - const uint8_t *packedsk: pointer to input serialized secret key
**************************************************/
static void unpack_sk(polyvec *sk, const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec_frombytes(sk, packedsk);
}

/*************************************************
* Name:        pack_ciphertext
*
* Description: Serialize the ciphertext as concatenation of the
*              compressed and serialized vector of polynomials b
*              and the compressed and serialized polynomial v
*
* Arguments:   uint8_t *r: pointer to the output serialized ciphertext
*              poly *pk: pointer to the input vector of polynomials b
*              poly *v: pointer to the input polynomial v
**************************************************/
static void pack_ciphertext(uint8_t r[KYBER_INDCPA_BYTES], polyvec *b, poly *v)
{
  polyvec_compress(r, b);
  poly_compress(r+KYBER_POLYVECCOMPRESSEDBYTES, v);
}

/*************************************************
* Name:        unpack_ciphertext
*
* Description: De-serialize and decompress ciphertext from a byte array;
*              approximate inverse of pack_ciphertext
*
* Arguments:   - polyvec *b: pointer to the output vector of polynomials b
*              - poly *v: pointer to the output polynomial v
*              - const uint8_t *c: pointer to the input serialized ciphertext
**************************************************/
static void unpack_ciphertext(polyvec *b, poly *v, const uint8_t c[KYBER_INDCPA_BYTES])
{
  polyvec_decompress(b, c);
  poly_decompress(v, c+KYBER_POLYVECCOMPRESSEDBYTES);
}

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r: pointer to output buffer
*              - unsigned int len: requested number of 16-bit integers (uniform mod q)
*              - const uint8_t *buf: pointer to input buffer (assumed to be uniformly random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint16_t val0, val1;

  ctr = pos = 0;
  while(ctr < len && pos + 3 <= buflen) {
    val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
    val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;

    if(val0 < KYBER_Q)
      r[ctr++] = val0;
    if(ctr < len && val1 < KYBER_Q)
      r[ctr++] = val1;
  }

  return ctr;
}

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

/*************************************************
* Name:        gen_matrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *a: pointer to ouptput matrix A
*              - const uint8_t *seed: pointer to input seed
*              - int transposed: boolean deciding whether A or A^T is generated
**************************************************/
#define GEN_MATRIX_NBLOCKS ((12*KYBER_N/8*(1 << 12)/KYBER_Q + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)
// Not static for benchmarking
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed)
{
  unsigned int ctr, i, j, k;
  unsigned int buflen, off;
  uint8_t buf[GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES+2];
  xof_state state;

  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      if(transposed)
        xof_absorb(&state, seed, i, j);
      else
        xof_absorb(&state, seed, j, i);

      xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
      buflen = GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES;
      ctr = rej_uniform(a[i].vec[j].coeffs, KYBER_N, buf, buflen);

      while(ctr < KYBER_N) {
        off = buflen % 3;
        for(k = 0; k < off; k++)
          buf[k] = buf[buflen - off + k];
        xof_squeezeblocks(buf + off, 1, &state);
        buflen = off + XOF_BLOCKBYTES;
        ctr += rej_uniform(a[i].vec[j].coeffs + ctr, KYBER_N - ctr, buf, buflen);
      }
    }
  }
}

/*************************************************
* Name:        indcpa_keypair
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
                              (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
**************************************************/
void indcpa_keypair(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  unsigned int i;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  uint8_t nonce = 0;
  polyvec a[KYBER_K], e, pkpv, skpv;

  randombytes(buf, KYBER_SYMBYTES);
  hash_g(buf, buf, KYBER_SYMBYTES);

  gen_a(a, publicseed);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(&skpv.vec[i], noiseseed, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(&e.vec[i], noiseseed, nonce++);

  polyvec_ntt(&skpv);
  polyvec_ntt(&e);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
    polyvec_basemul_acc_montgomery(&pkpv.vec[i], &a[i], &skpv);
    poly_tomont(&pkpv.vec[i]);
  }

  polyvec_add(&pkpv, &pkpv, &e);
  polyvec_reduce(&pkpv);

  pack_sk(sk, &skpv);
  pack_pk(pk, &pkpv, publicseed);
}

/*************************************************
* Name:        indcpa_enc
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *c: pointer to output ciphertext
*                            (of length KYBER_INDCPA_BYTES bytes)
*              - const uint8_t *m: pointer to input message
*                                  (of length KYBER_INDCPA_MSGBYTES bytes)
*              - const uint8_t *pk: pointer to input public key
*                                   (of length KYBER_INDCPA_PUBLICKEYBYTES)
*              - const uint8_t *coins: pointer to input random coins used as seed
*                                      (of length KYBER_SYMBYTES) to deterministically
*                                      generate all randomness
**************************************************/
void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES];
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, at[KYBER_K], b;
  poly v, k, epp;
 
  unpack_pk(&pkpv, seed, pk);
// for(int k=0;k<KYBER_SYMBYTES;k++)
//   	printf("%02X", seed[k]);
// printf("\n");
  poly_frommsg(&k, m);
  gen_at(at, seed);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(sp.vec+i, coins, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta2(ep.vec+i, coins, nonce++);
  poly_getnoise_eta2(&epp, coins, nonce++);
  polyvec_ntt(&sp);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++)
    polyvec_basemul_acc_montgomery(&b.vec[i], &at[i], &sp);
  polyvec_basemul_acc_montgomery(&v, &pkpv, &sp);
  polyvec_invntt_tomont(&b);
  poly_invntt_tomont(&v);

  polyvec_add(&b, &b, &ep);
  poly_add(&v, &v, &epp);
  poly_add(&v, &v, &k);
  polyvec_reduce(&b);
  poly_reduce(&v);
  pack_ciphertext(c, &b, &v);

	// for (int p = 0; p<KYBER_INDCPA_BYTES;p++)
	// 	printf("%02X", c[p]);
	// printf("\n");
}
 
void indcpa_enc_batch(uint8_t *cs[KYBER_INDCPA_BYTES],
                const uint8_t *ms[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[][2*KYBER_SYMBYTES],
                poly p[])
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES]; // rho
  uint8_t nonce = 0; // N
  polyvec sps[N_MSG], pkpv, eps[N_MSG], at[KYBER_K], b[N_MSG]; // pkpv t-hat
  poly vs[N_MSG], ks[N_MSG], epps[N_MSG];
 
  unpack_pk(&pkpv, seed, pk);
//   for(int t=0;t<KYBER_INDCPA_PUBLICKEYBYTES;t++)
//   	printf("%02X", pk[t]);
//   printf("\n");
  for (int j=0;j<N_MSG;j++)
	poly_frommsg(&ks[j], ms[j]);
    
  gen_at(at, seed);
    
  for (int j=0;j<N_MSG;j++){
    for(i=0;i<KYBER_K;i++)
      poly_getnoise_eta1(sps[j].vec+i, coins[j]+KYBER_SYMBYTES, nonce++); //sps[j] --> r[j]
    for(i=0;i<KYBER_K;i++)
      poly_getnoise_eta2(eps[j].vec+i, coins[j]+KYBER_SYMBYTES, nonce++); //eps[j] --> e1[j]
    poly_getnoise_eta2(&epps[j], coins[j]+KYBER_SYMBYTES, nonce++);       //epps[j] --> e2[j]
    polyvec_ntt(&sps[j]);
    nonce = 0;
  }
  // matrix-vector multiplication
  for (int j = 0; j<N_MSG-(N_MSG%BATCH_SIZE); j+=BATCH_SIZE)
    mult_p3(&b[j], at, &sps[j], p);
  // Remaining
  for (int j = N_MSG%BATCH_SIZE; j>0;j--){
	for(i=0;i<KYBER_K;i++){
    	polyvec_basemul_acc_montgomery(&b[N_MSG - j].vec[i], &at[i], &sps[N_MSG - j]);
	}
  } 
  
  for (int j=0;j<N_MSG;j++){
    polyvec_basemul_acc_montgomery(&vs[j], &pkpv, &sps[j]);
    polyvec_invntt_tomont(&b[j]);
    poly_invntt_tomont(&vs[j]);
    polyvec_add(&b[j], &b[j], &eps[j]);
    poly_add(&vs[j], &vs[j], &epps[j]);
    poly_add(&vs[j], &vs[j], &ks[j]);
    polyvec_reduce(&b[j]);
    poly_reduce(&vs[j]);
    pack_ciphertext(cs[j], &b[j], &vs[j]);
	// for (int p = 0; p<KYBER_INDCPA_BYTES;p++)
	// 	printf("%02X", cs[j][p]);
	// printf("\n");
  }
}


/*************************************************
* Name:        indcpa_dec
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *m: pointer to output decrypted message
*                            (of length KYBER_INDCPA_MSGBYTES)
*              - const uint8_t *c: pointer to input ciphertext
*                                  (of length KYBER_INDCPA_BYTES)
*              - const uint8_t *sk: pointer to input secret key
*                                   (of length KYBER_INDCPA_SECRETKEYBYTES)
**************************************************/
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec b, skpv;
  poly v, mp;

  unpack_ciphertext(&b, &v, c);
  unpack_sk(&skpv, sk);

  polyvec_ntt(&b);
  polyvec_basemul_acc_montgomery(&mp, &skpv, &b);
  poly_invntt_tomont(&mp);

  poly_sub(&mp, &v, &mp);
  poly_reduce(&mp);

  poly_tomsg(m, &mp);
}
