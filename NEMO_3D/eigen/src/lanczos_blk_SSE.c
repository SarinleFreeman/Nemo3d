#pragma message("INTEL_MODIFIED: lanczos first small loop modification is ON")
__declspec(align(16)) nml_dscalar a1=*alpha;
__declspec(align(16)) nml_dcscalar b2;
      b2.r=(*beta1).r;
      b2.i=(*beta1).i;
__declspec(align(16)) nml_dscalar a2=*alpha2;
__declspec(align(16)) nml_dcscalar b1;
      b1.r=(*beta2).r;
      b1.i=(*beta2).i;
__declspec(align(16)) nml_dscalar g1=*gamma1;
__declspec(align(16)) nml_dscalar g2=*gamma2;
      
	__m128d __a1, __a2, __b1, __b1_r, __b1_i, __b2, __b2_r, __b2_i, __g1;
	__m128d __q1_n, __q2_n, __q2_s, __q3_n, __q4_n, __q4_s, __r1_n;
	__m128d __b1_ir, __b2_ir;
	
	__declspec(align(16)) int imask[4] = {0,0,0, 0x80000000};
    __m128d mask = _mm_load_pd((double*) (&imask));

	__a1 = _mm_load_sd(&a1);
	__a1 = _mm_unpacklo_pd(__a1, __a1);

	__b1 = _mm_load_pd((double*) &b1.r);
	__b1_r = _mm_shuffle_pd(__b1, __b1, 0x3);
	__b1_i = _mm_shuffle_pd(__b1, __b1, 0x0);
	__b1_ir = _mm_xor_pd(__b1_i, mask);

	__g1 = _mm_load_sd(&g1);
	__g1 = _mm_unpacklo_pd(__g1, __g1);
	
	__b2 = _mm_load_pd((double*) &b2.r);
	__b2_r = _mm_shuffle_pd(__b2, __b2, 0x3);
	__b2_i = _mm_shuffle_pd(__b2, __b2, 0x0);
	__b2_ir = _mm_xor_pd(__b2_i, mask);
// prepare __a1, __b1_r, __b1_i, __b2_r, __b2_i, __g1
	for (j = 0; j < m; ++j) {
		//r1_n[j].r -= a1*q1_n[j].r + b1.r*q2_n[j].r - b1.i*q2_n[j].i + g1*q3_n[j].r + b2.r*q4_n[j].r + b2.i*q4_n[j].i;
		//r1_n[j].i -= a1*q1_n[j].i + b1.r*q2_n[j].i + b1.i*q2_n[j].r + g1*q3_n[j].i + b2.r*q4_n[j].i - b2.i*q4_n[j].r;		
		
		//a1*q1_n[j].r
		__q1_n = _mm_load_pd((double*)&q1_n[j]);
		__q1_n = _mm_mul_pd(__q1_n, __a1); 	
		//b1.r*q2_n[j].r
		__q2_n = _mm_load_pd((double*)&q2_n[j]);
		__q2_s = _mm_shuffle_pd(__q2_n,__q2_n,01);
		__q2_n = _mm_mul_pd(__q2_n, __b1_r);
		__q1_n = _mm_add_pd(__q1_n, __q2_n);
		//b1.i*q2_n[j].i
		__q2_s = _mm_mul_pd(__q2_s, __b1_ir);
		__q1_n = _mm_sub_pd(__q1_n, __q2_s);
		//g1*q3_n[j].r
		__q3_n = _mm_load_pd((double*)&q3_n[j]);
		__q3_n = _mm_mul_pd(__q3_n, __g1);	
		__q1_n = _mm_add_pd(__q1_n, __q3_n);
		//b2.r*q4_n[j].r
		__q4_n = _mm_load_pd((double*)&q4_n[j]);
		__q4_s = _mm_shuffle_pd(__q4_n,__q4_n,01);
		__q4_n = _mm_mul_pd(__q4_n, __b2_r);
		__q1_n = _mm_add_pd(__q1_n, __q4_n);
		//b2.i*q4_n[j].i
		__q4_s = _mm_mul_pd(__q4_s, __b2_ir);
		__q1_n = _mm_add_pd(__q1_n, __q4_s);
		// save result		r1_n[j].r -=
		__r1_n = _mm_load_pd((double*) &r1_n[j]);
		__r1_n = _mm_sub_pd(__r1_n, __q1_n);
		_mm_store_pd((double*)&r1_n[j], __r1_n);

      		r2_n[j].r -= b1.r*q1_n[j].r+b1.i*q1_n[j].i + a2*q2_n[j].r+g2*q4_n[j].r;
      		r2_n[j].i -= b1.r*q1_n[j].i-b1.i*q1_n[j].r + a2*q2_n[j].i+g2*q4_n[j].i;
		}
		
#pragma message("INTEL_MODIFIED: Vector copy (memcpy, lanczos_blk) optimization is ON")
	memcpy((void*)q3_n,(void*) q1_n, m*sizeof(nml_dscalar));
	memcpy((void*)q4_n,(void*) q2_n, m*sizeof(nml_dscalar));

