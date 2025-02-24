#pragma message ("INTEL_MODIFIED: Lanczos loop 2 vectorization is ON")
	
   // modified and tested by Hoon Ryu and Sunhee Lee
   // 2006.10.21

        nml_dscalar betap=0;
	__m128d __r_nj, __q_nj, __q_n1j, __al, __beta_n1, __betap;

	__al = _mm_load_sd((double*) alpha);
	__al = _mm_unpacklo_pd(__al, __al);

	__beta_n1 = _mm_load_sd((double *)&beta_n1);
	__beta_n1 = _mm_unpacklo_pd(__beta_n1, __beta_n1);

	for(j=0; j<m; ++j)
	{
		__q_nj = _mm_load_pd((double *)&q_n[j]);
		__r_nj = _mm_load_pd((double *)&r_n[j]);
		__q_n1j = _mm_load_pd((double *)&q_n1[j]);
		
		__q_nj = _mm_mul_pd(__q_nj, __al);
		__r_nj = _mm_sub_pd(__r_nj, __q_nj);
		
		__q_n1j = _mm_mul_pd(__q_n1j, __beta_n1);
		__r_nj = _mm_sub_pd(__r_nj, __q_n1j);

		__q_n1j = _mm_mul_pd(__r_nj, __r_nj);
		__q_nj = _mm_unpackhi_pd(__q_n1j, __q_n1j);
		
		nml_dscalar b1,b2;
		_mm_store_sd((double *)&b1, __q_nj);
		_mm_store_sd((double *)&b2, __q_n1j);
		betap += (b1+b2);
		
		__q_nj = _mm_load_pd((double *)&q_n[j]);
		_mm_store_pd((double *)&q_n1[j], __q_nj);
		
		_mm_store_pd((double *)&r_n[j], __r_nj);
	}

