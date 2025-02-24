#pragma message("INTEL_MODIFIED: Lanczos loop 1 vectorization is ON")
   __m128d __inv_beta_n1, __r_nj, __q_nj;
     nml_dscalar inv_beta_n1 = 1.0/beta_n1;
   __inv_beta_n1 = _mm_load_sd((double*)&inv_beta_n1);
   __inv_beta_n1 = _mm_unpacklo_pd(__inv_beta_n1, __inv_beta_n1);
   for (j = 0; j < m; ++j){
	__r_nj    = _mm_load_pd((double*)&r_n[j+1]);
	__q_nj    = _mm_mul_pd(__r_nj, __inv_beta_n1);
	_mm_store_pd((double*)&q_n[j], __q_nj);
	}
