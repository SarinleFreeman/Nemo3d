#pragma message ("INTEL_MODIFIED: Lanczos loop 2 vectorization is ON")
	memcpy(q_n1, q_n, m*sizeof(nml_dcscalar));
	__m128d __r_nj, __q_nj, __q_ns,  __al, __beta_n1, __betap;
	__al = _mm_load_sd((double*) alpha);
	__al = _mm_unpacklo_pd(__al, __al);

	__beta_n1 = _mm_load_sd((double*) &beta_n1);
	__beta_n1 = _mm_unpacklo_pd(__beta_n1, __beta_n1);
	__betap   = _mm_setzero_pd();

    for (j = 0; j < m; ++j) {
	  __q_nj = _mm_load_pd((double*) &q_n[j]);
	  __q_ns = _mm_shuffle_pd(__q_nj, __q_nj, 01);
	  __r_nj = _mm_load_pd((double*) &r_n[j]);
	  __q_nj = _mm_mul_pd(__q_nj,  __al);
	  __r_nj = _mm_sub_pd(__r_nj, __q_nj);
	  __q_ns = _mm_mul_pd(__q_ns, __beta_n1);
	  __r_nj = _mm_sub_pd(__r_nj, __q_ns);
	  _mm_store_pd((double*)&r_n[j], __r_nj);
	  __r_nj    = _mm_mul_pd(__r_nj, __r_nj);
      __betap = _mm_add_pd(__betap, __r_nj);
      } 
	
    __q_nj = _mm_unpackhi_pd(__betap, __betap); // use q_nj as temp
	__betap = _mm_add_sd(__betap, __q_nj);
	nml_dscalar betap;
	_mm_store_sd(&betap, __betap);
