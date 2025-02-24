#pragma message("INTEL_MODIFIED: Scale-loop modification (Vmem) is ON")
  __m128d __y, __x, __tmp;
//  int max_count = atom<d->Hzb.Natom * d->NBasisStates;
  int max_count = d->Hzb.Natom * d->NBasisStates;
  for (i=0; i<max_count; i++){
	__tmp = _mm_cvtps_pd(_mm_load_ss(&d->Hzb.Hdd[i]));
	//__tmp = _mm_shuffle_pd(__tmp, __tmp, 0);
	__tmp = _mm_unpacklo_pd(__tmp, __tmp);

	__x = _mm_load_pd((double*)(&x[i].r));
	__x = _mm_mul_pd(__x, __tmp);

	__y = _mm_load_pd((double*)(&y[i].r));
	__y = _mm_add_pd(__y, __x);

	_mm_store_pd((double*)(&y[i].r), __y);
    
  }
