#pragma message("INTEL_MODIFIED: Scale-loop modification (double) is ON")
_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
__m128d __y, __x, __tmp;
int max_count = d->Hzb.Natom * d->NBasisStates;
// #pragma unroll(4)
for (i = 0; i < max_count; i++) {
  __y = _mm_load_pd((double *)(&y[i].r));
#ifdef HamMemorySave // load as float
  __tmp = _mm_cvtps_pd(_mm_load_ss(&d->Hzb.Hdd[i]));
#else
  __tmp = _mm_load_sd(&d->Hzb.Hdd[i]);
#endif // Ham_Memory_save
  __x = _mm_load_pd((double *)(&x[i].r));
  //__tmp = _mm_shuffle_pd(__tmp, __tmp, 0);
  __tmp = _mm_unpacklo_pd(__tmp, __tmp);
  __x = _mm_mul_pd(__x, __tmp);
  __y = _mm_add_pd(__y, __x);

  _mm_store_pd((double *)(&y[i].r), __y);
}
