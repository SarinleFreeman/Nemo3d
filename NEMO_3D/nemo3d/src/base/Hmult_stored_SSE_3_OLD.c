#pragma message("INTEL_MODIFIED: Loop nest modification (double) is ON")
if (d->opt.ExecParam.Phys.MagneticFieldOn) {
  for (i = 0; i < Nb_2; i++) {
    int row = d->Hzb.indxHU_to_row[k] + i;
    complex *x_I = x + row;
    complex *y_I = y + row;
    complex *x_I2 = x + row + Nb_2;
    complex *y_I2 = y + row + Nb_2;
    for (j = 0; j < Nb_2; j++) {
      int col = d->Hzb.indxHU_to_col[k] + j;
      complex *x_J = x + col;
      complex *y_J = y + col;
      complex *x_J2 = x + col + Nb_2;
      complex *y_J2 = y + col + Nb_2;
      // inner if
      double Hu_indx_r = d->Hzb.Hu[indx] * d->Hzb.HphaseIn[k].r;
      double Hu_indx_i = d->Hzb.Hu[indx] * d->Hzb.HphaseIn[k].i;
      // do lower (upper?) triangular part
      y_I->r += Hu_indx_r * x_J->r - Hu_indx_i * x_J->i;
      y_I->i += Hu_indx_r * x_J->i + Hu_indx_i * x_J->r;
      y_I2->r += Hu_indx_r * x_J2->r - Hu_indx_i * x_J2->i;
      y_I2->i += Hu_indx_r * x_J2->i + Hu_indx_i * x_J2->r;
      // do upper (lower?) triangular part
      y_J->r += Hu_indx_r * x_I->r + Hu_indx_i * x_I->i;
      y_J->i += Hu_indx_r * x_I->i - Hu_indx_i * x_I->r;
      y_J2->r += Hu_indx_r * x_I2->r + Hu_indx_i * x_I2->i;
      y_J2->i += Hu_indx_r * x_I2->i - Hu_indx_i * x_I2->r;
      // end inner if
      indx++;
    }
  }
} else {
  // double tmp;
  __m128d tmp1, y_i, x_j, y_i2, x_j2, y_j, x_i, y_j2, x_i2, __mul;

  const int row = d->Hzb.indxHU_to_row[k];
  complex *x_I = x + row;
  complex *y_I = y + row;
  complex *x_I2 = x_I + Nb_2;
  complex *y_I2 = y_I + Nb_2;

  complex *x_J, *x_J2, *y_J, *y_J2;

  const int col = d->Hzb.indxHU_to_col[k];

  complex *x_J_start = x + col;
  complex *y_J_start = y + col;
  complex *x_J2_start = x + col + Nb_2;
  complex *y_J2_start = y + col + Nb_2;

  for (i = 0; i < Nb_2; i++) {

    x_J2 = x_J2_start;
    y_J = y_J_start;
    y_J2 = y_J2_start;

    x_i = _mm_load_pd((double *)(&x_I->r));
    y_i = _mm_load_pd((double *)(&y_I->r));
    x_i2 = _mm_load_pd((double *)(&x_I2->r));
    y_i2 = _mm_load_pd((double *)(&y_I2->r));

    for (complex *x_J = x + col; x_J < x + col + Nb_2; x_J++) {
      // for (j=0; j<Nb_2; j++) {

      // tmp = (double) d->Hzb.Hu[indx];

#ifdef HamMemorySave // load as float

      // tmp1 = _mm_cvtps_pd(_mm_load_ss((float*)&d->Hzb.Hu[indx]));
      tmp1 = _mm_cvtps_pd(_mm_load_ss((float *)__ptr));
#else
      // tmp1 = _mm_load_sd((double*)&d->Hzb.Hu[indx]);
      tmp1 = _mm_load_sd((double *)__ptr);
#endif // HamMemorySave

      x_j = _mm_load_pd((double *)(&x_J->r));

      y_j = _mm_load_pd((double *)(&y_J->r));

      tmp1 = _mm_unpacklo_pd(tmp1, tmp1);

      // y_J->r += tmp * x_I->r;
      // y_J->i += tmp * x_I->i;

      __mul = _mm_mul_pd(x_i, tmp1);
      y_j = _mm_add_pd(y_j, __mul);

      // y_I->r += tmp * x_J->r;
      // y_I->i += tmp * x_J->i;

      __mul = _mm_mul_pd(x_j, tmp1);
      y_i = _mm_add_pd(y_i, __mul);

      _mm_store_pd((double *)(&y_J->r), y_j);
      // y_I2->r += tmp * x_J2->r;
      // y_I2->i += tmp * x_J2->i;

      x_j2 = _mm_load_pd((double *)(&x_J2->r));
      __mul = _mm_mul_pd(x_j2, tmp1);
      y_i2 = _mm_add_pd(y_i2, __mul);

      // y_J2->r += tmp * x_I2->r;
      // y_J2->i += tmp * x_I2->i;
      y_j2 = _mm_load_pd((double *)(&y_J2->r));
      __mul = _mm_mul_pd(x_i2, tmp1);
      y_j2 = _mm_add_pd(y_j2, __mul);
      _mm_store_pd((double *)(&y_J2->r), y_j2);

      // indx++;
      __ptr++;
      // x_J++;
      y_J++;
      y_J2++;
      x_J2++;
    } // j
    y_J--;
    y_J2--;

    _mm_store_pd((double *)(&y_I->r), y_i);
    _mm_store_pd((double *)(&y_I2->r), y_i2);
    x_I++;
    y_I++;
    y_I2++;
    x_I2++;
  } // i

  y_I--;
  y_I2--;
  y_J--;
  y_J2--;
}
