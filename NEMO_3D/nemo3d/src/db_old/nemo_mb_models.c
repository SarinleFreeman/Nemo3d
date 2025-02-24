/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This application is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)

Written by:  Chris Bowen
             Gerhard Klimeck
             Fabiano Oyafuso
             Seungwon Lee
             Olga Lazarenkova
             Hook Hua

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/db_old/nemo_mb_models.c,v 1.1 2004/02/16 03:36:56 gekco Exp $
*****************************************************************************/
#ifdef DATABASE_CODE_FRAGMENT_TO_BE_FIXED_LATER

#include "nemo_mb_models.h"

/* Creates paramters for multiband materials used in creating Hamiltonian. */



#ifndef LIMIT_THEORY

/*
 * Boolean check for whether if the bandstructure model is independent,
 * non-coupled single bands. 
 */
Boolean isIndependentSingleBand(char *bandModel) {

  if (!bandModel) {
      n3d_warning("In isIndependentSingleBand():: Error - Band Model is NULL.\n");
      return FALSE;
  }

  if (
      !strcmp( bandModel, Coken_BM_Bands_1           ) ||
      !strcmp( bandModel, Coken_BM_Bands_1_ac        ) ||
      !strcmp( bandModel, Coken_BM_Bands_2           ) ||
      !strcmp( bandModel, Coken_BM_Bands_2_ac        ) ||
      !strcmp( bandModel, Coken_BM_Bands_10          ) ||
      !strcmp( bandModel, Coken_BM_Bands_10_so       ) ||
      !strcmp( bandModel, Coken_BM_Bands_10_so_s     ) ||
      !strcmp( bandModel, Coken_BM_Bands_10_so_s_ac  ) ||
      !strcmp( bandModel, Coken_BM_Bands_10_so_s_new ) ||
      !strcmp( bandModel, Coken_BM_Bands_10_2N       ) ||
      !strcmp( bandModel, Coken_BM_Bands_10_2N_s     ) ||
      !strcmp( bandModel, Coken_BM_Bands_18_s        ) ||
      !strcmp( bandModel, Coken_BM_Bands_18_2_s      ) ||
      !strcmp( bandModel, Coken_BM_Bands_20_s        ) ||
      !strcmp( bandModel, Coken_BM_Bands_20_2_s      ) ||
      !strcmp( bandModel, Coken_BM_Bands_28_s        ) ||
      !strcmp( bandModel, Coken_BM_Bands_kp8         )
      ) {
      return FALSE;
  }
  else {
      return TRUE;
  }

}



/*
 * Sets the semi material properties based on the given band model.
 */
int setMBModelParam(
                      char  *band_model,
                      int   *param_n,           /* Number of associated multiband parameters. */
		      char **model,
                      int   *bdim,             /* Block dimension of the bulk Hamiltonian. */
                      int   *band_inc,         /* Number of incoming bands. */
		      int   *nbr,              /* Number of neighbors. */
                      int   *abdim,            /* Block dimension for the RGF algorithm. */
                      int   *spin_orbit,       /* spin orbit coupling */
                      int   *spin_degeneracy,
                      int   *sigma_singular    /* singular tight binding boundary condition. */
                     ) {
   printf("Executing set MBModelParam\n\n\n\n\n@@@@@@@@@@@@@@@\n@@@@@@@\n@@@@@exit(-1);\n");
    if (!strcmp(band_model,Coken_BM_Bands_2)) {
        *param_n          = 3;
        if (model) { *model = n3d_copy_str("ecev"); }
        *bdim            = 2;
        *band_inc        = 1;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 0;
        *spin_degeneracy = 1;
        *sigma_singular  = 1;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_10)) {
        *param_n          = 15;
        if (model) { *model = n3d_copy_str("sp3s"); }
        *bdim            = 10;
        *band_inc        = 4;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 0;
        *spin_degeneracy = 1;
        *sigma_singular  = 1;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_10_so)) {
        *param_n          = 15;
        if (model) { *model = n3d_copy_str("sp3s"); }
        *bdim            = 10;
        *band_inc        = 4;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 1;
        *spin_degeneracy = 1;
        *sigma_singular  = 1;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_10_so_s)) {
        *param_n          = 15;
        if (model) { *model = n3d_copy_str("sp3s_s"); }
        *bdim            = 20;
        *band_inc        = 8;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 1;
        *spin_degeneracy = 2;   /* GEKCO change 11/10/97 */
        *sigma_singular  = 1;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_10_so_s_new)) {
        *param_n          = 15;
        if (model) { *model = n3d_copy_str("sp3s_s_n"); }
        *bdim            = 20;
        *band_inc        = 8;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 1;
        *spin_degeneracy = 2;   /* GEKCO change 11/10/97 */
        *sigma_singular  = 1;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_10_2N)) {
        *param_n          = 37;
        if (model) { *model = n3d_copy_str("sp3s_2n"); }
        *bdim            = 10;
        *band_inc        = 4;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 1;
        *spin_degeneracy = 1;
        *sigma_singular  = 0;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_10_2N_s)) {
        *param_n          = 37;
        if (model) { *model = n3d_copy_str("sp3s_2n_s"); }
        *bdim            = 20;
        *band_inc        = 8;
        *nbr             = 1;
        *abdim           = *bdim / *nbr; 
        *spin_orbit      = 1;
        *spin_degeneracy = 2;
        *sigma_singular  = 0;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_18_s)) {
        *param_n          = 35;
        if (model) { *model = n3d_copy_str("spd_s"); }
        *bdim            = 36;
        *band_inc        = 8;
        *nbr             = 1;
        *abdim           = *bdim / *nbr; 
        *spin_orbit      = 1;
        *spin_degeneracy = 2;
        *sigma_singular  = 1;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_18_2_s)) {
        *param_n          = 35;
        if (model) { *model = n3d_copy_str("spd2_s"); }
        *bdim            = 36;
        *band_inc        = 28;
        *nbr             = 1;
        *abdim           = *bdim / *nbr; 
        *spin_orbit      = 1;
        *spin_degeneracy = 2;
        *sigma_singular  = 1;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_20_s)) {
        *param_n          = 35;
        if (model) { *model = n3d_copy_str("spds_s"); }
        *bdim            = 40;
        *band_inc        = 8;
        *nbr             = 1;
        *abdim           = *bdim / *nbr; 
        *spin_orbit      = 1;
        *spin_degeneracy = 2;
        *sigma_singular  = 1;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_20_2_s)) {
        *param_n          = 35;
        if (model) { *model = n3d_copy_str("spds2_s"); }
        *bdim            = 40;
        *band_inc        = 28;
        *nbr             = 1;
        *abdim           = *bdim / *nbr; 
        *spin_orbit      = 1;
        *spin_degeneracy = 2;
        *sigma_singular  = 1;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_28_s)) {
        *param_n          = 38;
        if (model) { *model = n3d_copy_str("spdd_s"); }
        *bdim            = 56;
        *band_inc        = 8;
        *nbr             = 1;
        *abdim           = *bdim / *nbr; 
        *spin_orbit      = 1;
        *spin_degeneracy = 2;
        *sigma_singular  = 1;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_kp8)) {
        *param_n          = 18;
        if (model) { *model = n3d_copy_str("kp8"); }
        *bdim            = 16;
        *band_inc        = 4;
        *nbr             = 1;
        *abdim           = *bdim / *nbr; 
        *spin_orbit      = 1;
        *spin_degeneracy = 2;
        *sigma_singular  = 0;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_1)) {
        *param_n          = 5;
        if (model) { *model = n3d_copy_str("sb"); }
        *bdim            = 1;
        *band_inc        = 0;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 0;
        *spin_degeneracy = 1;
        *sigma_singular  = 0;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_1_ac)) {
        *param_n          = 5;
        if (model) { *model = n3d_copy_str("sb_ac"); }
        *bdim            = 1*(1+2*SideBands);
        *band_inc        = 0 + 1*SideBands;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 0;
        *spin_degeneracy = 1;
        *sigma_singular  = 0;
    }
    else if (!strcmp(band_model,Coken_BM_Bands_2_ac)) {
        *param_n          = 3;
        if (model) { *model = n3d_copy_str("ecev_ac"); }
        *bdim            = 2*(1+2*SideBands);
        *band_inc        = 1 + 2*SideBands;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 0;
        *spin_degeneracy = 1;
#ifdef  SPLIT_HAMILTONIAN_AC
        *sigma_singular  = 1;
#else
#ifdef SIMPLE_HAMILTONIAN_AC
        *sigma_singular  = 0;
#else
        "Should not compile with proper #defines";
#endif /* SIMPLE_HAMILTONIAN_AC */
#endif /*  SPLIT_HAMILTONIAN_AC */
    }
    else if (!strcmp(band_model,Coken_BM_Bands_10_so_s_ac)) {
        *param_n          = 15;
	if (model) { *model = n3d_copy_str("sp3s_s_ac"); }
        *bdim            = 20*(1+2*SideBands);
        *band_inc        = 8+20*SideBands;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 1;
        *spin_degeneracy = 2;   /* GEKCO change 11/10/97 */
        *sigma_singular  = 0;
    }
    else {
        *param_n          = 3;
        if (model) { *model = n3d_copy_str("sb"); }
        *bdim            = 1;
        *band_inc        = 0;
        *nbr             = 1;
        *abdim           = *bdim / *nbr;
        *spin_orbit      = 0;
        *spin_degeneracy = 1;
        *sigma_singular  = 0;
    }

    return 0;
}



/*
 * Adjusts the pointer to the multiband parameter lists based on the mutliband model.
 */
semiMaterial adjustMBParamPointer(semiMaterial sm) {

    if (
        !strcmp( sm->band_model, Coken_BM_Bands_1           ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_1_ac        ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_10          ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_10_so       ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_10_so_s     ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_10_so_s_ac  ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_10_so_s_new ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_10_2N       ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_10_2N_s     ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_18_s        ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_18_2_s      ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_20_s        ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_20_2_s      ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_28_s        ) ||
        !strcmp( sm->band_model, Coken_BM_Bands_kp8         )
       ) {
        sm->mb_p = sm->mbn_p;
    }
    else {
        sm->mb_p = sm->mb2_p;
    }

    if (sm->mb_new_3d) {
        sm->mb_p = sm->mb_pv;
    }

    return sm;
}



/*
 * Initialize the semi material structure.
 */
semiMaterial initMBSemiMaterial(semiMaterial sm) {

    *sm->Abort_result = setMBModelParam(
                                          sm->band_model,
                                          &sm->mb_param_n, 
	  			          &sm->mb_model,
		  		          &sm->mb_bdim,
				          &sm->mb_band_inc, 
				          &sm->mb_nbr,
				          &sm->mb_abdim,
				          &sm->mb_spin_orbit,
				          &sm->mb_spin_degeneracy,
				          &sm->mb_sigma_singular
				         );
  
    sm = adjustMBParamPointer(sm);

    return sm;
}


#endif /* LIMIT_THEORY */
#endif  /* DATABASE_CODE_FRAGMENT_TO_BE_FIXED_LATER  */
