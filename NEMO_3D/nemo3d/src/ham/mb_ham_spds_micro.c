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
$Header: /repo/nemo3d/src/ham/mb_ham_spds_micro.c,v 1.44 2009/04/27 22:22:46
nkharche Exp $
*****************************************************************************/

#include "mb_ham_spds_micro.h"

#define sqr3 1.73205080757          /* sqrt(3.0) */
#define sqr2 1.41421356237          /* sqrt(2.0) */
#define BohrMagneton 5.788381804e-5 /* in unit of eV/T */
#define MagneticFluxQuantum                                                    \
  2.06783372e3 /* (h/2e) in unit of Tesla*(nanometer)^2 */

double sp3d5s_Ebound = sp3d5s_Ebound_undefine;

/*
 * off diagonal coupling parameters
 * These parameters reference an array that
 * may include the strain induced modifications
 * according to a power law scaling with the distance.
 * Note that there are some minus signs built into these definitions
 */

#define sss vpar_e_const[pV_sss]
#define ststs vpar_e_const[pV_ststs]
#define stass vpar_e_const[pV_stass]
#define sasts vpar_e_const[pV_sasts]
#define sapcs vpar_e_const[pV_sapcs]
#define scpas (-vpar_e_const[pV_scpas])
#define stapcs vpar_e_const[pV_stapcs]
#define stcpas (-vpar_e_const[pV_stcpas])
#define sadcs vpar_e_const[pV_sadcs]
#define scdas vpar_e_const[pV_scdas]
#define stadcs vpar_e_const[pV_stadcs]
#define stcdas vpar_e_const[pV_stcdas]
#define pps vpar_e_const[pV_pps]
#define ppp vpar_e_const[pV_ppp]
#define padcs vpar_e_const[pV_padcs]
#define pcdas (-vpar_e_const[pV_pcdas])
#define padcp vpar_e_const[pV_padcp]
#define pcdap (-vpar_e_const[pV_pcdap])
#define dds vpar_e_const[pV_dds]
#define ddp vpar_e_const[pV_ddp]
#define ddd vpar_e_const[pV_ddd]

/*
   These definitions allow for the simple assembly of the
   sums of diagonal elements that are
   needed in the assembly of the on-diagonal
   corrections due to strain.
*/
#define pVpE_sa vpar_orig[pE_sa]
#define pVpE_pa vpar_orig[pE_pa]
#define pVpE_sc vpar_orig[pE_sc]
#define pVpE_pc vpar_orig[pE_pc]
#define pVpE_sta vpar_orig[pE_sta]
#define pVpE_stc vpar_orig[pE_stc]
#define pVpE_da vpar_orig[pE_da]
#define pVpE_dc vpar_orig[pE_dc]

#define Additional_Diagonal_constant
#ifdef Additional_Diagonal_constant
#define c_ss strain_vec[c_ind_ss]
#define c_stst strain_vec[c_ind_stst]
#define c_stas strain_vec[c_ind_stas]
#define c_sast strain_vec[c_ind_sast]
#define c_sapc strain_vec[c_ind_sapc]
#define c_pasc strain_vec[c_ind_pasc]
#define c_scpa c_pasc
#define c_stapc strain_vec[c_ind_stapc]
#define c_pastc strain_vec[c_ind_pastc]
#define c_stcpa c_pastc
#define c_sadc strain_vec[c_ind_sadc]
#define c_dasc strain_vec[c_ind_dasc]
#define c_scda c_dasc
#define c_stadc strain_vec[c_ind_stadc]
#define c_dastc strain_vec[c_ind_dastc]
#define c_stcda c_dastc
#define c_pp strain_vec[c_ind_pp]
#define c_padc strain_vec[c_ind_padc]
#define c_dapc strain_vec[c_ind_dapc]
#define c_pcda c_dapc
#define c_dd strain_vec[c_ind_dd]

#define pVpE_ss (pVpE_sa + pVpE_sc) * c_ss
#define pVpE_pp (pVpE_pa + pVpE_pc) * c_pp
#define pVpE_dd (pVpE_da + pVpE_dc) * c_dd
#define pVpE_stst (pVpE_sta + pVpE_stc) * c_stst

#define pVpE_sapc (pVpE_sa + pVpE_pc) * c_sapc
#define pVpE_sadc (pVpE_sa + pVpE_dc) * c_sadc
#define pVpE_sast (pVpE_sa + pVpE_stc) * c_sast

#define pVpE_pasc (pVpE_pa + pVpE_sc) * c_pasc
#define pVpE_padc (pVpE_pa + pVpE_dc) * c_padc
#define pVpE_pastc (pVpE_pa + pVpE_stc) * c_pastc
#define pVpE_scpa pVpE_pasc
#define pVpE_stcpa pVpE_pastc

#define pVpE_dasc (pVpE_da + pVpE_sc) * c_dasc
#define pVpE_dapc (pVpE_da + pVpE_pc) * c_dapc
#define pVpE_dastc (pVpE_da + pVpE_stc) * c_dastc
#define pVpE_scda pVpE_dasc
#define pVpE_pcda pVpE_dapc
#define pVpE_stcda pVpE_dastc

#define pVpE_stas (pVpE_sta + pVpE_sc) * c_stas
#define pVpE_stapc (pVpE_sta + pVpE_pc) * c_stapc
#define pVpE_stadc (pVpE_sta + pVpE_dc) * c_stadc

#else
#define pVpE_ss (pVpE_sa + pVpE_sc)
#define pVpE_pp (pVpE_pa + pVpE_pc)
#define pVpE_dd (pVpE_da + pVpE_dc)
#define pVpE_stst (pVpE_sta + pVpE_stc)

#define pVpE_sapc (pVpE_sa + pVpE_pc)
#define pVpE_sadc (pVpE_sa + pVpE_dc)
#define pVpE_sast (pVpE_sa + pVpE_stc)

#define pVpE_pasc (pVpE_pa + pVpE_sc)
#define pVpE_padc (pVpE_pa + pVpE_dc)
#define pVpE_pastc (pVpE_pa + pVpE_stc)
#define pVpE_scpa pVpE_pasc
#define pVpE_stcpa pVpE_pastc

#define pVpE_dasc (pVpE_da + pVpE_sc)
#define pVpE_dapc (pVpE_da + pVpE_pc)
#define pVpE_dastc (pVpE_da + pVpE_stc)
#define pVpE_scda pVpE_dasc
#define pVpE_pcda pVpE_dapc
#define pVpE_stcda pVpE_dastc

#define pVpE_stas (pVpE_sta + pVpE_sc)
#define pVpE_stapc (pVpE_sta + pVpE_pc)
#define pVpE_stadc (pVpE_sta + pVpE_dc)

#endif
/*
  These definitions are used for the fast access of the
  diagonal matrix element corrections due to strain.
  */
#define pVp_sss vpar_orig[pV_sss]
#define pVp_ststs vpar_orig[pV_ststs]
#define pVp_stass vpar_orig[pV_stass]
#define pVp_sasts vpar_orig[pV_sasts]
#define pVp_sapcs vpar_orig[pV_sapcs]
#define pVp_scpas vpar_orig[pV_scpas]
#define pVp_stapcs vpar_orig[pV_stapcs]
#define pVp_stcpas vpar_orig[pV_stcpas]
#define pVp_sadcs vpar_orig[pV_sadcs]
#define pVp_scdas vpar_orig[pV_scdas]
#define pVp_stadcs vpar_orig[pV_stadcs]
#define pVp_stcdas vpar_orig[pV_stcdas]
#define pVp_pps vpar_orig[pV_pps]
#define pVp_ppp vpar_orig[pV_ppp]
#define pVp_padcs vpar_orig[pV_padcs]
#define pVp_pcdas vpar_orig[pV_pcdas]
#define pVp_padcp vpar_orig[pV_padcp]
#define pVp_pcdap vpar_orig[pV_pcdap]
#define pVp_dds vpar_orig[pV_dds]
#define pVp_ddp vpar_orig[pV_ddp]
#define pVp_ddd vpar_orig[pV_ddd]

static int diag_map[20];

static int diag_map_unset = 1;
/* Establish a static vector used to map the parameter list
   into the diagonal elements needed for the Hamiltonian.  */
void diag_map_set(void) {
  diag_map_unset = 0;
  diag_map[STa] = pE_sta;
  diag_map[Sa] = pE_sa;
  diag_map[Pxa] = pE_pa;
  diag_map[Pya] = pE_pa;
  diag_map[Pza] = pE_pa;
  diag_map[Dxya] = pE_da;
  diag_map[Dyza] = pE_da;
  diag_map[Dzxa] = pE_da;
  diag_map[Dx2my2a] = pE_da;
  diag_map[Dz2a] = pE_da;
  diag_map[STc] = pE_stc;
  diag_map[Sc] = pE_sc;
  diag_map[Pxc] = pE_pc;
  diag_map[Pyc] = pE_pc;
  diag_map[Pzc] = pE_pc;
  diag_map[Dxyc] = pE_dc;
  diag_map[Dyzc] = pE_dc;
  diag_map[Dzxc] = pE_dc;
  diag_map[Dx2my2c] = pE_dc;
  diag_map[Dz2c] = pE_dc;
  return;
}

real E_const_get_old(rvectr vpar_e_const, int basis1_econst, int basis2_econst,
                     real l, real m, real n, real l2, real m2, real n2) {
  real Out_econst;

  Out_econst = 0.0;

#include "mb_ham_spds_E_const.h"

  return Out_econst;
}

#ifdef Additional_Diagonal_constant
real E_const_get(rvectr vpar_e_const, rvectr vpar_orig, rvectr strain_vec,
                 int basis1_econst, int basis2_econst, real l, real m, real n,
                 real l2, real m2, real n2, real l0, real m0, real n0, real l02,
                 real m02, real n02, real cstrain, real *out_diag_corr)
#else
real E_const_get(rvectr vpar_e_const, rvectr vpar_orig, int basis1_econst,
                 int basis2_econst, real l, real m, real n, real l2, real m2,
                 real n2, real l0, real m0, real n0, real l02, real m02,
                 real n02, real cstrain, real *out_diag_corr)
#endif
{
  real Out_econst;
  real r1, r2, r3;
  real Out_diag_corr = *out_diag_corr;

  Out_econst = 0.0;

  /* This must be a switch statement.
     Replacing this with get_*() macros also means rewriting ei() and eidk()
     below, which means expanding the loops in tbh_const_bulk(), tbh_1d_nemo(),
     and tbh_const_bulk_dk(), below.  It also probably implies rewriting
     spin_orbit_spds().
  */
#include "mb_ham_spds_inclVoffdEdiag.h"

  *out_diag_corr = Out_diag_corr;

  return Out_econst;
}
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

real diag(rvectr param, int basis) {
#define NEW_DIAG
#ifdef NEW_DIAG
  if (diag_map_unset)
    diag_map_set();
  return param[diag_map[basis]];
#else
  real E;

  switch (basis) {
  case STa:
    E = param[4];
    break;
  case Sa:
    E = param[0];
    break;
  case Pxa:
    E = param[1];
    break;
  case Pya:
    E = param[1];
    break;
  case Pza:
    E = param[1];
    break;
  case Dxya:
    E = param[6];
    break;
  case Dyza:
    E = param[6];
    break;
  case Dzxa:
    E = param[6];
    break;
  case Dx2my2a:
    E = param[6];
    break;
  case Dz2a:
    E = param[6];
    break;
  case STc:
    E = param[5];
    break;
  case Sc:
    E = param[2];
    break;
  case Pxc:
    E = param[3];
    break;
  case Pyc:
    E = param[3];
    break;
  case Pzc:
    E = param[3];
    break;
  case Dxyc:
    E = param[7];
    break;
  case Dyzc:
    E = param[7];
    break;
  case Dzxc:
    E = param[7];
    break;
  case Dx2my2c:
    E = param[7];
    break;
  case Dz2c:
    E = param[7];
  }

  return E;
#endif
}

complex spin_orbit_spds(int basis1, int basis2, int spin1, int spin2,
                        real delta_a, real delta_c) {
  complex E;

  E.r = E.i = 0.0;

  if (spin1 == 1 && spin2 == 0) {
    switch (basis1) {
    case Pxa:
      switch (basis2) {
      case Pza:
        E.r = delta_a;
        break;
      }
      break;
    case Pya:
      switch (basis2) {
      case Pza:
        E.i = -delta_a;
        break;
      }
      break;
    case Pza:
      switch (basis2) {
      case Pxa:
        E.r = -delta_a;
        break;
      case Pya:
        E.i = delta_a;
        break;
      }
      break;
    case Pxc:
      switch (basis2) {
      case Pzc:
        E.r = delta_c;
        break;
      }
      break;
    case Pyc:
      switch (basis2) {
      case Pzc:
        E.i = -delta_c;
        break;
      }
      break;
    case Pzc:
      switch (basis2) {
      case Pxc:
        E.r = -delta_c;
        break;
      case Pyc:
        E.i = delta_c;
        break;
      }
    }
  }
  if (spin1 == 0 && spin2 == 1) {
    switch (basis1) {
    case Pxa:
      switch (basis2) {
      case Pza:
        E.r = -delta_a;
        break;
      }
      break;
    case Pya:
      switch (basis2) {
      case Pza:
        E.i = -delta_a;
        break;
      }
      break;
    case Pza:
      switch (basis2) {
      case Pxa:
        E.r = delta_a;
        break;
      case Pya:
        E.i = delta_a;
        break;
      }
      break;
    case Pxc:
      switch (basis2) {
      case Pzc:
        E.r = -delta_c;
        break;
      }
      break;
    case Pyc:
      switch (basis2) {
      case Pzc:
        E.i = -delta_c;
        break;
      }
      break;
    case Pzc:
      switch (basis2) {
      case Pxc:
        E.r = delta_c;
        break;
      case Pyc:
        E.i = delta_c;
        break;
      }
    }
  }
  if (spin1 == 1 && spin2 == 1) {
    switch (basis1) {
    case Pxa:
      switch (basis2) {
      case Pya:
        E.i = -delta_a;
        break;
      }
      break;
    case Pya:
      switch (basis2) {
      case Pxa:
        E.i = delta_a;
        break;
      }
      break;
    case Pxc:
      switch (basis2) {
      case Pyc:
        E.i = -delta_c;
        break;
      }
      break;
    case Pyc:
      switch (basis2) {
      case Pxc:
        E.i = delta_c;
        break;
      }
    }
  }
  if (spin1 == 0 && spin2 == 0) {
    switch (basis1) {
    case Pxa:
      switch (basis2) {
      case Pya:
        E.i = delta_a;
        break;
      }
      break;
    case Pya:
      switch (basis2) {
      case Pxa:
        E.i = -delta_a;
        break;
      }
      break;
    case Pxc:
      switch (basis2) {
      case Pyc:
        E.i = delta_c;
        break;
      }
      break;
    case Pyc:
      switch (basis2) {
      case Pxc:
        E.i = -delta_c;
        break;
      }
    }
  }

  return E;
}

int h3d_offdiag(cmatrix h, rvectr nnv, ivectr basis0, ivectr basisn,
                ivectr spin, rvectr param, rvectr mb_strain,
                real unstrnd_cubic_cell_length, int nb2) {
  int i, j;
  static rvectr param_strained = NULL;
  rvectr param_keep = NULL;

  static real l, m, n, l2, m2, n2;
  if (param_strained == NULL) {
    param_strained = Rvectr(sp3d5s_nparam);
  } else {
    param_strained = check_rvectr(param_strained, sp3d5s_nparam);
  }

  param_keep = param_strained;

  E_const_prep(param, mb_strain, param_strained, unstrnd_cubic_cell_length,
               nnv[0], nnv[1], nnv[2], &l, &m, &n, &l2, &m2, &n2);

  if (mb_strain == NULL) {
    param_strained = param;
  }

  /* printf("\n");
  for ( i=0; i < sp3d5s_nparam; i++ )
  printf("\n%1.2f  %1.2f", param[i], param_strained[i]); */

  for (i = 0; i < nb2; i++)
    for (j = 0; j < nb2; j++)
      h[i][j].r = h[i][j].i = 0.0;

  /* Construct Hamiltonian */
  for (i = 0; i < nb2; i++) {
    for (j = 0; j < nb2; j++) {
      /*  Off-diagonal coupling */
      if (spin[i] == spin[j]) {
#define INLINE_E_CONST
#ifdef INLINE_E_CONST
#define basis1_econst basis0[i]
#define basis2_econst basisn[j]
#define vpar_e_const param_strained
#define Out_econst h[i][j].r
#define Out_diag_corr junk
#define vpar_orig param
#define vpar_orig param

/* #include "mb_ham_spds_inclVoffdEdiag.h" */
#include "mb_ham_spds_E_const.h"
#undef vpar_orig
#undef vpar_orig
#undef Out_diag_corr
#undef basis1_econst
#undef basis2_econst
#undef vpar_e_const
#undef Out_econst
#else
        h[i][j].r = E_const_get(param_strained, basis0[i], basisn[j], l, m, n,
                                l2, m2, n2);
#endif
      }
    }
  }

  if (param_strained != param_keep)
    param_strained = param_keep;

  return (0);
}

#include "mb_ham_spds_micro2.h" /* BM_Type enum */

int h3d_offdiag_diag_corr(cmatrix ho, cmatrix hd, rvectr nnv, rvectr nnv0,
                          ivectr basis_this, ivectr basis_nbr, ivectr spin,
                          rvectr param, rvectr mb_strain,
                          real unstrnd_cubic_cell_length, int Nbasis) {
  int i, j;
  static rvectr param_strained = NULL;
  rvectr param_keep = NULL;
  real r1, r2, r3;

  real cstrain = 0.0;

  static real l, m, n, l2, m2, n2;
  static real l0, m0, n0, l02, m02, n02;
#ifdef Additional_Diagonal_constant
  static rvectr mb_strain_temp = NULL;
  int set_temp = 0;
#endif
  register nml_double E;
  BM_Type BandModel = BM_1_s_nospin;

  if (param_strained == NULL) {
    param_strained = Rvectr(sp3d5s_nparam);
  } else {
    param_strained = check_rvectr(param_strained, sp3d5s_nparam);
  }

  if (mb_strain)
    cstrain = mb_strain[eta_const];

  param_keep = param_strained;

#ifdef Additional_Diagonal_constant
  if (!mb_strain) {
    if (mb_strain_temp == NULL) {
      mb_strain_temp = Rvectr(sp3d5s_nstrain);
    } else {
      mb_strain_temp = check_rvectr(mb_strain_temp, sp3d5s_nstrain);
    }
    mb_strain = mb_strain_temp;
    set_temp = 1;
  }
#endif

  /* returns (l,m,n) = (dx/dr, dy/dr, dz/dr), equilibrium
     distances (l0,m0,n0) and also their magnitude squared
     (l2,m2,n2) and (l02,m02,n02)
  Also, return strained parameters */
  strain_corr_prep(param, mb_strain, param_strained, unstrnd_cubic_cell_length,
                   nnv[0], nnv[1], nnv[2], nnv0[0], nnv0[1], nnv0[2], &l, &m,
                   &n, &l2, &m2, &n2, &l0, &m0, &n0, &l02, &m02, &n02);

  /* Determine the band model.  The values of Nbasis and spin[0] tested here
     are set in basis_init() in src/base/ham3d_final.c.  This is a hack; the
     proper way to do this is to simply pass the BandModel parameter of the
     current qd_struct to this function, which involves adding a parameter
     to every call to this function in the code. */
  switch (Nbasis) {
  case 10:
    switch (spin[0]) {
    case 1:
      BandModel = BM_10_sp3ss_spin;
      break;
    case -1:
      BandModel = BM_10_sp3d5ss_nospin;
      break;
    }
    break;
  case 20: /* sp3d5ss_spin */
    BandModel = BM_20_sp3d5ss_spin;
    break;
  }

  /*
   * Update the Hamiltonian.
   */

#define vpar_e_const param_strained
#define vpar_orig param
#define strain_vec mb_strain

#include "mb_ham_spds_VoffdEdiag_constructors.h"

  switch (BandModel) {
  case BM_10_sp3d5ss_nospin:
  case BM_20_sp3d5ss_spin:
    /* Explicitly list the calculations necessary for sp3d5s* spin up
       (0 <= row, col < 10) interactions. */
    if (basis_this[0] < STc) /* anion */
    {
      /* STa */
      E = hd[0][0].r;
      get_STa_STc(ho[0][0].r, E);
      get_STa_Sc(ho[1][0].r, E);
      get_STa_Pxc(ho[2][0].r, E);
      get_STa_Pyc(ho[3][0].r, E);
      get_STa_Pzc(ho[4][0].r, E);
      get_STa_Dxyc(ho[5][0].r, E);
      get_STa_Dyzc(ho[6][0].r, E);
      get_STa_Dzxc(ho[7][0].r, E);
      get_STa_Dx2my2c(ho[8][0].r, E);
      get_STa_Dz2c(ho[9][0].r, E);
      hd[0][0].r = E;
      /* Sa */
      E = hd[1][1].r;
      get_Sa_STc(ho[0][1].r, E);
      get_Sa_Sc(ho[1][1].r, E);
      get_Sa_Pxc(ho[2][1].r, E);
      get_Sa_Pyc(ho[3][1].r, E);
      get_Sa_Pzc(ho[4][1].r, E);
      get_Sa_Dxyc(ho[5][1].r, E);
      get_Sa_Dyzc(ho[6][1].r, E);
      get_Sa_Dzxc(ho[7][1].r, E);
      get_Sa_Dx2my2c(ho[8][1].r, E);
      get_Sa_Dz2c(ho[9][1].r, E);
      hd[1][1].r = E;
      /* Pxa */
      E = hd[2][2].r;
      get_Pxa_STc(ho[0][2].r, E);
      get_Pxa_Sc(ho[1][2].r, E);
      get_Pxa_Pxc(ho[2][2].r, E);
      get_Pxa_Pyc(ho[3][2].r, E);
      get_Pxa_Pzc(ho[4][2].r, E);
      get_Pxa_Dxyc(ho[5][2].r, E);
      get_Pxa_Dyzc(ho[6][2].r, E);
      get_Pxa_Dzxc(ho[7][2].r, E);
      get_Pxa_Dx2my2c(ho[8][2].r, E);
      get_Pxa_Dz2c(ho[9][2].r, E);
      hd[2][2].r = E;
      /* Pya */
      E = hd[3][3].r;
      get_Pya_STc(ho[0][3].r, E);
      get_Pya_Sc(ho[1][3].r, E);
      get_Pya_Pxc(ho[2][3].r, E);
      get_Pya_Pyc(ho[3][3].r, E);
      get_Pya_Pzc(ho[4][3].r, E);
      get_Pya_Dxyc(ho[5][3].r, E);
      get_Pya_Dyzc(ho[6][3].r, E);
      get_Pya_Dzxc(ho[7][3].r, E);
      get_Pya_Dx2my2c(ho[8][3].r, E);
      get_Pya_Dz2c(ho[9][3].r, E);
      hd[3][3].r = E;
      /* Pza */
      E = hd[4][4].r;
      get_Pza_STc(ho[0][4].r, E);
      get_Pza_Sc(ho[1][4].r, E);
      get_Pza_Pxc(ho[2][4].r, E);
      get_Pza_Pyc(ho[3][4].r, E);
      get_Pza_Pzc(ho[4][4].r, E);
      get_Pza_Dxyc(ho[5][4].r, E);
      get_Pza_Dyzc(ho[6][4].r, E);
      get_Pza_Dzxc(ho[7][4].r, E);
      get_Pza_Dx2my2c(ho[8][4].r, E);
      get_Pza_Dz2c(ho[9][4].r, E);
      hd[4][4].r = E;
      /* Dxya */
      E = hd[5][5].r;
      get_Dxya_STc(ho[0][5].r, E);
      get_Dxya_Sc(ho[1][5].r, E);
      get_Dxya_Pxc(ho[2][5].r, E);
      get_Dxya_Pyc(ho[3][5].r, E);
      get_Dxya_Pzc(ho[4][5].r, E);
      get_Dxya_Dxyc(ho[5][5].r, E);
      get_Dxya_Dyzc(ho[6][5].r, E);
      get_Dxya_Dzxc(ho[7][5].r, E);
      get_Dxya_Dx2my2c(ho[8][5].r, E);
      get_Dxya_Dz2c(ho[9][5].r, E);
      hd[5][5].r = E;
      /* Dyza */
      E = hd[6][6].r;
      get_Dyza_STc(ho[0][6].r, E);
      get_Dyza_Sc(ho[1][6].r, E);
      get_Dyza_Pxc(ho[2][6].r, E);
      get_Dyza_Pyc(ho[3][6].r, E);
      get_Dyza_Pzc(ho[4][6].r, E);
      get_Dyza_Dxyc(ho[5][6].r, E);
      get_Dyza_Dyzc(ho[6][6].r, E);
      get_Dyza_Dzxc(ho[7][6].r, E);
      get_Dyza_Dx2my2c(ho[8][6].r, E);
      get_Dyza_Dz2c(ho[9][6].r, E);
      hd[6][6].r = E;
      /* Dzxa */
      E = hd[7][7].r;
      get_Dzxa_STc(ho[0][7].r, E);
      get_Dzxa_Sc(ho[1][7].r, E);
      get_Dzxa_Pxc(ho[2][7].r, E);
      get_Dzxa_Pyc(ho[3][7].r, E);
      get_Dzxa_Pzc(ho[4][7].r, E);
      get_Dzxa_Dxyc(ho[5][7].r, E);
      get_Dzxa_Dyzc(ho[6][7].r, E);
      get_Dzxa_Dzxc(ho[7][7].r, E);
      get_Dzxa_Dx2my2c(ho[8][7].r, E);
      get_Dzxa_Dz2c(ho[9][7].r, E);
      hd[7][7].r = E;
      /* Dx2my2a */
      E = hd[8][8].r;
      get_Dx2my2a_STc(ho[0][8].r, E);
      get_Dx2my2a_Sc(ho[1][8].r, E);
      get_Dx2my2a_Pxc(ho[2][8].r, E);
      get_Dx2my2a_Pyc(ho[3][8].r, E);
      get_Dx2my2a_Pzc(ho[4][8].r, E);
      get_Dx2my2a_Dxyc(ho[5][8].r, E);
      get_Dx2my2a_Dyzc(ho[6][8].r, E);
      get_Dx2my2a_Dzxc(ho[7][8].r, E);
      get_Dx2my2a_Dx2my2c(ho[8][8].r, E);
      get_Dx2my2a_Dz2c(ho[9][8].r, E);
      hd[8][8].r = E;
      /* Dz2a */
      E = hd[9][9].r;
      get_Dz2a_STc(ho[0][9].r, E);
      get_Dz2a_Sc(ho[1][9].r, E);
      get_Dz2a_Pxc(ho[2][9].r, E);
      get_Dz2a_Pyc(ho[3][9].r, E);
      get_Dz2a_Pzc(ho[4][9].r, E);
      get_Dz2a_Dxyc(ho[5][9].r, E);
      get_Dz2a_Dyzc(ho[6][9].r, E);
      get_Dz2a_Dzxc(ho[7][9].r, E);
      get_Dz2a_Dx2my2c(ho[8][9].r, E);
      get_Dz2a_Dz2c(ho[9][9].r, E);
      hd[9][9].r = E;

      /* Repeate same calculations for spin down */
      if (BM_20_sp3d5ss_spin == BandModel) {
        /* STa */
        E = hd[10][10].r;
        get_STa_STc(ho[10][10].r, E);
        get_STa_Sc(ho[11][10].r, E);
        get_STa_Pxc(ho[12][10].r, E);
        get_STa_Pyc(ho[13][10].r, E);
        get_STa_Pzc(ho[14][10].r, E);
        get_STa_Dxyc(ho[15][10].r, E);
        get_STa_Dyzc(ho[16][10].r, E);
        get_STa_Dzxc(ho[17][10].r, E);
        get_STa_Dx2my2c(ho[18][10].r, E);
        get_STa_Dz2c(ho[19][10].r, E);
        hd[10][10].r = E;
        /* Sa */
        E = hd[11][11].r;
        get_Sa_STc(ho[10][11].r, E);
        get_Sa_Sc(ho[11][11].r, E);
        get_Sa_Pxc(ho[12][11].r, E);
        get_Sa_Pyc(ho[13][11].r, E);
        get_Sa_Pzc(ho[14][11].r, E);
        get_Sa_Dxyc(ho[15][11].r, E);
        get_Sa_Dyzc(ho[16][11].r, E);
        get_Sa_Dzxc(ho[17][11].r, E);
        get_Sa_Dx2my2c(ho[18][11].r, E);
        get_Sa_Dz2c(ho[19][11].r, E);
        hd[11][11].r = E;
        /* Pxa */
        E = hd[12][12].r;
        get_Pxa_STc(ho[10][12].r, E);
        get_Pxa_Sc(ho[11][12].r, E);
        get_Pxa_Pxc(ho[12][12].r, E);
        get_Pxa_Pyc(ho[13][12].r, E);
        get_Pxa_Pzc(ho[14][12].r, E);
        get_Pxa_Dxyc(ho[15][12].r, E);
        get_Pxa_Dyzc(ho[16][12].r, E);
        get_Pxa_Dzxc(ho[17][12].r, E);
        get_Pxa_Dx2my2c(ho[18][12].r, E);
        get_Pxa_Dz2c(ho[19][12].r, E);
        hd[12][12].r = E;
        /* Pya */
        E = hd[13][13].r;
        get_Pya_STc(ho[10][13].r, E);
        get_Pya_Sc(ho[11][13].r, E);
        get_Pya_Pxc(ho[12][13].r, E);
        get_Pya_Pyc(ho[13][13].r, E);
        get_Pya_Pzc(ho[14][13].r, E);
        get_Pya_Dxyc(ho[15][13].r, E);
        get_Pya_Dyzc(ho[16][13].r, E);
        get_Pya_Dzxc(ho[17][13].r, E);
        get_Pya_Dx2my2c(ho[18][13].r, E);
        get_Pya_Dz2c(ho[19][13].r, E);
        hd[13][13].r = E;
        /* Pza */
        E = hd[14][14].r;
        get_Pza_STc(ho[10][14].r, E);
        get_Pza_Sc(ho[11][14].r, E);
        get_Pza_Pxc(ho[12][14].r, E);
        get_Pza_Pyc(ho[13][14].r, E);
        get_Pza_Pzc(ho[14][14].r, E);
        get_Pza_Dxyc(ho[15][14].r, E);
        get_Pza_Dyzc(ho[16][14].r, E);
        get_Pza_Dzxc(ho[17][14].r, E);
        get_Pza_Dx2my2c(ho[18][14].r, E);
        get_Pza_Dz2c(ho[19][14].r, E);
        hd[14][14].r = E;
        /* Dxya */
        E = hd[15][15].r;
        get_Dxya_STc(ho[10][15].r, E);
        get_Dxya_Sc(ho[11][15].r, E);
        get_Dxya_Pxc(ho[12][15].r, E);
        get_Dxya_Pyc(ho[13][15].r, E);
        get_Dxya_Pzc(ho[14][15].r, E);
        get_Dxya_Dxyc(ho[15][15].r, E);
        get_Dxya_Dyzc(ho[16][15].r, E);
        get_Dxya_Dzxc(ho[17][15].r, E);
        get_Dxya_Dx2my2c(ho[18][15].r, E);
        get_Dxya_Dz2c(ho[19][15].r, E);
        hd[15][15].r = E;
        /* Dyza */
        E = hd[16][16].r;
        get_Dyza_STc(ho[10][16].r, E);
        get_Dyza_Sc(ho[11][16].r, E);
        get_Dyza_Pxc(ho[12][16].r, E);
        get_Dyza_Pyc(ho[13][16].r, E);
        get_Dyza_Pzc(ho[14][16].r, E);
        get_Dyza_Dxyc(ho[15][16].r, E);
        get_Dyza_Dyzc(ho[16][16].r, E);
        get_Dyza_Dzxc(ho[17][16].r, E);
        get_Dyza_Dx2my2c(ho[18][16].r, E);
        get_Dyza_Dz2c(ho[19][16].r, E);
        hd[16][16].r = E;
        /* Dzxa */
        E = hd[17][17].r;
        get_Dzxa_STc(ho[10][17].r, E);
        get_Dzxa_Sc(ho[11][17].r, E);
        get_Dzxa_Pxc(ho[12][17].r, E);
        get_Dzxa_Pyc(ho[13][17].r, E);
        get_Dzxa_Pzc(ho[14][17].r, E);
        get_Dzxa_Dxyc(ho[15][17].r, E);
        get_Dzxa_Dyzc(ho[16][17].r, E);
        get_Dzxa_Dzxc(ho[17][17].r, E);
        get_Dzxa_Dx2my2c(ho[18][17].r, E);
        get_Dzxa_Dz2c(ho[19][17].r, E);
        hd[17][17].r = E;
        /* Dx2my2a */
        E = hd[18][18].r;
        get_Dx2my2a_STc(ho[10][18].r, E);
        get_Dx2my2a_Sc(ho[11][18].r, E);
        get_Dx2my2a_Pxc(ho[12][18].r, E);
        get_Dx2my2a_Pyc(ho[13][18].r, E);
        get_Dx2my2a_Pzc(ho[14][18].r, E);
        get_Dx2my2a_Dxyc(ho[15][18].r, E);
        get_Dx2my2a_Dyzc(ho[16][18].r, E);
        get_Dx2my2a_Dzxc(ho[17][18].r, E);
        get_Dx2my2a_Dx2my2c(ho[18][18].r, E);
        get_Dx2my2a_Dz2c(ho[19][18].r, E);
        hd[18][18].r = E;
        /* Dz2a */
        E = hd[19][19].r;
        get_Dz2a_STc(ho[10][19].r, E);
        get_Dz2a_Sc(ho[11][19].r, E);
        get_Dz2a_Pxc(ho[12][19].r, E);
        get_Dz2a_Pyc(ho[13][19].r, E);
        get_Dz2a_Pzc(ho[14][19].r, E);
        get_Dz2a_Dxyc(ho[15][19].r, E);
        get_Dz2a_Dyzc(ho[16][19].r, E);
        get_Dz2a_Dzxc(ho[17][19].r, E);
        get_Dz2a_Dx2my2c(ho[18][19].r, E);
        get_Dz2a_Dz2c(ho[19][19].r, E);
        hd[19][19].r = E;
      } // if (BM_20_sp3d5ss_spin == BandModel)

    } else /* cation */
    {
      /* STc */
      E = hd[0][0].r;
      get_STc_STa(ho[0][0].r, E);
      get_STc_Sa(ho[1][0].r, E);
      get_STc_Pxa(ho[2][0].r, E);
      get_STc_Pya(ho[3][0].r, E);
      get_STc_Pza(ho[4][0].r, E);
      get_STc_Dxya(ho[5][0].r, E);
      get_STc_Dyza(ho[6][0].r, E);
      get_STc_Dzxa(ho[7][0].r, E);
      get_STc_Dx2my2a(ho[8][0].r, E);
      get_STc_Dz2a(ho[9][0].r, E);
      hd[0][0].r = E;
      /* Sc */
      E = hd[1][1].r;
      get_Sc_STa(ho[0][1].r, E);
      get_Sc_Sa(ho[1][1].r, E);
      get_Sc_Pxa(ho[2][1].r, E);
      get_Sc_Pya(ho[3][1].r, E);
      get_Sc_Pza(ho[4][1].r, E);
      get_Sc_Dxya(ho[5][1].r, E);
      get_Sc_Dyza(ho[6][1].r, E);
      get_Sc_Dzxa(ho[7][1].r, E);
      get_Sc_Dx2my2a(ho[8][1].r, E);
      get_Sc_Dz2a(ho[9][1].r, E);
      hd[1][1].r = E;
      /* Pxc */
      E = hd[2][2].r;
      get_Pxc_STa(ho[0][2].r, E);
      get_Pxc_Sa(ho[1][2].r, E);
      get_Pxc_Pxa(ho[2][2].r, E);
      get_Pxc_Pya(ho[3][2].r, E);
      get_Pxc_Pza(ho[4][2].r, E);
      get_Pxc_Dxya(ho[5][2].r, E);
      get_Pxc_Dyza(ho[6][2].r, E);
      get_Pxc_Dzxa(ho[7][2].r, E);
      get_Pxc_Dx2my2a(ho[8][2].r, E);
      get_Pxc_Dz2a(ho[9][2].r, E);
      hd[2][2].r = E;
      /* Pyc */
      E = hd[3][3].r;
      get_Pyc_STa(ho[0][3].r, E);
      get_Pyc_Sa(ho[1][3].r, E);
      get_Pyc_Pxa(ho[2][3].r, E);
      get_Pyc_Pya(ho[3][3].r, E);
      get_Pyc_Pza(ho[4][3].r, E);
      get_Pyc_Dxya(ho[5][3].r, E);
      get_Pyc_Dyza(ho[6][3].r, E);
      get_Pyc_Dzxa(ho[7][3].r, E);
      get_Pyc_Dx2my2a(ho[8][3].r, E);
      get_Pyc_Dz2a(ho[9][3].r, E);
      hd[3][3].r = E;
      /* Pzc */
      E = hd[4][4].r;
      get_Pzc_STa(ho[0][4].r, E);
      get_Pzc_Sa(ho[1][4].r, E);
      get_Pzc_Pxa(ho[2][4].r, E);
      get_Pzc_Pya(ho[3][4].r, E);
      get_Pzc_Pza(ho[4][4].r, E);
      get_Pzc_Dxya(ho[5][4].r, E);
      get_Pzc_Dyza(ho[6][4].r, E);
      get_Pzc_Dzxa(ho[7][4].r, E);
      get_Pzc_Dx2my2a(ho[8][4].r, E);
      get_Pzc_Dz2a(ho[9][4].r, E);
      hd[4][4].r = E;
      /* Dxyc */
      E = hd[5][5].r;
      get_Dxyc_STa(ho[0][5].r, E);
      get_Dxyc_Sa(ho[1][5].r, E);
      get_Dxyc_Pxa(ho[2][5].r, E);
      get_Dxyc_Pya(ho[3][5].r, E);
      get_Dxyc_Pza(ho[4][5].r, E);
      get_Dxyc_Dxya(ho[5][5].r, E);
      get_Dxyc_Dyza(ho[6][5].r, E);
      get_Dxyc_Dzxa(ho[7][5].r, E);
      get_Dxyc_Dx2my2a(ho[8][5].r, E);
      get_Dxyc_Dz2a(ho[9][5].r, E);
      hd[5][5].r = E;
      /* Dyzc */
      E = hd[6][6].r;
      get_Dyzc_STa(ho[0][6].r, E);
      get_Dyzc_Sa(ho[1][6].r, E);
      get_Dyzc_Pxa(ho[2][6].r, E);
      get_Dyzc_Pya(ho[3][6].r, E);
      get_Dyzc_Pza(ho[4][6].r, E);
      get_Dyzc_Dxya(ho[5][6].r, E);
      get_Dyzc_Dyza(ho[6][6].r, E);
      get_Dyzc_Dzxa(ho[7][6].r, E);
      get_Dyzc_Dx2my2a(ho[8][6].r, E);
      get_Dyzc_Dz2a(ho[9][6].r, E);
      hd[6][6].r = E;
      /* Dzxc */
      E = hd[7][7].r;
      get_Dzxc_STa(ho[0][7].r, E);
      get_Dzxc_Sa(ho[1][7].r, E);
      get_Dzxc_Pxa(ho[2][7].r, E);
      get_Dzxc_Pya(ho[3][7].r, E);
      get_Dzxc_Pza(ho[4][7].r, E);
      get_Dzxc_Dxya(ho[5][7].r, E);
      get_Dzxc_Dyza(ho[6][7].r, E);
      get_Dzxc_Dzxa(ho[7][7].r, E);
      get_Dzxc_Dx2my2a(ho[8][7].r, E);
      get_Dzxc_Dz2a(ho[9][7].r, E);
      hd[7][7].r = E;
      /* Dx2my2c */
      E = hd[8][8].r;
      get_Dx2my2c_STa(ho[0][8].r, E);
      get_Dx2my2c_Sa(ho[1][8].r, E);
      get_Dx2my2c_Pxa(ho[2][8].r, E);
      get_Dx2my2c_Pya(ho[3][8].r, E);
      get_Dx2my2c_Pza(ho[4][8].r, E);
      get_Dx2my2c_Dxya(ho[5][8].r, E);
      get_Dx2my2c_Dyza(ho[6][8].r, E);
      get_Dx2my2c_Dzxa(ho[7][8].r, E);
      get_Dx2my2c_Dx2my2a(ho[8][8].r, E);
      get_Dx2my2c_Dz2a(ho[9][8].r, E);
      hd[8][8].r = E;
      /* Dz2c */
      E = hd[9][9].r;
      get_Dz2c_STa(ho[0][9].r, E);
      get_Dz2c_Sa(ho[1][9].r, E);
      get_Dz2c_Pxa(ho[2][9].r, E);
      get_Dz2c_Pya(ho[3][9].r, E);
      get_Dz2c_Pza(ho[4][9].r, E);
      get_Dz2c_Dxya(ho[5][9].r, E);
      get_Dz2c_Dyza(ho[6][9].r, E);
      get_Dz2c_Dzxa(ho[7][9].r, E);
      get_Dz2c_Dx2my2a(ho[8][9].r, E);
      get_Dz2c_Dz2a(ho[9][9].r, E);
      hd[9][9].r = E;

      /* Repeate same calculations for spin down */
      if (BM_20_sp3d5ss_spin == BandModel) {
        /* STc */
        E = hd[10][10].r;
        get_STc_STa(ho[10][10].r, E);
        get_STc_Sa(ho[11][10].r, E);
        get_STc_Pxa(ho[12][10].r, E);
        get_STc_Pya(ho[13][10].r, E);
        get_STc_Pza(ho[14][10].r, E);
        get_STc_Dxya(ho[15][10].r, E);
        get_STc_Dyza(ho[16][10].r, E);
        get_STc_Dzxa(ho[17][10].r, E);
        get_STc_Dx2my2a(ho[18][10].r, E);
        get_STc_Dz2a(ho[19][10].r, E);
        hd[10][10].r = E;
        /* Sc */
        E = hd[11][11].r;
        get_Sc_STa(ho[10][11].r, E);
        get_Sc_Sa(ho[11][11].r, E);
        get_Sc_Pxa(ho[12][11].r, E);
        get_Sc_Pya(ho[13][11].r, E);
        get_Sc_Pza(ho[14][11].r, E);
        get_Sc_Dxya(ho[15][11].r, E);
        get_Sc_Dyza(ho[16][11].r, E);
        get_Sc_Dzxa(ho[17][11].r, E);
        get_Sc_Dx2my2a(ho[18][11].r, E);
        get_Sc_Dz2a(ho[19][11].r, E);
        hd[11][11].r = E;
        /* Pxc */
        E = hd[12][12].r;
        get_Pxc_STa(ho[10][12].r, E);
        get_Pxc_Sa(ho[11][12].r, E);
        get_Pxc_Pxa(ho[12][12].r, E);
        get_Pxc_Pya(ho[13][12].r, E);
        get_Pxc_Pza(ho[14][12].r, E);
        get_Pxc_Dxya(ho[15][12].r, E);
        get_Pxc_Dyza(ho[16][12].r, E);
        get_Pxc_Dzxa(ho[17][12].r, E);
        get_Pxc_Dx2my2a(ho[18][12].r, E);
        get_Pxc_Dz2a(ho[19][12].r, E);
        hd[12][12].r = E;
        /* Pyc */
        E = hd[13][13].r;
        get_Pyc_STa(ho[10][13].r, E);
        get_Pyc_Sa(ho[11][13].r, E);
        get_Pyc_Pxa(ho[12][13].r, E);
        get_Pyc_Pya(ho[13][13].r, E);
        get_Pyc_Pza(ho[14][13].r, E);
        get_Pyc_Dxya(ho[15][13].r, E);
        get_Pyc_Dyza(ho[16][13].r, E);
        get_Pyc_Dzxa(ho[17][13].r, E);
        get_Pyc_Dx2my2a(ho[18][13].r, E);
        get_Pyc_Dz2a(ho[19][13].r, E);
        hd[13][13].r = E;
        /* Pzc */
        E = hd[14][14].r;
        get_Pzc_STa(ho[10][14].r, E);
        get_Pzc_Sa(ho[11][14].r, E);
        get_Pzc_Pxa(ho[12][14].r, E);
        get_Pzc_Pya(ho[13][14].r, E);
        get_Pzc_Pza(ho[14][14].r, E);
        get_Pzc_Dxya(ho[15][14].r, E);
        get_Pzc_Dyza(ho[16][14].r, E);
        get_Pzc_Dzxa(ho[17][14].r, E);
        get_Pzc_Dx2my2a(ho[18][14].r, E);
        get_Pzc_Dz2a(ho[19][14].r, E);
        hd[14][14].r = E;
        /* Dxyc */
        E = hd[15][15].r;
        get_Dxyc_STa(ho[10][15].r, E);
        get_Dxyc_Sa(ho[11][15].r, E);
        get_Dxyc_Pxa(ho[12][15].r, E);
        get_Dxyc_Pya(ho[13][15].r, E);
        get_Dxyc_Pza(ho[14][15].r, E);
        get_Dxyc_Dxya(ho[15][15].r, E);
        get_Dxyc_Dyza(ho[16][15].r, E);
        get_Dxyc_Dzxa(ho[17][15].r, E);
        get_Dxyc_Dx2my2a(ho[18][15].r, E);
        get_Dxyc_Dz2a(ho[19][15].r, E);
        hd[15][15].r = E;
        /* Dyzc */
        E = hd[16][16].r;
        get_Dyzc_STa(ho[10][16].r, E);
        get_Dyzc_Sa(ho[11][16].r, E);
        get_Dyzc_Pxa(ho[12][16].r, E);
        get_Dyzc_Pya(ho[13][16].r, E);
        get_Dyzc_Pza(ho[14][16].r, E);
        get_Dyzc_Dxya(ho[15][16].r, E);
        get_Dyzc_Dyza(ho[16][16].r, E);
        get_Dyzc_Dzxa(ho[17][16].r, E);
        get_Dyzc_Dx2my2a(ho[18][16].r, E);
        get_Dyzc_Dz2a(ho[19][16].r, E);
        hd[16][16].r = E;
        /* Dzxc */
        E = hd[17][17].r;
        get_Dzxc_STa(ho[10][17].r, E);
        get_Dzxc_Sa(ho[11][17].r, E);
        get_Dzxc_Pxa(ho[12][17].r, E);
        get_Dzxc_Pya(ho[13][17].r, E);
        get_Dzxc_Pza(ho[14][17].r, E);
        get_Dzxc_Dxya(ho[15][17].r, E);
        get_Dzxc_Dyza(ho[16][17].r, E);
        get_Dzxc_Dzxa(ho[17][17].r, E);
        get_Dzxc_Dx2my2a(ho[18][17].r, E);
        get_Dzxc_Dz2a(ho[19][17].r, E);
        hd[17][17].r = E;
        /* Dx2my2c */
        E = hd[18][18].r;
        get_Dx2my2c_STa(ho[10][18].r, E);
        get_Dx2my2c_Sa(ho[11][18].r, E);
        get_Dx2my2c_Pxa(ho[12][18].r, E);
        get_Dx2my2c_Pya(ho[13][18].r, E);
        get_Dx2my2c_Pza(ho[14][18].r, E);
        get_Dx2my2c_Dxya(ho[15][18].r, E);
        get_Dx2my2c_Dyza(ho[16][18].r, E);
        get_Dx2my2c_Dzxa(ho[17][18].r, E);
        get_Dx2my2c_Dx2my2a(ho[18][18].r, E);
        get_Dx2my2c_Dz2a(ho[19][18].r, E);
        hd[18][18].r = E;
        /* Dz2c */
        E = hd[19][19].r;
        get_Dz2c_STa(ho[10][19].r, E);
        get_Dz2c_Sa(ho[11][19].r, E);
        get_Dz2c_Pxa(ho[12][19].r, E);
        get_Dz2c_Pya(ho[13][19].r, E);
        get_Dz2c_Pza(ho[14][19].r, E);
        get_Dz2c_Dxya(ho[15][19].r, E);
        get_Dz2c_Dyza(ho[16][19].r, E);
        get_Dz2c_Dzxa(ho[17][19].r, E);
        get_Dz2c_Dx2my2a(ho[18][19].r, E);
        get_Dz2c_Dz2a(ho[19][19].r, E);
        hd[19][19].r = E;
      } // if (BM_20_sp3d5ss_spin == BandModel)
    }
    /* Commented to retain spin spitting term in the Hamiltonian
          if (BM_20_sp3d5ss_spin == BandModel)
          {
             /* Now the calculations for spin up are finished.  The spin down
       calculations are identical, so copy the values we just calculated. */
    /*         for ( i=0; i < 10; i++ )
             {
                for ( j=0; j < 10; j++ )
                {
                   ho[i+10][j+10].r = ho[i][j].r;
                   hd[i+10][j+10].r = hd[i][j].r;
                }
             }
          }
       End Commented to retain spin spitting term in the Hamiltonian */
    break; /* end BM_10_sp3d5ss_nospin, BM_20_sp3d5ss_spin */

  case BM_10_sp3ss_spin:
    /* Explicitly list the calculations necessary for sp3s* spin up
       (0 <= row, col < 5) interactions. */
    if (basis_this[0] < STc) /* anion */
    {
      /* STa */
      E = hd[0][0].r;
      get_STa_STc(ho[0][0].r, E);
      get_STa_Sc(ho[1][0].r, E);
      get_STa_Pxc(ho[2][0].r, E);
      get_STa_Pyc(ho[3][0].r, E);
      get_STa_Pzc(ho[4][0].r, E);
      hd[0][0].r = E;
      /* Sa */
      E = hd[1][1].r;
      get_Sa_STc(ho[0][1].r, E);
      get_Sa_Sc(ho[1][1].r, E);
      get_Sa_Pxc(ho[2][1].r, E);
      get_Sa_Pyc(ho[3][1].r, E);
      get_Sa_Pzc(ho[4][1].r, E);
      hd[1][1].r = E;
      /* Pxa */
      E = hd[2][2].r;
      get_Pxa_STc(ho[0][2].r, E);
      get_Pxa_Sc(ho[1][2].r, E);
      get_Pxa_Pxc(ho[2][2].r, E);
      get_Pxa_Pyc(ho[3][2].r, E);
      get_Pxa_Pzc(ho[4][2].r, E);
      hd[2][2].r = E;
      /* Pya */
      E = hd[3][3].r;
      get_Pya_STc(ho[0][3].r, E);
      get_Pya_Sc(ho[1][3].r, E);
      get_Pya_Pxc(ho[2][3].r, E);
      get_Pya_Pyc(ho[3][3].r, E);
      get_Pya_Pzc(ho[4][3].r, E);
      hd[3][3].r = E;
      /* Pza */
      E = hd[4][4].r;
      get_Pza_STc(ho[0][4].r, E);
      get_Pza_Sc(ho[1][4].r, E);
      get_Pza_Pxc(ho[2][4].r, E);
      get_Pza_Pyc(ho[3][4].r, E);
      get_Pza_Pzc(ho[4][4].r, E);
      hd[4][4].r = E;
      /* Dxya */
      E = hd[5][5].r;
      get_Dxya_STc(ho[0][5].r, E);
      get_Dxya_Sc(ho[1][5].r, E);
      get_Dxya_Pxc(ho[2][5].r, E);
      get_Dxya_Pyc(ho[3][5].r, E);
      get_Dxya_Pzc(ho[4][5].r, E);
      hd[5][5].r = E;
      /* Dyza */
      E = hd[6][6].r;
      get_Dyza_STc(ho[0][6].r, E);
      get_Dyza_Sc(ho[1][6].r, E);
      get_Dyza_Pxc(ho[2][6].r, E);
      get_Dyza_Pyc(ho[3][6].r, E);
      get_Dyza_Pzc(ho[4][6].r, E);
      hd[6][6].r = E;
      /* Dzxa */
      E = hd[7][7].r;
      get_Dzxa_STc(ho[0][7].r, E);
      get_Dzxa_Sc(ho[1][7].r, E);
      get_Dzxa_Pxc(ho[2][7].r, E);
      get_Dzxa_Pyc(ho[3][7].r, E);
      get_Dzxa_Pzc(ho[4][7].r, E);
      hd[7][7].r = E;
      /* Dx2my2a */
      E = hd[8][8].r;
      get_Dx2my2a_STc(ho[0][8].r, E);
      get_Dx2my2a_Sc(ho[1][8].r, E);
      get_Dx2my2a_Pxc(ho[2][8].r, E);
      get_Dx2my2a_Pyc(ho[3][8].r, E);
      get_Dx2my2a_Pzc(ho[4][8].r, E);
      hd[8][8].r = E;
      /* Dz2a */
      E = hd[9][9].r;
      get_Dz2a_STc(ho[0][9].r, E);
      get_Dz2a_Sc(ho[1][9].r, E);
      get_Dz2a_Pxc(ho[2][9].r, E);
      get_Dz2a_Pyc(ho[3][9].r, E);
      get_Dz2a_Pzc(ho[4][9].r, E);
      hd[9][9].r = E;
    } else /* cation */
    {
      /* STc */
      E = hd[0][0].r;
      get_STc_STa(ho[0][0].r, E);
      get_STc_Sa(ho[1][0].r, E);
      get_STc_Pxa(ho[2][0].r, E);
      get_STc_Pya(ho[3][0].r, E);
      get_STc_Pza(ho[4][0].r, E);
      hd[0][0].r = E;
      /* Sc */
      E = hd[1][1].r;
      get_Sc_STa(ho[0][1].r, E);
      get_Sc_Sa(ho[1][1].r, E);
      get_Sc_Pxa(ho[2][1].r, E);
      get_Sc_Pya(ho[3][1].r, E);
      get_Sc_Pza(ho[4][1].r, E);
      hd[1][1].r = E;
      /* Pxc */
      E = hd[2][2].r;
      get_Pxc_STa(ho[0][2].r, E);
      get_Pxc_Sa(ho[1][2].r, E);
      get_Pxc_Pxa(ho[2][2].r, E);
      get_Pxc_Pya(ho[3][2].r, E);
      get_Pxc_Pza(ho[4][2].r, E);
      hd[2][2].r = E;
      /* Pyc */
      E = hd[3][3].r;
      get_Pyc_STa(ho[0][3].r, E);
      get_Pyc_Sa(ho[1][3].r, E);
      get_Pyc_Pxa(ho[2][3].r, E);
      get_Pyc_Pya(ho[3][3].r, E);
      get_Pyc_Pza(ho[4][3].r, E);
      hd[3][3].r = E;
      /* Pzc */
      E = hd[4][4].r;
      get_Pzc_STa(ho[0][4].r, E);
      get_Pzc_Sa(ho[1][4].r, E);
      get_Pzc_Pxa(ho[2][4].r, E);
      get_Pzc_Pya(ho[3][4].r, E);
      get_Pzc_Pza(ho[4][4].r, E);
      hd[4][4].r = E;
      /* Dxyc */
      E = hd[5][5].r;
      get_Dxyc_STa(ho[0][5].r, E);
      get_Dxyc_Sa(ho[1][5].r, E);
      get_Dxyc_Pxa(ho[2][5].r, E);
      get_Dxyc_Pya(ho[3][5].r, E);
      get_Dxyc_Pza(ho[4][5].r, E);
      hd[5][5].r = E;
      /* Dyzc */
      E = hd[6][6].r;
      get_Dyzc_STa(ho[0][6].r, E);
      get_Dyzc_Sa(ho[1][6].r, E);
      get_Dyzc_Pxa(ho[2][6].r, E);
      get_Dyzc_Pya(ho[3][6].r, E);
      get_Dyzc_Pza(ho[4][6].r, E);
      hd[6][6].r = E;
      /* Dzxc */
      E = hd[7][7].r;
      get_Dzxc_STa(ho[0][7].r, E);
      get_Dzxc_Sa(ho[1][7].r, E);
      get_Dzxc_Pxa(ho[2][7].r, E);
      get_Dzxc_Pya(ho[3][7].r, E);
      get_Dzxc_Pza(ho[4][7].r, E);
      hd[7][7].r = E;
      /* Dx2my2c */
      E = hd[8][8].r;
      get_Dx2my2c_STa(ho[0][8].r, E);
      get_Dx2my2c_Sa(ho[1][8].r, E);
      get_Dx2my2c_Pxa(ho[2][8].r, E);
      get_Dx2my2c_Pya(ho[3][8].r, E);
      get_Dx2my2c_Pza(ho[4][8].r, E);
      hd[8][8].r = E;
      /* Dz2c */
      E = hd[9][9].r;
      get_Dz2c_STa(ho[0][9].r, E);
      get_Dz2c_Sa(ho[1][9].r, E);
      get_Dz2c_Pxa(ho[2][9].r, E);
      get_Dz2c_Pya(ho[3][9].r, E);
      get_Dz2c_Pza(ho[4][9].r, E);
      hd[9][9].r = E;
    }
    /* Now the calculations for spin up are finished.  The spin down
       calculations are identical, so copy the values we just calculated. */
    for (i = 0; i < 5; i++) {
      for (j = 0; j < 5; j++) {
        ho[i + 5][j + 5].r = ho[i][j].r;
        hd[i + 5][j + 5].r = hd[i][j].r;
      }
    }
    break; /* end BM_10_sp3ss_spin */

  default: /* The most general case. */
#define basis1_econst basis_this[i]
#define basis2_econst basis_nbr[j]
#define Out_econst ho[j][i].r
#ifdef OLD_BROKEN_VERSION
#define Out_econst ho[i][j].r
#endif
#define Out_diag_corr hd[i][i].r
    for (i = 0; i < Nbasis; i++) {
      for (j = 0; j < Nbasis; j++) {
        /*  Off-diagonal coupling */
        if (spin[i] == spin[j]) {
          /* switch statement */
#include "mb_ham_spds_inclVoffdEdiag.h"
        }
      }
    }
#undef Out_diag_corr
#undef Out_econst
#undef basis2_econst
#undef basis1_econst
    break;
  } /* end switch (BandModel) */

#undef strain_vec
#undef vpar_orig
#undef vpar_e_const

  if (param_strained != param_keep) {
    param_strained = param_keep;
  }

#ifdef Additional_Diagonal_constant
  if (set_temp)
    mb_strain = NULL;
#endif

  return (0);
}

int h3d_diag(cmatrix h, ivectr basis0, ivectr spin, rvectr param, int Nbasis) {
  int i, j;
  complex cdum;

#define NEW_DIAG_replace
#ifdef NEW_DIAG_replace
  if (diag_map_unset)
    diag_map_set();
#endif

  for (i = 0; i < Nbasis; i++)
    for (j = 0; j < Nbasis; j++)
      h[i][j].r = h[i][j].i = 0.0;

  /* Construct Hamiltonian */
  for (i = 0; i < Nbasis; i++) {
    /* On-site Energies */
#ifdef NEW_DIAG_replace
    h[i][i].r = param[diag_map[basis0[i]]];
#else
    h[i][i].r = diag(param, basis0[i]);
#endif

    for (j = 0; j < Nbasis; j++) {
      /* Spin-Orbit interactions */
#define NEW_SPIN_ORBIT_FILL
#ifdef NEW_SPIN_ORBIT_FILL
#define Out_e_so cdum
#define e_so_basis1 basis0[i]
#define e_so_basis2 basis0[j]
#define e_so_spin1 spin[i]
#define e_so_spin2 spin[j]
#define e_so_delta_a param[pE_La]
#define e_so_delta_c param[pE_Lc]
#include "mb_ham_spds_so.h"
#undef Out_e_so
#undef e_so_basis1
#undef e_so_basis2
#undef e_so_spin1
#undef e_so_spin2
#undef e_so_delta_a
#undef e_so_delta_c
#else
      cdum = spin_orbit_spds(basis0[i], basis0[j], spin[i], spin[j],
                             param[pE_La], param[pE_Lc]);
#endif
      h[i][j].r += cdum.r;
      h[i][j].i += cdum.i;
    }
  }

  return (0);
}

int h3d_diag_surface_atom_shift(qd_struct d) {
  static int warnonce = 0;
  cmatrix h = d->hd;
  switch (d->BandModel) {
  case BM_10_sp3ss_spin:
    h[0][0].r += d->opt.ExecParam.Phys.BondShift_sstar;
    h[1][1].r += d->opt.ExecParam.Phys.BondShift_s;
    h[2][2].r += d->opt.ExecParam.Phys.BondShift_p;
    h[3][3].r += d->opt.ExecParam.Phys.BondShift_p;
    h[4][4].r += d->opt.ExecParam.Phys.BondShift_p;
    h[5][5].r += d->opt.ExecParam.Phys.BondShift_sstar;
    h[6][6].r += d->opt.ExecParam.Phys.BondShift_s;
    h[7][7].r += d->opt.ExecParam.Phys.BondShift_p;
    h[8][8].r += d->opt.ExecParam.Phys.BondShift_p;
    h[9][9].r += d->opt.ExecParam.Phys.BondShift_p;
    break;
  case BM_20_sp3d5ss_spin:
    h[0][0].r += d->opt.ExecParam.Phys.BondShift_sstar;
    h[1][1].r += d->opt.ExecParam.Phys.BondShift_s;
    h[2][2].r += d->opt.ExecParam.Phys.BondShift_p;
    h[3][3].r += d->opt.ExecParam.Phys.BondShift_p;
    h[4][4].r += d->opt.ExecParam.Phys.BondShift_p;
    h[5][5].r += d->opt.ExecParam.Phys.BondShift_d;
    h[6][6].r += d->opt.ExecParam.Phys.BondShift_d;
    h[7][7].r += d->opt.ExecParam.Phys.BondShift_d;
    h[8][8].r += d->opt.ExecParam.Phys.BondShift_d;
    h[9][9].r += d->opt.ExecParam.Phys.BondShift_d;
    h[10][10].r += d->opt.ExecParam.Phys.BondShift_sstar;
    h[11][11].r += d->opt.ExecParam.Phys.BondShift_s;
    h[12][12].r += d->opt.ExecParam.Phys.BondShift_p;
    h[13][13].r += d->opt.ExecParam.Phys.BondShift_p;
    h[14][14].r += d->opt.ExecParam.Phys.BondShift_p;
    h[15][15].r += d->opt.ExecParam.Phys.BondShift_d;
    h[16][16].r += d->opt.ExecParam.Phys.BondShift_d;
    h[17][17].r += d->opt.ExecParam.Phys.BondShift_d;
    h[18][18].r += d->opt.ExecParam.Phys.BondShift_d;
    h[19][19].r += d->opt.ExecParam.Phys.BondShift_d;
    break;
  case BM_1_s_nospin:
    if (warnonce == 0) {
      printf("Cannot shift orbital energies in Bands_1_s_nospin: "
             "h3d_diag_surface_atom_shift\n");
      warnonce = 1;
    }
    break;
  case BM_10_sp3d5ss_nospin:
    h[0][0].r += d->opt.ExecParam.Phys.BondShift_sstar;
    h[1][1].r += d->opt.ExecParam.Phys.BondShift_s;
    h[2][2].r += d->opt.ExecParam.Phys.BondShift_p;
    h[3][3].r += d->opt.ExecParam.Phys.BondShift_p;
    h[4][4].r += d->opt.ExecParam.Phys.BondShift_p;
    h[5][5].r += d->opt.ExecParam.Phys.BondShift_d;
    h[6][6].r += d->opt.ExecParam.Phys.BondShift_d;
    h[7][7].r += d->opt.ExecParam.Phys.BondShift_d;
    h[8][8].r += d->opt.ExecParam.Phys.BondShift_d;
    h[9][9].r += d->opt.ExecParam.Phys.BondShift_d;
    break;
  default:
    die("Unidentified bandstructure model: h3d_diag_surface_atom_shift");
    break;
  }

  return (0);
}

void h3d_diag_dangling_bond_shift(qd_struct d, int cindx1, int aindx1) {
  static int warnonce = 0;
  if (d->BandModel == BM_1_s_nospin) {
    if (warnonce == 0) {
      printf("Cannot shift dangling bonds in Bands_1_s_nospin: "
             "h3d_diag_dangling_bond_shift\n");
      warnonce = 1;
    }
    return;
  }

  if (d->geo.isCellWurtzite())
    h3d_diag_dangling_bond_shift_wurtzite(d, cindx1, aindx1);
  /* New ZB_110 unit cell */
  else if (d->geo.isCellZB_110() || d->geo.isCellZB_110_small() ||
           d->geo.isCellZB_111() ||
           d->geo.isCellUserDef()) /* Added for ZB_110 */
    h3d_diag_dangling_bond_shift_ZB_rot(d, cindx1, aindx1);
  /* End New ZB_110 unit cell */
  else
    h3d_diag_dangling_bond_shift_zincblende(d, cindx1, aindx1);
}

void h3d_diag_dangling_bond_shift_wurtzite(qd_struct d, int cindx1,
                                           int aindx1) {
  double shift = d->opt.ExecParam.Phys.DanglingBondEnergyShift;
  cmatrix h = d->hd;
  rmatrix H_sp3 = Rmatrix(4, 4);
  rmatrix T = Rmatrix(4, 4);

  H_sp3[0][0] = d->hd[1][1].r;
  H_sp3[1][1] = d->hd[2][2].r;
  H_sp3[2][2] = d->hd[3][3].r;
  H_sp3[3][3] = d->hd[4][4].r;

  rvectr v_dang = Rvectr(4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      T[i][j] = d->geo.SP3_hybrid(aindx1, i, j);
    }
  }

  int i1 = d->geo.cell__ijk[cindx1][0];
  int j1 = d->geo.cell__ijk[cindx1][1];
  int k1 = d->geo.cell__ijk[cindx1][2];

  for (int n = 0; n < d->geo.Neighbors(aindx1); n++) {
    int i2 = d->geo.NbrCell(aindx1, n, 0);
    int j2 = d->geo.NbrCell(aindx1, n, 1);
    int k2 = d->geo.NbrCell(aindx1, n, 2);
    int cindx2 = d->geo.ijk__cell[i1 + i2][j1 + j2][k1 + k2];
    int aindx2 = d->geo.NbrCell(aindx1, n, 3);
    if (cindx2 >= 0 && d->geo.AtomType[cindx2][aindx2])
      v_dang[n] = 0.0;

#ifdef ELIMINATE_SSMAP
    else if (d->geo.hasPeriodicity() && (cindx2 + d->geo.N_Cell + 1) >= 0) {
      int cindx3 = cindx2 + d->geo.N_Cell + 1;
#else  /* ELIMINATE_SSMAP */
    else if (d->geo.hasPeriodicity() && d->geo.ssmap[cindx1][aindx1][n] != -1) {
      int cindx3 = d->geo.ssmap[cindx1][aindx1][n];
#endif /* ELIMINATE_SSMAP */

      if (d->geo.AtomType[cindx3][aindx2])
        v_dang[n] = 0.0;
      else
        v_dang[n] = shift;
    } else
      v_dang[n] = shift;
  }

  // matrix-matrix multiplication
  rmatrix TH = Rmatrix(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      TH[i][j] = 0.0;
      for (int k = 0; k < 4; k++) {
        TH[i][j] += T[i][k] * H_sp3[k][j];
      }
    }
  }
  rmatrix THT_ = Rmatrix(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      THT_[i][j] = 0.0;
      for (int k = 0; k < 4; k++) {
        THT_[i][j] += TH[i][k] * T[j][k];
        // transform Hamiltonian to sp3 hybrid basis
      }
    }
    THT_[i][i] += v_dang[i];
    // shifting the energies of dangling bonds.
  }
  rmatrix T_THT_ = Rmatrix(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      T_THT_[i][j] = 0.0;
      for (int k = 0; k < 4; k++) {
        T_THT_[i][j] += T[k][i] * THT_[k][j];
      }
    }
  }
  rmatrix T_THT_T = Rmatrix(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      T_THT_T[i][j] = 0.0;
      for (int k = 0; k < 4; k++) {
        T_THT_T[i][j] += T_THT_[i][k] * T[k][j];
        // transform Hamiltonian back to s px,py,pz basis
      }
    }
  }

  // copying the new hamiltonian to d->hd
  switch (d->BandModel) {
  case BM_10_sp3ss_spin:
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        h[i + 1][j + 1].r = T_THT_T[i][j];
        h[i + 6][j + 6].r = T_THT_T[i][j];
      }
    }
    break;
  case BM_20_sp3d5ss_spin:
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        h[i + 1][j + 1].r = T_THT_T[i][j];
        h[i + 11][j + 11].r = T_THT_T[i][j];
      }
    }
    break;
  case BM_10_sp3d5ss_nospin:
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        h[i + 1][j + 1].r = T_THT_T[i][j];
      }
    }
    break;
  default:
    die("Illegal bandstructure model in "
        "h3d_diag_dangling_bond_shift_wurtzite\n");
  }

  rm_rmatrix(&T);
  rm_rmatrix(&H_sp3);
  rm_rmatrix(&TH);
  rm_rmatrix(&THT_);
  rm_rmatrix(&T_THT_);
  rm_rmatrix(&T_THT_T);
  rm_rvectr(&v_dang);
}

/* Added for ZB_110 and ZB_111 */
void h3d_diag_dangling_bond_shift_ZB_rot(qd_struct d, int cindx1, int aindx1) {
  double shift = d->opt.ExecParam.Phys.DanglingBondEnergyShift;
  cmatrix h = d->hd;
  rmatrix H_sp3 = Rmatrix(4, 4);
  rmatrix T = Rmatrix(4, 4);

  H_sp3[0][0] = d->hd[1][1].r;
  H_sp3[1][1] = d->hd[2][2].r;
  H_sp3[2][2] = d->hd[3][3].r;
  H_sp3[3][3] = d->hd[4][4].r;

  rmatrix T_axes = Rmatrix(4, 4);
  rvectr a_latt = Rvectr(3); /* Lattice constants of ZB_111 */
  /* Transformation matrix from ZB_* axes to x0,y0,z0 axes */
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      T_axes[i][j] = d->geo.get_T_rot(i, j);
      a_latt[j] = d->geo.get_a_latt_rot_norm(j);
    }

  int i1 = d->geo.cell__ijk[cindx1][0];
  int j1 = d->geo.cell__ijk[cindx1][1];
  int k1 = d->geo.cell__ijk[cindx1][2];

  rvectr v_dang = Rvectr(4);
  v_dang[0] = shift;
  v_dang[1] = shift;
  v_dang[2] = shift;
  v_dang[3] = shift;
  rvectr nnv = Rvectr(3);
  rvectr nnv_r = Rvectr(3);

  for (int n = 0; n < d->geo.Neighbors(aindx1); n++) {
    int i2 = d->geo.NbrCell(aindx1, n, 0);
    int j2 = d->geo.NbrCell(aindx1, n, 1);
    int k2 = d->geo.NbrCell(aindx1, n, 2);
    int cindx2 = d->geo.ijk__cell[i1 + i2][j1 + j2][k1 + k2];
    int aindx2 = d->geo.NbrCell(aindx1, n, 3);

    // approximate bond direction. not necessary to add displacement.
    nnv[0] = i2 + (d->geo.PositionInCell(aindx2, 0) -
                   d->geo.PositionInCell(aindx1, 0));
    nnv[1] = j2 + (d->geo.PositionInCell(aindx2, 1) -
                   d->geo.PositionInCell(aindx1, 1));
    nnv[2] = k2 + (d->geo.PositionInCell(aindx2, 2) -
                   d->geo.PositionInCell(aindx1, 2));

    /* Scale to new lattice constants */
    nnv[0] *= a_latt[0];
    nnv[1] *= a_latt[1];
    nnv[2] *= a_latt[2];

    /* nnv_r = transpose(T_axes)*(nnv.*a_latt) */
    for (int i = 0; i < 3; i++) {
      nnv_r[i] = 0.0;
      for (int j = 0; j < 3; j++) {
        nnv_r[i] += T_axes[j][i] * nnv[j];
      }
    }
    /* Assign vallues to T */
    T[n][0] = 0.5;
    for (int i = 0; i < 3; i++) {
      T[n][i + 1] = 0.5 * 4.0 * nnv_r[i];
    }

    if (cindx2 >= 0 && d->geo.AtomType[cindx2][aindx2])
      v_dang[n] = 0.0;

#ifdef ELIMINATE_SSMAP
    else if (d->geo
                 .hasPeriodicity() /*&& (cindx2 + d->geo.N_Cell + 1 ) >= 0*/) {
      int cindx3 = -9999;
      if ((cindx2 + d->geo.N_Cell + 1) >= 0)
        cindx3 = cindx2 + d->geo.N_Cell + 1;
      if (d->geo.is_tilted)
        cindx3 = d->geo.get_cindx_tilted(cindx1, aindx1, n);
#else  /* ELIMINATE_SSMAP */
    else if (d->geo.hasPeriodicity() && d->geo.ssmap[cindx1][aindx1][n] != -1) {
      int cindx3 = d->geo.ssmap[cindx1][aindx1][n];
#endif /* ELIMINATE_SSMAP */

      if (cindx3 >= 0)
        if (d->geo.AtomType[cindx3][aindx2])
          v_dang[n] = 0.0;
    }
  }

  // Transpose
  rmatrix T_ = Rmatrix(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      T_[i][j] = T[j][i];
    }
  }

  rmatrix TH = Rmatrix(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      TH[i][j] = 0.0;
      for (int k = 0; k < 4; k++) {
        TH[i][j] += T[i][k] * H_sp3[k][j];
      }
    }
  }

  rmatrix THT_ = Rmatrix(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      THT_[i][j] = 0.0;
      for (int k = 0; k < 4; k++) {
        THT_[i][j] += TH[i][k] * T_[k][j]; /*T[j][k]*/
        ;
        // transform Hamiltonian to sp3 hybrid basis
      }
    }
    THT_[i][i] += v_dang[i];
    // shifting the energies of dangling bonds.
  }

  rmatrix T_THT_ = Rmatrix(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      T_THT_[i][j] = 0.0;
      for (int k = 0; k < 4; k++) {
        T_THT_[i][j] += /*T[k][i]*/ T_[i][k] * THT_[k][j];
      }
    }
  }

  rmatrix T_THT_T = Rmatrix(4, 4);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      T_THT_T[i][j] = 0.0;
      for (int k = 0; k < 4; k++) {
        T_THT_T[i][j] += T_THT_[i][k] * T[k][j];
        // transform Hamiltonian back to s px,py,pz basis
      }
    }
  }

  // copying the new hamiltonian to d->hd
  switch (d->BandModel) {
  case BM_10_sp3ss_spin:
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        h[i + 1][j + 1].r = T_THT_T[i][j];
        h[i + 6][j + 6].r = T_THT_T[i][j];
      }
    }
    break;
  case BM_20_sp3d5ss_spin:
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        h[i + 1][j + 1].r = T_THT_T[i][j];
        h[i + 11][j + 11].r = T_THT_T[i][j];
      }
    }
    break;
  case BM_10_sp3d5ss_nospin:
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        h[i + 1][j + 1].r = T_THT_T[i][j];
      }
    }
    break;
  default:
    die("Illegal bandstructure model in "
        "h3d_diag_dangling_bond_shift_zincblende_rot\n");
  }

  rm_rmatrix(&T);
  rm_rmatrix(&H_sp3);
  rm_rmatrix(&T_);
  rm_rmatrix(&TH);
  rm_rmatrix(&THT_);
  rm_rmatrix(&T_THT_);
  rm_rmatrix(&T_THT_T);
  rm_rvectr(&v_dang);
  rm_rmatrix(&T_axes);
  rm_rvectr(&a_latt);
  rm_rvectr(&nnv);
  rm_rvectr(&nnv_r);
}
/* End Added for ZB_110 and ZB_111 */

void h3d_diag_dangling_bond_shift_zincblende(qd_struct d, int cindx1,
                                             int aindx1) {
  double shift = d->opt.ExecParam.Phys.DanglingBondEnergyShift;
  shift *= 0.25;

  cmatrix h = d->hd;
  int i1 = d->geo.cell__ijk[cindx1][0];
  int j1 = d->geo.cell__ijk[cindx1][1];
  int k1 = d->geo.cell__ijk[cindx1][2];

  rvectr nnv = Rvectr(3);

  rvectr *bond;
  bond = new rvectr[4];
  bond[0] = Rvectr(3);
  bond[1] = Rvectr(3);
  bond[2] = Rvectr(3);
  bond[3] = Rvectr(3);

  bond[0][0] = 0.25;
  bond[0][1] = 0.25;
  bond[0][2] = 0.25;
  bond[1][0] = -0.25;
  bond[1][1] = -0.25;
  bond[1][2] = 0.25;
  bond[2][0] = 0.25;
  bond[2][1] = -0.25;
  bond[2][2] = -0.25;
  bond[3][0] = -0.25;
  bond[3][1] = 0.25;
  bond[3][2] = -0.25;

  rvectr v_dang = Rvectr(4);
  v_dang[0] = shift;
  v_dang[1] = shift;
  v_dang[2] = shift;
  v_dang[3] = shift;

  double cosv = 0.0;
  for (int n = 0; n < d->geo.Neighbors(aindx1); n++) {
    int i2 = d->geo.NbrCell(aindx1, n, 0);
    int j2 = d->geo.NbrCell(aindx1, n, 1);
    int k2 = d->geo.NbrCell(aindx1, n, 2);
    int cindx2 = d->geo.ijk__cell[i1 + i2][j1 + j2][k1 + k2];
    int aindx2 = d->geo.NbrCell(aindx1, n, 3);
    // approximate bond direction. not necessary to add displacement.
    nnv[0] = i2 + (d->geo.PositionInCell(aindx2, 0) -
                   d->geo.PositionInCell(aindx1, 0));
    nnv[1] = j2 + (d->geo.PositionInCell(aindx2, 1) -
                   d->geo.PositionInCell(aindx1, 1));
    nnv[2] = k2 + (d->geo.PositionInCell(aindx2, 2) -
                   d->geo.PositionInCell(aindx1, 2));

    int ib = 0;
    cosv = vector_dot_vector(nnv, bond[ib], 3);
    while (fabs(cosv) < 0.15 && ib < 3)
      cosv = vector_dot_vector(nnv, bond[++ib], 3);

    if (cindx2 >= 0 && d->geo.AtomType[cindx2][aindx2]) {
      v_dang[ib] = 0.0;
    }

#ifdef ELIMINATE_SSMAP
    else if (d->geo
                 .hasPeriodicity() /*&& (cindx2 + d->geo.N_Cell + 1 ) >= 0*/) {
      int cindx3 = -9999;
      if ((cindx2 + d->geo.N_Cell + 1) >= 0)
        cindx3 = cindx2 + d->geo.N_Cell + 1;
      if (d->geo.is_tilted)
        cindx3 = d->geo.get_cindx_tilted(cindx1, aindx1, n);
#else  /* ELIMINATE_SSMAP */
    else if (d->geo.hasPeriodicity() && d->geo.ssmap[cindx1][aindx1][n] != -1) {
      int cindx3 = d->geo.ssmap[cindx1][aindx1][n];
#endif /* ELIMINATE_SSMAP */

      if (cindx3 >= 0)
        if (d->geo.AtomType[cindx3][aindx2])
          v_dang[ib] = 0.0;
    }
  }

  if (cosv > 0) {

    double vsum = v_dang[0] + v_dang[1] + v_dang[2] + v_dang[3];
    h[1][1].r += vsum;
    h[2][2].r += vsum;
    h[3][3].r += vsum;
    h[4][4].r += vsum;
    /* Bands_10_sp3d5ss_nospin:
       This bandmodel is spin-less.
       The modification of the Hamiltonian therefore only needs to be computed
       for a single block!
    */
    if (d->BandModel == BM_10_sp3ss_spin /* opposite spin component */) {
      h[6][6].r += vsum;
      h[7][7].r += vsum;
      h[8][8].r += vsum;
      h[9][9].r += vsum;
    } else if (d->BandModel ==
               BM_20_sp3d5ss_spin /* opposite spin component */) {
      h[11][11].r += vsum;
      h[12][12].r += vsum;
      h[13][13].r += vsum;
      h[14][14].r += vsum;
    }

    vsum = v_dang[0] - v_dang[1] + v_dang[2] - v_dang[3];
    h[1][2].r += vsum;
    h[2][1].r += vsum;
    h[3][4].r += vsum;
    h[4][3].r += vsum;
    if (d->BandModel == BM_10_sp3ss_spin /* opposite spin component */) {
      h[6][7].r += vsum;
      h[7][6].r += vsum;
      h[8][9].r += vsum;
      h[9][8].r += vsum;
    } else if (d->BandModel ==
               BM_20_sp3d5ss_spin /* opposite spin component */) {
      h[11][12].r += vsum;
      h[12][11].r += vsum;
      h[13][14].r += vsum;
      h[14][13].r += vsum;
    }

    vsum = v_dang[0] - v_dang[1] - v_dang[2] + v_dang[3];
    h[1][3].r += vsum;
    h[3][1].r += vsum;
    h[2][4].r += vsum;
    h[4][2].r += vsum;
    if (d->BandModel == BM_10_sp3ss_spin /* opposite spin component */) {
      h[6][8].r += vsum;
      h[8][6].r += vsum;
      h[7][9].r += vsum;
      h[9][7].r += vsum;
    } else if (d->BandModel ==
               BM_20_sp3d5ss_spin /* opposite spin component */) {
      h[11][13].r += vsum;
      h[13][11].r += vsum;
      h[12][14].r += vsum;
      h[14][12].r += vsum;
    }

    vsum = v_dang[0] + v_dang[1] - v_dang[2] - v_dang[3];
    h[1][4].r += vsum;
    h[4][1].r += vsum;
    h[2][3].r += vsum;
    h[3][2].r += vsum;
    if (d->BandModel == BM_10_sp3ss_spin /* opposite spin component */) {
      h[6][9].r += vsum;
      h[9][6].r += vsum;
      h[7][8].r += vsum;
      h[8][7].r += vsum;
    } else if (d->BandModel ==
               BM_20_sp3d5ss_spin /* opposite spin component */) {
      h[11][14].r += vsum;
      h[14][11].r += vsum;
      h[12][13].r += vsum;
      h[13][12].r += vsum;
    }
  }

  else {
    double vsum = v_dang[0] + v_dang[1] + v_dang[2] + v_dang[3];
    h[1][1].r += vsum;
    h[2][2].r += vsum;
    h[3][3].r += vsum;
    h[4][4].r += vsum;
    if (d->BandModel == BM_10_sp3ss_spin /* opposite spin component */) {
      h[6][6].r += vsum;
      h[7][7].r += vsum;
      h[8][8].r += vsum;
      h[9][9].r += vsum;
    } else if (d->BandModel ==
               BM_20_sp3d5ss_spin /* opposite spin component */) {
      h[11][11].r += vsum;
      h[12][12].r += vsum;
      h[13][13].r += vsum;
      h[14][14].r += vsum;
    }

    vsum = v_dang[0] - v_dang[1] + v_dang[2] - v_dang[3];
    h[1][2].r -= vsum;
    h[2][1].r -= vsum;
    h[3][4].r += vsum;
    h[4][3].r += vsum;
    if (d->BandModel == BM_10_sp3ss_spin /* opposite spin component */) {
      h[6][7].r -= vsum;
      h[7][6].r -= vsum;
      h[8][9].r += vsum;
      h[9][8].r += vsum;
    } else if (d->BandModel ==
               BM_20_sp3d5ss_spin /* opposite spin component */) {
      h[11][12].r -= vsum;
      h[12][11].r -= vsum;
      h[13][14].r += vsum;
      h[14][13].r += vsum;
    }

    vsum = v_dang[0] - v_dang[1] - v_dang[2] + v_dang[3];
    h[1][3].r -= vsum;
    h[3][1].r -= vsum;
    h[2][4].r += vsum;
    h[4][2].r += vsum;
    if (d->BandModel == BM_10_sp3ss_spin /* opposite spin component */) {
      h[6][8].r -= vsum;
      h[8][6].r -= vsum;
      h[7][9].r += vsum;
      h[9][7].r += vsum;
    } else if (d->BandModel ==
               BM_20_sp3d5ss_spin /* opposite spin component */) {
      h[11][13].r -= vsum;
      h[13][11].r -= vsum;
      h[12][14].r += vsum;
      h[14][12].r += vsum;
    }

    vsum = v_dang[0] + v_dang[1] - v_dang[2] - v_dang[3];
    h[1][4].r -= vsum;
    h[4][1].r -= vsum;
    h[2][3].r += vsum;
    h[3][2].r += vsum;
    if (d->BandModel == BM_10_sp3ss_spin /* opposite spin component */) {
      h[6][9].r -= vsum;
      h[9][6].r -= vsum;
      h[7][8].r += vsum;
      h[8][7].r += vsum;
    } else if (d->BandModel ==
               BM_20_sp3d5ss_spin /* opposite spin component */) {
      h[11][14].r -= vsum;
      h[14][11].r -= vsum;
      h[12][13].r += vsum;
      h[13][12].r += vsum;
    }
  }

  rm_rvectr(&nnv);
  rm_rvectr(&bond[0]);
  rm_rvectr(&bond[1]);
  rm_rvectr(&bond[2]);
  rm_rvectr(&bond[3]);
  rm_rvectr(&v_dang);
  delete[] bond;
}

int h3d_diag_pz_addition_for_wurtzite_crystal_splitting(qd_struct d,
                                                        double pz_add) {
  if (d->BandModel == BM_10_sp3ss_spin) {
    cmatrix h = d->hd;
    h[4][4].r += pz_add;
    h[9][9].r += pz_add;
    return (0);
  } else if (d->BandModel == BM_20_sp3d5ss_spin) {
    cmatrix h = d->hd;
    h[4][4].r += pz_add;
    h[14][14].r += pz_add;
    return (0);
  } else if (d->BandModel == BM_10_sp3d5ss_nospin) {
    cmatrix h = d->hd;
    h[4][4].r += pz_add;
    return (0);
  } else
    return (0);
}

double vector_dot_vector(rvectr vec1, rvectr vec2, int dim) {
  double dot_product = 0.0;
  for (int d = 0; d < dim; d++)
    dot_product += vec1[d] * vec2[d];
  return dot_product;
}

complex ei(rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length,
           rvectr k, rvectr d, rvectr d0, int ba, int bc, real *E_diag_corr) {
  int i, kdim;
  real arg, E;
  complex r;

  static rvectr param_strained = NULL;
  rvectr param_keep = NULL;
  static real l, m, n, l2, m2, n2;
  static real l0, m0, n0, l02, m02, n02;
  real cstrain = 0.0;
#ifdef Additional_Diagonal_constant
  static rvectr mb_strain_temp = NULL;
  int set_temp = 0;
#endif

  if (param_strained == NULL) {
    param_strained = Rvectr(nml_dv_extent(param));
  } else {
    param_strained = check_rvectr(param_strained, nml_dv_extent(param));
  }

  param_keep = param_strained;

  /* E_const_prep(param, mb_strain, param_strained, unstrnd_cubic_cell_length,
     d[0], d[1], d[2],&l, &m, &n, &l2, &m2, &n2); */
  strain_corr_prep(param, mb_strain, param_strained, unstrnd_cubic_cell_length,
                   d[0], d[1], d[2], d0[0], d0[1], d0[2], &l, &m, &n, &l2, &m2,
                   &n2, &l0, &m0, &n0, &l02, &m02, &n02);
  if (mb_strain == NULL)
    param_strained = param;

  kdim = nml_dv_extent(k);

  if (mb_strain)
    cstrain = mb_strain[eta_const];

#ifdef Additional_Diagonal_constant
  if (!mb_strain) {
    if (mb_strain_temp == NULL) {
      mb_strain_temp = Rvectr(sp3d5s_nstrain);
    } else {
      mb_strain_temp = check_rvectr(param_strained, nml_dv_extent(param));
    }
    mb_strain = mb_strain_temp;
    set_temp = 1;
  }
  E = E_const_get(param_strained, param, mb_strain, ba, bc, l, m, n, l2, m2, n2,
                  l0, m0, n0, l02, m02, n02, cstrain, E_diag_corr);
  if (set_temp)
    mb_strain = NULL;

#else
  E = E_const_get(param_strained, param, ba, bc, l, m, n, l2, m2, n2, l0, m0,
                  n0, l02, m02, n02, cstrain, E_diag_corr);
#endif

  arg = 0.0;

  if (ba > 9) {
    for (i = 0; i < kdim; i++)
      arg -= k[i] * d[i];
  } else {
    for (i = 0; i < kdim; i++)
      arg += k[i] * d[i];
  }

  r.r = E * cos(arg);
  r.i = E * sin(arg);

  if (param_strained != param_keep)
    param_strained = param_keep;

  return r;
}

complex eidk(rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length,
             rvectr k, rvectr d, rvectr d0, int ba, int bc, int dim1,
             real *E_diag_corr) {
  int i, kdim;
  real arg, E;
  complex r, iarg, eiarg;

  static rvectr param_strained = NULL;
  rvectr param_keep = NULL;
  static real l, m, n, l2, m2, n2;
  static real l0, m0, n0, l02, m02, n02;
  real cstrain = 0.0;
#ifdef Additional_Diagonal_constant
  static rvectr mb_strain_temp = NULL;
  int set_temp = 0;
#endif

  if (param_strained == NULL) {
    param_strained = Rvectr(nml_dv_extent(param));
  } else {
    param_strained = check_rvectr(param_strained, nml_dv_extent(param));
  }

  param_keep = param_strained;

  /* E_const_prep(param, mb_strain, param_strained, unstrnd_cubic_cell_length,
     d[0], d[1], d[2], &l, &m, &n, &l2, &m2, &n2); */
  strain_corr_prep(param, mb_strain, param_strained, unstrnd_cubic_cell_length,
                   d[0], d[1], d[2], d0[0], d0[1], d0[2], &l, &m, &n, &l2, &m2,
                   &n2, &l0, &m0, &n0, &l02, &m02, &n02);
  if (mb_strain == NULL)
    param_strained = param;

  if (mb_strain)
    cstrain = mb_strain[eta_const];

  kdim = nml_dv_extent(k);

#ifdef Additional_Diagonal_constant
  if (!mb_strain) {
    if (mb_strain_temp == NULL) {
      mb_strain_temp = Rvectr(sp3d5s_nstrain);
    } else {
      mb_strain_temp = check_rvectr(param_strained, nml_dv_extent(param));
    }
    mb_strain = mb_strain_temp;
    set_temp = 1;
  }
  E = E_const_get(param_strained, param, mb_strain, ba, bc, l, m, n, l2, m2, n2,
                  l0, m0, n0, l02, m02, n02, cstrain, E_diag_corr);
  if (set_temp)
    mb_strain = NULL;
#else
  E = E_const_get(param_strained, param, ba, bc, l, m, n, l2, m2, n2, l0, m0,
                  n0, l02, m02, n02, cstrain, E_diag_corr);
#endif

  /* arg = 0.0;
  for ( i=0; i < kdim; i++ )
  arg += k[i]*d[i]; */

  arg = iarg.r = 0.0;
  if (ba > 9) {
    iarg.i = -d[dim1];
    for (i = 0; i < kdim; i++)
      arg -= k[i] * d[i];
  } else {
    iarg.i = d[dim1];
    for (i = 0; i < kdim; i++)
      arg += k[i] * d[i];
  }

  eiarg = Complex(E * cos(arg), E * sin(arg));
  /* iarg.r = 0.0;
     iarg.i = d[dim1]; */

  r = cmul(iarg, eiarg);

  if (param_strained != param_keep)
    param_strained = param_keep;

  return r;
}

int tbh_const_bulk(cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis,
                   ivectr spin, rvectr param, rvectr mb_strain,
                   real unstrnd_cubic_cell_length, rvectr k) {
  int i, j, l, nn1, nbasis;
  complex cdum;

  nbasis = nml_iv_extent(basis);
  nn1 = rdim(nnv);

  /* Set matrix to zero */
  for (i = 0; i < nbasis; i++)
    for (j = 0; j < nbasis; j++)
      h[i][j].r = h[i][j].i = 0.0;

  /* Construct Hamiltonian */
  for (i = 0; i < nbasis; i++) {
    /* On-site Energies */
    h[i][i].r = diag(param, basis[i]);
    for (j = 0; j < nbasis; j++) {
      /* Nearest-neigbor interactions */
      /* if ( basis[i] == 1 && basis[j] == 1 ) j++; */
      if (spin[i] == spin[j]) {
        /* Sum over neigbors */
        for (l = 0; l < nn1; l++) {
          cdum = ei(param, mb_strain, unstrnd_cubic_cell_length, k, nnv[l],
                    nnv0[l], basis[i], basis[j], &(h[i][i].r));
          h[i][j].r += cdum.r;
          h[i][j].i += cdum.i;
        }
      }
      /* Spin-Orbit interactions */
      cdum = spin_orbit_spds(basis[i], basis[j], spin[i], spin[j], param[8],
                             param[9]);
      h[i][j].r += cdum.r;
      h[i][j].i += cdum.i;
    }
  }

  return (0);
}

int tbh_const_1d_nemo(cmatrix h, cmatrix hp1, rmatrix nnv, rmatrix nnv0,
                      ivectr basis, ivectr spin, rvectr param, rvectr mb_strain,
                      real unstrnd_cubic_cell_length, rvectr k) {
  int i, j, l, nn1, nbasis, nb;
  complex cdum;

  real StrainShift = 0.0;
#ifdef STRAIN_SHIFT_ENABLE_1D
  if (mb_strain) {
    if (sp3d5s_Ebound == sp3d5s_Ebound_undefine) {
      sp3d5s_Ebound = mb_strain[eta_cat_shift];
      printf("FIX:  currently use cation shift\n");
    }
    if (sp3d5s_Ebound != sp3d5s_Ebound_default) {
      printf("Only support a default value of %f as an overall\n"
             "energy shift for the strain corrections at this time.\n"
             "Adjust the strain parameters accordingly!!!\n",
             sp3d5s_Ebound_default);
      exit(-1);
    }
    StrainShift = sp3d5s_Ebound;
  } else {
    sp3d5s_Ebound = sp3d5s_Ebound_default;
    StrainShift = sp3d5s_Ebound;
  }

#endif

  nbasis = nml_iv_extent(basis);
  nn1 = rdim(nnv);
  nb = nbasis / 2;

  /* Set matrix to zero */
  for (i = 0; i < nbasis; i++)
    for (j = 0; j < nbasis; j++)
      h[i][j].r = h[i][j].i = hp1[i][j].r = hp1[i][j].i = 0.0;

  /* Construct Hamiltonian */
  for (i = 0; i < nbasis; i++) {
    /* On-site Energies */
    h[i][i].r = diag(param, basis[i]);
    for (j = 0; j < nbasis; j++) {
      /* Nearest-neigbor interactions */
      /* if ( basis[i] == 1 && basis[j] == 1 ) j++; */
      if (spin[i] == spin[j]) {
        /* Sum over neigbors */
        for (l = 0; l < nn1; l++) {
          cdum = ei(param, mb_strain, unstrnd_cubic_cell_length, k, nnv[l],
                    nnv0[l], basis[i], basis[j], &(h[i][i].r));
          if (nnv[l][2] > 0) {
            h[i][j].r += cdum.r;
            h[i][j].i += cdum.i;
          } else if (j >= nb && i < nb) {
            hp1[j][i].r -= cdum.r;
            hp1[j][i].i += cdum.i;
          }
        }
      }
      /* Spin-Orbit interactions */
      cdum = spin_orbit_spds(basis[i], basis[j], spin[i], spin[j], param[8],
                             param[9]);
      h[i][j].r += cdum.r;
      h[i][j].i += cdum.i;
    }
#ifdef STRAIN_SHIFT_ENABLE_1D
    h[i][i].r += StrainShift;
#endif
  }

  return (0);
}

int tbh_const_1d_nemo_old(cmatrix h, rvectr mb_strain,
                          real unstrnd_cubic_cell_length, cmatrix hp1,
                          rmatrix nnv, rmatrix nnv0, rvectr param, rvectr k) {
  int i, j, l, nn1, nbasis;
  complex cdum;

  nbasis = 20;
  nn1 = rdim(nnv);

  /* Set matrices to zero */
  for (i = 0; i < 2 * nbasis; i++)
    for (j = 0; j < 2 * nbasis; j++) {
      h[i][j].r = h[i][j].i = 0.0;
      hp1[i][j].r = hp1[i][j].i = 0.0;
    }

  /* Diagonal Elements */
  for (i = 0; i < nbasis; i++)
    h[i][i].r = diag(param, i);

  /* Off-Diagonal Elements */
  for (i = 0; i < nbasis / 2; i++)
    for (j = nbasis / 2; j < nbasis; j++)
      for (l = 0; l < nn1; l++) {
        cdum = ei(param, mb_strain, unstrnd_cubic_cell_length, k, nnv[l],
                  nnv0[l], i, j, &(h[i][i].r));
        if (nnv[l][2] > 0) {
          h[i][j].r += cdum.r;
          h[i][j].i += cdum.i;
        } else {
          hp1[j][i].r -= cdum.r;
          hp1[j][i].i += cdum.i;
        }
      }

  for (i = 0; i < nbasis / 2; i++)
    for (j = nbasis / 2; j < nbasis; j++)
      h[j][i] = cconj(h[i][j]);

  /* Add Spin Down */
  for (i = nbasis; i < 2 * nbasis; i++)
    for (j = nbasis; j < 2 * nbasis; j++) {
      h[i][j] = h[i - nbasis][j - nbasis];
      hp1[i][j] = hp1[i - nbasis][j - nbasis];
    }

  /* Spin Orbit */
  h[2][3].i = h[23][22].i = h[3][24].i = h[23][4].i = -param[8];
  h[3][2].i = h[22][23].i = h[24][3].i = h[4][23].i = param[8];
  h[4][22].r = h[22][4].r = -param[8];
  h[2][24].r = h[24][2].r = param[8];

  h[12][13].i = h[33][32].i = h[13][34].i = h[33][14].i = -param[9];
  h[13][12].i = h[32][33].i = h[34][13].i = h[14][33].i = param[9];
  h[14][32].r = h[32][14].r = -param[9];
  h[12][34].r = h[34][12].r = param[9];

  return (0);
}

int tbh_const_bulk_dk(cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis,
                      ivectr spin, rvectr param, rvectr mb_strain,
                      real unstrnd_cubic_cell_length, rvectr k, int dim1) {
  int i, j, l, nn1, nbasis;
  complex cdum;

  nbasis = nml_iv_extent(basis);
  nn1 = rdim(nnv);

  /* Set matrix to zero */
  for (i = 0; i < nbasis; i++)
    for (j = 0; j < nbasis; j++)
      h[i][j].r = h[i][j].i = 0.0;

  /* Construct Hamiltonian */
  for (i = 0; i < nbasis; i++) {
    for (j = 0; j < nbasis; j++) {
      /* Nearest-neigbor interactions */
      /* if ( basis[i] == 1 && basis[j] == 1 ) j++; */
      if (spin[i] == spin[j]) {
        /* Sum over neigbors */
        for (l = 0; l < nn1; l++) {
          cdum = eidk(param, mb_strain, unstrnd_cubic_cell_length, k, nnv[l],
                      nnv0[l], basis[i], basis[j], dim1, &(h[i][i].r));
          h[i][j].r += cdum.r;
          h[i][j].i += cdum.i;
        }
      }
    }
  }

  return (0);
}

int tbh_const_bulk_dk_old(cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis,
                          ivectr spin, rvectr param, rvectr mb_strain,
                          real unstrnd_cubic_cell_length, rvectr k, int dim1) {
  int i, j, l, nn1, nbasis;
  complex cdum;

  nbasis = nml_iv_extent(basis);
  nn1 = rdim(nnv);

  /* Set matrix to zero */
  for (i = 0; i < 2 * nbasis; i++)
    for (j = 0; j < 2 * nbasis; j++)
      h[i][j].r = h[i][j].i = 0.0;

  /* Off-Diagonal Elements */
  for (i = 0; i < nbasis / 2; i++)
    for (j = nbasis / 2; j < nbasis; j++)
      for (l = 0; l < nn1; l++) {
        cdum = eidk(param, mb_strain, unstrnd_cubic_cell_length, k, nnv[l],
                    nnv0[l], basis[i], basis[j], dim1, &(h[i][i].r));
        h[i][j].r += cdum.r;
        h[i][j].i += cdum.i;
      }

  for (i = 0; i < nbasis / 2; i++)
    for (j = nbasis / 2; j < nbasis; j++)
      h[j][i] = cconj(h[i][j]);

  /* Add Spin Down */
  for (i = nbasis; i < 2 * nbasis; i++)
    for (j = nbasis; j < 2 * nbasis; j++)
      h[i][j] = h[i - nbasis][j - nbasis];

  return (0);
}

real dE_dk(cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis, ivectr spin,
           rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length,
           rvectr k, int band, rvectr dE) {
  int i, j, ai;
  real e1, sum;
  complex csum, a;
  ivectr imap;
  cvectr eigv, x, y;
  cmatrix chi;

  imap = Ivectr(40);
  x = Cvectr(40);
  y = Cvectr(40);
  eigv = Cvectr(40);
  chi = Coperator(40);

  for (i = 0; i < 40; i++)
    imap[i] = i;

  tbh_const_bulk(h, nnv, nnv0, basis, spin, param, mb_strain,
                 unstrnd_cubic_cell_length, k);
  ceigvec_ns_full(eigv, chi, h);

  for (j = 1; j < 40; j++) {
    a = eigv[j];
    ai = imap[j];
    i = j - 1;
    while (i >= 0 && eigv[i].r > a.r) {
      imap[i + 1] = imap[i];
      eigv[i + 1] = eigv[i];
      i--;
    }
    imap[i + 1] = ai;
    eigv[i + 1] = a;
  }
  e1 = eigv[band].r;

  sum = 0.0;
  for (i = 0; i < 40; i++)
    x[i] = chi[i][imap[band]];
  for (i = 0; i < 40; i++)
    sum += x[i].r * x[i].r + x[i].i * x[i].i;
  sum = sqrt(sum);
  for (i = 0; i < 40; i++) {
    x[i].r /= sum;
    x[i].i /= sum;
  }

  for (j = 0; j < 3; j++) {
    csum.r = csum.i = 0.0;
    tbh_const_bulk_dk(h, nnv, nnv0, basis, spin, param, mb_strain,
                      unstrnd_cubic_cell_length, k, j);
    mul_cmatrix(y, h, x);
    for (i = 0; i < 40; i++)
      csum = cadd(csum, cmul(cconj(x[i]), y[i]));
    dE[j] = csum.r;
  }

  rm_ivectr(&imap);
  rm_cvectr(&eigv);
  rm_cvectr(&x);
  rm_cvectr(&y);
  rm_cmatrix(&chi);

  return e1;
}

real dE_dk_simple(cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis,
                  ivectr spin, rvectr param, rvectr mb_strain,
                  real unstrnd_cubic_cell_length, rvectr k, int band,
                  rvectr dE) {
  int i;
  real kdelta, e1;
  ivectr imap;
  cvectr eigv, x, y;
  cmatrix chi;

  kdelta = 2e-5;
  imap = Ivectr(40);
  x = Cvectr(40);
  y = Cvectr(40);
  eigv = Cvectr(40);
  chi = Coperator(40);

  for (i = 0; i < 40; i++)
    imap[i] = i;

  tbh_const_bulk(h, nnv, nnv0, basis, spin, param, mb_strain,
                 unstrnd_cubic_cell_length, k);
  ceigval_ns_full(eigv, h);
  e1 = eigv[band].r;

  for (i = 0; i < 3; i++) {
    k[i] -= kdelta;
    tbh_const_bulk(h, nnv, nnv0, basis, spin, param, mb_strain,
                   unstrnd_cubic_cell_length, k);
    k[i] += kdelta;
    ceigval_ns_full(eigv, h);
    dE[i] = (e1 - eigv[band].r) / kdelta;
  }

  rm_ivectr(&imap);
  rm_cvectr(&eigv);
  rm_cvectr(&x);
  rm_cvectr(&y);
  rm_cmatrix(&chi);

  return e1;
}

real mtensor(cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis, ivectr spin,
             rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length,
             rvectr k, int band, rmatrix mtens) {
  int i, j;
  real kdelta, E;
  rvectr dE1, dE2;

  kdelta = 1e-4;
  dE1 = Rvectr(3);
  dE2 = Rvectr(3);

  E = dE_dk(h, nnv, nnv0, basis, spin, param, mb_strain,
            unstrnd_cubic_cell_length, k, band, dE1);

  for (i = 0; i < 3; i++) {
    k[i] -= kdelta;
    dE_dk(h, nnv, nnv0, basis, spin, param, mb_strain,
          unstrnd_cubic_cell_length, k, band, dE2);
    k[i] += kdelta;
    for (j = 0; j < 3; j++)
      mtens[j][i] = (dE1[j] - dE2[j]) / kdelta;
  }

  /*
    for ( i=0; i < 3; i++ )
    for ( j=0; j < 3; j++ )
    mtens[i][j] = 2.0*HBAR_MASSFACTOR/mtens[i][j];
  */

  for (i = 0; i < 3; i++)
    mtens[i][i] = 2.0 * HBAR_MASSFACTOR / mtens[i][i];

  rm_rvectr(&dE1);
  rm_rvectr(&dE2);

  return E;
}

real Newton_micro(cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis,
                  ivectr spin, rvectr param, rvectr mb_strain,
                  real unstrnd_cubic_cell_length, rvectr k, int band,
                  rvectr mdE, rmatrix mtens) {
  int i, j;
  real kdelta, E;
  rvectr dE1, dE2;

  kdelta = 1e-4;
  dE1 = Rvectr(3);
  dE2 = Rvectr(3);

  E = dE_dk(h, nnv, nnv0, basis, spin, param, mb_strain,
            unstrnd_cubic_cell_length, k, band, dE1);

  for (i = 0; i < 3; i++) {
    mdE[i] = -1.0 * dE1[i];
    k[i] -= kdelta;
    dE_dk(h, nnv, nnv0, basis, spin, param, mb_strain,
          unstrnd_cubic_cell_length, k, band, dE2);
    k[i] += kdelta;
    for (j = 0; j < 3; j++)
      mtens[j][i] = (dE1[j] - dE2[j]) / kdelta;
  }

  rm_rvectr(&dE1);
  rm_rvectr(&dE2);

  return E;
}

/***********************************************************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

void E_const_prep(rvectr par_in, rvectr mb_strain, rvectr par_out,
                  real unstrnd_cubic_cell_length, real xd, real yd, real zd,
                  real *l, real *m, real *n, real *l2, real *m2, real *n2) {
  real den;

#if 0 
   /* Construct Direction Cosines and squares of unperturbed DCs */
    den = sqrt(xd0*xd0+yd0*yd0+zd0*zd0);
    *l0 = xd0/den;
    *m0 = yd0/den;
    *n0 = zd0/den;
    *l02 = *l0 * *l0;
    *m02 = *m0 * *m0;
    *n02 = *n0 * *n0;
#endif
  /* Construct Direction Cosines and squares of DCs */
  den = sqrt(xd * xd + yd * yd + zd * zd);
  *l = xd / den;
  *m = yd / den;
  *n = zd / den;
  *l2 = *l * *l;
  *m2 = *m * *m;
  *n2 = *n * *n;

#define STRAIN_SCALING
#ifdef STRAIN_SCALING
  /* Pass in a vector of exponents vexp[], order is indicated in Beltram */
  /* Pass in the lattice constant a in units of nm */
  /* Execute this only if vexp != NULL */
  if (mb_strain != NULL) {
    real brat, dummy;
    /* brat = equilibrium bond length divided by actual bond length */
    brat = (0.25 * sqr3 * unstrnd_cubic_cell_length) / den;
    /* Scale the two-center integrals */
    par_out[pV_sss] = par_in[pV_sss] * pow(brat, mb_strain[eta_sss]);
    par_out[pV_ststs] = par_in[pV_ststs] * pow(brat, mb_strain[eta_ststs]);
    dummy = pow(brat, mb_strain[eta_ssts]);
    par_out[pV_stass] = par_in[pV_stass] * dummy;
    par_out[pV_sasts] = par_in[pV_sasts] * dummy;
    dummy = pow(brat, mb_strain[eta_sps]);
    par_out[pV_sapcs] = par_in[pV_sapcs] * dummy;
    par_out[pV_scpas] = par_in[pV_scpas] * dummy;
    dummy = pow(brat, mb_strain[eta_stps]);
    par_out[pV_stapcs] = par_in[pV_stapcs] * dummy;
    par_out[pV_stcpas] = par_in[pV_stcpas] * dummy;
    dummy = pow(brat, mb_strain[eta_sds]);
    par_out[pV_sadcs] = par_in[pV_sadcs] * dummy;
    par_out[pV_scdas] = par_in[pV_scdas] * dummy;
    dummy = pow(brat, mb_strain[eta_stds]);
    par_out[pV_stadcs] = par_in[pV_stadcs] * dummy;
    par_out[pV_stcdas] = par_in[pV_stcdas] * dummy;
    par_out[pV_pps] = par_in[pV_pps] * pow(brat, mb_strain[eta_pps]);
    par_out[pV_ppp] = par_in[pV_ppp] * pow(brat, mb_strain[eta_ppp]);
    dummy = pow(brat, mb_strain[eta_pds]);
    par_out[pV_padcs] = par_in[pV_padcs] * dummy;
    par_out[pV_pcdas] = par_in[pV_pcdas] * dummy;
    dummy = pow(brat, mb_strain[eta_pdp]);
    par_out[pV_padcp] = par_in[pV_padcp] * dummy;
    par_out[pV_pcdap] = par_in[pV_pcdap] * dummy;
    par_out[pV_dds] = par_in[pV_dds] * pow(brat, mb_strain[eta_dds]);
    par_out[pV_ddp] = par_in[pV_ddp] * pow(brat, mb_strain[eta_ddp]);
    par_out[pV_ddd] = par_in[pV_ddd] * pow(brat, mb_strain[eta_ddd]);
  } else {
    par_out = par_in;
  }
#else
  par_out = par_in;
#endif

  return;
}

/*********************************************************************/

/* function provides the parameters for the off-site matrix element
   corrections due to strain.  */
void strain_corr_prep(rvectr par_in, rvectr mb_strain, rvectr par_scal_out,
                      real unstrnd_cubic_cell_length, real xd, real yd, real zd,
                      real xd0, real yd0, real zd0, real *l, real *m, real *n,
                      real *l2, real *m2, real *n2, real *l0, real *m0,
                      real *n0, real *l02, real *m02, real *n02) {
  real den;

  /* Construct Direction Cosines and squares of unperturbedDCs */
  den = sqrt(xd0 * xd0 + yd0 * yd0 + zd0 * zd0);
  *l0 = xd0 / den;
  *m0 = yd0 / den;
  *n0 = zd0 / den;
  *l02 = *l0 * *l0;
  *m02 = *m0 * *m0;
  *n02 = *n0 * *n0;

  /* Construct Direction Cosines and squares of DCs */
  den = sqrt(xd * xd + yd * yd + zd * zd);
  *l = xd / den;
  *m = yd / den;
  *n = zd / den;
  *l2 = *l * *l;
  *m2 = *m * *m;
  *n2 = *n * *n;

  /* Pass in a vector of exponents vexp[], order is indicated in Beltram */
  /* Pass in the lattice constant a in units of nm */
  /* Execute this only if vexp != NULL */
  if (mb_strain != NULL) {
    real brat, dummy, v, vp;

    /* brat = equilibrium bond length divided by actual bond length */
    brat = (0.25 * sqr3 * unstrnd_cubic_cell_length) / den;

    /* Scale the two-center integrals */
    dummy = pow(brat, mb_strain[eta_sss]);
    v = par_in[pV_sss];
    vp = v * dummy;
    par_scal_out[pV_sss] = vp;

    dummy = pow(brat, mb_strain[eta_ststs]);
    v = par_in[pV_ststs];
    vp = v * dummy;
    par_scal_out[pV_ststs] = vp;

    dummy = pow(brat, mb_strain[eta_ssts]);
    v = par_in[pV_stass];
    vp = v * dummy;
    par_scal_out[pV_stass] = vp;
    v = par_in[pV_sasts];
    vp = v * dummy;
    par_scal_out[pV_sasts] = vp;

    dummy = pow(brat, mb_strain[eta_sps]);
    v = par_in[pV_sapcs];
    vp = v * dummy;
    par_scal_out[pV_sapcs] = vp;
    v = par_in[pV_scpas];
    vp = v * dummy;
    par_scal_out[pV_scpas] = vp;

    dummy = pow(brat, mb_strain[eta_stps]);
    v = par_in[pV_stapcs];
    vp = v * dummy;
    par_scal_out[pV_stapcs] = vp;
    v = par_in[pV_stcpas];
    vp = v * dummy;
    par_scal_out[pV_stcpas] = vp;

    dummy = pow(brat, mb_strain[eta_sds]);
    v = par_in[pV_sadcs];
    vp = v * dummy;
    par_scal_out[pV_sadcs] = vp;
    v = par_in[pV_scdas];
    vp = v * dummy;
    par_scal_out[pV_scdas] = vp;

    dummy = pow(brat, mb_strain[eta_stds]);
    v = par_in[pV_stadcs];
    vp = v * dummy;
    par_scal_out[pV_stadcs] = vp;
    v = par_in[pV_stcdas];
    vp = v * dummy;
    par_scal_out[pV_stcdas] = vp;

    dummy = pow(brat, mb_strain[eta_pps]);
    v = par_in[pV_pps];
    vp = v * dummy;
    par_scal_out[pV_pps] = vp;

    dummy = pow(brat, mb_strain[eta_ppp]);
    v = par_in[pV_ppp];
    vp = v * dummy;
    par_scal_out[pV_ppp] = vp;

    dummy = pow(brat, mb_strain[eta_pds]);
    v = par_in[pV_padcs];
    vp = v * dummy;
    par_scal_out[pV_padcs] = vp;
    v = par_in[pV_pcdas];
    vp = v * dummy;
    par_scal_out[pV_pcdas] = vp;

    dummy = pow(brat, mb_strain[eta_pdp]);
    v = par_in[pV_padcp];
    vp = v * dummy;
    par_scal_out[pV_padcp] = vp;
    v = par_in[pV_pcdap];
    vp = v * dummy;
    par_scal_out[pV_pcdap] = vp;

    dummy = pow(brat, mb_strain[eta_dds]);
    v = par_in[pV_dds];
    vp = v * dummy;
    par_scal_out[pV_dds] = vp;

    dummy = pow(brat, mb_strain[eta_ddp]);
    v = par_in[pV_ddp];
    vp = v * dummy;
    par_scal_out[pV_ddp] = vp;

    dummy = pow(brat, mb_strain[eta_ddd]);
    v = par_in[pV_ddd];
    vp = v * dummy;
    par_scal_out[pV_ddd] = vp;
  } else {
    par_scal_out = par_in;
  }

  return;
}

/* This function converts nearest neighbor sp3s* parameters
   from the Nemo database into the two-center integral parameters
   used in run_3d.  Spin-orbit interactions are divided by 3,
   neighbor coupling is divided by 4 and the order is different */
rvectr param_sp3s_2_sp3d5s(rvectr param_in, rvectr param_out,
                           real unstrnd_cubic_cell_length, real delX, real delY,
                           real delZ) {
  int i, n = sp3d5s_nparam, n_short = sp3s_nparam;
  rvectr ptemp;

  ptemp = Rvectr(n_short);
  if (!param_out)
    param_out = Rvectr(n);

  real sqrt_3 = sqrt(3.0);

  for (i = 0; i < n_short; i++)
    ptemp[i] = param_in[i];

  /* Onsite energies are the same, so copy them. */
  for (i = 0; i < 6; i++)
    param_out[i] = param_in[i];
  for (i = 6; i < n; i++)
    param_out[i] = 0;

  /* spin-orbit */
  param_out[8] = ptemp[13] / 3;
  param_out[9] = ptemp[14] / 3;

  /* coupling */
  param_out[10] = ptemp[6] / 4;
  param_out[14] = ptemp[9] * sqrt_3 / 4;
  param_out[15] = ptemp[10] * sqrt_3 / 4;
  param_out[16] = ptemp[11] * sqrt_3 / 4;
  param_out[17] = ptemp[12] * sqrt_3 / 4;

  param_out[22] = (ptemp[7] + 2 * ptemp[8]) / 4;
  param_out[23] = (ptemp[7] - ptemp[8]) / 4;

  param_out[31] = unstrnd_cubic_cell_length;
  param_out[32] = delX;
  param_out[33] = delY;
  param_out[34] = delZ;

  rm_rvectr(&ptemp);

  return param_out;
}

complex vector_potential_phase(qd_struct d, int cindex, int aindex,
                               rvectr nnv) {

  rvectr R1 = Rvectr(3); // center atom position vector
  rvectr R2 = Rvectr(3); // neighbor atom position vector

  /* Following two lines were commented because of the following reason */
  /* If we change sign of nnv, Hamiltonian is not hermitian when more than one
     processors are used.  On every processor, upper triangular part of
     Hamiltonian is computed then Hermiticity of Hamiltonian is used for
     computing lower triangule part. When more than one processors are used
     off-diagonal blocks in upper triangular part reside on one processor and
     off-diagonal blocks in lower triangular part reside on neighboring
     processor. If we change the sign of nnv, then these blocks are equal. For
     Hamiltonian to be hermitian, these blocks should be complex conjugate of
     each other. With the following lines commented, they are complex conjugate
     and Hamiltonian is hermitian. I think this was a bug. */

  //  if ( d->geo.isCation(cindex,aindex) )
  //     nnv[0]*=-1.0; nnv[1]*=-1.0; nnv[2]*=-1.0; // bring back the original
  //     sign of nnv.
  // now nnv = R2-R1

  /* Shift origin to the center. Works only for Zincblende unit cells. */
  /* For other unit cells similar expression can be written. */
  double xminn = d->geo.lattice_x *
                 (d->geo.cell_xmin + d->geo.PositionInCell(6, 0) - 0.25);
  double xmaxx =
      d->geo.lattice_x * (d->geo.cell_xmax + d->geo.PositionInCell(0, 0));
  double yminn = d->geo.lattice_y *
                 (d->geo.cell_ymin + d->geo.PositionInCell(6, 1) - 0.25);
  double ymaxx =
      d->geo.lattice_y * (d->geo.cell_ymax + d->geo.PositionInCell(0, 1));
  double zminn = d->geo.lattice_z *
                 (d->geo.cell_zmin + d->geo.PositionInCell(7, 2) - 0.25);
  double zmaxx =
      d->geo.lattice_z * (d->geo.cell_zmax + d->geo.PositionInCell(0, 2));

  /* Origin in NEMO_3D co-ordinate system */
  double x0_n = (xminn + xmaxx) / 2;
  double y0_n = (yminn + ymaxx) / 2;
  double tz = d->geo.lattice_z * (d->geo.cell_zmax - d->geo.cell_zmin + 1 +
                                  d->geo.dz - 1 /* Not clear! */) -
              d->geo.adjThicknessZ;
  //  double z0_n = (x0_n/*/2*/)*(double(d->geo.dz)/double(d->geo.dx))*
  //     (d->geo.lattice_z/d->geo.lattice_x) + tz/2;
  double z0_n = (zminn + zmaxx) / 2;
  /* End: shift origin to the center */

  /* Real space co-ordinates of atom */
  R1[0] = d->geo.lattice_x *
              (d->geo.cell__ijk[cindex][0] + d->geo.PositionInCell(aindex, 0)) +
          d->geo.l_m_2_atomdsp[cindex][aindex][0];
  R1[1] = d->geo.lattice_y *
              (d->geo.cell__ijk[cindex][1] + d->geo.PositionInCell(aindex, 1)) +
          d->geo.l_m_2_atomdsp[cindex][aindex][1];
  R1[2] = d->geo.lattice_z *
              (d->geo.cell__ijk[cindex][2] + d->geo.PositionInCell(aindex, 2)) +
          d->geo.l_m_2_atomdsp[cindex][aindex][2];
  /* End: Real space co-ordinates of atom */

  /* Real space co-ordinates of neighbor */
  for (int i = 0; i < 3; i++) {
    /* Commented because following calculation is redundant.
        R1[i]= d->geo.lattice_x*
           (d->geo.cell__ijk[cindex][i]+d->geo.PositionInCell(aindex,i))
             + d->geo.l_m_2_atomdsp[cindex][aindex][i];
    */
    R2[i] = R1[i] + nnv[i];
  }
  /* End: Real space co-ordinates of neighbor */

  double Bx = d->opt.ExecParam.Phys.Bx;
  double By = d->opt.ExecParam.Phys.By;
  double Bz = d->opt.ExecParam.Phys.Bz;

  /* Rotate NEMO_3D co-ordinate system by tilt angle */
  double theta = 0.0; // atan(double(d->geo.dz)/double(d->geo.dx));
  /* if vector potential is defined in rotated gauge NEMO_3D doesn't converge on
     multiple processors.  Smaller tilt angle has negligible effect on rotations
     which can be ignored. */

  /* Origin in rotated co-ordinate system */

  /*  double x0 = x0_n*cos(theta) - z0_n*sin(theta);
    double y0 = y0_n;
    double z0 = x0_n*sin(theta) + z0_n*cos(theta);
  */

  /* Origin in rotated co-ordinate system */
  double x0 = x0_n;
  double y0 = y0_n;
  double z0 = z0_n;

  /* Components of vector potential */
  // B-field in z
  /*
    double A_0_R1 = -Bz*(R1[1]-y0);
    double A_1_R1 = Bz*(R1[0]*cos(theta)-R1[2]*sin(theta)-x0);
    double A_0_R2 = -Bz*(R2[1]-y0);
    double A_1_R2 = Bz*(R2[0]*cos(theta)-R2[2]*sin(theta)-x0);
  */

  // B-field in x

  /*
    double A_0_R1 = -Bx*(R1[2]-z0);
    double A_1_R1 = Bx*(R1[1]-y0);
    double A_0_R2 = -Bx*(R2[2]-z0);
    double A_1_R2 = Bx*(R2[1]-y0);
  */

  // B-field in y

  /*
    double A_0_R1 = By*(R1[2]-z0);
    double A_1_R1 = -By*(R1[0]-x0);
    double A_0_R2 = By*(R2[2]-z0);
    double A_1_R2 = -By*(R2[0]-x0);
  */

  // B-field in 3D (Bx, By, Bz)

  double A_0_R1 = By * (R1[2] - z0) - Bz * (R1[1] - y0);
  double A_0_R2 = By * (R2[2] - z0) - Bz * (R2[1] - y0);
  double A_1_R1 = Bz * (R1[0] - x0) - Bx * (R1[2] - z0);
  double A_1_R2 = Bz * (R2[0] - x0) - Bx * (R2[2] - z0);
  double A_2_R1 = Bx * (R1[1] - y0) - By * (R1[0] - x0);
  double A_2_R2 = Bx * (R2[1] - y0) - By * (R2[0] - x0);

  rvectr A = Rvectr(3);
  /* Choose appropriate gauge from following options */
  if (d->opt.ExecParam.Phys.Bgauge == "Symmetric") {
    /* Vector potential in symmetric gauge */

    // B-field in z

    /*     A[0] = 0.25*(A_0_R1 + A_0_R2);
         A[1] = 0.25*(A_1_R1 + A_1_R2);
         A[2] = 0.0;
    */

    // B-field in x
    /*     A[0] = 0.0; //0.25*(A_0_R1 + A_0_R2);
         A[1] = 0.25*(A_0_R1 + A_0_R2);
         A[2] = 0.25*(A_1_R1 + A_1_R2);
    */

    // B-field in 3D

    A[0] = 0.25 * (A_0_R1 + A_0_R2);
    A[1] = 0.25 * (A_1_R1 + A_1_R2);
    A[2] = 0.25 * (A_2_R1 + A_2_R2);

  } else if (d->opt.ExecParam.Phys.Bgauge == "Asymmetric_x") {
    /* Vector potential in asymmetric gauge */
    A[0] = 0.0;
    A[1] = 0.5 * (A_1_R1 + A_1_R2);
    A[2] = 0.0;
  } else if (d->opt.ExecParam.Phys.Bgauge == "Asymmetric_y") {
    /* Vector potential in asymmetric gauge */
    A[0] = 0.5 * (A_0_R1 + A_0_R2);
    A[1] = 0.0;
    A[2] = 0.0;
  } else
    die("Unknown magnetic field gauge");

  /* double MagneticFlux =
     ((R2[0]-R1[0])*cos(theta)+(R2[2]-R1[2])*sin(theta))*A[0]
                          + (R2[1]-R1[1])*A[1]
                          +
     ((R2[0]-R1[0])*sin(theta)+(R2[2]-R1[2])*cos(theta))*A[2];
  */

  double MagneticFlux =
      (R2[0] - R1[0]) * A[0] + (R2[1] - R1[1]) * A[1] + (R2[2] - R1[2]) * A[2];

  rm_rvectr(&R1);
  rm_rvectr(&R2);
  rm_rvectr(&A);

  double FluxRatio = (MagneticFlux / MagneticFluxQuantum);
  double Integral = FluxRatio * M_PI;
  complex PeirlsSubstitution;
  PeirlsSubstitution.r = cos(Integral);
  PeirlsSubstitution.i = -sin(Integral);
  return PeirlsSubstitution;
}

int apply_magnetic_field(qd_struct d, int cindex, int aindex) {
  // printf("\ntest: cindex = %d   aindex = %d\n",cindex,aindex);
  if (d->BandModel == BM_20_sp3d5ss_spin || d->BandModel == BM_10_sp3ss_spin) {
    int nb = d->NBasisStates / 2;

    double muBx = d->opt.ExecParam.Phys.Bx * BohrMagneton;
    double muBy = d->opt.ExecParam.Phys.By * BohrMagneton;
    double muBz = d->opt.ExecParam.Phys.Bz * BohrMagneton;

    for (int i = 0; i < nb; i++) {
      d->hd[i][i].r += muBz;
      d->hd[i + nb][i + nb].r -= muBz;
      d->hd[i][i + nb].r += muBx;
      d->hd[i][i + nb].i -= muBy;
      d->hd[i + nb][i].r += muBx;
      d->hd[i + nb][i].i += muBy;
    }
  }

  return (0);
}

int apply_magnetic_field(qd_struct d, int cindex, int aindex, rvectr nnv) {
  if (d->BandModel == BM_20_sp3d5ss_spin || d->BandModel == BM_10_sp3ss_spin) {
    complex phase = vector_potential_phase(d, cindex, aindex, nnv);
    int nb = d->NBasisStates / 2;
    for (int irow = 0; irow < nb; irow++) {
      for (int icol = 0; icol < nb; icol++) {
        /* Commented/added for fixing Magnetic field */
        /*            d->ho[irow][icol].i = d->ho[irow][icol].r*phase.i;
                    d->ho[irow][icol].r = d->ho[irow][icol].r*phase.r;
                    d->ho[irow+nb][icol+nb].i =
           d->ho[irow+nb][icol+nb].r*phase.i; d->ho[irow+nb][icol+nb].r =
           d->ho[irow+nb][icol+nb].r*phase.r;
        */
        double ho_i =
            d->ho[irow][icol].i * phase.r + d->ho[irow][icol].r * phase.i;
        double ho_r =
            d->ho[irow][icol].r * phase.r - d->ho[irow][icol].i * phase.i;
        d->ho[irow][icol].i = ho_i;
        d->ho[irow][icol].r = ho_r;
        ho_i = d->ho[irow + nb][icol + nb].i * phase.r +
               d->ho[irow + nb][icol + nb].r * phase.i;
        ho_r = d->ho[irow + nb][icol + nb].r * phase.r -
               d->ho[irow + nb][icol + nb].i * phase.i;
        d->ho[irow + nb][icol + nb].i = ho_i;
        d->ho[irow + nb][icol + nb].r = ho_r;
        /* End Commented/added for fixing Magnetic field */
      }
    }
  }

  else if (d->BandModel == BM_1_s_nospin ||
           d->BandModel == BM_10_sp3d5ss_nospin) {
    complex phase = vector_potential_phase(d, cindex, aindex, nnv);
    int nb = d->NBasisStates;
    for (int irow = 0; irow < nb; irow++) {
      for (int icol = 0; icol < nb; icol++) {
        /* Commented/added for fixing Magnetic field */
        /*
                    d->ho[irow][icol].i = d->ho[irow][icol].r*phase.i;
                    d->ho[irow][icol].r = d->ho[irow][icol].r*phase.r;
        */
        double ho_i =
            d->ho[irow][icol].i * phase.r + d->ho[irow][icol].r * phase.i;
        double ho_r =
            d->ho[irow][icol].r * phase.r - d->ho[irow][icol].i * phase.i;
        d->ho[irow][icol].i = ho_i;
        d->ho[irow][icol].r = ho_r;
        /* End Commented/added for fixing Magnetic field */
      }
    }
  }

  else {
    die("apply_magnetic_field is not implemented for this band model");
  }

  return (0);
}

int compute_strain_tensor(real *exx, real *eyy, real *ezz, real *exy, real *exz,
                          real *eyz, qd_struct d, int Zc, int Za) {
  bool has_periodicity = d->geo.hasPeriodicity();
  int b, m, n, p;
  real dl;
  ivectr cell_nbr = NULL;
  real dR[4][3], dR0[4][3];

  int map_mm_to_n[8][8];
  int map_mb_to_m[8][4]; /* map (cation m, block label of nbr) to m of nbr */
  int map_n_to_l[4];

  real a_lattice[3];
  a_lattice[0] = d->geo.lattice_x;
  a_lattice[1] = d->geo.lattice_y;
  a_lattice[2] = d->geo.lattice_z;

  /* need map from nbr indx (n) to atom indx (m) */
  for (m = 0; m < d->geo.AtomsPerCellMax(); m++) {
    for (n = 0; n < d->geo.Neighbors(m); n++) {
      int m_nbr = d->geo.NbrCell(m, n, 3);
      map_mm_to_n[m][m_nbr] = n;
    }
  }

  for (m = 0; m < 4; m++) {
    map_mb_to_m[m][0] = 4 + m;
    map_mb_to_m[m][1] = 4 + (7 - m) % 4;
    map_mb_to_m[m][2] = 4 + (2 + m) % 4;
    map_mb_to_m[m][3] = 4 + (5 - m) % 4;
  }
  for (m = 4; m < 8; m++) {
    map_mb_to_m[m][0] = (m - 2) % 4;
    map_mb_to_m[m][1] = (9 - m) % 4;
    map_mb_to_m[m][2] = m % 4;
    map_mb_to_m[m][3] = (7 - m) % 4;
  }

  int flag;
  int i = d->geo.cell__ijk[Zc][0];
  int j = d->geo.cell__ijk[Zc][1];
  int k = d->geo.cell__ijk[Zc][2];

  if (!d->geo.AtomType[Zc][Za])
    return -1;

  /* check that (l,m) has all its neighbors */
  flag = 0;
  for (n = 0; n < d->geo.Neighbors(Za); n++) {
    int i_nbr = i + d->geo.NbrCell(Za, n, 0);
    int j_nbr = j + d->geo.NbrCell(Za, n, 1);
    int k_nbr = k + d->geo.NbrCell(Za, n, 2);

    int aindx = d->geo.NbrCell(Za, n, 3);
    int cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

#ifdef ELIMINATE_SSMAP
    if (cindx < 0 && has_periodicity && (cindx + d->geo.N_Cell + 1) >= 0) {
      cindx = cindx + d->geo.N_Cell + 1;
    }
    if (d->geo.is_tilted)
      cindx = d->geo.get_cindx_tilted(Zc, Za, n);
#else  /* ELIMINATE_SSMAP */
    /* Commented for ssmap. Addition for tilt */
    if (/*cindx < 0 &&*/ has_periodicity && d->geo.ssmap[Zc][Za][n] != -1) {
      cindx = d->geo.ssmap[Zc][Za][n];
    }
#endif /* ELIMINATE_SSMAP */

    map_n_to_l[n] = cindx;

    if (cindx < 0 || !d->geo.AtomType[cindx][aindx]) {
      flag = 1;
      break;
    }
  }

  if (flag)
    return -1; // incomplete set of neighbors found

#ifdef DISPL_SINGLE_PRECISION
  fvectr disp_this = d->geo.l_m_2_atomdsp[Zc][Za];
#else  /* DISPL_SINGLE_PRECISION */
  rvectr disp_this = d->geo.l_m_2_atomdsp[Zc][Za];
#endif /* DISPL_SINGLE_PRECISION */

  for (b = 0; b < d->geo.Neighbors(Za); b++) {
    int aindx = map_mb_to_m[Za][b];
    n = map_mm_to_n[Za][aindx];
    int cindx = map_n_to_l[n];

    real alat = d->strain.amat[(int)d->geo.AtomType[Zc][Za]]
                              [(int)d->geo.AtomType[cindx][aindx]];
#ifdef DISPL_SINGLE_PRECISION
    fvectr disp_nbr = d->geo.l_m_2_atomdsp[cindx][aindx];
#else  /* DISPL_SINGLE_PRECISION */
    rvectr disp_nbr = d->geo.l_m_2_atomdsp[cindx][aindx];
#endif /* DISPL_SINGLE_PRECISION */

    rvectr rAtom_nbr = d->geo.PositionInCell(aindx);
    rvectr rAtom_this = d->geo.PositionInCell(Za);
    cell_nbr = d->geo.NbrCell(Za, n);

    p = 0;
    for (p = 0; p < 3; p++) {
      dl = cell_nbr[p] + rAtom_nbr[p] - rAtom_this[p];
      dR[b][p] = a_lattice[p] * dl + (disp_nbr[p] - disp_this[p]);
      dR0[b][p] = alat * dl;
    }
  }

  rvectr a = Rvectr(9);
  rvectr a0 = Rvectr(9);

  for (int i = 0; i < 3; i++) {
    a[i * 3 + 0] = dR[1][i] + dR[2][i] - dR[3][i] - dR[0][i];
    a[i * 3 + 1] = dR[2][i] + dR[3][i] - dR[1][i] - dR[0][i];
    a[i * 3 + 2] = dR[1][i] + dR[3][i] - dR[2][i] - dR[0][i];

    a0[i * 3 + 0] = dR0[1][i] + dR0[2][i] - dR0[3][i] - dR0[0][i];
    a0[i * 3 + 1] = dR0[2][i] + dR0[3][i] - dR0[1][i] - dR0[0][i];
    a0[i * 3 + 2] = dR0[1][i] + dR0[3][i] - dR0[2][i] - dR0[0][i];
  }

  if (!d->geo.isCation(Zc, Za)) {
    for (int i = 0; i < 3; i++) {
      a[i * 3 + 2] = -a[i * 3 + 2];
      a0[i * 3 + 2] = -a0[i * 3 + 2];
    }
  }

  int info;
  ivectr ipiv = Ivectr(9);

  //  strain_tensor = a*inv(a0) - 1;
  int order = 3;
  dgesv(&order, &order, &a0[0], &order, &ipiv[0], &a[0], &order, &info);
  *exx = a[0] - 1.0;
  *eyy = a[4] - 1.0;
  *ezz = a[8] - 1.0;
  *exy = a[1];
  *exz = a[2];
  *eyz = a[5];

  /*
  printf("[status=%d] %d  %d %e %e %e %e %e %e\n",
         info, Zc, Za, *exx, *eyy, *ezz, *exy, *exz, *eyz);
  */
  rm_ivectr(&ipiv);
  rm_rvectr(&a);
  rm_rvectr(&a0);

  return 0;
}

int on_site_d_orbital_shift(qd_struct d, int Zc, int Za, real E_d, real b_d,
                            real *Exy, real *Exz, real *Eyz) {
  double exx, eyy, ezz, exy, exz, eyz;

  if (compute_strain_tensor(&exx, &eyy, &ezz, &exy, &exz, &eyz, d, Zc, Za) ==
      -1)
    return -1;

  double eps_zz_xx = ezz - exx;

  *Exy = E_d * (1.0 + 2.0 * b_d * eps_zz_xx);
  *Exz = E_d * (1.0 - b_d * eps_zz_xx);
  *Eyz = *Exz;

  return 0;
}

rvectr eigval_h_full(rvectr d, coperator a, int *error_code) {
  fprintf(stderr, "*** Error: eigval_h_full is unimplemented. ***\n");
  *error_code = 1;
  return d;
}
