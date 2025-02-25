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
$Header: /repo/nemo3d/src/base/StrainParam.c,v 1.12 2005/01/14 20:24:38 marek
Exp $
*****************************************************************************/

#include "StrainParam.h"
#include "mat_def.h"
#include "run3d_mpi.h"
#include "util_messages.h"

#define how_close_to_maximum_of_distorsion_energy 0.001
// #define PR_WARNING_ON_VFFC_TRUNCATION

/**
 * Get the effective value of alpha scaled with strain.
 * @param m_atom_type is the type of the m-th atom
 * @param n_atom_type is the type of the n-th atom
 * @param d2 is the square of the relaxed m-n2 bond of the V doublet
 * @param r2 is the square of the strained m-n2 bond of the V doublet
 * @param anh is an integer 0/1/2 specifying the choice of the anharmonicity
 * corrections
 * @param pr is aboolean value, specifying either to print the wornings at
 * truncation
 * @return the effective value of alphaa scaled with strain
 */
real StrainParam::alpha_eff(int m_atom_type, // type of the m-th atom
                            int n_atom_type, // type of the n-th atom
                            real r2,         // square of the real bond length
                            real d2, // square of the equilibrium bond length
                            int anh, // 1-use acoustical anharmonicity
                                     // corrections 2-use optical anharmonicity
                                     // corrections 3-use linear interpolation
                            bool pr // 1/0=do/not print if alpha is truncated
                            )
// Get the effective value of alpha scaled with strain
{
  real alpha;
  real delt = 1 - how_close_to_maximum_of_distorsion_energy;
  real A;
  if (anh == 1) {
    A = this->aa_VFF[m_atom_type][n_atom_type];
  } else if (anh == 2) {
    A = this->ao_VFF[m_atom_type][n_atom_type];
  } else {
    if (pr && !mpi_n3d_id) {
      printf("\nOption anh=%d has not been implemented yet!!!\n"
             "Use anharmonicity corrections for better strain distribution "
             "instead (anh=0).",
             anh);
    }
    A = this->aa_VFF[m_atom_type][n_atom_type];
  }
  // Check if close to the maximum
  if (r2 < (d2 * (1 + delt * 2 / (3 * A))))
    alpha = (this->alpha[m_atom_type][n_atom_type] * (1 - A * (r2 / d2 - 1)));
  else { // Truncate near the maximum
    alpha = (this->alpha[m_atom_type][n_atom_type] * (1 - delt / 1.5));
#ifdef PR_WARNING_ON_VFFC_TRUNCATION
    if (pr && !mpi_n3d_id) {
      printf("!!!%d-%d bond : sqrt(r2/d2)=%g/%g=%g : alpha=%g\n", m_atom_type,
             n_atom_type, sqrt(r2), sqrt(d2), sqrt(r2 / d2), alpha);
    }
#endif // PR_WARNING_ON_VFFC_TRUNCATION
  } // Truncate near the maximum

  return (alpha);
}

/**
 * Get the effective value of beta scaled with strain.
 * @param m_atom_type is the type of the m-th atom
 * @param n1_atom_type is the type of the n1-th atom
 * @param n2_atom_type is the type of the n2-th atom
 * @param d1 is the vector of the relaxed m-n1 bond of the V doublet
 * @param d2 is the vector of the relaxed m-n2 bond of the V doublet
 * @param r1 is the vector of the strained m-n1 bond of the V doublet
 * @param r2 is the vector of the strained m-n2 bond of the V doublet
 * @param anh is an integer 0/1/2 specifying the choice of the anharmonicity
 * corrections
 * @param pr is aboolean value, specifying either to print the wornings at
 * truncation
 * @return the effective value of beta scaled with strain
 */
real StrainParam::beta_eff(
    int m_atom_type,      // type of the m-th atom
    int n1_atom_type,     // type of the n1-th atom
    int n2_atom_type,     // type of the n2-th atom
    rvectr d1, rvectr d2, // vectors of the relaxed bonds of the V doublet
    rvectr r1, rvectr r2, // vectors of the strained bonds of the V doublet
    int anh,              // 1-use acoustical anharmonicity corrections
    // 2-use optical anharmonicity corrections
    // 3-use linear interpolation
    bool pr // 1/0=do/not print if beta is truncated
    )
// Get the effective value of beta scaled with strain
{
  real d10 = sqrt(d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]);
  real d20 = sqrt(d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]);
  real rr = r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2];
  real r10 = sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
  real r20 = sqrt(r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]);
  real rd = r10 * r20 / (d10 * d20);
  real tt = rr / (r10 * r20); //=cos(theta)
  // cos(theta0)=-1/3 in ideal diamond and zinc-blende crystal
  real be0 = sqrt(this->beta[m_atom_type][n1_atom_type] *
                  this->beta[m_atom_type][n1_atom_type]);
  real beta = be0;
  real B, C;
  if (anh == 1) {
    B = sqrt(this->ba_VFF[m_atom_type][n1_atom_type] *
             this->ba_VFF[m_atom_type][n2_atom_type]);
    C = sqrt(this->ca_VFF[m_atom_type][n1_atom_type] *
             this->ca_VFF[m_atom_type][n2_atom_type]);
  } else if (anh == 2) {
    B = sqrt(this->bo_VFF[m_atom_type][n1_atom_type] *
             this->bo_VFF[m_atom_type][n2_atom_type]);
    C = sqrt(this->co_VFF[m_atom_type][n1_atom_type] *
             this->co_VFF[m_atom_type][n2_atom_type]);
  } else {
    if (pr && !mpi_n3d_id) {
      printf("\nOption anh=%d has not been implemented yet!!!\n"
             "Use anharmonicity corrections for better strain distribution "
             "instead (anh=0).",
             anh);
    }
    B = sqrt(this->ba_VFF[m_atom_type][n1_atom_type] *
             this->ba_VFF[m_atom_type][n2_atom_type]);
    C = sqrt(this->ca_VFF[m_atom_type][n1_atom_type] *
             this->ca_VFF[m_atom_type][n2_atom_type]);
  }
  //////////////

  real delt = 1 + how_close_to_maximum_of_distorsion_energy;
  real delh = 1 - how_close_to_maximum_of_distorsion_energy;
  bool flag = true;

  // Check if close to max shear strain
  if (tt > delt / 3.0 * (2 * (1 / B - 1 / 3.0) - 1 / (3 * rd))) {
    // truncate near the maximum of the angle dependence
    beta =
        be0 *
        ((1.0 - 2 * delt / 3.0) + B / 3.0 * (delt / 3.0 * (1 / rd + 2) - 1)) *
        (1 - C * (rd - 1));
    flag = false;
#ifdef PR_WARNING_ON_VFFC_TRUNCATION
    if (pr && !mpi_n3d_id) {
      printf("  !!!r1/d1=%g/%g=%g r2/d2=%g/%g=%g cos(tt)=%g => !!angle!!: "
             "beta=%g\n",
             r10, d10, r10 / d10, r20, d20, r20 / d20, tt, beta);
    }
#endif // PR_WARNING_ON_VFFC_TRUNCATION
  }

  // Ckeck if close to max hydrostatic strain
  if (rd > (1 + delh * 2 / (3 * C))) {
    // truncate near the maximum of the bond length dependence
    real be = be0 * (1 - B * (tt + 1 / 3.0)) * (1 - delh * 2 / 3.0);
    if (flag || (!flag && (be > beta))) {
      beta = be;
#ifdef PR_WARNING_ON_VFFC_TRUNCATION
      if (pr && !mpi_n3d_id) {
        printf("  !!!r1/d1=%g/%g=%g r2/d2=%g/%g=%g cos(tt)=%g => !!length!!: "
               "beta=%g\n",
               r10, d10, r10 / d10, r20, d20, r20 / d20, tt, beta);
      }
#endif // PR_WARNING_ON_VFFC_TRUNCATION
    }
    flag = false;
  }

  // Scale beta with strain
  if (flag) {
    beta = be0 * (1 - B * (tt + 1 / 3.0)) * (1 - C * (rd - 1));
  }

  return (beta);
}

/*  Set up data structures amat, alpha, beta */
void StrainParam::SetMatParm(vector<MaterialHandle> sMatList) {
  extern int str_to_atomID(const char *);

  this->amat = Rmatrix(EL_TOTAL, EL_TOTAL);
  this->alpha = Rmatrix(EL_TOTAL, EL_TOTAL);
  this->beta = Rmatrix(EL_TOTAL, EL_TOTAL);
  this->aa_VFF = Rmatrix(
      EL_TOTAL,
      EL_TOTAL); // Coefficients of alpha strain dependence at low energy
  this->ba_VFF =
      Rmatrix(EL_TOTAL,
              EL_TOTAL); // Coefficients of beta strain dependence at low energy
  this->ca_VFF =
      Rmatrix(EL_TOTAL,
              EL_TOTAL); // Coefficients of beta strain dependence at low energy
  this->ao_VFF = Rmatrix(
      EL_TOTAL,
      EL_TOTAL); // Coefficients of alpha strain dependence at high energy
  this->bo_VFF = Rmatrix(
      EL_TOTAL,
      EL_TOTAL); // Coefficients of beta strain dependence at high energy
  this->co_VFF = Rmatrix(
      EL_TOTAL,
      EL_TOTAL); // Coefficients of beta strain dependence at high energy

  for (vector<MaterialHandle>::iterator sm = sMatList.begin();
       sm != sMatList.end(); sm++) {

    const Material_struct *const smh = sm->hndl;
    if (!smh)
      continue;

    int anion = str_to_atomID(sm->anion.c_str());
    int cation = str_to_atomID(sm->cation.c_str());

    if (smh->unstrnd_cubic_cell_length <= 0.0 &&
        mpi_n3d_id == mpi_n3d_masterid) {
      n3d_warning("Bad unstrained cubic cell length=%g for the material %s\n"
                  "Will run into trouble in a strain calculation\n",
                  smh->unstrnd_cubic_cell_length, smh->material.c_str());
    }

    if (smh->strain_alpha <= 0.0 && mpi_n3d_id == mpi_n3d_masterid) {
      n3d_warning("Bad alpha strain value=%g for the material %s\n"
                  "Will run into trouble in a strain calculation\n",
                  smh->strain_alpha, smh->material.c_str());
    }

    if (smh->strain_beta <= 0.0 && mpi_n3d_id == mpi_n3d_masterid) {
      n3d_warning("Bad beta strain value=%g for the material %s\n"
                  "Will run into trouble in a strain calculation\n",
                  smh->strain_beta, smh->material.c_str());
    }

    this->alpha[anion][cation] = this->alpha[cation][anion] = smh->strain_alpha;
    this->beta[anion][cation] = this->beta[cation][anion] = smh->strain_beta;
    this->aa_VFF[anion][cation] = this->aa_VFF[cation][anion] =
        smh->strain_VFF[0];
    this->ba_VFF[anion][cation] = this->ba_VFF[cation][anion] =
        smh->strain_VFF[1];
    this->ca_VFF[anion][cation] = this->ca_VFF[cation][anion] =
        smh->strain_VFF[2];
    this->ao_VFF[anion][cation] = this->ao_VFF[cation][anion] =
        smh->strain_VFF[3];
    this->bo_VFF[anion][cation] = this->bo_VFF[cation][anion] =
        smh->strain_VFF[4];
    this->co_VFF[anion][cation] = this->co_VFF[cation][anion] =
        smh->strain_VFF[5];
    this->amat[anion][cation] = this->amat[cation][anion] =
        smh->unstrnd_cubic_cell_length;
  }
}

void StrainParam::Deallocate() {
  if (this->amat) {
    rm_rmatrix(&this->amat);
    this->amat = NULL;
  }
  if (this->alpha) {
    rm_rmatrix(&this->alpha);
    this->alpha = NULL;
  }
  if (this->beta) {
    rm_rmatrix(&this->beta);
    this->beta = NULL;
  }
  if (this->aa_VFF) {
    rm_rmatrix(&this->aa_VFF);
    this->aa_VFF = NULL;
  }
  if (this->ba_VFF) {
    rm_rmatrix(&this->ba_VFF);
    this->ba_VFF = NULL;
  }
  if (this->ca_VFF) {
    rm_rmatrix(&this->ca_VFF);
    this->ca_VFF = NULL;
  }
  if (this->ao_VFF) {
    rm_rmatrix(&this->ao_VFF);
    this->ao_VFF = NULL;
  }
  if (this->bo_VFF) {
    rm_rmatrix(&this->bo_VFF);
    this->bo_VFF = NULL;
  }
  if (this->co_VFF) {
    rm_rmatrix(&this->co_VFF);
    this->co_VFF = NULL;
  }
  if (this->lm_sndL) {
    rm_ivectr(&this->lm_sndL);
    this->lm_sndL = NULL;
  }
  if (this->lm_sndR) {
    rm_ivectr(&this->lm_sndR);
    this->lm_sndR = NULL;
  }
  if (this->lm_rcvL) {
    rm_ivectr(&this->lm_rcvL);
    this->lm_rcvL = NULL;
  }
  if (this->lm_rcvR) {
    rm_ivectr(&this->lm_rcvR);
    this->lm_rcvR = NULL;
  }
  if (this->rcvL_lm) {
    rm_imatrix(&this->rcvL_lm);
    this->rcvL_lm = NULL;
  }
  if (this->rcvR_lm) {
    rm_imatrix(&this->rcvR_lm);
    this->rcvR_lm = NULL;
  }
}

void StrainParam::SetOptions(const top_struct &opt) {
  const Keating_struct &keating = opt.ExecParam.Strain.Keating;

  if (keating.MinWrtLatt == Keating_struct::x) {
    MinWrtX = true;
    MinWrtY = false;
    MinWrtZ = false;
  } else if (keating.MinWrtLatt == Keating_struct::y) {
    MinWrtX = false;
    MinWrtY = true;
    MinWrtZ = false;
  } else if (keating.MinWrtLatt == Keating_struct::z) {
    MinWrtX = false;
    MinWrtY = false;
    MinWrtZ = true;
  } else if (keating.MinWrtLatt == Keating_struct::xy) {
    MinWrtX = true;
    MinWrtY = true;
    MinWrtZ = false;
  } else if (keating.MinWrtLatt == Keating_struct::xz) {
    MinWrtX = true;
    MinWrtY = false;
    MinWrtZ = true;
  } else if (keating.MinWrtLatt == Keating_struct::yz) {
    MinWrtX = false;
    MinWrtY = true;
    MinWrtZ = true;
  } else if (keating.MinWrtLatt == Keating_struct::xyz) {
    MinWrtX = true;
    MinWrtY = true;
    MinWrtZ = true;
  } else {
    MinWrtX = false;
    MinWrtY = false;
    MinWrtZ = false;
  }

  // we should not minimize the lattice constant in a direction
  // in which we do not have periodic boundary conditions.
  if (keating.Periodicity_strain == Keating_struct::SameAsGlobal) {
    const enum ElGeo_struct::PeriodicityType &bc =
        opt.ExecParam.ElCalc.ElGeo.Periodicity;

    MinWrtX =
        MinWrtX &&
        (bc == ElGeo_struct::Periodic_x || bc == ElGeo_struct::Periodic_xy ||
         bc == ElGeo_struct::Periodic_xz || bc == ElGeo_struct::Periodic_xyz);

    MinWrtY =
        MinWrtY &&
        (bc == ElGeo_struct::Periodic_y || bc == ElGeo_struct::Periodic_xy ||
         bc == ElGeo_struct::Periodic_yz || bc == ElGeo_struct::Periodic_xyz);

    MinWrtZ =
        MinWrtZ &&
        (bc == ElGeo_struct::Periodic_z || bc == ElGeo_struct::Periodic_yz ||
         bc == ElGeo_struct::Periodic_xz || bc == ElGeo_struct::Periodic_xyz);

  } else {
    const enum Keating_struct::Periodicity_strainType &bc =
        keating.Periodicity_strain;

    MinWrtX = MinWrtX && (bc == Keating_struct::Periodic_x ||
                          bc == Keating_struct::Periodic_xy ||
                          bc == Keating_struct::Periodic_xz ||
                          bc == Keating_struct::Periodic_xyz);

    MinWrtY = MinWrtY && (bc == Keating_struct::Periodic_y ||
                          bc == Keating_struct::Periodic_xy ||
                          bc == Keating_struct::Periodic_yz ||
                          bc == Keating_struct::Periodic_xyz);

    MinWrtZ = MinWrtZ && (bc == Keating_struct::Periodic_z ||
                          bc == Keating_struct::Periodic_yz ||
                          bc == Keating_struct::Periodic_xz ||
                          bc == Keating_struct::Periodic_xyz);
  }

  SCALE_FCTR = keating.SCALE_FCTR;
}
