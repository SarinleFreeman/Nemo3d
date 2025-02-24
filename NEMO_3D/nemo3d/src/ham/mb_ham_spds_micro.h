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
$Header: /repo/nemo3d/src/ham/mb_ham_spds_micro.h,v 1.16 2007/06/12 20:15:22 nkharche Exp $
*****************************************************************************/

#ifndef MB_HAM_SPDS_MICRO_H
#define MB_HAM_SPDS_MICRO_H 1


#define SIDEBANDS 1

#undef SIMPLESTAMP
#undef  SPLIT_HAMILTONIAN_AC
#define SIMPLE_HAMILTONIAN_AC

#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include "realtype.h"
#include <cvector.h>
#include "cvector_lib.h"
#include <rvector.h>
#include <fvector.h>
#include <ivector.h>
#include <complex.h>
#include <matrix.h> /* used for cdim() and vdim() definitions */
#include <cmatrix.h>
#include "cmatrix_lib.h"
#include <rmatrix.h>
#include "constants_nemo.h"
#include "qd_struct.h"
//#include <stream.h>  // this does not seem to be used
#include "mat_def.h"

void ceigval_ns_full( cvectr eigv, cmatrix h );
void ceigvec_ns_full( cvectr eigv , cmatrix z, cmatrix h ); 
rvectr eigval_h_full (rvectr d, coperator a, int *error_code);


#define STRAIN_SHIFT_ENABLE			       
#define STRAIN_SHIFT_ENABLE_1D			       

#define sp3s_nparam 15
#define sp3d5s_nparam  35
#define sp3d5s_nparam_copy_size  10
#define Additional_Diagonal_constant
#ifdef Additional_Diagonal_constant
#define sp3d5s_nstrain_extra 16
#define sp3d5s_nstrain 33
#else
#define sp3d5s_nstrain_extra 0
#define sp3d5s_nstrain 33
#endif

#define sp3d5s_nstrain_old 15
#define sp3d5s_Ebound_default 27.0
#define sp3d5s_Ebound_undefine -100.0
#define sp3d5s_strain_ideal 2.0
#define sp3d5s_c_strain_ideal 1.0

extern double  sp3d5s_Ebound;


#define STa      0
#define Sa       1
#define Pxa      2
#define Pya      3
#define Pza      4
#define Dxya     5
#define Dyza     6
#define Dzxa     7
#define Dx2my2a  8
#define Dz2a     9
#define STc     10
#define Sc      11
#define Pxc     12
#define Pyc     13
#define Pzc     14
#define Dxyc    15
#define Dyzc    16
#define Dzxc    17
#define Dx2my2c 18
#define Dz2c    19

/* normal order in the parameter list (35 members) for sp3d5s* band model
    0	pAs_s_s 
    1	pAs_p_p 
    2	pIn_s_s 	
    3	pIn_p_p 
    4	pAs_ss_ss 
    5	pIn_ss_ss 
    6	pAs_d_d 
    7	pIn_d_d 
    8   pAs_so_so 
    9	pIn_so_so 
   10   pInAs_s_s
   11	pInAs_ss_ss
   12	pInAs_ssa_sc
   13	pInAs_sa_ssc
   14	pInAs_sa_pc 
   15	pInAs_sc_pa 
   16	pInAs_ssa_pc 
   17	pInAs_ssc_pa 
   18	pInAs_sa_dc 
   19	pInAs_sc_da 
   20	pInAs_ssa_dc 
   21	pInAs_ssc_da 
   22	pInAs_p_p_sigm 
   23	pInAs_p_p_pi 
   24	pInAs_pa_dc_sigm 
   25	pInAs_pc_da_sigm 
   26	pInAs_pa_dc_pi 
   27	pInAs_pc_da_pi 
   28	pInAs_d_d_sigm 
   29	pInAs_d_d_pi 
   30	pInAs_d_d_del 
   31	0.60583         lattice constant 
   32	1.0             bond distortion in x
   33	1.0             bond distortion in y
   34	1.0             bond distortion in z
*/

/* 
   These macos provide the symbolic
   mapping from an orbital description / name
   into the tight binding parameter vector.
*/

/* normal order in the parameter list (15 members) for sp3s* band model
   Notations are defined in J. Phys. Chem. Solids Vol 44, pp 365 (1983)
   The parameters will be converted to the notation used in sp3d5s* band model

    0   E(s,a)  = pAs_s_s   
    1   E(p,a)  = pAs_p_p 
    2   E(s,c)  = pIn_s_s         
    3   E(p,c)  = pIn_p_p 
    4   E(s*,a) = pAs_ss_ss 
    5   E(s*,c) = pIn_ss_ss
    6   V(s, s) = 4 * pInAs_s_s 
    7   V(x,x)  = (4/3) * (pInAs_p_p_sigm + 2* pInAs_p_p_pi)
    8   V(x,y)  = (4/3) * (pInAs_p_p_sigma - pInAs_p_p_pi)
    9   E(sa,pc)  = (4/sqrt(3)) * pInAs_sa_pc
   10   E(sc,pa)  = (4/sqrt(3)) * pInAs_sc_pa 
   11   E(s*a,pc) = (4/sqrt(3)) * pInAs_ssa_pc  
   12   E(s*c,pa) = -(4/sqrt(3)) * pInAs_ssc_pa 
   13   Spin-Orbit coupling for anion  = 3 * pAs_so_so
   14   Spin-orbit coupling for cation = 3* pIn_so_so
*/

#define    pE_sa       0
#define    pE_pa       1
#define    pE_sc       2
#define    pE_pc       3
#define    pE_sta      4
#define    pE_stc      5
#define    pE_da       6
#define    pE_dc       7
#define    pE_La       8
#define    pE_Lc       9

#define    pV_sss     10
#define    pV_ststs   11
#define    pV_stass   12
#define    pV_sasts   13
#define    pV_sapcs   14
#define    pV_scpas   15
#define    pV_stapcs  16
#define    pV_stcpas  17
#define    pV_sadcs   18
#define    pV_scdas   19
#define    pV_stadcs  20
#define    pV_stcdas  21
#define    pV_pps     22
#define    pV_ppp     23
#define    pV_padcs   24
#define    pV_pcdas   25
#define    pV_padcp   26
#define    pV_pcdap   27
#define    pV_dds     28
#define    pV_ddp     29
#define    pV_ddd     30

#define    pV_unstr_latt 31
#define    pV_frac_x     32
#define    pV_frac_y     33
#define    pV_frac_z     34

/* 
   Macros for the strain induced off-diagonal parameter scaling
   mapping from an orbital interaction description to a position
   in the parameter vector.
*/
#define    eta_ssts     0
#define    eta_ststs    1
#define    eta_sss      2
#define    eta_sps      3
#define    eta_pps      4
#define    eta_ppp      5
#define    eta_sds      6
#define    eta_stps     7
#define    eta_pds      8
#define    eta_pdp      9
#define    eta_const   10
#define    eta_dds     11
#define    eta_ddp     12
#define    eta_ddd     13
#define    eta_stds    14
#define    eta_cat_shift 15
#define    eta_an_shift  16

#define    c_ind_ss     17  
#define    c_ind_stst	18
#define    c_ind_stas	19
#define    c_ind_sast   20
#define    c_ind_sapc	21
#define    c_ind_pasc   22
#define    c_ind_scpa	c_ind_pasc
#define    c_ind_stapc	23
#define    c_ind_pastc	24
#define    c_ind_stcpa	c_ind_pastc
#define    c_ind_sadc	25
#define    c_ind_dasc	26
#define    c_ind_scda	c_ind_dasc	
#define    c_ind_stadc	27
#define    c_ind_dastc	28
#define    c_ind_stcda	c_ind_dastc
#define    c_ind_pp	29
#define    c_ind_padc	30
#define    c_ind_dapc	31
#define    c_ind_pcda	c_ind_dapc
#define    c_ind_dd	32


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef FORTRAN_UNDERSCORE
#define dgesv   dgesv_
#endif

#ifdef FORTRAN_UNDERSCORE2
#define dgesv   dgesv__
#endif

#ifdef FORTRAN_ALLCAPS
#define dgesv   DGESV
#endif
   void dgesv(int*, int*, double*, int*, int*, double*, int*, int*);
#ifdef __cplusplus
}
#endif /* __cplusplus */


void diag_map_set(void);

int h3d_offdiag( cmatrix h, rvectr nnv, ivectr basis0, ivectr basisn, ivectr spin,
rvectr param, rvectr mb_strain, real unstrnd_lattice_const, int nb2 );

int h3d_offdiag_diag_corr( cmatrix ho, cmatrix hd, rvectr nnv, rvectr nnv0, ivectr
basis0, ivectr basisn, ivectr spin, rvectr param, rvectr mb_strain, real
unstrnd_lattice_const, int nb2 );

int h3d_diag( cmatrix h, ivectr basis0, ivectr spin, rvectr param, int nb2 );

int h3d_diag_strain_corr( cmatrix ho, cmatrix hd, rvectr nnv, ivectr basis0, ivectr
basisn, ivectr spin, rvectr param, rvectr mb_strain, real unstrnd_lattice_const , int nb2 );

int h3d_diag_pz_addition_for_wurtzite_crystal_splitting ( qd_struct d, double pz_addition);

int h3d_diag_surface_atom_shift ( qd_struct d );

void h3d_diag_dangling_bond_shift ( qd_struct d, int cindx, int aindx );
void h3d_diag_dangling_bond_shift_zincblende ( qd_struct d, int cindx, int aindx );
void h3d_diag_dangling_bond_shift_wurtzite ( qd_struct d, int cindx, int aindx );
void h3d_diag_dangling_bond_shift_ZB_rot ( qd_struct d, int cindx, int aindx ); /* Added for ZB_110 and ZB_111 */

double vector_dot_vector( rvectr, rvectr, int);

complex spin_orbit_spds( int basis1, int basis2, int spin1, int spin2, real delta_a, real delta_c );

/* Returns the diagonal element of spds hamiltonian */
real diag( rvectr param, int basis );

/* Returns off diagonal element between state ba and bc.  Vector
between anion and cation is input as [xd,yd,zd] */

/* Gets off diagonal elements and multiplys by appropriate phase
   factor */
complex ei( rvectr v, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k, rvectr d, rvectr d0,  int ba, int bc, real *E_diag_corr );

complex eidk( rvectr v, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k, rvectr d, rvectr d0,  int ba, int bc, int dim1, real *E_diag_corr );

/* Nemo interface for spds hamiltonian */
int tbh_const_1d_nemo( cmatrix h, cmatrix hp1, rmatrix nnv, rmatrix nnv0, ivectr basis, ivectr spin, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k );

int tbh_const_bulk( cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis, ivectr spin, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k );

int tbh_const_bulk_dk( cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis, ivectr spin, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k, int dim1 );

real mtensor( cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis, ivectr spin, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k, int band, rmatrix mtens);

real Newton_micro( cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis, ivectr spin, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k, int band, rvectr mdE, rmatrix mtens);

/* This one constructs the bulk Hamiltonian using the 3d hamiltonian calls */
int tbh_const_bulk_3d( cmatrix h, rmatrix nnv, ivectr basis, ivectr spin, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k );

real dE_dk( cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis, ivectr spin, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k, int band, rvectr dE );

void E_const_prep( rvectr par_in, rvectr mb_strain, rvectr par_out, real unstrnd_cubic_cell_length, real xd, real yd, real zd, 
		   real *l, real *m, real *n, real *l2, real *m2, real *n2 );

#ifdef Additional_Diagonal_constant
real E_const_get( rvectr vpar_e_const, rvectr vpar_orig, rvectr strain_vec, int basis1_econst, int basis2_econst, 
		  real l, real m, real n, real l2, real m2, real n2, 
		  real l0, real m0, real n0, real l02, real m02, real n02,
		  real cstrain, real *out_diag_corr);
#else
real E_const_get( rvectr vpar_e_const, rvectr vpar_orig, int basis1_econst, int basis2_econst, 
		  real l, real m, real n, real l2, real m2, real n2, 
		  real l0, real m0, real n0, real l02, real m02, real n02,
		  real cstrain, real *out_diag_corr);
#endif

void strain_corr_prep( rvectr par_in, rvectr mb_strain, rvectr par_scal_out, real unstrnd_cubic_cell_length, 
		       real xd, real yd, real zd, 
		       real xd0, real yd0, real zd0, 
		       real *l, real *m, real *n, real *l2, real *m2, real *n2 ,
		       real *l0, real *m0, real *n0, real *l02, real *m02, real *n02 
		       );


int tbh_const_1d_nemo_old( cmatrix h, rvectr mb_strain, real unstrnd_cubic_cell_length, cmatrix hp1, rmatrix nnv, rmatrix nnv0, rvectr param, rvectr k );


int tbh_const_bulk_dk_old( cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis, ivectr spin, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k, int dim1 );
real dE_dk_simple( cmatrix h, rmatrix nnv, rmatrix nnv0, ivectr basis, ivectr spin, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, rvectr k, int band, rvectr dE );
rvectr param_sp3s_2_sp3d5s( rvectr param_in, rvectr param_out, 
                            real unstrnd_lattice_const, 
                            real delX, real delY, real delZ );

complex vector_potential_phase(qd_struct d, int cindex, int aindex, rvectr nnv);
int     apply_magnetic_field(qd_struct d, int cindex, int aindex); 
int     apply_magnetic_field(qd_struct d, int cindex, int aindex, rvectr nnv); 

int compute_strain_tensor(real* exx, real* eyy, real* ezz, real* exy,
                          real* exz, real* eyz,
                          qd_struct d, int Zc, int Za);

int on_site_d_orbital_shift(qd_struct d, int Zc, int Za, real E_d, real b_d,
                            real* Exy, real* Exz, real* Eyz);
#endif
