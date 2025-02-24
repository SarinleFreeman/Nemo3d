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
$Header: /repo/nemo3d/src/base/qd_struct.h,v 1.20 2007/03/08 06:48:18 lee509 Exp $ 
*****************************************************************************/


#ifndef QD_STRUCT_STRUCT_H
#define QD_STRUCT_STRUCT_H 1

#include <i2tensor.h>
#include <i3tensor.h>
#include <r2tensor.h>
#include <r3tensor.h>
#include <rvector.h>
#include <cvector.h>
#include <rmatrix.h>
#include <cmatrix.h>
#include <imatrix.h>
#include <r2tensor.h>
#include "realtype.h"
#include "ivectori3matrix.h"
   /* #include "mat_struct.h" */

#include "ivectorvector.h"
#include "cmatrixmatrix.h"
#include "rmatrixrvector.h"

#include "top_struct.H"
#include "MaterialHandle.h"

#include "geo_struct.h"
#include "waveFunction.h"

#include "HamZincBlende.h"
#include "HamZincBlendeComplex.h"

#include "phon_struct.h"

#include "StrainParam.h"
#include "mb_ham_spds_micro2.h"

#define MEMORYCHECK
#ifdef MPI3d
#define MEM_ID 1
#else
#define MEM_ID 0
#endif

#define MAX_MAT_USED 20

#define DELIMITERS " \n\t"

#define E_LOOSE_BOND_SHIFT  20.0

#undef HAMILTONIAN_SCALE
#define HAMILTONIAN_SCALE_VALUE 1.0



/*  additional properties of semi_materials:
    matid3d
    anion_id
    cation_id
    anion
    cation
*/


/*! \class QD_struct
 *  \brief Represents the entire structure to be simulated.

 *  This class encapsulates data relevant to the entire structure to be simulated by NEMO3D. It also has resultant output data that come from simulation. */

typedef class QD_struct {
 public:
   char      *inputfile;
   /*!< name of input deck */
   top_struct opt;
   /*!< this has several information extracted from input deck.
    * @see top_struct::top_struct() */
   
   vector<MaterialHandle> sMatList;
   /*!< 1D array of data type <materialHandle>. This is the list of materials avaiable in a specific simulation domain. */ 
   int Abort_result;
   /*!< Flag to chec whether simulation was forced to quit or not. Set to 1 if process is forced to quit. */
   geo_struct geo;
   /*!< geometry infomation. 
    * @see geo_struct::geo_struct() */

   imatrix Atom_to_lm;
   /*!< 2D integer matrix. Atom_to_lm[k][0] returns index of a cell this atom belongs. 
    * Atom_to_lm[k][1] returns relative index of this atom in a cell of index l. 
    * (Note : apart from its absolute atomic index, every atom has another relative index 
    * in a unit cell that starts from 1 to geo.AtomperCellMax())
    * @see geo_struct::AtomperCellMax() 
    */
   
   ivectr surfaceAtoms;
   /*!< list of atoms on surface. surfaceAtoms[k] returns atomic index of k-th surface atom. */
   vector<waveFunction> WF;
   /*!< This matrix has raw eigenvalues. This data will be additionally scaled by 
    * HAMILTONIAN_SCALE_VALUE, which produces actual eigenvalues that we see from the result. */
   int neigv;
   /*!< Number of eigenvalues */
   real shift, eta;
   rvectr mstar, te;
   imatrix tindx;
   rmatrix wf, opmat;
   r3tensor opmat_pxyz;
   cmatrix wfc;
   /*!< complex matrix that contains calculated wavefunction values. 
    * wfc[ev][geo.l_m_2_offset_ham_tot[l][m] + n] indicates wave function of 
    * (nth basis of the atom whose cell index = l and atomic index = m) for ev-th eigen value
    * @see geo_struct::l_m_2_offset_ham_tot
    */
   
   cvectr eigv;
   /*!< complex vecgtor that contains calculated eigenvectors.*/
   
   i3tensor mat;
 
   /* Bandstructure stuff.  Will seperate this into a substructure. */
   int NBasisStates;
   /*!< number of basis per atom. It will be set according to the selected band model. */
   
   int  n_atom_tot, n_ham_tot, npair;
   BM_Type BandModel;
  
   real kxL, kyL, kzL;
   /*!< increment of k value into x, y and z direction. 
    * initial, final k values and number of steps between init/final values are extraced 
    * from input deck. 
    */
   
   ivectr Basis_Spin, BasisAnion, BasisCation;
   rvectr param;
   cmatrix hd, ho;
    
   r3tensor parmat;
   
   // temporary kludge to add electrostatic potential
   r2tensor phi;
   

   HamZincBlende Hzb;
   HamZincBlendeComplex HzbComplex;

   #if 1 
   // sunnylee
   // Cubic 1s structure (sparse/full store)
   rvectr CubicSparse1s;
   bool is_CubicSparse1sSet;
   int CubicMaxNbr; // value will be maxnbr+1 (includes self energy)
   #endif

   /* granularity of physical domain over which electronic structure 
      is to be computed */
   int _primitiveCell_granularity;
   /*!< Set to 1 if electronic structure has granularity on a primitive cell */
   
   int _bandstruct_on_full_domain;
   /*!< Set to 1 if electronic calculation is being done in full domain. */ 
   
   int _bandstruct_on_local_domain; 
   /*!< Set to 1 if electronic calculation is being done in local domain. */
   
   /* Hamiltonian storage stuff */
   cmatrixmatrix Ham;
   ivectr Ham_neighbor;
   imatrix Ham_offset_x, Ham_offset_y, Ham_yc;
   int Ham_use_hermiticity;
   imatrix Ham_transp;
   imatrix Ham_atom_nbr_2_m;
   imatrix Ham_atom_nbr_2_l;
   imatrix l_m_2_atom;
   /*!< l_m_2_atom[l][m] returns absolute atomic index of an atom whose cell index = l 
    * and relative atomic index in a unit cell is m.*/
    
   /* Phonon stuff */
   rvectr Amass;//Atomic masses
   phon_struct phon;
   
   StrainParam strain;

   rmatrixrvectr strnmat;

   /* lanczos stuff */
   cvectr rc, r0, q0, q_m10, r0_start, cvec_aux;
   
   /* matrix multiplier stuff */
   cmatrix mm;

   /* Parallel stuff */
   int maxseg;
   /*!< largest size of vector segment. */
    
   int nproc;
   /*!< # of processes */
   
   ivectr cell_s;
   /*!< starting index of a unitcell segment in each process. cell_s[n] has the starting index of a unitcell segment 
    * allocated in n-th process. */
   
   ivectr cell_ln;
   /*!< size of a unitcell segment (# of unitcell, equivalently) in each process. cell_ln[n] has # of unitcells 
    * allocated in n-th process. */
   ivectr seg_s;
   /*!< staring index of a vector segment in each process. seg_s[n] has the starting index of a vector segment 
    * allocated in n-th process. */
   ivectr seg_ln;
   /*!< size (length) of a vector segment in each process. seg_ln[n] has the length of a vector segment
    * allocated in n-th process. */
   ivectr nvmap;
   imatrix commsize, commorder;
#define	NEW_VECMAP
#ifdef NEW_VECMAP
   ivectrvectr vecmap;
#else
   imatrix vecmap;
#endif
   /* Rayleigh Quotient stuff */
   cvectr drq;

   QD_struct();
   void setTopStruct(const char* f_xml, const char* f_in);
   void checkInput();
   void set_sMatList();
   void cleanup();
} *qd_struct;



#endif


