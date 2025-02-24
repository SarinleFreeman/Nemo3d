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
             Olga L. Lazarenkova
             Hook Hua

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/base/StrainParam.h,v 1.8 2005/01/14 20:24:38 marek Exp $
*****************************************************************************/

#ifndef _STRAIN_PARAM
#define _STRAIN_PARAM

#include "realtype.h"
#include <ivector.h>
#include <rvector.h>
#include <rmatrix.h>
#include <imatrix.h>

#include "top_struct.H"
#include "MaterialHandle.h"

class StrainComm {
 public:
   int N;       // size of raw data in units of sizeof(int)
   void* data;  // points to raw data 
   int* indx;    // points to index data 
   real* rdata;  // points to real data

   StrainComm() {;} // dummy constructor
};


class StrainParam {
 public:
   // shoule we minimize the lattice constant wrt to x, y, z?
   bool MinWrtX, MinWrtY, MinWrtZ;

   int debugp;
   int first_or_second; 

   double SCALE_FCTR;

   rvectr dp;
   rmatrix amat, alpha, beta;
   rmatrix aa_VFF,ba_VFF,ca_VFF;//Coefficients of the strain dependence of the VFF constants in the Keating model. Used for the strain calculations and the low-energy acousticl phonons
   rmatrix ao_VFF,bo_VFF,co_VFF;//Coefficients of the strain dependence of the VFF constants in the Keating model. Used for optical phonons.

   /* mapping stuff for strain */
   ivectr lm_sndL, lm_sndR, lm_rcvL, lm_rcvR;
   int cellLMin, cellLMax, atomLMin, atomLMax, cellRMin, cellRMax, atomRMin, atomRMax ;
   imatrix rcvL_lm, rcvR_lm ;
   StrainComm sndL, sndR, rcvL, rcvR;
   

   StrainParam() {;} // dummy constructor
   void SetMatParm(vector<MaterialHandle> sMatList);
   void SetOptions(const top_struct& opt);

   void Deallocate();
   ~StrainParam() { Deallocate(); }

   bool doPeriodMinimization() const { return MinWrtX || MinWrtY || MinWrtZ; }
   //Get the effective values of alpha and beta if they were schaled with strain
   /*   real parabola(real Anh,//anharmonicity parameter
                 real Dnh,//second-order anharmonicity parameter
                 real arg//argument corresponding to the strain
                 );
   */
   real alpha_eff(int m_atom_type,//type of the m-th atom
                  int n_atom_type,//type of the n-th atom
                  real r2,//square of the real bond length
                  real d2,//square of the equilibrium bond length
                  int anh,//1-use acoustical anharmonicity corrections
                  //2-use optical anharmonicity corrections
                  //3-use linear interpolation
                  bool pr//true for printing truncated alpha
                  );
   real beta_eff(int m_atom_type,//type of the m-th atom
                 int n1_atom_type,//type of the n1-th atom
                 int n2_atom_type,//type of the n2-th atom
                 rvectr d1,rvectr d2,//vectors of the relaxed bonds of the V doublet
                 rvectr r1,rvectr r2,//vectors of the strained bonds of the V doublet
                 int anh,//1-use acoustical anharmonicity corrections
                 //2-use optical anharmonicity corrections
                 //3-use linear interpolation
                 bool pr//true for printing truncated beta
                 );

};

#endif
