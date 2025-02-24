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
$Header: /repo/nemo3d/src/base/mat3d.c,v 1.15 2005/11/12 19:14:41 gekco Exp $ 
*****************************************************************************/

#include "mat3d.h"



int trivial_mat_init( qd_struct d )
{
    d->mstar = Rvectr(10);

    d->mstar[0] = 1; /* Free mass */
    d->mstar[1] = .067; /* GaAs */

    return(0);
}


#ifdef NEED_unit_cell
/* Does not appear to be called anywhere!!!  */
int unit_cell(Dev_struct::CrystalStructType structure, int site, int anion, int cation, int cation1, real ratio1, int cation2, real ratio2)
{
    int element;

    if ( ratio1 != 0 || ratio2 != 0 )
    {
        /* Create a new cation out of the mixture using VCA */
    }
    
    switch(structure)
    {
      case(Zincblende):
        switch(site)
	{
	  case(0):
	    element = cation;
	    break;
	  case(1):
	    element = cation;
	    break;
	  case(2):
	    element = cation;
	    break;
	  case(3):
	    element = cation;
	    break;
	  case(4):
	    element = anion;
	    break;
	  case(5):
	    element = anion;
	    break;
	  case(6):
	    element = anion;
	    break;
	  case(7):
	    element = anion;
	}
	break;
      case(Cubic):
	switch(site)
	{
	  case(0):
	    element = anion;
	    break;
	default:
	    die("Should not reach this code in unit_cell\n");
	    break;
	}
    default:
	die("Should not reach this code in unit_cell 2\n");
	    break;
    }   
    return(0);
}
#endif /* NEED_Unit_cell */


void strnmat_init(qd_struct d, const vector<double>& mb_strain, 
                  int cation, int anion)
{
   // strnmat should not be needed if Hamiltonian is position-independent
   if (!d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian)  return;
   
   rvectr tmp_rvectr = Rvectr(mb_strain.size());
   for (unsigned int i=0; i<mb_strain.size(); i++) {
      tmp_rvectr[i] = mb_strain[i];
   }

   add_elem_rmatrixrvectr(tmp_rvectr, cation, anion, d->strnmat);

   cout_master << "   Initializing strain parameters for ("
               << atomID_to_str(cation) << atomID_to_str(anion) << ")\n";

   if (!d->strnmat[anion][cation]) {
      cout_master << "   Copying strain parameters for (" 
                  << atomID_to_str(anion) << atomID_to_str(cation)
                  << ") from (" 
                  << atomID_to_str(cation) << atomID_to_str(anion)
                  << ")\n"; 
      rvectr tmp_rvectr2 = Rvectr(mb_strain.size());
      for (unsigned int i=0; i<mb_strain.size(); i++) {
         tmp_rvectr2[i] = mb_strain[i];
      }
      add_elem_rmatrixrvectr(tmp_rvectr2, anion, cation, d->strnmat);
   }
}


void param_init( qd_struct d )
{
   // masterPrint("\nsp3d5s_nparam=%d\n",sp3d5s_nparam);
   d->strnmat = Rmatrixrvectr(EL_TOTAL-1,EL_TOTAL);
   d->param = Rvectr(sp3d5s_nparam);
   d->parmat = R3tensor(0,EL_TOTAL-1,0,EL_TOTAL-1,0,sp3d5s_nparam-1);

  /* Loop through the material list and initialize the paramter matrix 
     for the known and used materials */
   for (vector<MaterialHandle>::iterator sm = d->sMatList.begin();
        sm != d->sMatList.end(); sm++) {

      Material_struct* smh = sm->hndl;
      if (!smh) continue;

      int anion  = str_to_atomID(sm->anion.c_str());
      int cation = str_to_atomID(sm->cation.c_str());

      real StrainShift_cation=0.0, StrainShift_anion=0.0;

#ifdef STRAIN_SHIFT_ENABLE
      if (smh->mb_strain.size()!=0 && 
          d->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian) {

         if (sp3d5s_Ebound==sp3d5s_Ebound_undefine) {
            cout_master << "   WARNING:  using cation energy shift as default\n";
            sp3d5s_Ebound=smh->mb_strain[eta_cat_shift];
         }
         
         StrainShift_cation = smh->mb_strain[eta_cat_shift];
         StrainShift_anion = smh->mb_strain[eta_an_shift];
      }
#endif
       
      if (d->BandModel==BM_20_sp3d5ss_spin      ||
	  d->BandModel==BM_10_sp3d5ss_nospin ){
         strnmat_init(d, smh->mb_strain, cation, anion);

	 for (int i=0; i<sp3d5s_nparam; i++) {
	    d->parmat[cation][anion][i] = smh->mb_pv[i];
         }
         
         d->parmat[cation][anion][pE_sa]  += smh->Ev_offset-StrainShift_anion;
         d->parmat[cation][anion][pE_pa]  += smh->Ev_offset-StrainShift_anion;
         d->parmat[cation][anion][pE_sta] += smh->Ev_offset-StrainShift_anion;
         d->parmat[cation][anion][pE_da]  += smh->Ev_offset-StrainShift_anion;

         d->parmat[cation][anion][pE_sc]  += smh->Ev_offset-StrainShift_cation;
         d->parmat[cation][anion][pE_pc]  += smh->Ev_offset-StrainShift_cation;
         d->parmat[cation][anion][pE_stc] += smh->Ev_offset-StrainShift_cation;
         d->parmat[cation][anion][pE_dc]  += smh->Ev_offset-StrainShift_cation;
      }
      else if (d->BandModel==BM_10_sp3ss_spin) {
         strnmat_init(d, smh->mb_strain, cation, anion);

	 for (int i=0; i<15; i++)
	    d->parmat[cation][anion][i] = smh->mb_pv[i];

	 param_sp3s_2_sp3d5s(d->parmat[cation][anion],
                             d->parmat[cation][anion],
                             smh->unstrnd_cubic_cell_length,
                             1.0,1.0,1.0);

         d->parmat[cation][anion][pE_sa]  += smh->Ev_offset-StrainShift_anion;
         d->parmat[cation][anion][pE_pa]  += smh->Ev_offset-StrainShift_anion;
         d->parmat[cation][anion][pE_sta] += smh->Ev_offset-StrainShift_anion;

         d->parmat[cation][anion][pE_sc]  += smh->Ev_offset-StrainShift_cation;
         d->parmat[cation][anion][pE_pc]  += smh->Ev_offset-StrainShift_cation;
         d->parmat[cation][anion][pE_stc] += smh->Ev_offset-StrainShift_cation;
      }
      else if (d->BandModel==BM_1_s_nospin) {
	 /* In the single band case we do not perform an electronic 
            strain calculation  and the filed strnmat is not filled */

	 for (int i=0; i<5; i++)
	    d->parmat[cation][anion][i] = smh->mb_pv[i];

	/* The 1-D NEMO parameter set for Bands_1_s_nospin consists of 
	   5 paramters:
	   param[0] = Eg
	   param[1] = mass_left
	   param[2] = mass_center
	   param[3] = mass_right
	   param[4] = lattice constant = a

	   Modify the parameter list that was obtained from the
	   1-D code to work somewhat with the 3-D code.
	   For Zincblende structures param[2] will be mapped to the 
           diagonal element of the Hamiltonian.
	   For Cubic structures param[0] will be mapped to the 
           diagonal element of the Hamiltonian.
	   param[10] will be mapped to the off-diagonal element of the
	   Hamiltonian.
	   For the diagonal element we would like to have 
	   d = HBAR_MASSFACTOR / (a*a) * 2 / m* + Ec
	   For the off diagonal element we would like to have:
	   s = HBAR_MASSFACTOR / (a*a) / m*
	*/
         switch(d->opt.Dev.CrystalStruct) {
         case Dev_struct::Cubic:	
	    /*triple the bandwidth value.
              note that there is effectively three folded bands...
              see notes in notebook #4 pg32. */
	    d->parmat[anion][anion][0] = 
               d->parmat[cation][cation][0] = 
               d->parmat[cation][anion][0] = 
               d->parmat[cation][anion][2] = 
               3.0 * HBAR_MASSFACTOR / (smh->mb_pv[4] * smh->mb_pv[4]) * 
               2.0 / smh->mb_pv[2] + smh->mb_pv[0] + smh->Ev_offset;
	    d->parmat[cation][anion][10] = 
               d->parmat[anion][anion][10] = 
               d->parmat[cation][cation][10] = 
               - HBAR_MASSFACTOR / (smh->mb_pv[4] * smh->mb_pv[4]) / smh->mb_pv[2] ;
	    break;
         case Dev_struct::Zincblende:
	    /* double the effective mass value 
               double the bandwidth 
               => modify hopping but not the diagonal ....*/
	    d->parmat[anion][anion][0] = 
               d->parmat[cation][cation][0] = 
               d->parmat[cation][anion][0] = 
               d->parmat[cation][anion][2] = 
               2.0* HBAR_MASSFACTOR / (smh->mb_pv[4] * smh->mb_pv[4]) * 
               2.0 / (2.0 * smh->mb_pv[2]) + smh->mb_pv[0] + smh->Ev_offset;
	    d->parmat[cation][anion][10] = 
               d->parmat[anion][anion][10] = 
               d->parmat[cation][cation][10] = 
               - HBAR_MASSFACTOR / (smh->mb_pv[4] * smh->mb_pv[4]) / 
               (2.0 * smh->mb_pv[2]) ;
	    break;
	    
         default:
	    d->parmat[anion][anion][0] = 
               d->parmat[cation][cation][0] = 
               d->parmat[cation][anion][0] = 
               d->parmat[cation][anion][2] = 
               HBAR_MASSFACTOR / (smh->mb_pv[4] * smh->mb_pv[4]) * 
               2.0 / smh->mb_pv[2] + smh->mb_pv[0] + smh->Ev_offset;
	    d->parmat[cation][anion][10] = 
               d->parmat[anion][anion][10] = 
               d->parmat[cation][cation][10] = 
               - HBAR_MASSFACTOR / (smh->mb_pv[4] * smh->mb_pv[4]) / smh->mb_pv[2] ;
	    break;
         }
      }
      else{
	die("Unimplemented bandstructure model in mat_param_init.\n");
      }

   }
   
   /* for the case of covalent bonding (eg. SiGe) we need to allow each
      atom type to be either cation or anion.  Therefore, we need to define
      material parameters that describe both A(cation)B(anion) and
      B(cation)A(anion).  For the case of polar bonding (e.g. GaAs), each
      atom type can only be either cation or anion so that separated GaAs
      and AsGa are not needed. Since the indexing order of the TB parameters
      depends on whether an atom is a cation or anion, we need to impose
      an additional constraint -- d->parmat *must* be accessed as 
      d->parmat[cation][anion].  For now, we perform the follow kludge:
      we check whether the lattice constant for AB is the same as for BA.  
      If AB (BA) is non-zero while BA (AB) is zero, then
      we copy the material parameters for AB (BA) to BA (AB).  If they
      are both non-zero we leave the material parameters untouched.

      Note that a more sensible way to proceed would have been to order
      parmat such that the TB matrix elements depended on 
      (AtomType_this,AtomType_nbr) instead of on (cation,anion).
   */
   for (int i=0; i<EL_TOTAL; i++) {
      for (int j=i+1; j<EL_TOTAL; j++) {
         if (d->parmat[i][j][pV_unstr_latt]==d->parmat[j][i][pV_unstr_latt])
            continue;
         // lattice constants differ; 
         // one set of material parameters has probably not yet been defined
         if ( d->parmat[i][j][pV_unstr_latt]==0.0) {
            cp_rarray2rarray(d->parmat[j][i],d->parmat[i][j],sp3d5s_nparam);
         } 
         else if ( d->parmat[j][i][pV_unstr_latt]==0.0) {
            cp_rarray2rarray(d->parmat[i][j],d->parmat[j][i],sp3d5s_nparam);
         } 
         else
            die("Lattice constants of AB and BA differ!");
      }
   }
}
