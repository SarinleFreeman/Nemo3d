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
$Header: /repo/nemo3d/src/base/phon_Keating.c,v 1.17 2007/08/10 15:03:33 baeh Exp $
*****************************************************************************/

/*
 * Calculation of phonon dispersion within Keating model
 * =====================================================
 * The dynamical matrix DM_{iD,jD} to be solved for eigenvalues
 * has maximum states_in_row=51 nonzero elements
 * with xyz block for every bond connected to the m-th atom in l-th unit cell
 * i=3*(l*d->N_MaxAtomCell+m)+xi, where xi counts the components [012] of the displacement vector.
 */

#include "phon_Keating.h"

#define A1DM
#define A2DM
#define BXXDM
#define BX1X2DM

#define DISREGARD_SMALL_ELEMENTS
//#ifdef PR_ZEROS_SUBSTITUTION

#define PR_DISP
//#define PR_DM_MATLAB
//#define PR_DMMAP
//#define PR_nbr_on_proc
//#define PR_PHON_PARAMETERS
//#define PIC_DM

#define TALK_A_LOT 100


/**
 * Returns index of the element in the Dynamical Matrix.
 * @param d the qd_struct that ...
 *
 *
 *
 * @return the index of the element in the Dynamical Matrix
 * @see d->geo.AtomsPerCellMax()
 */
inline int  N_DM(qd_struct d,
		 int l,//-th unit cell
		 int m,//-th atom
		 int xi//-th component
		 )
{
  return(3*(l*d->geo.AtomsPerCellMax()+m)+xi);
}




void Map_DM(qd_struct d)
     /*Mapping of the indexies of the dynamical matrix in
       d->phon.DMmap[d->geo.N_Cell][d->geo.cellgeom.AtomsPerCellMax()]*/
{
  int ind_i=0,ind_f=0;//indexies in the Dynamical matrix
  int n_shifted=0;//Count the number of shifted lines
  int n_rem=0;//Count the number of removed lines
  int n_free=0;//Current number of free indexies in DM
  int n_wait=0;//Current number of lines waiting to be shifted inside DM
  
  imatrix WaitInd;/*Temporary storage for indexies in the DM
                    [# of removed atoms][0..5]:
                    [n_shifted][0]=0/1 <=> waits/not
                    [n_shifted][1]=waiting cell
                    [n_shifted][2]=waiting atom
                    [n_shifted][3]=waiting projection
                    [n_shifted][4]=ind_i -- just to check
                    */
  ivectr FreeInd;//Indexies in DM that would correspond to the removed atoms
  int size_of_the_matrix=(3*(d->geo.N_Cell*d->geo.cellgeom.AtomsPerCellMax()-
                             d->n_atom_tot)>0 ?
                          3*(d->geo.N_Cell*d->geo.cellgeom.AtomsPerCellMax()-
                             d->n_atom_tot):
                          1);

  WaitInd=Imatrix(size_of_the_matrix,5);
  FreeInd=Ivectr(size_of_the_matrix);
  
  for(int cell=0; cell<d->geo.N_Cell; cell++){
    for(int atom=0; atom<d->geo.cellgeom.AtomsPerCellMax(); atom++){
      for(int proj=0; proj<3; proj++){
        ind_i=N_DM(d,cell,atom,proj);
        if(!d->geo.StrainHasPeriodicity()){//Necessary for closed boundary conditions only
          if(d->geo.AtomType[cell][atom]){
            //if atom has not been removed
            if(ind_i>=3*d->n_atom_tot){
              //Check if there are indexies avilable
              if(n_free>0){
                for(int j=0; j<n_rem; j++){
                  if(FreeInd[j]>0){
                    ind_f=FreeInd[j];
                    //Shift the stack of free indexies
                    for(int i=j; i<n_free; i++){
                      FreeInd[i]=FreeInd[i+1];
                    }//Shift the stack of free indexies
                    n_free--;
                    d->phon.DMmap[cell][atom][proj]=ind_f;
                    break;
                  }
                }//Search for a free index in DM
              }//There are positions available
              else{//No positions available => stand in the waiting list
                n_wait++;
                if(WaitInd[n_wait][0]){
                  printf("\nError with waiting list!!\n");
                  printf("n_wait=%d WaitInd[n_wait][0]=%d\n",n_wait,WaitInd[n_wait][0]);
                  printf("WaitInd[%d][1](cell)=%d ",n_wait,WaitInd[n_wait][1]);
                  printf("WaitInd[%d][2](atom)=%d ",n_wait,WaitInd[n_wait][2]);
                  printf("WaitInd[%d][3](proj)=%d\n",n_wait,WaitInd[n_wait][3]);
                  
                }
                WaitInd[n_wait][0]=1;
                WaitInd[n_wait][1]=cell;
                WaitInd[n_wait][2]=atom;
                WaitInd[n_wait][3]=proj;
                WaitInd[n_wait][4]=ind_i;
              }//No positions availabe
              n_shifted++;
            }//ind_i>3*d->n_atom_tot
            else{
              //Check if there are indexies avilable
              if(n_free>0){
                for(int j=0; j<n_rem; j++){
                  if(FreeInd[j]>0){
                    if(FreeInd[j]<ind_i){//Free index in DM is smaller then ind_i
                      ind_f=FreeInd[j];
                      //Shift the stack of free indexies
                      for(int i=j; i<n_free; i++){
                        FreeInd[i]=FreeInd[i+1];
                      }//Shift the stack of free indexies
                      
                      FreeInd[n_free]=ind_i;
                      d->phon.DMmap[cell][atom][proj]=ind_f;
                      break;
                    }//Free index in DM is smaller then ind_i
                    else{//Free index in DM is lager then ind_i
                      d->phon.DMmap[cell][atom][proj]=ind_i;
                    }//(else)Free index in DM is lager then ind_i
                  }//If n_free>0
                }//Search for a free index in DM
              }//There are positions available
              else{//No positions available => Assign to ind_i
                d->phon.DMmap[cell][atom][proj]=ind_i;
              }//No free positions available => Assigned to ind_i
            }//ind_i<=3*d->n_atom_tot
          }//if atom has not been removed
          else{//if atom has been removed
            d->phon.DMmap[cell][atom][proj]=-1;
            n_rem++;
            if(ind_i<3*d->n_atom_tot){
              //There is an available index in DM => Check in the waiting list
              if(n_wait>0){
                d->phon.DMmap
                  [WaitInd[n_wait][1]]
                  [WaitInd[n_wait][2]]
                  [WaitInd[n_wait][3]]=ind_i;
                WaitInd[n_wait][0]=0;
                WaitInd[n_wait][1]=0;
                WaitInd[n_wait][2]=0;
                WaitInd[n_wait][3]=0;
                WaitInd[n_wait][4]=0;
                n_wait--;
              }//n_wait>0
              else if(n_wait==0){
                FreeInd[n_rem]=ind_i;
                n_free++;
              }//n_wait=0
              else{
                printf("\nERROR: n_wait<0!! Continue anyway...\n");
              }//n_wait<0
            }//free index is within DM range
          }//if atom has been removed
        }//if boundary conditions are closed
        else{//if boundary conditions are 3D periodic
          d->phon.DMmap[cell][atom][proj]=ind_i;
        }
      }//xyz components of the displacement
    }//atoms
  }//cells
  //Printout the map:
#ifdef PR_DMMAP
  if(!mpi_n3d_id)
    Map_print(d);
#endif//PR_DMMAP
  //Deallocate memory
  rm_imatrix(&WaitInd);
  rm_ivectr(&FreeInd);
}//end of void Map_DM


void initial_assignments(qd_struct d,
                         int n_qph,//point in the reciprocal space
                         int l,//cell
                         int m,//atom
                         cvectr phase,//[n_max] phase associated with n-th neighbor
                         rvectr alpha,//[n_max] bond length VFF constant associated with n-th neighbor of m-th atom
                         rvectr beta,//[n_max] bond bending VFF constant associated with n-th neighbor of m-th atom
                         rvectr sqrtMn,//[n_max] square root of the atomic mass associated with n-th neighbor
                         rmatrix nnv,//[n_max] strained bond length between n-th neighbor and m-th atom
                         rmatrix nnv0,//[n_max] realxed bond length between n-th neighbor and m-th atom
                         rvectr dssq,//[n_max] strained bond length squared between n-th neighbor and m-th atom
                         rvectr d0sq//[n_max] realxed bond length squared between n-th neighbor and m-th atom
                         )
//Assign a set of arrays associated to the n-th neighbor of m-th atom in l-th unit cell
//n_max=d->geo.Neighbors(m)
{
  rvectr qR;//q_ph.R=q_ph_x*R_x+q_ph_y*R_y+q_ph_z*R_z to find the phase
  int cindx,aindx;//cell and atom indexes for the n-th neighbor
  rvectr D_ijk;//vector of the periods in 3 directions
  real r1;
  int AtomType_cindx_aindx;// Type of the neighbor atom
  real a_lattice[3];
  int AtomType_l_m=(int) d->geo.AtomType[l][m];//Type of the m-th atom
  
  qR=Rvectr(d->geo.Neighbors(m));
  if(d->geo.StrainHasPeriodicity()){
    //Assign the periods of the structure in all 3 directions
    D_ijk = Rvectr(3);
    D_ijk[0]=d->geo.lattice_x*((d->geo.cell_xmax-d->geo.cell_xmin)+1);
    D_ijk[1]=d->geo.lattice_y*((d->geo.cell_ymax-d->geo.cell_ymin)+1);
    D_ijk[2]=d->geo.lattice_z*((d->geo.cell_zmax-d->geo.cell_zmin)+1);
  }
  a_lattice[0]=d->geo.lattice_x;
  a_lattice[1]=d->geo.lattice_y;
  a_lattice[2]=d->geo.lattice_z;

  //loop through n neighbors of m-th atom
  for (int n = 0; n < d->geo.Neighbors(m); n++ ) {
    qR[n]=0.0;	    
    /* Get relative position of unit cell where this neighbor lives. */ 
    int i_nbr = d->geo.cell__ijk[l][0] + d->geo.NbrCell(m,n,0);
    int j_nbr = d->geo.cell__ijk[l][1] + d->geo.NbrCell(m,n,1);
    int k_nbr = d->geo.cell__ijk[l][2] + d->geo.NbrCell(m,n,2);
    
    /* Get the unit cell and atomic index for this neighbor */
    cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

#ifdef ELIMINATE_SSMAP
    if (cindx < 0 && d->geo.StrainHasPeriodicity() && (cindx + d->geo.N_Cell + 1 ) >=0){
      //This neighbor comes through the boundary with periodical conditions =>
      //1-change the index of the atom unit cell
      cindx = cindx + d->geo.N_Cell + 1 ;
#else /* ELIMINATE_SSMAP */
    if (cindx < 0 && d->geo.StrainHasPeriodicity() && d->geo.ssmap[l][m][n]!=-1){
      //This neighbor comes through the boundary with periodical conditions =>
      //1-change the index of the atom unit cell
      cindx = d->geo.ssmap[l][m][n];
#endif /* ELIMINATE_SSMAP */

      //2-check phase for all xyz projections
      for(int p=0;p<3;p++){	 
        if(d->geo.NbrCell(m,n,p)!=0 && d->phon.E_ph[n_qph][p]!=0){
          qR[n]+=d->geo.NbrCell(m,n,p)*D_ijk[p]*d->phon.E_ph[n_qph][p];
        }
      }//check phase for all xyz projections
      phase[n].r=cos(qR[n]);
      phase[n].i=sin(qR[n]);
    }//This neighbor comes through the boundary with periodical conditions
    else{//This neighbor is far from the boundary of the simulation domain
      phase[n].r=1.0;
      phase[n].i=0.0;
    }
    //Atomic index of the neighbor atom
    aindx = d->geo.NbrCell(m,n,3);
    
    /* Check if neighbor's unit cell is within simulation domain. */
    if ( cindx >= 0 ) {
      AtomType_cindx_aindx = (int) d->geo.AtomType[cindx][aindx];
      /* Check if the neighbor is within the simulation domain */
      if (AtomType_cindx_aindx) {
        //Neighbor is within simulation domain =>
        sqrtMn[n]=sqrt(d->Amass[AtomType_cindx_aindx]);//atomic mass square root
        if (sqrtMn[n]==0.0){
          die("Encountered a zero sqrtMn[n]\n");
        }

#ifdef DISPL_SINGLE_PRECISION
        nml_float* rv1 = d->geo.l_m_2_atomdsp[cindx][aindx];
        nml_float* rv2 = d->geo.l_m_2_atomdsp[l][m];
#else /* DISPL_SINGLE_PRECISION */
        real* rv1 = d->geo.l_m_2_atomdsp[cindx][aindx];
        real* rv2 = d->geo.l_m_2_atomdsp[l][m];
#endif /* DISPL_SINGLE_PRECISION */
        real* rv3 = d->geo.PositionInCell(aindx);
        real* rv4 = d->geo.PositionInCell(m);
        int* iv1 = d->geo.NbrCell(m,n);
        
        real alat = d->strain.amat[AtomType_l_m][AtomType_cindx_aindx];
        
        /* Get relative real-space location for this neighbor.  There are two
           terms inside the brackets.  The first gives the relative location
           of the neighbor's unit cell.  The second gives the relative location
           of the neighbor's atom within the unit cell.  The term disp is the
           atom's displacement due to strain */
        for (int p=0; p < 3; p++, rv1++, rv2++, rv3++, rv4++, iv1++ ) {
          real disp = *rv1 - *rv2;
          r1 = (double)(*iv1) + (*rv3 - *rv4);
          nnv[n][p] = (a_lattice[p] * r1 + disp);
          nnv0[n][p] = alat*r1;
        }
        if (nnv0[n][0]==0.0 && nnv0[n][1]==0.0 && nnv0[n][2]==0.0) {
          printf("WARNING:  nnv0 vanishes!!! on proc %d -- Try to analyse:\n", mpi_n3d_id); 
          printf("W   (l,m,n) = (%d,%d,%d)\n", l, m, n);
          printf("W   (l_nbr,m_nbr) = (%d,%d)\n", cindx, aindx);
          printf("W   alat = (%g %g)\n", alat, d->strain.amat[AtomType_l_m][AtomType_cindx_aindx]);
          printf("W   (AtomType_l_m,AtomType_cindx_aindx) = (%d,%d)\n",
                 AtomType_l_m, AtomType_cindx_aindx);
          printf("W   d->nnmap[m][n]=(%d,%d,%d)\n",
                 d->geo.NbrCell(m,n,0),
                 d->geo.NbrCell(m,n,1),
                 d->geo.NbrCell(m,n,2));
          printf("W CellAtom_2_ijk[aindx]=(%e %e %e)\n", 
                 d->geo.PositionInCell(aindx,0), 
                 d->geo.PositionInCell(aindx,1), 
                 d->geo.PositionInCell(aindx,2));
          printf("W CellAtom_2_ijk[m]=(%e %e %e)\n", 
                 d->geo.PositionInCell(m,0),
                 d->geo.PositionInCell(m,1), 
                 d->geo.PositionInCell(m,2));
        }//nnv0 vanishes
                /* Bond length distortion */
        dssq[n] =  nnv[n][0]* nnv[n][0] +  nnv[n][1]* nnv[n][1] +  nnv[n][2]* nnv[n][2];
        if (dssq[n]==0.0) printf("ERROR:  dssq[%d] vanishes!!!\n",n), fflush(stdout);
        d0sq[n] =  nnv0[n][0]* nnv0[n][0] +  nnv0[n][1]* nnv0[n][1] +  nnv0[n][2]* nnv0[n][2];
        if (d0sq[n]==0.0) printf("ERROR:  d0sq[%d] vanishes!!!\n",n), fflush(stdout);

        //Valence field force constants
        if(d->opt.ExecParam.Strain.StrainModel==Strain_struct::VFF_keating_strained) {
          int anh;
          if(((d->opt.ExecParam.PhonCalc.anh==PhonCalc_struct::Default)
              &&(d->opt.ExecParam.PhonCalc.PS==PhonCalc_struct::Lower))
             ||(d->opt.ExecParam.PhonCalc.anh==PhonCalc_struct::AcousticalSet)){
            anh=1;
          }
          else if (((d->opt.ExecParam.PhonCalc.anh==PhonCalc_struct::Default)
                    &&(d->opt.ExecParam.PhonCalc.PS==PhonCalc_struct::Higher))
                   ||(d->opt.ExecParam.PhonCalc.anh==PhonCalc_struct::OpticalSet)){
            anh=2;
          }
          else{
            if(!mpi_n3d_id){
              printf("\nThe requested option for the optimization of the anharmonicity corrections \n"
                     "has not been implemented yet!!!\n"
                     "Use anharmonicity corrections for better strain distribution instead.\n"
                     );
            }
            anh=1;
          }
          alpha[n] = d->strain.alpha_eff(AtomType_l_m,AtomType_cindx_aindx,dssq[n],d0sq[n],anh,true);
          /*  printf("alpha_eff[%d]=%g   alpha_eff[%d]=%g\n",
                 n,alpha[n],
                 n,d->strain.alpha[AtomType_l_m][AtomType_cindx_aindx]);*/
        }
        else if(d->opt.ExecParam.Strain.StrainModel==Strain_struct::VFF_keating) 
          alpha[n] = d->strain.alpha[AtomType_l_m][AtomType_cindx_aindx];
        beta[n]  =  d->strain.beta[AtomType_l_m][AtomType_cindx_aindx];
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        /*printf("\nd->strain.a_VFF[%d][%d]=%g  d->strain.b_VFF[%d][%d]=%g\n",
               AtomType_l_m,AtomType_cindx_aindx,
               d->strain.a_VFF[AtomType_l_m][AtomType_cindx_aindx],
               AtomType_l_m,AtomType_cindx_aindx,
               d->strain.b_VFF[AtomType_l_m][AtomType_cindx_aindx]);
        */

        
      }//Neighbor is within simulation domain
    }//Neighbor exists (cindx!=0)
  }//loop through n neighbors of m-th atom

  rm_rvectr(&qR);
  if(d->geo.StrainHasPeriodicity()){
    //Remove the vector of the periods of the structure in all 3 directions
    rm_rvectr(&D_ijk);
  }

}//end void initial_assignments

void Add_element(qd_struct d,
                 int n_row,int n_col,
                 real tail_r,real tail_i
                 )
/*Adds tail to the Dynamical matrix d->phon.DM[row][i-1] and 
  the corresponding column indexes to d->phon.Ind_DM[row][i].
  Here 0<i<=d->phon.states_in_row */
{
    int i;
    int IndDM_nrow=d->phon.Ind_DM[n_row][0];
    for(i=1; i<=IndDM_nrow; i++){
      /*search if the element present */
      if (d->phon.Ind_DM[n_row][i]==n_col){
        break;
      }
    }
    if ((i>d->phon.states_in_row)||(i>d->phon.N_cols[mpi_n3d_id])){
      printf("Too few positions reserved:\n");
    }
    else if (i>d->phon.Ind_DM[n_row][0])
      {//New element inserted
        (d->phon.Ind_DM[n_row][0])++;
        d->phon.Ind_DM[n_row][i]=n_col;
      }
    //Add tail to the element of DM
    d->phon.DM[n_row][i-1].r+=tail_r;
    d->phon.DM[n_row][i-1].i+=tail_i;
}//end of void Add_element

void Add(qd_struct d,
         int row, int col,
         real tail_r,real tail_i
         )
     /*Decide either add the element to DM on the current processor,
       what DM indexies n_row and n_col should be
       and call Add_element*/
{
  int YZblock=((d->geo.cell_ymax-d->geo.cell_ymin+1)
               *(d->geo.cell_zmax-d->geo.cell_zmin+1)
               *3*d->geo.cellgeom.AtomsPerCellMax());//# of DM rows/cols in 1 u.c. slab in x-direction
  int DM_total_size=3*d->n_atom_tot;//Total # of DM rows/cols
  int n_col,n_row;//Indexies in the DM block stored on the processor

  //column
  if(col>=d->phon.seg_col_s[mpi_n3d_id] && col<=d->phon.seg_col_s[mpi_n3d_id]+d->phon.seg_col_ln[mpi_n3d_id]-1){
    //Column is within the range on the processor -- keep
    if(row<d->phon.min_row[mpi_n3d_id]){ // discard unless this is the
                                     // last processor under periodic-X b.c.
      if(d->geo.StrainIsPeriodicX()){
        if(mpi_n3d_id==mpi_n3d_numprocs-1){
          if(row>=YZblock){ // should not happen
            printf("row=%d >= %d=YZblock!!! => There must be a mistake on %d processor!!!!",
                   row, YZblock,mpi_n3d_id);fflush(stdout);
          }
          else if (row<0){  // should not happen
            printf("row=%d < 0!!! => There must be a mistake on %d processor!!!!",
                   row, mpi_n3d_id);fflush(stdout);
          }
          else{//add to the upper block on the last proc under periodic b.c.
            n_col=col-d->phon.seg_col_s[mpi_n3d_id];
            n_row=row+d->phon.max_row[mpi_n3d_id]-d->phon.min_row[mpi_n3d_id]+1;
            Add_element(d,n_row,n_col,tail_r,tail_i);
            if(n_row>d->phon.N_rows[mpi_n3d_id]){// should not happen
              printf("n_row=%d > %d=d->phon.N_rows[%d] => "
                     "!!!!!!!!!There must be a mistake on %d processor!!!!!!!!",
                     n_row,d->phon.N_rows[mpi_n3d_id],mpi_n3d_id,mpi_n3d_id);fflush(stdout);
            }
            if(n_row<=d->phon.max_row[mpi_n3d_id]-d->phon.min_row[mpi_n3d_id]){// should not happen
              printf("n_row=%d <= %d=(d->phon.max_row[%d]=%d)-(d->phon.min_row[%d]=%d) => "
                     "!!!!!!!!!There is OVERWRITE on %d processor!!!!!!!!",
                     n_row,d->phon.max_row[mpi_n3d_id]-d->phon.min_row[mpi_n3d_id],
                     mpi_n3d_id,d->phon.max_row[mpi_n3d_id],mpi_n3d_id,d->phon.min_row[mpi_n3d_id],mpi_n3d_id);
              fflush(stdout);
            }
          }//element is within the upper block on the last processor
        }//we are on the last processor
      }//strain is periodic
    }//row is less then d->phon.min_row[mpi_n3d_id]
    else if (row > d->phon.max_row[mpi_n3d_id]){ // discard unless this is the
                                             // first processor under periodic-X b.c.
      if(d->geo.StrainIsPeriodicX()){
        if(!mpi_n3d_id){
          if(row>=DM_total_size){// should not happen
            printf("row=%d >= %d=DM_total_size => "
                   "!!!!!There must be a mistake on %d processor!!!!",
                   row,DM_total_size,mpi_n3d_id);fflush(stdout);
          }
          else if (row<DM_total_size-YZblock){// should not happen
            printf("row=%d < %d=(DM_total_size=%d)-(YZblock=%d) => "
                   "!!!!!There must be a mistake on %d processor!!!!",
                   row,DM_total_size-YZblock,DM_total_size,YZblock, mpi_n3d_id);
            fflush(stdout);
          }
          else{// Ok: lower left corner of DM 
            n_col=col-d->phon.seg_col_s[mpi_n3d_id];
            n_row=row-(DM_total_size-d->phon.N_rows[mpi_n3d_id]);
            Add_element(d,n_row,n_col,tail_r,tail_i);
            if(n_row>d->phon.N_rows[mpi_n3d_id]){// should not happen
              printf("n_row=%d > %d=d->phon.N_rows[%d] => "
                     "!!!!!!!!!There must be a mistake on %d processor!!!!!!!!",
                     n_row,d->phon.N_rows[mpi_n3d_id],mpi_n3d_id,mpi_n3d_id);fflush(stdout);
            }
            if(n_row<=d->phon.max_row[mpi_n3d_id]-d->phon.min_row[mpi_n3d_id]){// should not happen
              printf("n_row=%d <= %d=(d->phon.max_row[%d]=%d)-(d->phon.min_row[%d]=%d) => "
                     "!!!!!!!!!There is OVERWRITE on %d processor!!!!!!!!",
                     n_row,d->phon.max_row[mpi_n3d_id]-d->phon.min_row[mpi_n3d_id],
                     mpi_n3d_id,d->phon.max_row[mpi_n3d_id],mpi_n3d_id,d->phon.min_row[mpi_n3d_id],mpi_n3d_id);
              fflush(stdout);
            }
          }//element is within lower block on the first processor
        }//we are on the first processor
      }//Strain is periodic
    }//row > d->phon.max_row[mpi_n3d_id]
    else{ // regular case
      n_col=col-d->phon.seg_col_s[mpi_n3d_id];
      n_row=row-d->phon.min_row[mpi_n3d_id];
      Add_element(d,n_row,n_col,tail_r,tail_i);
      if(n_row>=d->phon.N_rows[mpi_n3d_id]){// should not happen
        printf("n_row=%d > %d=d->phon.N_rows[%d] => "
               "!!!!!!!!!There must be a mistake on %d processor!!!!!!!!",
               n_row,d->phon.N_rows[mpi_n3d_id],mpi_n3d_id,mpi_n3d_id);fflush(stdout);
      }
      if(n_row<0){// should not happen
        printf("n_row=%d < 0 => "
               "!!!!!!!!!There must be a mistake on %d processor!!!!!!!!",
               n_row,mpi_n3d_id);fflush(stdout);
      }
    }//d->phon.max_row[mpi_n3d_id] >= row >= d->phon.min_row[mpi_n3d_id]
  }//column is within the right range
}//end of void Add


void Add_Herm(qd_struct d,
	      int row,int col,
	      real tail_r,real tail_i
	      )
/*Add tail and tail* to "hermitian" parts of DM*/
{
  Add(d,row,col,tail_r,tail_i);
  Add(d,col,row,tail_r,-tail_i);
}

real DM_A1(
           real alpha,//bond length distorsion force constant array of N_Neighbor[m] elements
           real dssq,//square of the length of the strained bond
           real d0sq,//square of the length of the relaxed bond
           real nnv//component of the strained bond
           )
     /*Computes A1:
       A1=d^2 E_alpha/du_x(n)du_x(n)=d^2 E_alpha/du_x(m)du_x(m)
       -A1=d^2 E_alpha/du_x(n)du_x(m)=d^2 E_alpha/du_x(m)du_x(n)
     */
{
  return(4*0.375*alpha*(dssq-d0sq+2*nnv*nnv)/d0sq);
}

real DM_A2(
           real alpha,//bond length distorsion force constant array of N_Neighbor[m] elements
           real d0sq,//square of the length of the relaxed bond
           real x_ij,real y_ij//components of the strained bond
           )
     /*Computes A2:
       A2=d^2 E_alpha/du_x(n)du_y(n)=d^2 E_alpha/du_x(m)du_y(m)
       -A2=d^2 E_alpha/du_x(n)du_y(m)=d^2 E_alpha/du_x(m)du_y(n)
     */
{
  return(8*0.375*alpha*(x_ij*y_ij)/d0sq);
}

real doublet(real beta_eff,//effective beta=sqrt(beta1*beta2), which may be scaled with strain
             //real beta1,real beta2,//bond bending distorsion force constant
             real d0sq1,real d0sq2//square of the length of the relaxed bond
             )
     /*Computes coefficient
       0.375*sqrt(beta1*beta2)/sqrt(d0sq1*d0sq2)
       for the given V-shaped bond doublet*/
{
  return(2*0.375*beta_eff/*sqrt(beta1*beta2)*//sqrt(d0sq1*d0sq2));
}

real DM_D(
          rvectr d1,rvectr d2,//vectors of the relaxed bonds of the V doublet
          rvectr r1,rvectr r2//vectors of the strained bonds of the V doublet
          )
     /*Computes D:
       Term for bond bending distorsion
       to be added in different elements of DM
       with the same displacement projection
     */
{
  return((r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2]-d1[0]*d2[0]-d1[1]*d2[1]-d1[2]*d2[2]));
}

void Add_A1(qd_struct d,
            real A1,//element to be added
            int mDM,//row
            int nDM,//column
            real sqrtMm,//square root of the atomic mass of m-th atom
            real sqrtMn,//square root of the atomic mass ofthe n-th neighbor of m-th atom
            real phase_r,//real part of the phase
            real phase_i//imaginary part of the phase
            )
     /*Add elements A1 to Dynamical Matrix*/
{
  Add(d,mDM,mDM,A1/(sqrtMm*sqrtMm),0.0);
  Add(d,nDM,nDM,A1/(sqrtMn*sqrtMn),0.0);
  Add_Herm(d,mDM,nDM,
           -A1/(sqrtMn*sqrtMm)*phase_r,
           -A1/(sqrtMn*sqrtMm)*phase_i);
}


void Add_A2(qd_struct d,
            real A2,//element to be added
            int l,//cell
            int m,//atom
            int xim,//displacement component of the m-th atom
            int cindx,//cell of the n-th neighbor of the m-th atom
            int aindx,//atom of the n-th neighbor of the m-th atom
            int xin,//displacement component of the n-th neighbor of the m-th atom
            real sqrtMm,//square root of the atomic mass of m-th atom
            real sqrtMn,//square root of the atomic mass ofthe n-th neighbor of m-th atom
            real phase_r,//real part of the phase
            real phase_i//imaginary part of the phase
            )
     /*Add elements A2 to Dynamical Matrix*/
{
  int mDM=d->phon.DMmap[l][m][xim];
  int m1DM=d->phon.DMmap[l][m][xin];
  int nDM=d->phon.DMmap[cindx][aindx][xin];
  int n1DM=d->phon.DMmap[cindx][aindx][xim];
  Add_Herm(d,mDM,m1DM,A2/(sqrtMm*sqrtMm),0.0);
  Add_Herm(d,nDM,n1DM,A2/(sqrtMn*sqrtMn),0.0);
  Add_Herm(d,mDM,nDM,
           -A2/(sqrtMn*sqrtMm)*phase_r,
           -A2/(sqrtMn*sqrtMm)*phase_i);
  Add_Herm(d,m1DM,n1DM,
           -A2/(sqrtMn*sqrtMm)*phase_r,
           -A2/(sqrtMn*sqrtMm)*phase_i);
}

void Add_BXX(qd_struct d,
             real V,//factor for the V doublet
             real D,//=(r_mn1*r_mn2-d_mn1*d_mn2)
             real xn1,real xn2,//components of the bond length between m-th atom and n1(n2) neighbor
             int cell_m,int atom_m,//cell and atom indexes of the vertex atom
             int cell_n1,int atom_n1,//cell and atom indexes of the n1-th neighbor atom
             int cell_n2,int atom_n2,//cell and atom indexes of the n2-th neighbor atom
             int xi,//projection of the displacement
             real sqrtMm,//square root of the atomic mass of m-th atom
             real sqrtMn1,//square root of the atomic mass of the n1-th neighbor of m-th atom
             real sqrtMn2,//square root of the atomic mass of the n2-th neighbor of m-th atom
             real phase_n1_r,//real part of the phase of the n1-th neighbor of m-th atom
             real phase_n1_i,//imaginary part of the phase of the n1-th neighbor of m-th atom
             real phase_n2_r,//real part of the phase of the n2-th neighbor of m-th atom
             real phase_n2_i//imaginary part of the phase of the n2-th neighbor of m-th atom            
             )
     /*Add elements B to Dynamical Matrix
       when projections of the displacements are the same*/
{
  Add(d,d->phon.DMmap[cell_m][atom_m][xi],d->phon.DMmap[cell_m][atom_m][xi],
      V*(2*D+(xn1+xn2)*(xn1+xn2))/(sqrtMm*sqrtMm),0.0);
  Add(d,d->phon.DMmap[cell_n1][atom_n1][xi],d->phon.DMmap[cell_n1][atom_n1][xi],
      V*xn2*xn2/(sqrtMn1*sqrtMn1),0.0);
  Add(d,d->phon.DMmap[cell_n2][atom_n2][xi],d->phon.DMmap[cell_n2][atom_n2][xi],
      V*xn1*xn1/(sqrtMn2*sqrtMn2),0.0);
  Add_Herm(d,d->phon.DMmap[cell_n1][atom_n1][xi],d->phon.DMmap[cell_n2][atom_n2][xi],
           V*(D+xn1*xn2)/(sqrtMn1*sqrtMn2)*(phase_n1_r*phase_n2_r+phase_n1_i*phase_n2_i),
           V*(D+xn1*xn2)/(sqrtMn1*sqrtMn2)*(phase_n1_r*phase_n2_i-phase_n1_i*phase_n2_r));
  Add_Herm(d,d->phon.DMmap[cell_m][atom_m][xi],d->phon.DMmap[cell_n1][atom_n1][xi],
           -V*(D+xn2*(xn1+xn2))/(sqrtMm*sqrtMn1)*phase_n1_r,
           -V*(D+xn2*(xn1+xn2))/(sqrtMm*sqrtMn1)*phase_n1_i);
  Add_Herm(d,d->phon.DMmap[cell_m][atom_m][xi],d->phon.DMmap[cell_n2][atom_n2][xi],
           -V*(D+xn1*(xn1+xn2))/(sqrtMm*sqrtMn2)*phase_n2_r,
           -V*(D+xn1*(xn1+xn2))/(sqrtMm*sqrtMn2)*phase_n2_i);
  
}

void Add_BX1X2(qd_struct d,
             real V,//factor for the V doublet
             real x1n1,real x1n2,real x2n1,real x2n2,//components of the bond length between m-th atom and n1(n2) neighbor
             int cell_m,int atom_m,//cell and atom indexes of the vertex atom
             int cell_n1,int atom_n1,//cell and atom indexes of the n1-th neighbor atom
             int cell_n2,int atom_n2,//cell and atom indexes of the n2-th neighbor atom
             int xi1,int xi2,//projections of the displacement
             real sqrtMm,//square root of the atomic mass of m-th atom
             real sqrtMn1,//square root of the atomic mass of the n1-th neighbor of m-th atom
             real sqrtMn2,//square root of the atomic mass of the n2-th neighbor of m-th atom
             real phase_n1_r,//real part of the phase of the n1-th neighbor of m-th atom
             real phase_n1_i,//imaginary part of the phase of the n1-th neighbor of m-th atom
             real phase_n2_r,//real part of the phase of the n2-th neighbor of m-th atom
             real phase_n2_i//imaginary part of the phase of the n2-th neighbor of m-th atom            
             )
     /*Add elements B to Dynamical Matrix
       when projections of the displacements are different*/
{
  Add_Herm(d,d->phon.DMmap[cell_m][atom_m][xi1],d->phon.DMmap[cell_m][atom_m][xi2],
           V*((x1n1+x1n2)*(x2n1+x2n2))/(sqrtMm*sqrtMm),0.0);
  Add_Herm(d,d->phon.DMmap[cell_n1][atom_n1][xi1],d->phon.DMmap[cell_n1][atom_n1][xi2],
           V*x1n2*x2n2/(sqrtMn1*sqrtMn1),0.0);
  Add_Herm(d,d->phon.DMmap[cell_n2][atom_n2][xi1],d->phon.DMmap[cell_n2][atom_n2][xi2],
           V*x1n1*x2n1/(sqrtMn2*sqrtMn2),0.0);
  Add_Herm(d,d->phon.DMmap[cell_n1][atom_n1][xi1],d->phon.DMmap[cell_n2][atom_n2][xi2],
           V*x1n2*x2n1/(sqrtMn1*sqrtMn2)*(phase_n1_r*phase_n2_r+phase_n1_i*phase_n2_i),
           V*x1n2*x2n1/(sqrtMn1*sqrtMn2)*(phase_n1_r*phase_n2_i-phase_n1_i*phase_n2_r));
  Add_Herm(d,d->phon.DMmap[cell_n1][atom_n1][xi2],d->phon.DMmap[cell_n2][atom_n2][xi1],
           V*x1n1*x2n2/(sqrtMn1*sqrtMn2)*(phase_n1_r*phase_n2_r+phase_n1_i*phase_n2_i),
           V*x1n1*x2n2/(sqrtMn1*sqrtMn2)*(phase_n1_r*phase_n2_i-phase_n1_i*phase_n2_r));
  Add_Herm(d,d->phon.DMmap[cell_m][atom_m][xi1],d->phon.DMmap[cell_n1][atom_n1][xi2],
           -V*(x1n1+x1n2)*x2n2/(sqrtMm*sqrtMn1)*phase_n1_r,
           -V*(x1n1+x1n2)*x2n2/(sqrtMm*sqrtMn1)*phase_n1_i);
  Add_Herm(d,d->phon.DMmap[cell_m][atom_m][xi2],d->phon.DMmap[cell_n1][atom_n1][xi1],
           -V*(x2n1+x2n2)*x1n2/(sqrtMm*sqrtMn1)*phase_n1_r,
           -V*(x2n1+x2n2)*x1n2/(sqrtMm*sqrtMn1)*phase_n1_i);
  Add_Herm(d,d->phon.DMmap[cell_m][atom_m][xi1],d->phon.DMmap[cell_n2][atom_n2][xi2],
           -V*(x1n1+x1n2)*x2n1/(sqrtMm*sqrtMn2)*phase_n2_r,
           -V*(x1n1+x1n2)*x2n1/(sqrtMm*sqrtMn2)*phase_n2_i);
  Add_Herm(d,d->phon.DMmap[cell_m][atom_m][xi2],d->phon.DMmap[cell_n2][atom_n2][xi1],
           -V*(x2n1+x2n2)*x1n1/(sqrtMm*sqrtMn2)*phase_n2_r,
           -V*(x2n1+x2n2)*x1n1/(sqrtMm*sqrtMn2)*phase_n2_i);

}

#define swapInt(g, h)  { int y=(g); (g)=(h); (h)=y; }
#define swapReal(g, h)  { real y=(g); (g)=(h); (h)=y; }

void comb_DM(ivectr index,//[d->phon.Ind_DM[n_row][i]+1] vector d->phon.Ind_DM[n_row][i]
                    cvectr DMrow,//[d->phon.Ind_DM[n_row][i]] vector d->phon.DM[n_row][i-1]
                    int beg, int end//lower and upper bounds of the arrays
                    )
     /*Modification of quick_sort to
      sort d->phon.Ind_DM[n_row][i]
      in ascending order.
      Also complex elements of DM
      d->phon.DM[n_row][i-1] are swapped.*/
{
   int p, i;
   
   if (beg>=end) return;
   
   swapInt(index[beg], index[(beg+end)/2]);
   swapReal(DMrow[beg-1].r, DMrow[((beg+end)/2)-1].r);
   swapReal(DMrow[beg-1].i, DMrow[((beg+end)/2)-1].i);
   
   p = beg;
   for (i = beg+1; i<=end; i++)
      if(index[i] < index[beg]){
	 ++p;
	 swapInt(index[i], index[p]);
         swapReal(DMrow[i-1].r, DMrow[p-1].r);
         swapReal(DMrow[i-1].i, DMrow[p-1].i);
      }
   
   swapInt(index[beg], index[p]);
   swapReal(DMrow[beg-1].r, DMrow[p-1].r);
   swapReal(DMrow[beg-1].i, DMrow[p-1].i);
   
   comb_DM(index,DMrow, beg, p-1);
   comb_DM(index,DMrow, p+1, end);
}

void DM_zeros(qd_struct d)
     //Assignes zeros to elements with value less then zero
{
  real zero=1e-10;
  for(int row=0;row<d->phon.N_rows[mpi_n3d_id];row++){
    for(int i=1;i<=d->phon.Ind_DM[row][0];i++){
      if(fabs(d->phon.DM[row][i-1].r)<zero){
        if(d->phon.DM[row][i-1].r)
#ifdef PR_ZEROS_SUBSTITUTION
          printf("\nDM[%d][%d].r=%g<%g => DM[%d][%d].r=0",
                 row,d->phon.Ind_DM[row][i],d->phon.DM[row][i-1].r,
                 zero,row,d->phon.Ind_DM[row][i]);
#endif// PR_ZEROS_SUBSTITUTION
        d->phon.DM[row][i-1].r=0;
      }
      if(fabs(d->phon.DM[row][i-1].i)<zero){
        if(d->phon.DM[row][i-1].i)
#ifdef PR_ZEROS_SUBSTITUTION
          printf("\nDM[%d][%d].i=%g<%g => DM[%d][%d].r=0",
                 row,d->phon.Ind_DM[row][i],d->phon.DM[row][i-1].i,
                 zero,row,d->phon.Ind_DM[row][i]);
#endif// PR_ZEROS_SUBSTITUTION
        d->phon.DM[row][i-1].i=0;
      }
    }
  }
}

void DM_Keating(qd_struct d,
		int n_qph, //-th point in the reciprocal space
                int cell0,//initial cell index
                int cellF//final cell index
		)
     // Generates valence-field-force part of dynamical matrix
{
  cvectr phase;//=exp[i(q_ph.R)] -- phase factor for periodic boundary conditions
  real sqrtMm;//atomic masses square roots of atom m
  rvectr sqrtMn;//arrays of atomic masses square roots neighbor atoms
  int nDM,mDM;//index of the atom in the Dinamical Matrix
  real A1,A2,D;/*dummies for terms to be added to different elements of Dynamical Matrix*/
  int cnbr2,anbr2;//cell and atom of the second neighbor
  real Vcoeff;//coefficient for V doublet
  const int maxNbr = d->geo.NeighborsMax();
  int i, j, k;//indexies of the atom's u.c.
  int i_nbr,j_nbr,k_nbr;//indexies of the neighbor's u.c.
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   int AtomType_l_m;
   int aindx, cindx;

   bool has_periodicity = d->geo.StrainHasPeriodicity();

   
   //Memory allocation and initial assignments -- begin
   real Beta;
   rvectr  alpha = Rvectr(maxNbr);
   rvectr  beta = Rvectr(maxNbr);
   rmatrix nnv = Rmatrix(maxNbr,3);
   rmatrix nnv0 = Rmatrix(maxNbr,3);
   rvectr  dssq = Rvectr(maxNbr);//strained bonds squared between m-th atom and its neighbors
   rvectr  d0sq = Rvectr(maxNbr);//relaxed bonds squared between m-th atom and its neighbors

   phase = Cvectr(maxNbr);//Phase m->n
   sqrtMn = Rvectr(maxNbr);//Square roots of mass of neighbor atoms

   //Memory allocation and initial assignments -- end

   // MAIN LOOPS
   /* Loop through all unit cells associated with this processor */
   for (int l=cell0; l <= cellF; l++ ){
     /* These are real space indexes for the unit cell */
     i = d->geo.cell__ijk[l][0];
     j = d->geo.cell__ijk[l][1];
     k = d->geo.cell__ijk[l][2];
     
     /* Loop through all the atoms in the unit cell within simulation domain */
     
     for (int m=0; m < d->geo.AtomsPerCellMax(); m++ ) {
       AtomType_l_m = (int) d->geo.AtomType[l][m];
       
       /*  Check if atom is in simulation domain */ 
       if ( !AtomType_l_m ) continue;
       
       sqrtMm=sqrt(d->Amass[AtomType_l_m]);//atomic mass square root
       if (sqrtMm==0.0){//Should not happen
         die("Encountered a zero sqrtMm\n");
       }
       
       /* Loop through the atom's nearest neighbors */
       //First loop -- to make all the assignments
       initial_assignments(d,n_qph,l,m,phase,alpha,beta,sqrtMn,nnv,nnv0,dssq,d0sq);
       for(int xim=0; xim<3; xim++){
         mDM=d->phon.DMmap[l][m][xim];	     
         for (int n=0; n < d->geo.Neighbors(m); n++ ) {
           /* Get relative position of unit cell where this neighbor lives. */ 
           i_nbr = i + d->geo.NbrCell(m,n,0);
           j_nbr = j + d->geo.NbrCell(m,n,1);
           k_nbr = k + d->geo.NbrCell(m,n,2);
           /* Get the unit cell and atomic index for this neighbor */

#ifdef ELIMINATE_SSMAP
           cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];
           if (cindx < 0 && has_periodicity && (cindx + d->geo.N_Cell + 1 ) >= 0) {
             cindx = cindx + d->geo.N_Cell + 1 ;
           }
#else /* ELIMINATE_SSMAP */
           cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];
           if (cindx < 0 && has_periodicity && d->geo.ssmap[l][m][n]!=-1) {
             cindx = d->geo.ssmap[l][m][n];
           }
#endif /* ELIMINATE_SSMAP */


           aindx = d->geo.NbrCell(m,n,3);
           
           /* Check if neighbor's unit cell is within simulation domain. */
           if ( cindx >= 0 ) {
             if ( d->geo.AtomType[cindx][aindx] ){
               //////////////////////////////////////////////////////////////////
               //                                                              //
               // Call functions to fill the elements in the dynamical matrix  //
               //                                                              //
               //////////////////////////////////////////////////////////////////
               
               /* Loop through  components of the displacement vector of n-th neighbor of atom m*/
               for (int xin=xim; xin<3; xin++){
                 if((nDM=d->phon.DMmap[cindx][aindx][xin])>mDM){
                   // *** make sure no bond connects to itself ***//
                   
                   //If this bond has not been considered yet
                   ///////////////////////////////////
                   // Compute bond stretching terms //
                   ///////////////////////////////////
                   if(xim==xin){//xx
#ifdef A1DM
                     A1=DM_A1(alpha[n],dssq[n],d0sq[n],nnv[n][xim]);
                     Add_A1(d,A1,mDM,nDM,sqrtMm,sqrtMn[n],phase[n].r,phase[n].i);
#endif/*A1DM*/
                   }//xx
                   else{//xy
#ifdef A2DM
                     A2=DM_A2(alpha[n],d0sq[n],nnv[n][xin],nnv[n][xim]);
                     Add_A2(d,A2,l,m,xim,cindx,aindx,xin,sqrtMm,sqrtMn[n],phase[n].r,phase[n].i);
#endif/*A2DM*/
                   }//xy
                 }//if this bond has not been considered yet (i.e., nDM>mDM)
                 /*Second-nearest neighbors loop*/
                 for(int nbr2=n+1; nbr2<d->geo.Neighbors(m); nbr2++){
                   cnbr2=d->geo.ijk__cell
                     [i+d->geo.NbrCell(m,nbr2,0)]
                     [j+d->geo.NbrCell(m,nbr2,1)]
                     [k+d->geo.NbrCell(m,nbr2,2)];
#ifdef ELIMINATE_SSMAP
                   if(cnbr2<0 && has_periodicity
                      && (cnbr2 + d->geo.N_Cell + 1 ) >= 0){//cell index
                     cnbr2= cnbr2 + d->geo.N_Cell + 1 ;
                   }
#else /* ELIMINATE_SSMAP */
                   if(cnbr2<0 && has_periodicity
                      && d->geo.ssmap[l][m][nbr2]!=-1){//cell index
                     cnbr2=d->geo.ssmap[l][m][nbr2];
                   }
#endif /* ELIMINATE_SSMAP */
                   //atom index
                   anbr2=d->geo.NbrCell(m,nbr2,3);
                   //Check if this second neihgbor is within simulation domain
                   if(cnbr2>=0){//nbr2 is within simulation domain
                     if(d->geo.AtomType[cnbr2][anbr2]){//and nbr2 has not been removed
                       ////////////////////////////////
                       // Compute bond bending terms //
                       ////////////////////////////////
                       Beta=0.0;
                       if(d->opt.ExecParam.Strain.StrainModel==Strain_struct::VFF_keating_strained) {
                         /*  cout << "\nInput data for getting Beta:"
                              << "\n Type of "
                              << m << "-th atom: " << (int) d->geo.AtomType[l][m]
                              << "\n Type of the "
                              << n <<"-th neighbor: "<< (int) d->geo.AtomType[cindx][aindx]
                              << "\n Type of the "
                              << nbr2 << "-th neighbor: " << (int) d->geo.AtomType[cnbr2][anbr2]
                              << "\n Realxed bond to "
                              << n <<"-th neighbor:  ("
                              << nnv0[n][0] << "," << nnv0[n][1] << "," <<nnv0[n][2] <<")"
                              << "\n Realxed bond to "
                              << nbr2 <<"-th neighbor:  ("
                              << nnv0[nbr2][0] << "," << nnv0[nbr2][1] << "," <<nnv0[nbr2][2] <<")"
                              << "\n Strained bond to "
                              << n <<"-th neighbor:  ("
                              << nnv[n][0] << "," << nnv[n][1] << "," <<nnv[n][2] <<")"
                              << "\n Strained bond to "
                              << nbr2 <<"-th neighbor:  ("
                              << nnv[nbr2][0] << "," << nnv[nbr2][1] << "," <<nnv[nbr2][2] <<")"
                              << endl;*/
                         int anh;
                         if(((d->opt.ExecParam.PhonCalc.anh==PhonCalc_struct::Default)
                             &&(d->opt.ExecParam.PhonCalc.PS==PhonCalc_struct::Lower))
                            ||(d->opt.ExecParam.PhonCalc.anh==PhonCalc_struct::AcousticalSet)){
                           anh=1;
                         }
                         else if (((d->opt.ExecParam.PhonCalc.anh==PhonCalc_struct::Default)
                                   &&(d->opt.ExecParam.PhonCalc.PS==PhonCalc_struct::Higher))
                                  ||(d->opt.ExecParam.PhonCalc.anh==PhonCalc_struct::OpticalSet)){
                           anh=2;
                         }
                         else{
                           if(!mpi_n3d_id){
                             printf("\nThe requested option for the optimization of the anharmonicity corrections \n"
                                    "has not been implemented yet!!!\n"
                                    "Use anharmonicity corrections for better strain distribution instead.\n"
                                    );
                           }
                           anh=1;
                         }
                         Beta=d->strain.beta_eff((int) d->geo.AtomType[l][m],//type of the m-th atom
                                                 (int) d->geo.AtomType[cindx][aindx],//type of the n1-th atom
                                                 (int) d->geo.AtomType[cnbr2][anbr2],//type of the n2-th atom
                                                 nnv0[n],nnv0[nbr2],//vectors of the relaxed bonds of the V doublet
                                                 nnv[n],nnv[nbr2],//vectors of the strained bonds of the V doublet
                                                 anh,//variant of the anharmonicity corrections optimization
                                                 true);
                         /*printf("%d(%d)<-%d->%d(%d) : beta_eff=%g   beta=%g\n",
                                n,aindx,m,nbr2,anbr2,
                                Beta,
                                sqrt(beta[n]*beta[nbr2]));*/
                       }
                       else if(d->opt.ExecParam.Strain.StrainModel==Strain_struct::VFF_keating)
                         Beta=sqrt(beta[n]*beta[nbr2]);
                       
                       Vcoeff=doublet(Beta,d0sq[n],d0sq[nbr2]);
                       if(xim==xin){//x1=x2
#ifdef BXXDM
                         D=DM_D(nnv0[n],nnv0[nbr2],nnv[n],nnv[nbr2]);
                         Add_BXX(d,
                                 Vcoeff,
                                 D,nnv[n][xim],nnv[nbr2][xim],
                                 l,m,
                                 cindx,aindx,
                                 cnbr2,anbr2,
                                 xim,
                                 sqrtMm,sqrtMn[n],sqrtMn[nbr2],
                                 phase[n].r,phase[n].i,
                                 phase[nbr2].r,phase[nbr2].i);
#endif/*BXXDM*/
                       }//x1=x2
                       else{//x1<x2
#ifdef BX1X2DM
                         Add_BX1X2(d,Vcoeff,
                                   nnv[n][xim],nnv[nbr2][xim],nnv[n][xin],nnv[nbr2][xin],
                                   l,m,cindx,aindx,cnbr2,anbr2,xim,xin,
                                   sqrtMm,sqrtMn[n],sqrtMn[nbr2],
                                   phase[n].r,phase[n].i,phase[nbr2].r,phase[nbr2].i);
#endif/*BX1X2DM*/
                       }//x1<x2
                     }//and nbr2 has not been removed
                   }//nbr2 is within simulation domain
                 }//Second-neighbor nbr loop
                 
               }//xin loop through components of the displacement vector of n-th neighbor of m-th atom
             }//if n-th neighbor of m-th  atom has not been removed
           }//if cindx within simulation domain
         }//loop through n neighbors of m-th atom
       }//loop through xim projections of m-th atom displacement
     }//loop through m atoms of l-th unit cell
   }//loop through l unit cells on the processor
   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DISREGARD_SMALL_ELEMENTS
   //Assigne zeros to elements with value less then "zero"
   DM_zeros(d);
#endif//DISREGARD_SMALL_ELEMENTS
   //Sort elements of DM in ascending order of their indexes
   for(int row=0;row<d->phon.N_rows[mpi_n3d_id];row++){
     comb_DM(d->phon.Ind_DM[row],
             d->phon.DM[row],
             1,d->phon.Ind_DM[row][0]//lower and upper bounds of the arrays
             );
   }//sort elements of DM in ascending order of their indexes
   
   rm_rmatrix(&nnv);
   rm_rmatrix(&nnv0);
   rm_rvectr(&alpha);
   rm_rvectr(&beta);
   rm_rvectr(&dssq);
   rm_rvectr(&d0sq);
   rm_cvectr(&phase);
   rm_rvectr(&sqrtMn);
 
}//END of DM_Keating


/**
    Computes phonon dispersion for 
    d->Nq_ph points in the reciprocal space starting from
    d->q_ph_min to d->q_ph_min.
*/
void phon_disp(qd_struct d)

{

  real q_ph_max,q_ph_min;//maximal and minimal values of the phonon wave vector
  
#ifdef TALK_A_LOT
  masterPrint("I am in phon_disp function now...\n");
#endif /* TALK_A_LOT*/

  d->phon.ph_s_alloc(d->geo.N_Cell,
                     d->geo.cellgeom.AtomsPerCellMax(),
		     d->opt.ExecParam.PhonCalc.N_branches,
		     d->opt.ExecParam.PhonCalc.Nq_ph,
                     d->nproc);
  /* If we are on edge, processor to past edge is NULL (-1) unless we
     have periodicity in x.*/
   d->phon.mpi_left[mpi_n3d_id] = ( mpi_n3d_id 
                                ? mpi_n3d_id-1 
                                : (d->geo.StrainIsPeriodicX() ? mpi_n3d_numprocs-1 : -1) );
   d->phon.mpi_right[mpi_n3d_id] = ( mpi_n3d_id < mpi_n3d_numprocs-1 
                                 ? mpi_n3d_id+1 
                                 : (d->geo.StrainIsPeriodicX() ? 0 : -1) );
   int YZslab=((d->geo.cell_ymax-d->geo.cell_ymin+1)
               *(d->geo.cell_zmax-d->geo.cell_zmin+1));//# of unit cells in 1 u.c. slab in x-direction
   int rem_slabs=(d->geo.cell_xmax-d->geo.cell_xmin+1)%mpi_n3d_numprocs;/*# of processors that have
                                                                      1 extra XYslab.
                                                                      Those slabs are added to the
                                                                      last processors.*/
   int cell_ln = YZslab*((d->geo.cell_xmax-d->geo.cell_xmin+1)/mpi_n3d_numprocs
                         +((mpi_n3d_id<(mpi_n3d_numprocs-rem_slabs))
                           ? 0 : 1));
   int cell_s = (mpi_n3d_id*YZslab*((d->geo.cell_xmax-d->geo.cell_xmin+1)/mpi_n3d_numprocs)+
                 ((mpi_n3d_numprocs-mpi_n3d_id<rem_slabs) ? YZslab*(rem_slabs-(mpi_n3d_numprocs-mpi_n3d_id)):0));
   int cell0 = cell_s;
   int cellF = cell0 + cell_ln-1;
   d->phon.seg_col_s[mpi_n3d_id]=3*cell_s*d->geo.cellgeom.AtomsPerCellMax();
   d->phon.seg_col_ln[mpi_n3d_id]=3*cell_ln*d->geo.cellgeom.AtomsPerCellMax();
   if(mpi_n3d_id)
     cell0 -= YZslab;
   if(mpi_n3d_id!=mpi_n3d_numprocs-1)
     cellF += YZslab;
   d->phon.szCommL[mpi_n3d_id] = ((d->phon.mpi_left[mpi_n3d_id]>=0)&&(d->phon.mpi_left[mpi_n3d_id]!=mpi_n3d_id) ?
                              3*YZslab*d->geo.cellgeom.AtomsPerCellMax() : 0);//Number of rows
                                                                              //in the buffer
                                                                              //to be sent to the
                                                                              //left processor
   d->phon.szCommR[mpi_n3d_id] = ((d->phon.mpi_right[mpi_n3d_id]>=0)&&(d->phon.mpi_right[mpi_n3d_id]!=mpi_n3d_id) ?
                              3*YZslab*d->geo.cellgeom.AtomsPerCellMax() : 0);//Number of rows
                                                                              //in the buffer
                                                                              //to be sent to the
                                                                              //right processor

   d->phon.N_cols[mpi_n3d_id]=3*cell_ln*d->geo.cellgeom.AtomsPerCellMax();//Number of columns of DM
                                                                      //stored on the processor
   d->phon.N_rows[mpi_n3d_id]=3*d->geo.cellgeom.AtomsPerCellMax()*
     (cellF-cell0+1+((d->geo.StrainIsPeriodicX() && (mpi_n3d_numprocs!=1) &&
                      ((!mpi_n3d_id) || (mpi_n3d_id == mpi_n3d_numprocs-1)))
                     ? ((d->geo.cell_ymax-d->geo.cell_ymin+1)
                        *(d->geo.cell_zmax-d->geo.cell_zmin+1)) : 0));//Number of rows of DM
                                                                      //stored on the processor
   if(d->phon.N_rows[mpi_n3d_id]>3*d->n_atom_tot)
     d->phon.N_rows[mpi_n3d_id]=3*d->n_atom_tot;
   for(int i=0;i<mpi_n3d_numprocs;i++){
     MPI_Bcast(&d->phon.szCommL[i], 1, MPI_INT, i, MPI_COMM_WORLD);
     MPI_Bcast(&d->phon.szCommR[i], 1, MPI_INT, i, MPI_COMM_WORLD);
     MPI_Bcast(&d->phon.N_rows[i], 1, MPI_INT, i, MPI_COMM_WORLD);
     MPI_Bcast(&d->phon.N_cols[i], 1, MPI_INT, i, MPI_COMM_WORLD);
   }
   //Check if the matrix is folded at periodic X b.c.
   if(mpi_n3d_numprocs==2 && (d->phon.N_cols[mpi_n3d_numprocs-1]+
                          d->phon.szCommL[mpi_n3d_numprocs-1]+
                          d->phon.szCommR[mpi_n3d_numprocs-1]>
                          d->phon.N_rows[mpi_n3d_numprocs-1])){
     //Matrix is "folded" => Update the size of communication buffers
     //                   => Reduce number of communication
     d->phon.szCommR[mpi_n3d_numprocs-1]=0;
     if(d->phon.N_cols[0]+d->phon.szCommL[0]+d->phon.szCommR[0]<=d->phon.N_rows[0]){
       d->phon.szCommR[d->phon.mpi_right[mpi_n3d_numprocs-1]]+=
         d->phon.szCommL[d->phon.mpi_right[mpi_n3d_numprocs-1]];
     }
     d->phon.szCommL[d->phon.mpi_right[mpi_n3d_numprocs-1]]=0;
     d->phon.mpi_left[d->phon.mpi_right[mpi_n3d_numprocs-1]]=-1;
     d->phon.mpi_right[mpi_n3d_numprocs-1]=-1;
   }//communication update
   if(d->phon.N_cols[mpi_n3d_id]+d->phon.szCommL[mpi_n3d_id]+d->phon.szCommR[mpi_n3d_id]>d->phon.N_rows[mpi_n3d_id]){
     //Should not happen
     printf("\n!!!!![(d->phon.N_cols[%d]=%d)+(d->phon.szCommL[%d]=%d)+"
            "(d->phon.szCommR[%d]=%d)]=%d > %d=d->phon.N_rows[%d]!!!!!\n",
            mpi_n3d_id,d->phon.N_cols[mpi_n3d_id],
            mpi_n3d_id,d->phon.szCommL[mpi_n3d_id],
            mpi_n3d_id,d->phon.szCommR[mpi_n3d_id],
            d->phon.N_cols[mpi_n3d_id]+d->phon.szCommL[mpi_n3d_id]+d->phon.szCommR[mpi_n3d_id],
            d->phon.N_rows[mpi_n3d_id],mpi_n3d_id);
     fflush(stdout);
     }
   ///////////////////////////////////////////////////////////////////
   //Allocate the Dynamical matrix and the buffers for communication//
   ///////////////////////////////////////////////////////////////////
   d->phon.DM_allocate(d->phon.N_rows[mpi_n3d_id],//Number of rows
                       (d->phon.N_cols[mpi_n3d_id]<d->phon.states_in_row ?
                        d->phon.N_cols[mpi_n3d_id] : d->phon.states_in_row)//Max number of elements
                                                                       //per row on
                                                                       //mpi_n3d_id proc
                       );
   
#ifdef TALK_A_LOT
   masterPrint("%d points of %d branches of phonon dispersion to be computed within Keating model",
               d->opt.ExecParam.PhonCalc.Nq_ph,d->opt.ExecParam.PhonCalc.N_branches);
#endif /*TALK_A_LOT*/
  ///////////////////////////////////////////////////////////////////
  //      Mapping of the indexies of the dynamical matrix in       //
  //d->phon.DMmap[d->geo.N_Cell][d->geo.cellgeom.AtomsPerCellMax()]//
  ///////////////////////////////////////////////////////////////////
  Map_DM(d);
  ///////////////////////////////////////////////////////////////////

  d->phon.min_row[mpi_n3d_id]=d->phon.DMmap[cell0][0][0];
  d->phon.max_row[mpi_n3d_id]=d->phon.DMmap[cellF][d->geo.cellgeom.AtomsPerCellMax()-1][2];

#ifdef PR_PHON_PARAMETERS
  d->phon.Print_parameters();
#endif//PR_PHON_PARAMETERS
   real alat=d->strain.amat[(int) d->geo.AtomType[0][0]][(int) d->geo.AtomType[0][4]];
   //Measure of hydrostatic strain
   real V_over_V0=((d->geo.lattice_x*((d->geo.cell_xmax-d->geo.cell_xmin)+1)
                    *d->geo.lattice_y*((d->geo.cell_ymax-d->geo.cell_ymin)+1)
                    *d->geo.lattice_z*((d->geo.cell_zmax-d->geo.cell_zmin)+1))/
                   (alat*((d->geo.cell_xmax-d->geo.cell_xmin)+1)
                    *alat*((d->geo.cell_ymax-d->geo.cell_ymin)+1)
                    *alat*((d->geo.cell_zmax-d->geo.cell_zmin)+1)));
   //Measure of biaxial strain
   real a_over_a0=(d->geo.lattice_x/alat);
   //Change for output:
   V_over_V0=a_over_a0;
  /*Assignment of qx,qy,qz in the first 3 columns of 
    d->phon.Eph[Nq_ph][N_branches+3]*/ 
#define BZ_BOUNDARY
  real dq;
  for (int i=0;i<3;i++){
  //Assignment of increments
#ifdef BZ_BOUNDARY
    real D;
      if(i==0)
        D=d->geo.lattice_x*((d->geo.cell_xmax-d->geo.cell_xmin)+1);
      if(i==1)
        D=d->geo.lattice_y*((d->geo.cell_ymax-d->geo.cell_ymin)+1);
      if(i==2)
        D=d->geo.lattice_z*((d->geo.cell_zmax-d->geo.cell_zmin)+1);

      q_ph_max=Pi/D*d->opt.ExecParam.PhonCalc.q_ph_max[i];
      q_ph_min=Pi/D*d->opt.ExecParam.PhonCalc.q_ph_min[i];
#endif// BZ_BOUNDARY
#ifndef BZ_BOUNDARY
    q_ph_max=d->opt.ExecParam.PhonCalc.q_ph_max[i];
    q_ph_min=d->opt.ExecParam.PhonCalc.q_ph_min[i];
#endif// BZ_BOUNDARY
    if(d->opt.ExecParam.PhonCalc.Nq_ph!=1){
      dq=(q_ph_max-q_ph_min)
        /(d->opt.ExecParam.PhonCalc.Nq_ph-1);
    }
    else
      dq=0.0;
    if(dq<1e-10)
      dq=0;
    /*Assignment of qx,qy,qz in the first 3 columns of 
      d->phon.Eph[Nq_ph][N_branches+3]*/
    real t=d->opt.ExecParam.PhonCalc.q_ph_min[i];
    for(int j=0;j<d->opt.ExecParam.PhonCalc.Nq_ph;j++){
      d->phon.E_ph[j][i]=t;
      t+=dq;
    }//Assignment of the wave vector values
  }//loop over xyz components
  
  ////////////////////////////////////////////////////////
  //Loop over all requested points in the Brillouin zone//
  ////////////////////////////////////////////////////////
  for(int n_qph=0;n_qph<d->opt.ExecParam.PhonCalc.Nq_ph;n_qph++){
    masterPrint(">>>>>>>>>>> %d : q_ph=( %g, %g, %g )<<<<<<<<<<<<",
                n_qph,
                d->phon.E_ph[n_qph][0],
                d->phon.E_ph[n_qph][1],
                d->phon.E_ph[n_qph][2]);
    //Clean the Dynamical matrix before rebuilding
    d->phon.DM_reallocate(d->phon.N_rows[mpi_n3d_id],
                          (d->phon.N_cols[mpi_n3d_id]<d->phon.states_in_row ?
                           d->phon.N_cols[mpi_n3d_id] : d->phon.states_in_row)//Max number of elements
                                                                          //per row on
                                                                          //mpi_n3d_id proc
                          );

    //////////////////////////////////////////////////////////////
    //Build valence-field-force elements of the Dynamical Matrix//
    //////////////////////////////////////////////////////////////

    //Main range
    DM_Keating(d,n_qph,cell0,cellF);
    //Lower block on the first processor at periodic X b.c.
    if((!mpi_n3d_id)&&d->geo.StrainIsPeriodicX() && (cellF<d->geo.N_Cell-1)){
      DM_Keating(d,n_qph,d->geo.N_Cell-YZslab,d->geo.N_Cell-1);
    }
    //Upper block on the last processor at periodic X b.c.
    if((mpi_n3d_id==mpi_n3d_numprocs-1)&&d->geo.StrainIsPeriodicX() && (cell0!=0)){
      DM_Keating(d,n_qph,0,YZslab-1);
    }
    ///////////////////////////////////////////////////////////
    
#ifdef PIC_DM
    d->phon.Picture_DM((3*d->geo.cellgeom.AtomsPerCellMax()*d->geo.N_Cell));
#endif//PIC_DM

#ifdef PR_DM_MATLAB
    masterPrint("\nDynamical matrix before going to parpack...\n");
    int proc_to_print=0;
    DM_print_matlab(d,proc_to_print);
#endif//PR_DM_MATLAB

    //////////////////////////////////////
    // Dynamical matrix diagonalization //
    //////////////////////////////////////
#ifndef NO_PARPACK
    phon_parpack(d,n_qph);
#else
    phon_arpack(d,n_qph);
#endif /* NO_PARPACK */
    //////////////////////////////////////
    if(d->opt.ExecParam.PhonCalc.PhonOut.PhonEq && !mpi_n3d_id){
      char* filename = n3d_strdup_n(d->inputfile);
      n3d_FileTypeSet(&filename, "PhonDisp", TRUE);
      
      print_phonon_dispersion(filename,d,V_over_V0);
      str_free(filename);
    }
  }//Loop over all requested points in the Brillouin zone
  
  if(d->opt.ExecParam.PhonCalc.PhonOut.PhonEq && !mpi_n3d_id){
    char* filename = n3d_strdup_n(d->inputfile);
    n3d_FileTypeSet(&filename, "PhonDisp", TRUE);

    print_phonon_dispersion(filename,d,V_over_V0);
    str_free(filename);
  }
#ifdef PR_DISP
  if(!mpi_n3d_id){
  real q;
  printf("\n%d branches of phonon dispersion in eV:\n",d->opt.ExecParam.PhonCalc.N_branches);
  printf("\n# a_over_a0  q_x,nm^-1     q_y,nm^-1    q_z,nm^-1  |    q,nm^-1  |");
  for(int n=0;n<d->opt.ExecParam.PhonCalc.N_branches;n++)
    printf("   %d-th branch",n);
  printf("\n");
  for(int n=0;n<d->opt.ExecParam.PhonCalc.Nq_ph;n++){
    q=sqrt(d->phon.E_ph[n][0]*d->phon.E_ph[n][0]+
           d->phon.E_ph[n][1]*d->phon.E_ph[n][1]+
           d->phon.E_ph[n][2]*d->phon.E_ph[n][2]);
    printf("\n%d  %g  %g  %g  %g   %g    ",
           n,V_over_V0,d->phon.E_ph[n][0],d->phon.E_ph[n][1],d->phon.E_ph[n][2],q);

    for(int i=0;i<d->opt.ExecParam.PhonCalc.N_branches;i++)
      printf("   %g",d->phon.E_ph[n][i+3]);
  }
  printf("\n");
  }
#endif/*PR_DISP*/
      
#ifdef PR_nbr_on_proc
  if(mpi_n3d_id==0)
    nbr_on_proc_print(d,cell0,cellF);
#endif//PR_nbr_on_proc

  d->phon.ph_s_clear();

#ifdef TALK_A_LOT
  masterPrint("I am leaving phon_disp now...");
#endif /*TALK_A_LOT*/
      
}//end of void phon_disp
