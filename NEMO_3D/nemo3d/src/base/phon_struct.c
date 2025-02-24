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
$Header: /repo/nemo3d/src/base/phon_struct.c,v 1.3 2003/10/08 16:17:07 hook Exp $ 
*****************************************************************************/


/* Structure to define all computed phonon stuff*/

#include "phon_struct.h"

void phon_struct::ph_s_alloc(int N_uc,//number of unit cells
                             int N_a,//number of atoms per unit cell
			     int N_branches,//number of phonon branches
			     int Nq_ph,//number of points in the reciprocal space
                             int nproc//number of processors
			     )
//Allocates memory for phonon-related matrixes
{
  int nbr_uc=1;/*All the interactions for the standart Keating model with
                 second-neares neighbors interaction limit
                 are within one neighbor unit cell.*/
  this->nbr_uc=nbr_uc;/*Number of neighbor unit cells in x direction, whose
                        atoms interact with the boundary atoms of the given unit cell*/
  int states_in_row=51;//=51 Redefine for Coulomb interaction!!!
  this->states_in_row=states_in_row;/*Maximum value on nonzero elements per row in Keating's model.
				      Change later to defined in the code*/
  this->DMmap=I3matrix(N_uc,N_a,3);//Map of indexies in the dynamical matrix

  this->seg_col_s=Ivectr(nproc);//Array of the initial indexies of the DM on the processor
  this->seg_col_ln=Ivectr(nproc);//Array of the number of elements of the DM on the processor
  this->min_row=Ivectr(nproc);/* minimal index of the rows of DM
                                 to be stored on every processor*/
  this->max_row=Ivectr(nproc);/* maximal index of the rows of DM
                                 to be stored on every processor*/
  this->szCommL=Ivectr(nproc);/*Number of rows in the buffer
                                to be exchanged with the left processor*/
  this->szCommR=Ivectr(nproc);/*Number of rows in the buffer
                                to be exchanged with the right processor*/
  this->N_cols=Ivectr(nproc);/*Number of columns of DM stored on every processor*/
  this->N_rows=Ivectr(nproc);/*Number of rows of DM stored on every processor*/
  this->mpi_left=Ivectr(nproc);/*Left processors mpi_n3d_id  for communication*/
  this->mpi_right=Ivectr(nproc);/*Right processors mpi_n3d_id  for communication*/

  
  this->E_ph=Rmatrix(Nq_ph,N_branches+3);/*[Nq_ph][N_branches+3] matrix of phonon energy 
					  dispersion branches.
					  The first 3 columns contain 
					  q_ph_x,q_ph_y,q_ph_z components*/
}

void phon_struct::comm_vectors_allocate()
     //Allocates memory for parts of the eigenvector that participate in communication
{
    this->wL=Cvectr(szCommL[mpi_n3d_id]);//part of the vector to be transfered to Left processor
    this->wR=Cvectr(szCommR[mpi_n3d_id]);//part of the vector to be transfered to Right processor
    this->vL=Cvectr(szCommR[mpi_left[mpi_n3d_id]]);//part of the vector to be received from Left processor
    this->vR=Cvectr(szCommL[mpi_right[mpi_n3d_id]]);//part of the vector to be received from Right processor

}

void phon_struct::comm_vectors_deallocate()
     //Deallocates memory for parts of the eigenvector that participate in communication
{
  if (this->wL){
    rm_cvectr(&this->wL);
    this->wL=NULL;
  }
  if (this->wR){
    rm_cvectr(&this->wR);
    this->wR=NULL;
  }
  if (this->vL){
    rm_cvectr(&this->vL);
    this->vL=NULL;
  }
  if (this->vR){
    rm_cvectr(&this->vR);
    this->vR=NULL;
  }  
}

void phon_struct::DM_allocate(
                              int N_rows,//total number of rows ON THE PROCESSOR
                              int DM_states_in_row//max number of elements per row on mpi_n3d_id proc
                              )
     //Allocate the Dynamical matrix and buffers for exchange with neighbor processors
{
  this->DM=Cmatrix(N_rows,states_in_row);//dynamical matrix
  this->Ind_DM=Imatrix(N_rows,states_in_row+1);/*indexies of nonzero elements 
							 in m-th row of DM.
							 [i][0] element contains the 
							 number of nonzero elements in 
							 i-th row*/
}

void phon_struct::DM_reallocate(
                                int N_rows,//total number of rows ON THE PROCESSOR
                                int DM_states_in_row//max number of elements per row on mpi_n3d_id proc
                                )
     //Refresh the Dynamical matrix and buffers for exchange with neighbor processors
{
  if (this->DM){
    rm_cmatrix(&this->DM);
    this->DM=Cmatrix(N_rows,DM_states_in_row);
  }
  else
    printf("\n!!!!!!!!!There is an error in DM_reallocate!!!!!!!!!\n");
}

void phon_struct::ph_s_dealloc()
//Deallocates memory for phonon-related matrixes
{
  //Map of indexies
  if (this->DMmap){
    rm_i3matrix(&this->DMmap);
    this->DMmap=NULL;
  }
  //Indexies
  if (this->seg_col_s){
    rm_ivectr(&this->seg_col_s);
    this->seg_col_s=NULL;
  }
  if (this->seg_col_ln){
    rm_ivectr(&this->seg_col_ln);
    this->seg_col_ln=NULL;
  }
  if (this->min_row){
    rm_ivectr(&this->min_row);
    this->min_row=NULL;
  }
  if (this->max_row){
    rm_ivectr(&this->max_row);
    this->max_row=NULL;
  }
  if (this->szCommL){
    rm_ivectr(&this->szCommL);
    this->szCommL=NULL;
  }
  if (this->szCommR){
    rm_ivectr(&this->szCommR);
    this->szCommR=NULL;
  }
  if (this->N_cols){
    rm_ivectr(&this->N_cols);
    this->N_cols=NULL;
  }
  if (this->N_rows){
    rm_ivectr(&this->N_rows);
    this->N_rows=NULL;
  }
  if (this->mpi_left){
    rm_ivectr(&this->mpi_left);
    this->mpi_left=NULL;
  }
  if (this->mpi_right){
    rm_ivectr(&this->mpi_right);
    this->mpi_right=NULL;
  }
  //Dynamical matrix and buffers for communication
  if (this->DM){
    rm_cmatrix(&this->DM);
    this->DM=NULL;
  }
  if (this->Ind_DM){
    rm_imatrix(&this->Ind_DM);
    this->Ind_DM=NULL;
  }
  //phonon energy spectrum
  if (this->E_ph){
    rm_rmatrix(&this->E_ph);
    this->E_ph=NULL;
  }
}

void phon_struct::ph_s_clear()
     //Clears all elements of phon_struct
{
  this->nbr_uc=0;
  this->states_in_row=0;
  this->ph_s_dealloc();
}

void phon_struct::Print_DM()
     //Prints out the Dynamical matrix stored on the processor
{
  printf("\n");
  for(int i=0; i<N_rows[mpi_n3d_id]; i++){
    for(int j=0;j<Ind_DM[i][0];j++){
      printf("mpi_n3d_id=%d  i=%d <=> row=i+min_row=%d; j=%d <=> Ind_DM=%d <=> col=Ind_DM+seg_col_s=%d : %g+i%g\n",
             mpi_n3d_id,i,(i+min_row[mpi_n3d_id]),j,Ind_DM[i][j+1],(Ind_DM[i][j+1]+seg_col_s[mpi_n3d_id]),
             DM[i][j].r,DM[i][j].i);
      fflush(stdout);
    }
  }
}

void phon_struct::Print_parameters()
     //Prints all the parameters on all processors
{
  printf("\n");
  fflush(stdout);
  printf("mpi_n3d_id=%d : pointer to DM -- %p; to Ind_DM --%p;\n"
         "      to DMmap -- %p; to E_ph -- %p\n",
         mpi_n3d_id,DM,Ind_DM,DMmap,E_ph);
  fflush(stdout);
  printf("mpi_n3d_id=%d : N_cols=%d\n     Index of the initial column stored on the proc seg_col_s=%d\n"
         "                                            # of nonzero elements per row seg_col_ln=%d\n",
         mpi_n3d_id,N_cols[mpi_n3d_id],seg_col_s[mpi_n3d_id],seg_col_ln[mpi_n3d_id]);
  fflush(stdout);
  printf("mpi_n3d_id=%d : N_rows=%d\n     Index of the initial row stored on the proc min_row=%d\n"
         "                               Index of the last row stored on the proc max_row=%d\n",
         mpi_n3d_id,N_rows[mpi_n3d_id],min_row[mpi_n3d_id],max_row[mpi_n3d_id]);
  fflush(stdout);
  printf("mpi_n3d_id=%d : Communicates with mpi_left=%d and mpi_right=%d processors\n"
         "            It sends them szCommL=%d and szCommR=%d elements of the eigenvector\n",
         mpi_n3d_id,mpi_left[mpi_n3d_id],mpi_right[mpi_n3d_id],szCommL[mpi_n3d_id],szCommR[mpi_n3d_id]);
  fflush(stdout);
  
}

void phon_struct::Picture_DM(int DM_total_size)
     /*Picture the Dynamical matrix stored on the processor
       "0" - zero, "X" -- nonzero element,
       every element is followed by (mpi_n3d_id) and
       pre... by <row,col> in the whole DM*/
{
  ostringstream DM_ostr;
  string el_str;
  bool flag;
  printf("\n");
  DM_ostr << "mpi_n3d_id=" << mpi_n3d_id<<"\n";
  DM_ostr << "     ";
  for(int col=0;col<seg_col_ln[mpi_n3d_id];col++){
    if((col+seg_col_s[mpi_n3d_id])%10){
      if(((col+seg_col_s[mpi_n3d_id])%10)!=5 &&((col+seg_col_s[mpi_n3d_id])%10)!=3)
        DM_ostr << " ";
      else if(((col+seg_col_s[mpi_n3d_id])%10)==3 && (col+seg_col_s[mpi_n3d_id])<=100)
        DM_ostr << " ";
      else if((col+seg_col_s[mpi_n3d_id])==5)
        DM_ostr << " ";
    }
    else
      DM_ostr << col+seg_col_s[mpi_n3d_id];
  }
  DM_ostr << "\n";
  DM_ostr << "  ";
  for(int col=0;col<seg_col_ln[mpi_n3d_id];col++){
    if((col+seg_col_s[mpi_n3d_id])%10)
      DM_ostr << " ";
    else
      DM_ostr << "|";
  }
  for(int row=0; row<N_rows[mpi_n3d_id]; row++){
    if(row<(max_row[mpi_n3d_id]-min_row[mpi_n3d_id]+1)){
      DM_ostr<<"\n"<< setw(3) << row+min_row[mpi_n3d_id]<<" ";
    }
    else if(mpi_n3d_numprocs!=1 &&!mpi_n3d_id){
      DM_ostr<<"\n"<< setw(3)
             << DM_total_size-(N_rows[mpi_n3d_id]-row) <<" ";
    }
    else if(mpi_n3d_numprocs!=1 && mpi_n3d_id==mpi_n3d_numprocs-1){
      DM_ostr<<"\n"<< setw(3) << row+min_row[mpi_n3d_id]-max_row[mpi_n3d_id]-1<<" ";
    }
    for(int col=0;col<seg_col_ln[mpi_n3d_id];col++){
      flag=true;
      if(Ind_DM[row][0]){
        for(int j=1;j<=Ind_DM[row][0];j++){
          if(Ind_DM[row][j]==col){
#ifndef PR_VALUES
            DM_ostr<<"X";
#endif// PR_VALUES
#ifdef PR_VALUES
            DM_ostr<< (int)DM[row][j-1].r;
            if((int)DM[row][j-1].r<10){
              if(DM[row][j-1].r<0)
                DM_ostr<< " ";
              else
                DM_ostr<< "  ";
            }
            else if((int)DM[row][j-1].r<100){
              if(DM[row][j-1].r>0){
                DM_ostr<< " ";
              }
            }
#endif// PR_VALUES
            flag=false;
            break;
          }
          /*else{
            el_str=".";
            //  DM_ostr<<".";
            }*/
        }//check if equals to 0 or not
        if(flag)
         DM_ostr<<".";
      }//if there are nonzero elements in the row
      else{
        //        el_str=".";
        DM_ostr<<".";
      }
      // DM_ostr<<el_str;
    }//col
  }//row
  DM_ostr<<"\n";
  printf("%s",DM_ostr.str().c_str());fflush(stdout);
}
