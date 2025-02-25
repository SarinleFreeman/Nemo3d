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
$Header: /repo/nemo3d/src/base/phon_output.c,v 1.6 2005/01/14 20:24:39 marek Exp
$
*****************************************************************************/

/* Print-outs for phon_struct*/

#include "phon_output.h"

void print_phonon_dispersion(
    ostream &output_file, qd_struct d,
    real a_over_a0 // Measure of strain (in x-direction)
    )
// Prints phonon dispersion in file_name
{
  real q;

  output_file << "\n%" << d->opt.ExecParam.PhonCalc.Nq_ph << " points of "
              << d->opt.ExecParam.PhonCalc.N_branches
              << " branches of phonon dispersion in meV:\n"
              << "%(If line is filled with zeros, it means, that it has not "
                 "been computed yet.)"
              << endl;
  output_file
      << "\n% a_x/a0   q_x,nm^-1     q_y,nm^-1    q_z,nm^-1  |    q,nm^-1  |";
  for (int n = 0; n < d->opt.ExecParam.PhonCalc.N_branches; n++)
    output_file << "   " << n << "-br";
  output_file << "\n";

  for (int n = 0; n < d->opt.ExecParam.PhonCalc.Nq_ph; n++) {
    q = sqrt(d->phon.E_ph[n][0] * d->phon.E_ph[n][0] +
             d->phon.E_ph[n][1] * d->phon.E_ph[n][1] +
             d->phon.E_ph[n][2] * d->phon.E_ph[n][2]);
    output_file << "\n " << a_over_a0 << "  " << d->phon.E_ph[n][0] << "  "
                << d->phon.E_ph[n][1] << "  " << d->phon.E_ph[n][2] << "  "
                << q;
    for (int i = 0; i < d->opt.ExecParam.PhonCalc.N_branches; i++)
      output_file << "  " << 1000 * d->phon.E_ph[n][i + 3];
  }
  output_file << "\n";
}

/**
 * Print the phonon energy  into the file run_name.PhonDisp
 * @param file_name is the name of the file (run_name.PhonDisp)
 * @param d is the qd_struct structure
 * @param a_over_a0 is the ratio between strained and relaxed lattice constant
 * in x-direction
 * @return the file ofthe phonon energy run_name.PhonDisp
 */
void print_phonon_dispersion(
    string file_name, qd_struct d,
    real a_over_a0 // Measure of strain (in x-direction)
) {
  if (file_name != "") {
    ofstream output_file(file_name.c_str());
    print_phonon_dispersion(output_file, d, a_over_a0 /*Delete me!!!*/);
    output_file.close();
  } else
    print_phonon_dispersion(cout, d, a_over_a0);
}

void DM_print_matlab(qd_struct d, int proc_to_print)
// Prints DM for Matlab
{
  real element_r, element_i;
  if (mpi_n3d_id == proc_to_print) {
    printf(
        "\nReal part of DM for matlab\n"
        "Number of rows on %d-th processor=N_rows[mpi_n3d_id]=%d and columns "
        "d->phon.N_cols[mpi_n3d_id]=%d, d->phon.seg_col_ln[mpi_n3d_id]=%d\n",
        mpi_n3d_id, d->phon.N_rows[mpi_n3d_id], d->phon.N_cols[mpi_n3d_id],
        d->phon.seg_col_ln[mpi_n3d_id]);
    for (int row = 0; row < d->phon.N_rows[mpi_n3d_id]; row++) {
      if (row <
          (d->phon.max_row[mpi_n3d_id] - d->phon.min_row[mpi_n3d_id] + 1)) {
        cout << "\n"
             << mpi_n3d_id << " " << setw(3)
             << row + d->phon.min_row[mpi_n3d_id] << " ";
      } else if (mpi_n3d_numprocs != 1 && !mpi_n3d_id) {
        cout << "\n"
             << mpi_n3d_id << " " << setw(3)
             << ((3 * d->geo.cellgeom.AtomsPerCellMax() * d->geo.N_Cell) -
                 (d->phon.N_rows[mpi_n3d_id] - row))
             << " ";
      } else if (mpi_n3d_numprocs != 1 && mpi_n3d_id == mpi_n3d_numprocs - 1) {
        cout << "\n"
             << mpi_n3d_id << " " << setw(3)
             << row + d->phon.min_row[mpi_n3d_id] -
                    d->phon.max_row[mpi_n3d_id] - 1
             << " ";
      }
      for (int col = 0; col < d->phon.seg_col_ln[mpi_n3d_id]; col++) {
        element_r = 0.0;
        element_i = 0.0;
        for (int j = 1; j <= d->phon.Ind_DM[row][0]; j++) {
          if (d->phon.Ind_DM[row][j] == col) {
            element_r = d->phon.DM[row][j - 1].r;
            element_i = d->phon.DM[row][j - 1].i;
            break;
          } else {
            element_r = 0.0;
            element_i = 0.0;
          }
        }
        printf(" %g", element_r);
      }
      printf("\n");
    }
  }
} // end of void DM_print_matlab

void Map_print(qd_struct d) {
  for (int cell = 0; cell < d->geo.N_Cell; cell++) {
    printf("\nmpi_n3d_id=%d : %d-th unit cell:\n", mpi_n3d_id, cell);
    fflush(stdout);
    for (int atom = 0; atom < d->geo.cellgeom.AtomsPerCellMax(); atom++) {
      for (int proj = 0; proj < 3; proj++) {
        printf(" d->phon.DMmap[%d][%d][%d]=%d d->geo.AtomType[%d][%d]=%d\n",
               cell, atom, proj, d->phon.DMmap[cell][atom][proj], cell, atom,
               (int)d->geo.AtomType[cell][atom]);
        fflush(stdout);
      } // proj
    } // atom
  } // cell
}

void nbr_on_proc_print(qd_struct d, int cell0, int cellF)
//
{
  int l, m, n, nbr2, cindx, aindx, cnbr2, anbr2;
  bool has_periodicity = d->geo.StrainHasPeriodicity();
  for (l = cell0; l <= cellF; l++) {
    printf("\nmpi_n3d_id=%d : %d-th unit cell:\n", mpi_n3d_id, l);
    fflush(stdout);
    int i = d->geo.cell__ijk[l][0];
    int j = d->geo.cell__ijk[l][1];
    int k = d->geo.cell__ijk[l][2];
    for (m = 0; m < d->geo.AtomsPerCellMax(); m++) {
      for (n = 0; n < d->geo.Neighbors(m); n++) {
        /* Get relative position of unit cell where this neighbor lives. */
        int i_nbr = i + d->geo.NbrCell(m, n, 0);
        int j_nbr = j + d->geo.NbrCell(m, n, 1);
        int k_nbr = k + d->geo.NbrCell(m, n, 2);
        /* Get the unit cell and atomic index for this neighbor */
        cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

#ifdef ELIMINATE_SSMAP
        if (cindx < 0 && has_periodicity && (cindx + d->geo.N_Cell + 1) >= 0) {
          cindx = cindx + d->geo.N_Cell + 1;
        }
#else  /* ELIMINATE_SSMAP */
        if (cindx < 0 && has_periodicity && d->geo.ssmap[l][m][n] != -1) {
          cindx = d->geo.ssmap[l][m][n];
        }
#endif /* ELIMINATE_SSMAP */

        aindx = d->geo.NbrCell(m, n, 3);
        printf("\n l=%d m=%d    n=%d => cindx =%d aindx=%d", l, m, n, cindx,
               aindx);
        fflush(stdout);
        for (nbr2 = n + 1; nbr2 < d->geo.Neighbors(m); nbr2++) {
          cnbr2 = d->geo.ijk__cell[i + d->geo.NbrCell(m, nbr2, 0)]
                                  [j + d->geo.NbrCell(m, nbr2, 1)]
                                  [k + d->geo.NbrCell(m, nbr2, 2)];

#ifdef ELIMINATE_SSMAP
          if (cnbr2 < 0 && has_periodicity &&
              (cnbr2 + d->geo.N_Cell + 1) >= 0) {
            cnbr2 = cnbr2 + d->geo.N_Cell + 1;
          }
#else  /* ELIMINATE_SSMAP */
          if (cnbr2 < 0 && has_periodicity && d->geo.ssmap[l][m][nbr2] != -1) {
            cnbr2 = d->geo.ssmap[l][m][nbr2];
          }
#endif /* ELIMINATE_SSMAP */

          // atom index
          anbr2 = d->geo.NbrCell(m, nbr2, 3);
          printf("\n         nbr2=%d => cnbr2 =%d anbr2=%d", nbr2, cnbr2,
                 anbr2);
          fflush(stdout);
        } // second-nearest neighbors loop
      } // nearest neighbors
    } // atom
  } // cell
} // end of void nbr_on_proc_print
#if 0
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
#endif // PR_VALUES
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
#endif // PR_VALUES
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

#endif // 0
