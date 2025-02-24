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
$Header: /repo/nemo3d/src/postprocessing/Coulomb.h,v 1.1 2004/08/26 21:28:07 swlee Exp $
*****************************************************************************/

#ifndef __coulomb_h_
#define __coulomb_h_
#include <fstream>
#include <vector>
#include "qd_struct.h"
#include "io_utils.h"
#include "cvector_lib.h"

class Coulomb{

public:
 static const double CONSTANT1, CONSTANT2;
	
 Coulomb(qd_struct din, int e, int h, int mpiID)
  : d(din), e_max(e), h_max(h), my_id(mpiID) 
  {
    g_min = d->cell_s[my_id]*8;
    g_max = d->cell_ln[my_id]*8+g_min;
    o_max = d->NBasisStates/2;
    o_total = o_max*o_max;
    numAtoms = d->n_atom_tot;
    _H_coul = Cvectr(e_max*h_max);
    _H_coul = Cvectr(e_max*h_max);

   MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);
   if(mpi_numprocs==1) g_max = numAtoms;

  }

  ~Coulomb() 
  { rm_cvectr(&_H_coul); rm_cvectr(&_H_exch);}  

  inline complex &cf_e(int n,int g, int s, int o)
  { return _cf_e[((n*numAtoms+g)*2+s)*o_max+o];}

  inline complex &cf_h(int n, int g, int s, int o)
  { return _cf_h[((n*numAtoms+g)*2+s)*o_max+o];}

  inline complex &H_coul(int n1,int n2)
  { return _H_coul[n1*h_max+n2];}

  inline complex &H_exch(int n1, int n2)
  { return _H_exch[n1*h_max+n2];}

  inline int &nn(int g1,int nadj)
  { return _nn[g1*4+nadj];}
  
  inline double &lattice(int g1, int dim)
  { return _lattice[g1*3+dim];}	  
  
  inline complex* get_H_coul(void)
    { return _H_coul;}
  inline complex* get_H_exch(void)
    { return _H_exch;}

  inline double &coul_on(int o1, int o2, int at)
    { return _coul_on[(o1*o_max+o2)*num_atom_type+at];}

  inline double &coul_on_avg(int o1, int o2)
    { return _coul_on[o1*o_max+o2];}

  inline double &exch_on(int o1, int o2, int at)
    { return _exch_on[(o1*o_max+o2)*num_atom_type+at];}

  inline double &exch_off(int o1, int o2, int pt)
    { return _exch_off[(o1*o_max+o2)*num_atom_pair+pt];}

  inline complex multiply(complex A, complex B)//conj(A)*B
    {
      complex result;
      result.r=A.r*B.r+A.i*B.i;
      result.i=A.r*B.i-A.i*B.r;
      return result;
    }

  void assign_table(ifstream&);
  double w_coul(double,int,int,int);
  double w_exch(int,int,int,int);
  double distance_sq(int, int);
  void compute_coulomb_exchange();
  void compute_exchange_vs_dot_separation(int,int);
  void compute_overlap_vs_dot_separation(int,int);
  complex compute_distant_exchange_interaction(int,int);
  void get_data(complex*, complex*, int*, int*, double*);
  double coulomb_cell(int, int);
  int  pair(int, int);
  void print_phase(void);
  int find_cell_overlap(int, int, int);
  bool AreCellsClose(int, int);
  
private:
  qd_struct d;
  int e_max, h_max, g_min, g_max, o_max, o_total;
  int numAtoms, num_atom_type, num_atom_pair, my_id;
  double dielect, inv_dielect;
  int* _atom;
  int* _nn;
  complex* _cf_e;
  complex* _cf_h;
  double* _lattice;
  double* _pair;
  vector<int>    atom_pair;
  vector<double> _coul_on;
  vector<double> _exch_on;
  vector<double> _exch_off;
  vector<double> _coul_on_avg;
  complex* _H_coul;
  complex* _H_exch;
  int mpi_numprocs;
};  

#endif
  
