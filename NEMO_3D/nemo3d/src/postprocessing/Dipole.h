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
$Header: /repo/nemo3d/src/postprocessing/Dipole.h,v 1.1 2004/08/26 21:28:07 swlee Exp $
*****************************************************************************/

#ifndef __dipole_h_
#define __dipole_h_
#include <vector>
#include <fstream>
#include "qd_struct.h"
#include "io_utils.h"
#include "cvector_lib.h"

class Dipole{

public:
  
  Dipole(qd_struct din, int e, int h, int mpiID)
  : d(din), e_max(e), h_max(h), my_id(mpiID) 
  {
    g_min = d->cell_s[my_id]*8;
    g_max = d->cell_ln[my_id]*8+g_min;
    o_max = d->NBasisStates/2;
    o_total = o_max*o_max;
    numAtoms = d->n_atom_tot;   
    _H_dipole_X = Cvectr(e_max*h_max);
    _H_dipole_Y = Cvectr(e_max*h_max);
    _H_dipole_Z = Cvectr(e_max*h_max);

   MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);
   if(mpi_numprocs==1) g_max = numAtoms;

  }
  
  ~Dipole()
  { 
    rm_cvectr(&_H_dipole_X);
    rm_cvectr(&_H_dipole_Y);
    rm_cvectr(&_H_dipole_Z);
  }
  
  void assign_table(ifstream&);
  int  pair(int, int);
  void compute_dipole();
  void compute_dipole_approximation();
  void compute_dipole_with_spin_suppressed();
  void get_data(complex*, complex*, int*, int*, int*, double*);
  void get_data(complex*, complex*, double*);
  complex dot_product(complex, complex);
  double amplitude_e(int, int, int);
  double amplitude_h(int, int, int);

  inline complex& H_dipole_X(int e, int h) 
    { return _H_dipole_X[e*h_max+h];}
  inline complex& H_dipole_Y(int e, int h) 
    { return _H_dipole_Y[e*h_max+h];}
  inline complex& H_dipole_Z(int e, int h) 
    { return _H_dipole_Z[e*h_max+h];}
  inline complex& cf_e(int e, int g, int s, int o) 
    { return _cf_e[((e*numAtoms+g)*2+s)*o_max+o];}
  inline complex& cf_h(int h, int g, int s, int o)
    { return _cf_h[((h*numAtoms+g)*2+s)*o_max+o];}
  inline int& nn(int iat, int nadj)
    { return _nn[iat*4+nadj];}
  inline double& lattice(int iat, int dim)
    { return _lattice[iat*3+dim];}	  
  inline complex* get_H_dipole_X(void)
    { return _H_dipole_X;}
  inline complex* get_H_dipole_Y(void)
    { return _H_dipole_Y;}
  inline complex* get_H_dipole_Z(void)
    { return _H_dipole_Z;}
  inline double& dipole_on(int ob, int at)
    { return _dipole_on[ob*num_atom_type+at];}   
  inline double& dipole_off(int ob1, int ob2, int pt)
    { return _dipole_off[(ob1*o_max+ob2)*num_atom_pair+pt];}

private:
 qd_struct d;
 int e_max, h_max, g_min, g_max, o_max, o_total;
 int numAtoms, num_atom_type, num_atom_pair, my_id;
 complex* _H_dipole_X;
 complex* _H_dipole_Y;
 complex* _H_dipole_Z;
 complex* _cf_e;
 complex* _cf_h;
 int* _atom;
 int* _shape;
 int* _nn;
 double* _lattice;
 vector<int>    atom_pair;
 vector<double> _dipole_on;
 vector<double> _dipole_off;
 int mpi_numprocs;
};  

#endif
