/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D
PostProcessing package.
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
  Seungwon Lee (Seungwon.Lee@jpl.nasa.gov)

Written by:  Seungwon Lee
*****************************************************************************/


#ifndef __dipole_h_
#define __dipole_h_
#include <vector>
#include <fstream>
#include <complex>
using namespace std;



class Dipole{

public:
  typedef complex<double> Complex_t;
  
  Dipole(int e, int h, int na, int g1, int g2, int o)
  : e_max(e), h_max(h), numAtoms(na), g_min(g1), g_max(g2), o_max(o)
  {
    o_total = o_max*o_max;
    _H_dipole_X = new complex<double>[e_max*h_max];
    _H_dipole_Y = new complex<double>[e_max*h_max];
    _H_dipole_Z = new complex<double>[e_max*h_max];
  }
  
  ~Dipole()
  {
    delete [] _H_dipole_X;
    delete [] _H_dipole_Y;
    delete [] _H_dipole_Z;
  }

  void assign_table(ifstream&);
  int  pair(int, int);
  void compute_dipole();
  void get_data(Complex_t*, Complex_t*, int*, int*, double*);
  inline complex<double>& H_dipole_X(int e, int h) 
    { return _H_dipole_X[e*h_max+h];}
  inline complex<double>& H_dipole_Y(int e, int h) 
    { return _H_dipole_Y[e*h_max+h];}
  inline complex<double>& H_dipole_Z(int e, int h) 
    { return _H_dipole_Z[e*h_max+h];}
  inline Complex_t& cf_e(int e, int g, int s, int o) 
    { return _cf_e[((e*numAtoms+g)*2+s)*o_max+o];}
  inline Complex_t& cf_h(int h, int g, int s, int o)
    { return _cf_h[((h*numAtoms+g)*2+s)*o_max+o];}
  inline int& nn(int iat, int nadj)
    { return _nn[iat*4+nadj];}
  inline double& lattice(int iat, int dim)
    { return _lattice[iat*3+dim];}	  
  inline complex<double>* get_H_dipole_X(void)
    { return _H_dipole_X;}
  inline complex<double>* get_H_dipole_Y(void)
    { return _H_dipole_Y;}
  inline complex<double>* get_H_dipole_Z(void)
    { return _H_dipole_Z;}
  inline double& dipole_on(int ob, int at)
    { return _dipole_on[ob*num_atom_type+at];}   
  inline double& dipole_off(int ob1, int ob2, int pt)
    { return _dipole_off[(ob1*o_max+ob2)*num_atom_pair+pt];}

private:
 int e_max, h_max, numAtoms, g_min, g_max, o_max, o_total;
 int num_atom_type, num_atom_pair;
 complex<double>* _H_dipole_X;
 complex<double>* _H_dipole_Y;
 complex<double>* _H_dipole_Z;
 Complex_t* _cf_e;
 Complex_t* _cf_h;
 int* _atom;
 int* _nn;
 double* _lattice;
 vector<int>    atom_pair;
 vector<double> _dipole_on;
 vector<double> _dipole_off;

};  

#endif
