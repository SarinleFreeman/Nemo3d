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


#ifndef __coulomb_h_
#define __coulomb_h_
#define MATRIX_SIZE 21

#include "nml_dcvector.h"
#include <fstream>
#include <complex>

#include <vector>
using namespace std;

class Coulomb{

public:
	
 typedef complex<double> Complex_t;


 Coulomb(int na, int g1, int g2, int g3, int g4, int o)
  : numAtoms(na), g_min1(g1), g_max1(g2), g_min2(g3), g_max2(g4), o_max(o)
  {
    
    o_total=o_max*o_max;
    _H_coul = new complex<double>[1];
    _H_exch = new complex<double>[1];

    for (int i = 1; i < MATRIX_SIZE; i++) {

    	coul_Matrix[i] = 0.0;
        exch_Matrix[i] = 0.0;

    }
    	
  }

  ~Coulomb() 
  {delete [] _H_coul; delete [] _H_exch;}  

  inline Complex_t &cf_e(int n,int g, int s, int o)
  { return _cf_e[((n*(g_max1 - g_min1)+g)*2+s)*o_max+o];}

  //inline Complex_t &cf_h(int n, int g, int s, int o)
  //{ return _cf_h[((n*numAtoms+g)*2+s)*o_max+o];}

  inline complex<double> &H_coul(void)
  { return _H_coul[0];}

  inline complex<double> &H_exch(void)
  { return _H_exch[0];}

  inline int &nn(int g1,int nadj)
  { return _nn[g1*4+nadj];}
  
  inline double &lattice(int g1, int dim)
  { return _lattice[g1*3+dim];}	  
  
  inline complex<double>* get_H_coul(void)
    { return _H_coul;}
  inline complex<double>* get_H_exch(void)
    { return _H_exch;}

  inline double &coul_on(int o1, int o2, int at)
    { return _coul_on[(o1*o_max+o2)*num_atom_type+at];}

  inline double &coul_on_avg(int o1, int o2)
    { return _coul_on_avg[o1*o_max+o2];}

  inline double &coul_on_rec_sq(int o1, int o2)
    { return _coul_on_rec_sq[o1*o_max+o2];}

  inline double &exch_on(int o1, int o2, int at)
    { return _exch_on[(o1*o_max+o2)*num_atom_type+at];}

  inline double &exch_off(int o1, int o2, int pt)
    { return _exch_off[(o1*o_max+o2)*num_atom_pair+pt];}

  inline double &coul(int o1, int o2){
  	return _coul[o1*o_max+o2];
  }

  inline double &exch(int o1, int o2){
        return _exch[o1*o_max+o2];
  }
  
  inline complex<double>* get_coul_Matrix(void) {
        return &coul_Matrix[0];
  }

  inline complex<double>* get_exch_Matrix(void) {
        return &exch_Matrix[0];
  }

  void assign_table1();
  void assign_table(ifstream&);
  double w_coul(double,int,int,int);
  double w_exch(int,int,int,int);
  double distance_sq(int, int);
  void compute_coulomb_exchange2(double, const char*, const int*, const int);
  void compute_coulomb_exchange4(double, const char*, const int*, const int);
  void hartree_matrix(double, const char*, const int*, int);
  void compute_one_component(double, complex<double>*, complex<double>*, complex<double>*, complex<double>*, complex<double>*, int, int, int, int, int); 
  //void compute_exchange(int, int, int, int, nml_dcvector*);	
  void compute_exchange(int, int, int, int, int);
  void get_data(Complex_t*, int*, int*, double*);
  void compute_exchange_cutoff(int, int, int, int, int, double);
  void print_data(void);
  int  pair(int, int);
  void getImpurityInfo(int); 
  void dotProduct(int, int, Complex_t*);
  void compute_HL_exchange(int m, int n, int p, int q);
  void spinPolarization(int, Complex_t*, Complex_t*);
  void dotProduct_Slater(Complex_t*, Complex_t*, Complex_t*, Complex_t*);
  void coulomb_Slater(int, int, Complex_t*, Complex_t*, Complex_t*, Complex_t*);
  void compute_coulomb_exchange4_no_dielectric_correction(double, const char*, const int*, const int);
	 
private:
  int numAtoms, g_min1, g_max1, g_min2, g_max2, o_max;
  int o_total;
  int num_atom_type, num_atom_pair;
  double dielect, inv_dielect;
  int* _atom;
  int* _nn;
  Complex_t* _cf_e;
  //Complex_t* _cf_h;
  double* _lattice;
  double* _pair;
  vector<int>    atom_pair;
  vector<double> _coul_on;
  vector<double> _exch_on;
  vector<double> _exch_off;
  vector<double> _coul_on_avg;
  vector<double> _coul_on_rec_sq;
  complex<double>* _H_coul;
  complex<double>* _H_exch;
  complex<double> coul_Matrix[MATRIX_SIZE];
  complex<double> exch_Matrix[MATRIX_SIZE];
  vector<double> _coul;
  vector<double> _exch;
};  

#endif
  
