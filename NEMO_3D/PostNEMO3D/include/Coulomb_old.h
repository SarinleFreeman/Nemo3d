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
#include <fstream>
#include <complex>
#include <vector>
using namespace std;


class Coulomb{

public:
	
 typedef complex<double> Complex_t;

 Coulomb(int e, int h, int na, int g1, int g2, int o)
  : e_max(e), h_max(h), numAtoms(na), g_min(g1), g_max(g2), o_max(o)
  {
    n_total=e_max*h_max;
    o_total=o_max*o_max;
    _H_coul = new complex<double>[n_total];
    _H_exch = new complex<double>[n_total];
  }

  ~Coulomb() 
  {delete [] _H_coul; delete [] _H_exch;}  

  inline Complex_t &cf_e(int n,int g, int s, int o)
  { return _cf_e[((n*numAtoms+g)*2+s)*o_max+o];}

  inline Complex_t &cf_h(int n, int g, int s, int o)
  { return _cf_h[((n*numAtoms+g)*2+s)*o_max+o];}

  inline complex<double> &H_coul(int n1,int n2)
  { return _H_coul[n1*h_max+n2];}

  inline complex<double> &H_exch(int n1, int n2)
  { return _H_exch[n1*h_max+n2];}

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
    { return _coul_on[o1*o_max+o2];}

  inline double &exch_on(int o1, int o2, int at)
    { return _exch_on[(o1*o_max+o2)*num_atom_type+at];}

  inline double &exch_off(int o1, int o2, int pt)
    { return _exch_off[(o1*o_max+o2)*num_atom_pair+pt];}

  void assign_table(ifstream&);
  double w_coul(double,int,int,int);
  double w_exch(int,int,int,int);
  double distance_sq(int, int);
  void compute_coulomb_exchange();
  void get_data(Complex_t*, Complex_t*, int*, int*, double*);
  void print_data(void);
  int  pair(int, int);

private:
  int e_max, h_max, numAtoms, g_min, g_max, o_max;
  int n_total, e_total, h_total, o_total;
  int num_atom_type, num_atom_pair;
  double dielect, inv_dielect;
  int* _atom;
  int* _nn;
  Complex_t* _cf_e;
  Complex_t* _cf_h;
  double* _lattice;
  double* _pair;
  vector<int>    atom_pair;
  vector<double> _coul_on;
  vector<double> _exch_on;
  vector<double> _exch_off;
  vector<double> _coul_on_avg;
  complex<double>* _H_coul;
  complex<double>* _H_exch;

};  

#endif
  
