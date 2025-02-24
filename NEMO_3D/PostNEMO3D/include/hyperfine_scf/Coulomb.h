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
#include "nml_dcvector.h"
#include <fstream>
#include <complex>

#include <vector>
using namespace std;

class Coulomb{

public:
	
 typedef complex<double> Complex_t;


 Coulomb(int na, int g1, int g2, int o)
  : numAtoms(na), g_min(g1), g_max(g2), o_max(o)
  {
    
    o_total=o_max*o_max;
    _H_coul = new complex<double>[1];
    _H_exch = new complex<double>[1];
    	
  }

  ~Coulomb() 
  {delete [] _H_coul; delete [] _H_exch;}  

  inline Complex_t &cf_e(int n,int g, int s, int o)
  { return _cf_e[((n*numAtoms+g)*2+s)*o_max+o];}

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

  void assign_table1();
  void assign_table(ifstream&);
  double w_coul(double,int,int,int);
  double w_exch(int,int,int,int);
  double distance_sq(int, int);
  void compute_coulomb_exchange();
  //void compute_exchange(int, int, int, int, nml_dcvector*);	
  void compute_exchange(int, int, int, int, int);
  void get_data(Complex_t*, int*, int*, double*);
  void compute_exchange_cutoff(int, int, int, int, int, double);
  void print_data(void);
  int  pair(int, int);
  void getImpurityInfo(int); 
  void dotProduct(int, int, Complex_t*);
  void compute_HL_exchange(int m, int n, int p, int q);
  void oneDSlice(int, int, int, int, double*, double*); 
  void AHF(double*, double*, double*, double*);
  void AHF_wf(double*, double*);
  void getSCF(double*, double*, double*);
  void getSCF_image(double*, double*, double*);
  void getSCF_CI(double*, double*);
  void AHF_shell(double*, double*, double*);
  void x2(Complex_t* , Complex_t* , Complex_t* , Complex_t* , Complex_t* , Complex_t* , Complex_t*, double, double, double );
  void moment(Complex_t* , Complex_t* , Complex_t* , double xP, double yP, double zP);
  //void dipole_moment(Complex_t* , Complex_t* , Complex_t*);
  void r2(Complex_t* , Complex_t* , Complex_t* , double xP, double yP, double zP);
  void getSCF_new(double*, double*, double*, double*, double*, double*);
  void density_Chandra(double*, double*);
  void HF(int, int*, double*, double*);

private:
  int numAtoms, g_min, g_max, o_max;
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
  vector<double> _coul;
  vector<double> _exch;
};  

#endif
  
