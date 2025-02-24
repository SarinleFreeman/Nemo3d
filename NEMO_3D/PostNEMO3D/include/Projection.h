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


#ifndef __projection_h_
#define __projection_h_
#include <vector>
#include <fstream>
#include <complex>
using namespace std;

class Projection{

public:
  typedef complex<double> Complex_t;
  
  Projection(int e, int h, int na, int g1, int g2, int o)
  : e_max(e), h_max(h), numAtoms(na), g_min(g1), g_max(g2), o_max(o)
  {
    _projection = new double[(e_max+h_max)*8];
    _coef       = new Complex_t[2*o_max];
  }
  
  ~Projection()
  {
    delete [] _projection;
    delete [] _coef;
  }

  void compute_projection();
  void get_data(Complex_t*, Complex_t*, int*, double*);
  inline Complex_t& cf_e(int e, int g, int s, int o) 
    { return _cf_e[((e*numAtoms+g)*2+s)*o_max+o];}
  inline Complex_t& cf_h(int h, int g, int s, int o)
    { return _cf_h[((h*numAtoms+g)*2+s)*o_max+o];}
  inline Complex_t& coef(int s, int o)
    { return _coef[s*o_max+o];}
  inline double& lattice(int iat, int dim)
    { return _lattice[iat*3+dim];}	  
  inline double& projection(int n, int orbital)
    { return _projection[n*8+orbital];}
  inline double* get_projection()
    { return _projection;}

//projection(n,0): total amplitude of s* orbital in envelope functions
//projection(n,1): total amplitude of s  orbital in envelope functions
//projection(n,2): total amplitude of heavy-hole Bloch function (J=3/2, Jz=+/-3/2)
//projection(n,3): total amplitude of light-hole Bloch function (J=3/2, Jz=+/-1/2)
//projection(n,4): total amplitude of spin-split Bloch function (J=1/2, Jz=+/-1/2)
//projection(n,5): total amplitude of d orbital  in envelope functions 

private:
 int e_max, h_max, numAtoms, g_min, g_max, o_max;
 int num_atom_type, num_atom_pair;
 Complex_t* _cf_e;
 Complex_t* _cf_h;
 Complex_t* _coef;
 int* _atom;
 int* _nn;
 double* _lattice;
 double* _projection;
 
};  

#endif
