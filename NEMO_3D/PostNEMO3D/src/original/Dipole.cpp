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


#include "Dipole.h"
#include "SimpleParser.h"


void Dipole::assign_table(ifstream& table)
{

 vector<string> slist;
  int p = 0;
  do {
    int argc = getwords(slist, table);
    if(argc==0) continue;
    if(slist[0] == "num_atom_type") {
       num_atom_type = atoi(slist[1].c_str());
       _dipole_on  = vector<double>(num_atom_type*6);
    } else if(slist[0] == "num_atom_pair") {
       num_atom_pair = atoi(slist[1].c_str());
       atom_pair = vector<int>(2*num_atom_pair);
       _dipole_off = vector<double>(num_atom_pair*o_total);
    } else if(slist[0] == "atom_pair") {
       atom_pair[p*2] = atoi(slist[1].c_str());
       atom_pair[p*2+1] = atoi(slist[2].c_str());
       p++;
    } else if(slist[0] == "dipole_onsite_s_Z_pz"){
      for(int at=0; at<num_atom_type; at++)
        dipole_on(0,at) = atof(slist[at+1].c_str());
    } else if(slist[0] == "dipole_onsite_ss_Z_pz"){
      for(int at=0; at<num_atom_type; at++)
        dipole_on(1,at) = atof(slist[at+1].c_str());
    } else if(slist[0] == "dipole_onsite_px_Z_dzx"){
      for(int at=0; at<num_atom_type; at++)
        dipole_on(2,at) = atof(slist[at+1].c_str());
    } else if(slist[0] == "dipole_onsite_px_X_dxx_yy"){
      for(int at=0; at<num_atom_type; at++)
        dipole_on(3,at) = atof(slist[at+1].c_str());
    } else if(slist[0] == "dipole_onsite_px_X_d3zz_rr"){
      for(int at=0; at<num_atom_type; at++)
        dipole_on(4,at) = atof(slist[at+1].c_str());
    } else if(slist[0] == "dipole_onsite_pz_Z_d3zz_rr"){
      for(int at=0; at<num_atom_type; at++)
        dipole_on(5,at) = atof(slist[at+1].c_str());
    }
  } while (!table.eof());

}  

int Dipole::pair(int g1, int g2)
{
  int pair_index=-1;
  for(int p=0; p<num_atom_pair; p++) {
    int a1 = _atom[g1]; int a2 = _atom[g2];
    if(a1==atom_pair[p*num_atom_pair] && a2==atom_pair[p*num_atom_pair+1]){
       pair_index = p;
       break;
    }
  }
  if(pair_index==-1) {
    cout<<"This atom pair is not registered. Fix file dipole.table."<<endl;
    exit(0);
  }
  return pair_index;
}

void Dipole::get_data(Complex_t* coef_e, Complex_t* coef_h,
                      int* atomid, int* nn, double* lat)
{
   _cf_e = coef_e;
   _cf_h = coef_h;
   _atom = atomid;
   _nn   = nn;
   _lattice = lat;
}

void Dipole::compute_dipole(void)
{
  for(int n1=0;n1<e_max;n1++){
    for(int n2=0;n2<h_max;n2++){
      Complex_t X1=0.0, X2=0.0, X3=0.0, X4=0.0, X5=0.0, X6=0.0; 
      Complex_t Y1=0.0, Y2=0.0, Y3=0.0, Y4=0.0, Y5=0.0, Y6=0.0; 
      Complex_t Z1=0.0, Z2=0.0, Z3=0.0, Z4=0.0, Z5=0.0; 
      for(int g1=g_min;g1<g_max;g1++){
            Complex_t x1=0.0, x2=0.0, x3=0.0, x4=0.0, x5=0.0;
            Complex_t y1=0.0, y2=0.0, y3=0.0, y4=0.0, y5=0.0;
            Complex_t z1=0.0, z2=0.0, z3=0.0, z4=0.0;
            Complex_t overlap=0.0;
            for(int s=0; s<2; s++){

              x1+=cf_e(n1,g1,s,1)*conj(cf_h(n2,g1,s,2));
              x1+=cf_e(n1,g1,s,2)*conj(cf_h(n2,g1,s,1));
              x2+=cf_e(n1,g1,s,0)*conj(cf_h(n2,g1,s,2));
              x2+=cf_e(n1,g1,s,2)*conj(cf_h(n2,g1,s,0));
              x3+=cf_e(n1,g1,s,3)*conj(cf_h(n2,g1,s,5));
              x3+=cf_e(n1,g1,s,5)*conj(cf_h(n2,g1,s,3));
              x3+=cf_e(n1,g1,s,4)*conj(cf_h(n2,g1,s,7));
              x3+=cf_e(n1,g1,s,7)*conj(cf_h(n2,g1,s,4));
              x4+=cf_e(n1,g1,s,2)*conj(cf_h(n2,g1,s,8));
              x4+=cf_e(n1,g1,s,8)*conj(cf_h(n2,g1,s,2));
              x5+=cf_e(n1,g1,s,2)*conj(cf_h(n2,g1,s,9));
              x5+=cf_e(n1,g1,s,9)*conj(cf_h(n2,g1,s,2));

              y1+=cf_e(n1,g1,s,1)*conj(cf_h(n2,g1,s,3));
              y1+=cf_e(n1,g1,s,3)*conj(cf_h(n2,g1,s,1));
              y2+=cf_e(n1,g1,s,0)*conj(cf_h(n2,g1,s,3));
              y2+=cf_e(n1,g1,s,3)*conj(cf_h(n2,g1,s,0));
              y3+=cf_e(n1,g1,s,2)*conj(cf_h(n2,g1,s,5));
              y3+=cf_e(n1,g1,s,5)*conj(cf_h(n2,g1,s,2));
              y3+=cf_e(n1,g1,s,4)*conj(cf_h(n2,g1,s,6));
              y3+=cf_e(n1,g1,s,6)*conj(cf_h(n2,g1,s,4));
              y4+=cf_e(n1,g1,s,3)*conj(cf_h(n2,g1,s,8));
              y4+=cf_e(n1,g1,s,8)*conj(cf_h(n2,g1,s,3));
              y5+=cf_e(n1,g1,s,3)*conj(cf_h(n2,g1,s,9));
              y5+=cf_e(n1,g1,s,9)*conj(cf_h(n2,g1,s,3));

              z1+=cf_e(n1,g1,s,1)*conj(cf_h(n2,g1,s,4));
              z1+=cf_e(n1,g1,s,4)*conj(cf_h(n2,g1,s,1));
              z2+=cf_e(n1,g1,s,0)*conj(cf_h(n2,g1,s,4));
              z2+=cf_e(n1,g1,s,4)*conj(cf_h(n2,g1,s,0));
              z3+=cf_e(n1,g1,s,2)*conj(cf_h(n2,g1,s,7));
              z3+=cf_e(n1,g1,s,7)*conj(cf_h(n2,g1,s,2));
              z3+=cf_e(n1,g1,s,3)*conj(cf_h(n2,g1,s,6));
              z3+=cf_e(n1,g1,s,6)*conj(cf_h(n2,g1,s,3));
              z4+=cf_e(n1,g1,s,4)*conj(cf_h(n2,g1,s,9));
              z4+=cf_e(n1,g1,s,9)*conj(cf_h(n2,g1,s,4));

              for(int orbital=0; orbital<o_max; orbital++)
                overlap+=cf_e(n1,g1,s,orbital)*conj(cf_h(n2,g1,s,orbital));

            }

	    X1+=x1*dipole_on(0, _atom[g1]);
	    X2+=x2*dipole_on(1, _atom[g1]);
	    X3+=x3*dipole_on(2, _atom[g1]);
            X4+=x4*dipole_on(3, _atom[g1]);
            X5+=x5*dipole_on(4, _atom[g1]);
            X6+=overlap*lattice(g1,0);	    

	    Y1+=y1*dipole_on(0, _atom[g1]);
	    Y2+=y2*dipole_on(1, _atom[g1]);
	    Y3+=y3*dipole_on(2, _atom[g1]);
            Y4+=y4*dipole_on(3, _atom[g1]);
            Y5+=y5*dipole_on(4, _atom[g1]);
            Y6+=overlap*lattice(g1,1);	    

	    Z1+=z1*dipole_on(0, _atom[g1]);
	    Z2+=z2*dipole_on(1, _atom[g1]);
	    Z3+=z3*dipole_on(2, _atom[g1]);
            Z4+=z4*dipole_on(5, _atom[g1]);
            Z5+=overlap*lattice(g1,2);	    

      }  
      H_dipole_X(n1,n2) = X1+X2+X3+X4+X5+X6; 
      H_dipole_Y(n1,n2) = Y1+Y2+Y3+Y4+Y5+Y6;
      H_dipole_Z(n1,n2) = Z1+Z2+Z3+Z4+Z5; 
    }
  }
}

