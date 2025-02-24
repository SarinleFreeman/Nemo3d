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
$Header: /repo/nemo3d/src/postprocessing/Dipole.cpp,v 1.1 2004/08/26 21:28:07 swlee Exp $
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

void Dipole::get_data(complex* coef_e, complex* coef_h,
                      int* atomid, int* shapeid, int* nn, double* lat)
{
   _cf_e = coef_e;
   _cf_h = coef_h;
   _atom = atomid;
   _shape = shapeid;
   _nn   = nn;
   _lattice = lat;
}

void Dipole::get_data(complex* coef_e, complex* coef_h, double* lat)
{
   _cf_e = coef_e;
   _cf_h = coef_h;
   _lattice = lat;
}

complex Dipole::dot_product(complex A, complex B){
    complex overlap;
    overlap.r=A.r*B.r+A.i*B.i;
    overlap.i=A.r*B.i-A.i*B.r;
    return overlap;
}

double Dipole::amplitude_e(int n, int g, int o) {
    complex UpSpin = cf_e(n, g, 0, o);
    complex DownSpin = cf_e(n, g, 1, o);
    double amp;
    amp  = UpSpin.r*UpSpin.r+UpSpin.i*UpSpin.i;
    amp += DownSpin.r*DownSpin.r+DownSpin.i*DownSpin.i;
    return sqrt(amp);
}

double Dipole::amplitude_h(int n, int g, int o) {
    complex UpSpin = cf_h(n, g, 0, o);
    complex DownSpin = cf_h(n, g, 1, o);
    double amp;
    amp  = UpSpin.r*UpSpin.r+UpSpin.i*UpSpin.i;
    amp += DownSpin.r*DownSpin.r+DownSpin.i*DownSpin.i;
    return sqrt(amp);
}

void Dipole::compute_dipole_approximation(void)
{

  for(int n1=0;n1<e_max;n1++){ //electron state
    for(int n2=0;n2<h_max;n2++){ //hole state
      complex X, Y, Z;
      X.r=0.0; X.i=0.0; Y.r=0.0; Y.i=0.0; Z.r=0.0; Z.i=0.0;

      for(int atom=g_min;atom<g_max;atom++){ //atom index
         complex overlap; overlap.r=0.0; overlap.i=0.0;
         for(int spin=0; spin<2; spin++){ //spin index
            for(int orbital=0; orbital<o_max; orbital++){ //orbital index
                complex tmp = dot_product(cf_e(n1,atom,spin ,orbital), cf_h(n2,atom, spin, orbital));
                overlap.r+=tmp.r; overlap.i+=tmp.i; //overlap between two wave functions
            }
         }

         X.r+=overlap.r*lattice(atom,0); // x coordinate of atomic position
         X.i+=overlap.i*lattice(atom,0); // x coordinate of atomic position
         Y.r+=overlap.r*lattice(atom,1); // y coordinate of atomic position
         Y.i+=overlap.i*lattice(atom,1); // y coordinate of atomic position
         Z.r+=overlap.r*lattice(atom,2); // z coordinate of atomic position
         Z.i+=overlap.i*lattice(atom,2); // z coordinate of atomic position

      } 
      cout_master<<"n1, n2, <n1|X|n2>: "<<n1<<" "<<n2<<" "<<"("<<X.r<<", "<<X.i<<")"<<endl;
      cout_master<<"n1, n2, <n1|Y|n2>: "<<n1<<" "<<n2<<" "<<"("<<Y.r<<", "<<Y.i<<")"<<endl;
      cout_master<<"n1, n2, <n1|Z|n2>: "<<n1<<" "<<n2<<" "<<"("<<Y.r<<", "<<Y.i<<")"<<endl;
    }
  }

}

void Dipole::compute_dipole(void)
{
  for(int n1=0;n1<e_max;n1++){
    for(int n2=0;n2<h_max;n2++){
      complex X1, X2, X3, X4, X5, X6; 
      complex Y1, Y2, Y3, Y4, Y5, Y6; 
      complex Z1, Z2, Z3, Z4, Z5; 
      X1.r=0.0; X1.i=0.0; 
      X2.r=0.0; X2.i=0.0; 
      X3.r=0.0; X3.i=0.0; 
      X4.r=0.0; X4.i=0.0; 
      X5.r=0.0; X5.i=0.0; 
      X6.r=0.0; X6.i=0.0; 
      Y1.r=0.0; Y1.i=0.0; 
      Y2.r=0.0; Y2.i=0.0; 
      Y3.r=0.0; Y3.i=0.0; 
      Y4.r=0.0; Y4.i=0.0; 
      Y5.r=0.0; Y5.i=0.0; 
      Y6.r=0.0; Y6.i=0.0; 
      Z1.r=0.0; Z1.i=0.0; 
      Z2.r=0.0; Z2.i=0.0; 
      Z3.r=0.0; Z3.i=0.0; 
      Z4.r=0.0; Z4.i=0.0; 
      Z5.r=0.0; Z5.i=0.0; 
      for(int g1=g_min;g1<g_max;g1++){
            complex x1, x2, x3, x4, x5;
            complex y1, y2, y3, y4, y5;
            complex z1, z2, z3, z4;
            complex overlap;
            x1.r=0.0; x1.i=0.0; 
            x2.r=0.0; x2.i=0.0; 
            x3.r=0.0; x3.i=0.0; 
            x4.r=0.0; x4.i=0.0; 
            x5.r=0.0; x5.i=0.0; 
            y1.r=0.0; y1.i=0.0; 
            y2.r=0.0; y2.i=0.0; 
            y3.r=0.0; y3.i=0.0; 
            y4.r=0.0; y4.i=0.0; 
            y5.r=0.0; y5.i=0.0; 
            z1.r=0.0; z1.i=0.0; 
            z2.r=0.0; z2.i=0.0; 
            z3.r=0.0; z3.i=0.0; 
            z4.r=0.0; z4.i=0.0; 
            overlap.r=0.0; overlap.i=0.0;

            for(int s=0; s<2; s++){
              for(int orbital=0; orbital<o_max; orbital++){
                complex tmp = dot_product(cf_e(n1,g1,s,orbital), cf_h(n2,g1,s,orbital));  
                overlap.r+=tmp.r; overlap.i+=tmp.i;
              }  

              complex tmp1, tmp2, tmp3, tmp4;

              tmp1 = dot_product(cf_e(n1,g1,s,1), cf_h(n2,g1,s,2));
              tmp2 = dot_product(cf_e(n1,g1,s,2), cf_h(n2,g1,s,1));
              x1.r+=tmp1.r+tmp2.r; x1.i+=tmp1.i+tmp2.i;
              
              tmp1 = dot_product(cf_e(n1,g1,s,0), cf_h(n2,g1,s,2));
              tmp2 = dot_product(cf_e(n1,g1,s,2), cf_h(n2,g1,s,0));
              x2.r+=tmp1.r+tmp2.r; x2.i+=tmp1.i+tmp2.i;
              
              tmp1 = dot_product(cf_e(n1,g1,s,3), cf_h(n2,g1,s,5));
              tmp2 = dot_product(cf_e(n1,g1,s,5), cf_h(n2,g1,s,3));
              tmp3 = dot_product(cf_e(n1,g1,s,4), cf_h(n2,g1,s,7));
              tmp4 = dot_product(cf_e(n1,g1,s,7), cf_h(n2,g1,s,4));
              x3.r+=tmp1.r+tmp2.r+tmp3.r+tmp4.r; 
              x3.i+=tmp1.i+tmp2.i+tmp3.i+tmp4.i; 
              
              tmp1 = dot_product(cf_e(n1,g1,s,2), cf_h(n2,g1,s,8));
              tmp2 = dot_product(cf_e(n1,g1,s,8), cf_h(n2,g1,s,2));
              x4.r+=tmp1.r+tmp2.r; x4.i+=tmp1.i+tmp2.i;
              
              tmp1 = dot_product(cf_e(n1,g1,s,2), cf_h(n2,g1,s,9));
              tmp2 = dot_product(cf_e(n1,g1,s,9), cf_h(n2,g1,s,2));
              x5.r+=tmp1.r+tmp2.r; x5.i+=tmp1.i+tmp2.i;
              
              tmp1 = dot_product(cf_e(n1,g1,s,1), cf_h(n2,g1,s,3));
              tmp2 = dot_product(cf_e(n1,g1,s,3), cf_h(n2,g1,s,1));
              y1.r+=tmp1.r+tmp2.r; y1.i+=tmp1.i+tmp2.i;
              
              tmp1 = dot_product(cf_e(n1,g1,s,0), cf_h(n2,g1,s,3));
              tmp2 = dot_product(cf_e(n1,g1,s,3), cf_h(n2,g1,s,0));
              y2.r+=tmp1.r+tmp2.r; y2.i+=tmp1.i+tmp2.i;
              
              tmp1 = dot_product(cf_e(n1,g1,s,2), cf_h(n2,g1,s,5));
              tmp2 = dot_product(cf_e(n1,g1,s,5), cf_h(n2,g1,s,2));
              tmp3 = dot_product(cf_e(n1,g1,s,4), cf_h(n2,g1,s,6));
              tmp4 = dot_product(cf_e(n1,g1,s,6), cf_h(n2,g1,s,4));
              y3.r+=tmp1.r+tmp2.r+tmp3.r+tmp4.r; 
              y3.i+=tmp1.i+tmp2.i+tmp3.i+tmp4.i; 
              
              tmp1 = dot_product(cf_e(n1,g1,s,3), cf_h(n2,g1,s,8));
              tmp2 = dot_product(cf_e(n1,g1,s,8), cf_h(n2,g1,s,3));
              y4.r+=tmp1.r+tmp2.r; y4.i+=tmp1.i+tmp2.i;
              
              tmp1 = dot_product(cf_e(n1,g1,s,3), cf_h(n2,g1,s,9));
              tmp2 = dot_product(cf_e(n1,g1,s,9), cf_h(n2,g1,s,3));
              y5.r+=tmp1.r+tmp2.r; y5.i+=tmp1.i+tmp2.i;
              
              tmp1 = dot_product(cf_e(n1,g1,s,1), cf_h(n2,g1,s,4));
              tmp2 = dot_product(cf_e(n1,g1,s,4), cf_h(n2,g1,s,1));
              z1.r+=tmp1.r+tmp2.r; z1.i+=tmp1.i+tmp2.i;
              
              tmp1 = dot_product(cf_e(n1,g1,s,0), cf_h(n2,g1,s,4));
              tmp2 = dot_product(cf_e(n1,g1,s,4), cf_h(n2,g1,s,0));
              z2.r+=tmp1.r+tmp2.r; z2.i+=tmp1.i+tmp2.i;
              
              tmp1 = dot_product(cf_e(n1,g1,s,2), cf_h(n2,g1,s,7));
              tmp2 = dot_product(cf_e(n1,g1,s,7), cf_h(n2,g1,s,2));
              tmp3 = dot_product(cf_e(n1,g1,s,3), cf_h(n2,g1,s,6));
              tmp4 = dot_product(cf_e(n1,g1,s,6), cf_h(n2,g1,s,3));
              z3.r+=tmp1.r+tmp2.r+tmp3.r+tmp4.r; 
              z3.i+=tmp1.i+tmp2.i+tmp3.i+tmp4.i; 
              
              tmp1 = dot_product(cf_e(n1,g1,s,4), cf_h(n2,g1,s,9));
              tmp2 = dot_product(cf_e(n1,g1,s,9), cf_h(n2,g1,s,4));
              z4.r+=tmp1.r+tmp2.r; z4.i+=tmp1.i+tmp2.i;
              
            }

	    X1.r+=x1.r*dipole_on(0, _atom[g1]);
	    X1.i+=x1.i*dipole_on(0, _atom[g1]);
	    X2.r+=x2.r*dipole_on(1, _atom[g1]);
	    X2.i+=x2.i*dipole_on(1, _atom[g1]);
	    X3.r+=x3.r*dipole_on(2, _atom[g1]);
	    X3.i+=x3.i*dipole_on(2, _atom[g1]);
            X4.r+=x4.r*dipole_on(3, _atom[g1]);
            X4.i+=x4.i*dipole_on(3, _atom[g1]);
            X5.r+=x5.r*dipole_on(4, _atom[g1]);
            X5.i+=x5.i*dipole_on(4, _atom[g1]);
            X6.r+=overlap.r*lattice(g1,0);	    
            X6.i+=overlap.i*lattice(g1,0);	    

	    Y1.r+=y1.r*dipole_on(0, _atom[g1]);
	    Y1.i+=y1.i*dipole_on(0, _atom[g1]);
	    Y2.r+=y2.r*dipole_on(1, _atom[g1]);
	    Y2.i+=y2.i*dipole_on(1, _atom[g1]);
	    Y3.r+=y3.r*dipole_on(2, _atom[g1]);
	    Y3.i+=y3.i*dipole_on(2, _atom[g1]);
            Y4.r+=y4.r*dipole_on(3, _atom[g1]);
            Y4.i+=y4.i*dipole_on(3, _atom[g1]);
            Y5.r+=y5.r*dipole_on(4, _atom[g1]);
            Y5.i+=y5.i*dipole_on(4, _atom[g1]);
            Y6.r+=overlap.r*lattice(g1,1);	    
            Y6.i+=overlap.i*lattice(g1,1);	    

	    Z1.r+=z1.r*dipole_on(0, _atom[g1]);
	    Z1.i+=z1.i*dipole_on(0, _atom[g1]);
	    Z2.r+=z2.r*dipole_on(1, _atom[g1]);
	    Z2.i+=z2.i*dipole_on(1, _atom[g1]);
	    Z3.r+=z3.r*dipole_on(2, _atom[g1]);
	    Z3.i+=z3.i*dipole_on(2, _atom[g1]);
            Z4.r+=z4.r*dipole_on(5, _atom[g1]);
            Z4.i+=z4.i*dipole_on(5, _atom[g1]);
            Z5.r+=overlap.r*lattice(g1,2);	    
            Z5.i+=overlap.i*lattice(g1,2);	    

      }  
      H_dipole_X(n1,n2).r = X1.r+X2.r+X3.r+X4.r+X5.r+X6.r; 
      H_dipole_X(n1,n2).i = X1.i+X2.i+X3.i+X4.i+X5.i+X6.i; 
      H_dipole_Y(n1,n2).r = Y1.r+Y2.r+Y3.r+Y4.r+Y5.r+Y6.r;
      H_dipole_Y(n1,n2).i = Y1.i+Y2.i+Y3.i+Y4.i+Y5.i+Y6.i;
      H_dipole_Z(n1,n2).r = Z1.r+Z2.r+Z3.r+Z4.r+Z5.r; 
      H_dipole_Z(n1,n2).i = Z1.i+Z2.i+Z3.i+Z4.i+Z5.i; 

      double XX = H_dipole_X(n1,n2).r *H_dipole_X(n1,n2).r;
      XX += H_dipole_X(n1,n2).i *H_dipole_X(n1,n2).i;
      double YY = H_dipole_Y(n1,n2).r *H_dipole_Y(n1,n2).r;
      YY += H_dipole_Y(n1,n2).i *H_dipole_Y(n1,n2).i;
      double ZZ = H_dipole_Z(n1,n2).r *H_dipole_Z(n1,n2).r;
      ZZ += H_dipole_Z(n1,n2).i *H_dipole_Z(n1,n2).i;
     
      double T_XX, T_YY, T_ZZ; 
      MPI_Allreduce(&XX, &T_XX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
      MPI_Allreduce(&YY, &T_YY, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
      MPI_Allreduce(&ZZ, &T_ZZ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

      cout_master<<"n1, n2 electron+hole, X Oscillator Strength: "<<n1<<" "<<n2<<" "<<T_XX<<endl;
      cout_master<<"n1, n2 electron+hole, Y Oscillator Strength: "<<n1<<" "<<n2<<" "<<T_YY<<endl;
      cout_master<<"n1, n2 electron+hole, Z Oscillator Strength: "<<n1<<" "<<n2<<" "<<T_ZZ<<endl<<endl;
    }
  }
}

void Dipole::compute_dipole_with_spin_suppressed(void)
{
  for(int n1=0;n1<e_max;n1++){
    for(int n2=0;n2<h_max;n2++){
      double X1, X2, X3, X4, X5, X6; 
      double Y1, Y2, Y3, Y4, Y5, Y6; 
      double Z1, Z2, Z3, Z4, Z5; 
      double Overlap=0.0;
      X1=0.0; X2=0.0; X3=0.0; X4=0.0; X5=0.0; X6=0.0;
      Y1=0.0; Y2=0.0; Y3=0.0; Y4=0.0; Y5=0.0; Y6=0.0;
      Z1=0.0; Z2=0.0; Z3=0.0; Z4=0.0; Z5=0.0;

      for(int g1=g_min;g1<g_max;g1++){

            double x1, x2, x3, x4, x5;
            double y1, y2, y3, y4, y5;
            double z1, z2, z3, z4;
            double overlap;
            double tmp1, tmp2, tmp3, tmp4;
            overlap=0.0;

            for(int orbital=0; orbital<o_max; orbital++)
               overlap += amplitude_e(n1,g1,orbital)*amplitude_h(n2,g1,orbital);  

            Overlap+=overlap;
            tmp1 = amplitude_e(n1,g1,1)*amplitude_h(n2,g1,2);
            tmp2 = amplitude_e(n1,g1,2)*amplitude_h(n2,g1,1);
            x1=tmp1+tmp2;
              
            tmp1 = amplitude_e(n1,g1,0)*amplitude_h(n2,g1,2);
            tmp2 = amplitude_e(n1,g1,2)*amplitude_h(n2,g1,0);
            x2=tmp1+tmp2; 
              
            tmp1 = amplitude_e(n1,g1,3)*amplitude_h(n2,g1,5);
            tmp2 = amplitude_e(n1,g1,5)*amplitude_h(n2,g1,3);
            tmp3 = amplitude_e(n1,g1,4)*amplitude_h(n2,g1,7);
            tmp4 = amplitude_e(n1,g1,7)*amplitude_h(n2,g1,4);
            x3=tmp1+tmp2+tmp3+tmp4; 
              
            tmp1 = amplitude_e(n1,g1,2)*amplitude_h(n2,g1,8);
            tmp2 = amplitude_e(n1,g1,8)*amplitude_h(n2,g1,2);
            x4=tmp1+tmp2; 
              
            tmp1 = amplitude_e(n1,g1,2)*amplitude_h(n2,g1,9);
            tmp2 = amplitude_e(n1,g1,9)*amplitude_h(n2,g1,2);
            x5=tmp1+tmp2;

            tmp1 = amplitude_e(n1,g1,1)*amplitude_h(n2,g1,3);
            tmp2 = amplitude_e(n1,g1,3)*amplitude_h(n2,g1,1);
            y1=tmp1+tmp2;
              
            tmp1 = amplitude_e(n1,g1,0)*amplitude_h(n2,g1,3);
            tmp2 = amplitude_e(n1,g1,3)*amplitude_h(n2,g1,0);
            y2=tmp1+tmp2; 
              
            tmp1 = amplitude_e(n1,g1,2)*amplitude_h(n2,g1,5);
            tmp2 = amplitude_e(n1,g1,5)*amplitude_h(n2,g1,2);
            tmp3 = amplitude_e(n1,g1,4)*amplitude_h(n2,g1,6);
            tmp4 = amplitude_e(n1,g1,6)*amplitude_h(n2,g1,4);
            y3=tmp1+tmp2+tmp3+tmp4; 
              
            tmp1 = amplitude_e(n1,g1,3)*amplitude_h(n2,g1,8);
            tmp2 = amplitude_e(n1,g1,8)*amplitude_h(n2,g1,3);
            y4=tmp1+tmp2; 
              
            tmp1 = amplitude_e(n1,g1,3)*amplitude_h(n2,g1,9);
            tmp2 = amplitude_e(n1,g1,9)*amplitude_h(n2,g1,3);
            y5=tmp1+tmp2;

            tmp1 = amplitude_e(n1,g1,1)*amplitude_h(n2,g1,4);
            tmp2 = amplitude_e(n1,g1,4)*amplitude_h(n2,g1,1);
            z1=tmp1+tmp2;
              
            tmp1 = amplitude_e(n1,g1,0)*amplitude_h(n2,g1,4);
            tmp2 = amplitude_e(n1,g1,4)*amplitude_h(n2,g1,0);
            z2=tmp1+tmp2; 
              
            tmp1 = amplitude_e(n1,g1,2)*amplitude_h(n2,g1,7);
            tmp2 = amplitude_e(n1,g1,7)*amplitude_h(n2,g1,2);
            tmp3 = amplitude_e(n1,g1,3)*amplitude_h(n2,g1,6);
            tmp4 = amplitude_e(n1,g1,6)*amplitude_h(n2,g1,3);
            z3=tmp1+tmp2+tmp3+tmp4; 
              
            tmp1 = amplitude_e(n1,g1,4)*amplitude_h(n2,g1,9);
            tmp2 = amplitude_e(n1,g1,9)*amplitude_h(n2,g1,4);
            z4=tmp1+tmp2; 

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
  
      H_dipole_X(n1,n2).r = X1+X2+X3+X4+X5+X6; 
      H_dipole_Y(n1,n2).r = Y1+Y2+Y3+Y4+Y5+Y6;
      H_dipole_Z(n1,n2).r = Z1+Z2+Z3+Z4+Z5; 

      double XX = H_dipole_X(n1,n2).r *H_dipole_X(n1,n2).r;
      double YY = H_dipole_Y(n1,n2).r *H_dipole_Y(n1,n2).r;
      double ZZ = H_dipole_Z(n1,n2).r *H_dipole_Z(n1,n2).r;

      double T_XX, T_YY, T_ZZ;
      double T_Overlap, T_X6, T_Y6, T_Z5;

      MPI_Allreduce(&XX, &T_XX, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&YY, &T_YY, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&ZZ, &T_ZZ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&X6, &T_X6, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&Y6, &T_Y6, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&Z5, &T_Z5, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&Overlap, &T_Overlap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      cout_master<<"n1, n2 electron+hole, X Oscillator Strength Approx: "<<n1<<" "<<n2<<" "<<T_XX<<endl;
      cout_master<<"n1, n2 electron+hole, Y Oscillator Strength Approx: "<<n1<<" "<<n2<<" "<<T_YY<<endl;
      cout_master<<"n1, n2 electron+hole, Z Oscillator Strength Approx: "<<n1<<" "<<n2<<" "<<T_ZZ<<endl;
      cout_master<<"n1, n2 electron+hole, X Envelope Oscillator: "<<n1<<" "<<n2<<" "<<T_X6<<endl;
      cout_master<<"n1, n2 electron+hole, Y Envelope Oscillator: "<<n1<<" "<<n2<<" "<<T_Y6<<endl;
      cout_master<<"n1, n2 electron+hole, Z Envelope Oscillator: "<<n1<<" "<<n2<<" "<<T_Z5<<endl;
      cout_master<<"n1, n2 electron+hole, Overlap Approx: "<<n1<<" "<<n2<<" "<<T_Overlap<<endl<<endl;
    }
  }
}

