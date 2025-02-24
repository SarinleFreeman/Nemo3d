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


#include "Coulomb.h"
#include "SimpleParser.h"
#include <complex>
#include <sstream>
#include <fstream>
#include <time.h>
#include <math.h>
using namespace std;
using std::string;

//#define VERBAL
#undef VERBAL


void Coulomb::assign_table1(){

num_atom_type = 1 ;
_coul_on  = vector<double>(o_total*num_atom_type);
_exch_on  = vector<double>(o_total*num_atom_type);
_coul_on_avg  = vector<double>(o_total);
_coul_on_rec_sq  = vector<double>(o_total);
num_atom_pair = 1;
atom_pair = vector<int>(2*num_atom_pair);
_exch_off = vector<double>(o_total*num_atom_pair);
dielect = 11.9;
inv_dielect = 1.0/dielect;

atom_pair[0]=0;
atom_pair[1]=0;
int at = 0;

double tmp=0.0;

//coulomb_onsite_ss_ss
tmp=2.29;
coul_on(0,0,at)=tmp;

//coulomb_onsite_ss_s
tmp=2.86;
coul_on(0,1,at)=tmp;
coul_on(1,0,at)=tmp;

//coulomb_onsite_ss_p
tmp=2.86;
coul_on(0,2,at)=tmp; coul_on(0,3,at)=tmp; coul_on(0,4,at)=tmp;
coul_on(2,0,at)=tmp; coul_on(3,0,at)=tmp; coul_on(4,0,at)=tmp;

//coulomb_onsite_ss_d
tmp=2.65;
coul_on(0,5,at)=tmp; coul_on(0,6,at)=tmp; coul_on(0,7,at)=tmp;
coul_on(0,8,at)=tmp; coul_on(0,9,at)=tmp;
coul_on(5,0,at)=tmp; coul_on(6,0,at)=tmp; coul_on(7,0,at)=tmp;
coul_on(8,0,at)=tmp; coul_on(9,0,at)=tmp;

//coulomb_onsite_s_s
tmp=9.69;
coul_on(1,1,at)=tmp;

//coulomb_onsite_s_p
tmp=9.68;
coul_on(1,2,at)=tmp; coul_on(1,3,at)=tmp; coul_on(1,4,at)=tmp;
coul_on(2,1,at)=tmp; coul_on(3,1,at)=tmp; coul_on(4,1,at)=tmp;

//coulomb_onsite_s_d
tmp=4.46;
coul_on(1,5,at)=tmp; coul_on(1,6,at)=tmp; coul_on(1,7,at)=tmp;
coul_on(1,8,at)=tmp; coul_on(1,9,at)=tmp;
coul_on(5,1,at)=tmp; coul_on(6,1,at)=tmp; coul_on(7,1,at)=tmp;
coul_on(8,1,at)=tmp; coul_on(9,1,at)=tmp;

//coulomb_onsite_p_p_homo
tmp=10.49;
coul_on(2,2,at)=tmp; coul_on(3,3,at)=tmp; coul_on(4,4,at)=tmp;

//coulomb_onsite_p_p_hetero
tmp=9.28;
coul_on(2,3,at)=tmp; coul_on(2,4,at)=tmp; coul_on(3,2,at)=tmp;
coul_on(3,4,at)=tmp; coul_on(4,2,at)=tmp; coul_on(4,3,at)=tmp;

//coulomb_onsite_p_d
tmp=4.53;
for(int p=2; p<5; p++) {
    for(int d=5; d<10; d++) {
         coul_on(p,d,at)=tmp; coul_on(d,p,at)=tmp;
    }
}

//coulomb_onsite_d_d
tmp=3.77;
for(int d1=5; d1<10; d1++)
   for(int d2=5; d2<10; d2++)
      coul_on(d1,d2,at)=tmp;

//exchange_onsite_ss_s
tmp=0.0;
exch_on(0,1,at)=tmp; exch_on(1,0,at)=tmp;

//exchange_onsite_ss_p
tmp=9.69;
exch_on(0,2,at)=tmp; exch_on(0,3,at)=tmp; exch_on(0,4,at)=tmp;
exch_on(2,0,at)=tmp; exch_on(3,0,at)=tmp; exch_on(4,0,at)=tmp;

//exchange_onsite_ss_d
tmp=0.0;
exch_on(0,5,at)=tmp; exch_on(0,6,at)=tmp; exch_on(0,7,at)=tmp;
exch_on(0,8,at)=tmp; exch_on(0,9,at)=tmp;
exch_on(5,0,at)=tmp; exch_on(6,0,at)=tmp; exch_on(7,0,at)=tmp;
exch_on(8,0,at)=tmp; exch_on(9,0,at)=tmp;

//exchange_onsite_s_p
tmp=2.28;
exch_on(1,2,at)=tmp; exch_on(1,3,at)=tmp; exch_on(1,4,at)=tmp;
exch_on(2,1,at)=tmp; exch_on(3,1,at)=tmp; exch_on(4,1,at)=tmp;

//exchange_onsite_s_d
tmp=0.0;
exch_on(1,5,at)=tmp; exch_on(1,6,at)=tmp; exch_on(1,7,at)=tmp;
exch_on(1,8,at)=tmp; exch_on(1,9,at)=tmp;
exch_on(5,1,at)=tmp; exch_on(6,1,at)=tmp; exch_on(7,1,at)=tmp;
exch_on(8,1,at)=tmp; exch_on(9,1,at)=tmp;

//exchange_onsite_px_py
tmp =0.62;
exch_on(2,3,at)=tmp; exch_on(2,4,at)=tmp;
exch_on(3,2,at)=tmp; exch_on(3,4,at)=tmp;
exch_on(4,2,at)=tmp; exch_on(4,3,at)=tmp;

//exchange_onsite_p_d
tmp=9.69;
for(int p=2; p<5; p++) {
   for(int d=5; d<10; d++) {
       exch_on(p,d,at)=tmp; exch_on(d,p,at)=tmp;
   }
}

//exchange_onsite_d_d
tmp = 0.0;
for(int d1=5; d1<10; d1++)
   for(int d2=5; d2<10; d2++)
        exch_on(d1,d2,at)=tmp;


//coul_avg not required
  for(int m=0;m<o_total;m++) {
    for(int iatom=0; iatom<num_atom_type; iatom++){
      _coul_on_avg[m]+=_coul_on[m*num_atom_type+iatom];
        }
    _coul_on_avg[m]/=num_atom_type;
    _coul_on_rec_sq[m]=1.0/(_coul_on_avg[m] * _coul_on_avg[m]);
  }

}



void Coulomb::assign_table(ifstream& table)
{

  vector<string> slist;
  int p = 0;
  do {
    int argc = getwords(slist, table);
    if(argc==0) continue;
    if(slist[0] == "num_atom_type") {

       num_atom_type = atoi(slist[1].c_str());
       _coul_on  = vector<double>(o_total*num_atom_type);
       _exch_on  = vector<double>(o_total*num_atom_type);
       _coul_on_avg  = vector<double>(o_total);
       _coul_on_rec_sq  = vector<double>(o_total);

    } else if(slist[0] == "num_atom_pair") {
       num_atom_pair = atoi(slist[1].c_str());
       atom_pair = vector<int>(2*num_atom_pair);
       _exch_off = vector<double>(o_total*num_atom_pair); 

    } else if(slist[0] == "dielectric_constant") {
       dielect = atof(slist[1].c_str());
       inv_dielect = 1.0/dielect;

    } else if(slist[0] == "atom_pair") {
       atom_pair[p*2] = atoi(slist[1].c_str()); 
       atom_pair[p*2+1] = atoi(slist[2].c_str()); 
       p++;
    } else if(slist[0] == "coulomb_onsite_ss_ss"){
      for(int at=0; at<num_atom_type; at++) 
        coul_on(0,0,at) = atof(slist[at+1].c_str());
	
    } else if(slist[0] == "coulomb_onsite_ss_s"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        coul_on(0,1,at)=tmp; coul_on(1,0,at)=tmp;
      }

    } else if(slist[0] == "coulomb_onsite_ss_p"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	coul_on(0,2,at)=tmp; coul_on(0,3,at)=tmp; coul_on(0,4,at)=tmp;       
	coul_on(2,0,at)=tmp; coul_on(3,0,at)=tmp; coul_on(4,0,at)=tmp;       
      }
    } else if(slist[0] == "coulomb_onsite_ss_d"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	coul_on(0,5,at)=tmp; coul_on(0,6,at)=tmp; coul_on(0,7,at)=tmp;       
	coul_on(0,8,at)=tmp; coul_on(0,9,at)=tmp; 
	coul_on(5,0,at)=tmp; coul_on(6,0,at)=tmp; coul_on(7,0,at)=tmp;       
	coul_on(8,0,at)=tmp; coul_on(9,0,at)=tmp; 
      }
    } else if(slist[0] == "coulomb_onsite_s_s"){
      for(int at=0; at<num_atom_type; at++) 
        coul_on(1,1,at) = atof(slist[at+1].c_str());
    } else if(slist[0] == "coulomb_onsite_s_p"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	coul_on(1,2,at)=tmp; coul_on(1,3,at)=tmp; coul_on(1,4,at)=tmp;       
	coul_on(2,1,at)=tmp; coul_on(3,1,at)=tmp; coul_on(4,1,at)=tmp;       
      }
    } else if(slist[0] == "coulomb_onsite_s_d"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	coul_on(1,5,at)=tmp; coul_on(1,6,at)=tmp; coul_on(1,7,at)=tmp;       
	coul_on(1,8,at)=tmp; coul_on(1,9,at)=tmp; 
	coul_on(5,1,at)=tmp; coul_on(6,1,at)=tmp; coul_on(7,1,at)=tmp;       
	coul_on(8,1,at)=tmp; coul_on(9,1,at)=tmp; 
      }
    } else if(slist[0] == "coulomb_onsite_p_p_homo"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        coul_on(2,2,at)=tmp; coul_on(3,3,at)=tmp; coul_on(4,4,at)=tmp;
      }
    } else if(slist[0] == "coulomb_onsite_p_p_hetero"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        coul_on(2,3,at)=tmp; coul_on(2,4,at)=tmp; coul_on(3,2,at)=tmp;
        coul_on(3,4,at)=tmp; coul_on(4,2,at)=tmp; coul_on(4,3,at)=tmp;
      }
    } else if(slist[0] == "coulomb_onsite_p_d"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        for(int p=2; p<5; p++) {
          for(int d=5; d<10; d++) {
             coul_on(p,d,at)=tmp; coul_on(d,p,at)=tmp;
          }
        }      
      }
    } else if(slist[0] == "coulomb_onsite_d_d"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        for(int d1=5; d1<10; d1++) 
          for(int d2=5; d2<10; d2++) 
             coul_on(d1,d2,at)=tmp;
		 
      }
    } else if(slist[0] == "exchange_onsite_ss_s"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        exch_on(0,1,at)=tmp; exch_on(1,0,at)=tmp;
      }
    } else if(slist[0] == "exchange_onsite_ss_p"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	exch_on(0,2,at)=tmp; exch_on(0,3,at)=tmp; exch_on(0,4,at)=tmp;       
	exch_on(2,0,at)=tmp; exch_on(3,0,at)=tmp; exch_on(4,0,at)=tmp;       
      }
    } else if(slist[0] == "exchange_onsite_ss_d"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	exch_on(0,5,at)=tmp; exch_on(0,6,at)=tmp; exch_on(0,7,at)=tmp;       
	exch_on(0,8,at)=tmp; exch_on(0,9,at)=tmp; 
	exch_on(5,0,at)=tmp; exch_on(6,0,at)=tmp; exch_on(7,0,at)=tmp;       
	exch_on(8,0,at)=tmp; exch_on(9,0,at)=tmp; 
      }
    } else if(slist[0] == "exchange_onsite_s_p"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	exch_on(1,2,at)=tmp; exch_on(1,3,at)=tmp; exch_on(1,4,at)=tmp;       
	exch_on(2,1,at)=tmp; exch_on(3,1,at)=tmp; exch_on(4,1,at)=tmp;       
      }
    } else if(slist[0] == "exchange_onsite_s_d"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	exch_on(1,5,at)=tmp; exch_on(1,6,at)=tmp; exch_on(1,7,at)=tmp;       
	exch_on(1,8,at)=tmp; exch_on(1,9,at)=tmp; 
	exch_on(5,1,at)=tmp; exch_on(6,1,at)=tmp; exch_on(7,1,at)=tmp;       
	exch_on(8,1,at)=tmp; exch_on(9,1,at)=tmp; 
      }
    } else if(slist[0] == "exchange_onsite_px_py"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        exch_on(2,3,at)=tmp; exch_on(2,4,at)=tmp;
        exch_on(3,2,at)=tmp; exch_on(3,4,at)=tmp;
        exch_on(4,2,at)=tmp; exch_on(4,3,at)=tmp; 
      }
    } else if(slist[0] == "exchange_onsite_p_d"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        for(int p=2; p<5; p++) {
          for(int d=5; d<10; d++) {
             exch_on(p,d,at)=tmp; exch_on(d,p,at)=tmp;
          }
        }      
      }
    } else if(slist[0] == "exchange_onsite_d_d"){
      for(int at=0; at<num_atom_type; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        for(int d1=5; d1<10; d1++) 
          for(int d2=5; d2<10; d2++) 
             exch_on(d1,d2,at)=tmp; 
      }
    } else if(slist[0] == "exchange_offsite_ss_ss"){
      for(int at=0; at<num_atom_pair; at++) 
        exch_off(0,0,at) = atof(slist[at+1].c_str());
    } else if(slist[0] == "exchange_offsite_ss_s"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        exch_off(0,1,at)=tmp; exch_on(1,0,at)=tmp;
      }
    } else if(slist[0] == "exchange_offsite_ss_p"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	exch_off(0,2,at)=tmp; exch_off(0,3,at)=tmp; exch_off(0,4,at)=tmp;       
      }
    } else if(slist[0] == "exchange_offsite_ss_d"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	exch_off(0,5,at)=tmp; exch_off(0,6,at)=tmp; exch_off(0,7,at)=tmp;       
	exch_off(0,8,at)=tmp; exch_off(0,9,at)=tmp; 
      }
    } else if(slist[0] == "exchange_offsite_s_ss"){
      for(int at=0; at<num_atom_pair; at++) 
        exch_off(1,0,at) = atof(slist[at+1].c_str());
    } else if(slist[0] == "exchange_offsite_s_s"){
      for(int at=0; at<num_atom_pair; at++) 
        exch_off(1,1,at) = atof(slist[at+1].c_str());
    } else if(slist[0] == "exchange_offsite_s_p"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	exch_off(1,2,at)=tmp; exch_off(1,3,at)=tmp; exch_off(1,4,at)=tmp;       
      }
    } else if(slist[0] == "exchange_offsite_s_d"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	exch_off(1,5,at)=tmp; exch_off(1,6,at)=tmp; exch_off(1,7,at)=tmp;       
	exch_off(1,8,at)=tmp; exch_off(1,9,at)=tmp; 
      }
    } else if(slist[0] == "exchange_offsite_p_ss"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	exch_off(2,0,at)=tmp; exch_off(3,0,at)=tmp; exch_off(4,0,at)=tmp;       
      }
    } else if(slist[0] == "exchange_offsite_p_s"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
	exch_off(2,1,at)=tmp; exch_off(3,1,at)=tmp; exch_off(4,1,at)=tmp;       
      }
    } else if(slist[0] == "exchange_offsite_p_p"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        exch_off(2,2,at)=tmp; exch_off(2,3,at)=tmp; exch_off(2,4,at)=tmp;
        exch_off(3,2,at)=tmp; exch_off(3,3,at)=tmp; exch_off(3,4,at)=tmp;
        exch_off(4,2,at)=tmp; exch_off(4,3,at)=tmp; exch_off(4,4,at)=tmp;
      }
    } else if(slist[0] == "exchange_offsite_p_d"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        for(int p=2; p<5; p++) 
          for(int d=5; d<10; d++) 
             exch_off(p,d,at)=tmp; 
      }
    } else if(slist[0] == "exchange_offsite_d_ss"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
          for(int d=5; d<10; d++) 
             exch_off(d,0,at)=tmp; 
      }
    } else if(slist[0] == "exchange_offsite_d_s"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
          for(int d=5; d<10; d++) 
             exch_off(d,1,at)=tmp; 
      }
    } else if(slist[0] == "exchange_offsite_d_p"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        for(int d=5; d<10; d++) 
          for(int p=2; p<5; p++) 
             exch_off(d,p,at)=tmp; 
      }
    } else if(slist[0] == "exchange_offsite_d_d"){
      for(int at=0; at<num_atom_pair; at++) {
        double tmp = atof(slist[at+1].c_str()); 
        for(int d1=5; d1<10; d1++) 
          for(int d2=5; d2<10; d2++) 
             exch_off(d1,d2,at)=tmp;
      }
    } 
  } while (!table.eof());   

//coul_avg not required
  for(int m=0;m<o_total;m++) {
    for(int iatom=0; iatom<num_atom_type; iatom++){ 
      _coul_on_avg[m]+=_coul_on[m*num_atom_type+iatom];
	}
    _coul_on_avg[m]/=num_atom_type;
    _coul_on_rec_sq[m]=1.0/(_coul_on_avg[m] * _coul_on_avg[m]);
  }

/*for(int m=0;m<o_total;m++) {
	cout<<"avg:m="<<m<<_coul_on_avg[m]<<endl;
	cout<<"rec_sq:m="<<m<<_coul_on_rec_sq[m]<<endl;
}*/

}


//always si-si
int Coulomb::pair(int g1, int g2)
{
  int pair_index=-1;
  for(int p=0; p<num_atom_pair; p++) {
    int a1 = _atom[g1]; int a2 = _atom[g2];
    if(a1==atom_pair[p*2] && a2==atom_pair[p*2+1]){
       pair_index = p;
       break;
    }
  }  
  if(pair_index==-1) {
    cout<<"This atom pair is not registered. Fix file coulomb.table."<<endl;
    exit(0);
  }
  return pair_index;
}


void Coulomb::get_data(Complex_t* coef_e, int* atomid, int* nn, double* lat) 
{
   _cf_e = coef_e;
  // _cf_h = coef_h;
   _atom = atomid;
   _nn   = nn;  
   _lattice = lat; 
}

void Coulomb::print_data(void)
{
   cout<<"electron wavefunction"<<endl;
   for(int g=g_min; g<g_max; g++) 
     for(int s=0; s<2; s++) 
         for(int o=0; o<o_max; o++)
            cout<<cf_e(0,g,s,o)<<endl; 
  /* cout<<endl<<"hole wavefunction"<<endl;
   for(int g=g_min; g<g_max; g++) 
     for(int s=0; s<2; s++) 
         for(int o=0; o<o_max; o++)
            cout<<cf_h(0,g,s,o)<<endl; */
   cout<<endl<<"atomic id"<<endl;
   for(int g=g_min; g<g_max; g++) 
      cout<<_atom[g]<<endl;
   cout<<endl<<"lattice point"<<endl;
   for(int g=g_min; g<g_max; g++)
      cout<<lattice(g,0)<<" "<<lattice(g,1)<<" "<<lattice(g,2)<<endl;
   cout<<endl<<"nearest neighbors"<<endl;
   for(int g=g_min; g<g_max; g++)
      cout<<nn(g,0)<<" "<<nn(g,1)<<" "<<nn(g,2)<<" "<<nn(g,3)<<endl;
}

void Coulomb::oneDSlice(int m, int n, int p, int wf, double *r_init, double *dr){
double x=0.0;
double y=0.0;
double z=0.0;
wf--;
if((m==0) && (n==0) && (p==0)){
        cout<<"error: gradient all zero"<<endl;
}else{

if(m!=0){
  for(int i=0;i<numAtoms;i++){
        x=lattice(i,0);
        y=r_init[1]+n*(x-r_init[0])/(double)m;
        z=r_init[2]+p*(x-r_init[0])/(double)m;
        if((lattice(i,1)<=y+dr[1])&&(lattice(i,1)>=y-dr[1])&&(lattice(i,2)<=z+dr[2])&&(lattice(i,2)>=z-dr[2])){
        double sum=0.0;
        Complex_t *psi1=&cf_e(wf,i,0,0);
        for(int j=0;j<2*o_max;j++){
                sum+=(*psi1).real()*(*psi1).real()+(*psi1).imag()*(*psi1).imag();
                psi1++;
        }
        printf("%e %e %e %e\n",lattice(i,0), lattice(i,1), lattice(i,2), sum);
        }
  }
}else if(n!=0){
        for(int i=0;i<numAtoms;i++){
        y=lattice(i,1);
        x=r_init[0]+m*(y-r_init[1])/(double)n;
        z=r_init[2]+p*(y-r_init[1])/(double)n;
        if((lattice(i,0)<=x+dr[0])&&(lattice(i,0)>=x-dr[0])&&(lattice(i,2)<=z+dr[2])&&(lattice(i,2)>=z-dr[2])){
        double sum=0.0;
        Complex_t *psi1=&cf_e(wf,i,0,0);
        for(int j=0;j<2*o_max;j++){
                sum+=(*psi1).real()*(*psi1).real()+(*psi1).imag()*(*psi1).imag();
                psi1++;
        }
        //printf("%e %e %e %e\n",lattice(i,0), lattice(i,1), lattice(i,2), sum);
        printf("%e %e\n", lattice(i,1), sum);
        }
 }
}else{
        for(int i=0;i<numAtoms;i++){
        z=lattice(i,3);
        x=r_init[0]+m*(z-r_init[2])/(double)p;
        y=r_init[1]+n*(z-r_init[2])/(double)p;
        if((lattice(i,0)<=x+dr[0])&&(lattice(i,0)>=x-dr[0])&&(lattice(i,1)<=y+dr[1])&&(lattice(i,1)>=y-dr[1])){
        double sum=0.0;
        Complex_t *psi1=&cf_e(wf,i,0,0);
        for(int j=0;j<2*o_max;j++){
                sum+=(*psi1).real()*(*psi1).real()+(*psi1).imag()*(*psi1).imag();
                psi1++;
        }
        printf("%e %e %e %e\n",lattice(i,0), lattice(i,1), lattice(i,2), sum);
        //printf("%e %e\n", lattice(i,1), sum);
        }
  }
}
}// end of outer if
}

/*
void Coulomb::get_LS(){

 cout<<"break point 5"<<endl;
    cvectr X_WF = Cvectr(d->seg_ln[my_id]);
    cvectr Y_WF = Cvectr(d->seg_ln[my_id]);
    cvectr Z_WF = Cvectr(d->seg_ln[my_id]);
    cvectr WF = Cvectr(d->seg_ln[my_id]);
    complex cf;
    int n=0;
    for(int iat=g_min; iat<g_max; iat++){
      double X, Y, Z;
      X=lattice(iat,0)-X0;
      Y=lattice(iat,1)-Y0;
      Z=lattice(iat,2)-Z0;
      for(int spin=0; spin<2; spin++){
        for(int iorb=0; iorb<o_max; iorb++){
          if(carrier=='e') cf = cf_e(n1,iat,spin,iorb);
          else             cf = cf_h(n1,iat,spin,iorb);
          X_WF[n].r = cf.r*X;
          X_WF[n].i = cf.i*X;
          Y_WF[n].r = cf.r*Y;
          Y_WF[n].i = cf.i*Y;
          Z_WF[n].r = cf.r*Z;
          Z_WF[n].i = cf.i*Z;
          WF[n].r = cf.r;
          WF[n].i = cf.i;
          n++;
        }
      }
    }
   if(direction=='x'){
      cvectr H_Y_WF = Cvectr(d->seg_ln[my_id]);
      cvectr H_Z_WF = Cvectr(d->seg_ln[my_id]);
      Hmult_spds_complete(H_Y_WF, d, 0.0, Y_WF);
      Hmult_spds_complete(H_Z_WF, d, 0.0, Z_WF);
      complex Y_H_Z = vect_dot_vect(Y_WF, H_Z_WF, d);
      complex Z_H_Y = vect_dot_vect(Z_WF, H_Y_WF, d);
      complex global_L;
      global_L.r=13.22*(Z_H_Y.i-Y_H_Z.i);
      global_L.i=13.22*(Y_H_Z.r-Z_H_Y.r);
      rm_cvectr(&X_WF); rm_cvectr(&Y_WF); rm_cvectr(&Z_WF); rm_cvectr(&WF);
      rm_cvectr(&H_Y_WF); rm_cvectr(&H_Z_WF);
      return global_L;
    }
    else if(direction=='y'){
      cvectr H_X_WF = Cvectr(d->seg_ln[my_id]);
      cvectr H_Z_WF = Cvectr(d->seg_ln[my_id]);
      Hmult_spds_complete(H_X_WF, d, 0.0, X_WF);
      Hmult_spds_complete(H_Z_WF, d, 0.0, Z_WF);
      complex Z_H_X = vect_dot_vect(Z_WF, H_X_WF, d);
      complex X_H_Z = vect_dot_vect(X_WF, H_Z_WF, d);
      complex global_L;
      global_L.r=13.22*(X_H_Z.i-Z_H_X.i);
      global_L.i=13.22*(Z_H_X.r-X_H_Z.r);
      rm_cvectr(&X_WF); rm_cvectr(&Y_WF); rm_cvectr(&Z_WF); rm_cvectr(&WF);
      rm_cvectr(&H_X_WF); rm_cvectr(&H_Z_WF);
      return global_L;
    }
  else {
         cout<<"break point 7"<<endl;

      cvectr H_X_WF = Cvectr(d->seg_ln[my_id]);
      cvectr H_Y_WF = Cvectr(d->seg_ln[my_id]);
      Hmult_spds_complete(H_X_WF, d, 0.0, X_WF);
      Hmult_spds_complete(H_Y_WF, d, 0.0, Y_WF);
         cout<<"break point 8"<<endl;

      complex X_H_Y = vect_dot_vect(X_WF, H_Y_WF, d);
      complex Y_H_X = vect_dot_vect(Y_WF, H_X_WF, d);
      complex global_L;
      global_L.r=13.22*(Y_H_X.i-X_H_Y.i);
      global_L.i=13.22*(X_H_Y.r-Y_H_X.r);
         cout<<"break point 9"<<endl;

      rm_cvectr(&X_WF); rm_cvectr(&Y_WF); rm_cvectr(&Z_WF); rm_cvectr(&WF);
      rm_cvectr(&H_X_WF); rm_cvectr(&H_Y_WF);
         cout<<"break point 10"<<endl;

      return global_L;
    }


  complex local_s = Complex(0.0,0.0);
	
  
    for(int at=g_min; at<g_max; at++){
      for(int iorb=0; iorb<o_max; iorb++){
        complex c1_up, c1_down, c2_up, c2_down;
        if(carrier=='e') {
          c1_up=cf_e(n1,at,0,iorb);
          c2_up=cf_e(n2,at,0,iorb);
          c1_down=cf_e(n1,at,1,iorb);
          c2_down=cf_e(n2,at,1,iorb);
        }
        else {
          c1_up=cf_h(n1,at,0,iorb);
          c2_up=cf_h(n2,at,0,iorb);
          c1_down=cf_h(n1,at,1,iorb);
          c2_down=cf_h(n2,at,1,iorb);
        }
        if(direction=='x'){
          local_s.r +=c1_up.r*c2_down.r+c1_up.i*c2_down.i
                   +c1_down.r*c2_up.r+c1_down.i*c2_up.i;
          local_s.i +=c1_up.r*c2_down.i-c1_up.i*c2_down.r
                   +c1_down.r*c2_up.i-c1_down.i*c2_up.r;
        }
        else if(direction=='y'){
          local_s.r +=c1_up.r*c2_down.i-c1_up.i*c2_down.r
                   +c1_down.i*c2_up.r-c1_down.r*c2_up.i;
          local_s.i +=c1_down.r*c2_up.r+c1_down.i*c2_up.i
                   -c1_up.r*c2_down.r-c1_up.i*c2_down.i;
        }
        else {
          local_s.r +=c1_up.r*c2_up.r+c1_up.i*c2_up.i
                   -c1_down.r*c2_down.r-c1_down.i*c2_down.i;
          local_s.i +=c1_up.r*c2_up.i-c1_up.i*c2_up.r
                    -c1_down.r*c2_down.i+c1_down.i*c2_down.r;
        }

*/

void Coulomb::density_Chandra(double *nr, double *wf){

double x_coor=0.0;
double y_coor=0.0;
double z_coor=0.0;
double r12=0.0;
double r1=0.0;
double r2=0.0;
double x_coor_prime=0.0;
double y_coor_prime=0.0;
double z_coor_prime=0.0;
double ND=1.0;
double alpha=1.075;
double beta=0.478;
double gamma=0.312;
double a0=0.543095;
double donor_x=28*a0;
double donor_y=28*a0;
double donor_z=28*a0;
double A=0.0;
double B=0.0;
double AA=0.0;
double BB=0.0;
double gamma2=gamma*gamma;
for(int g=g_min;g<g_max;g++){
        x_coor=lattice(g,0);
        y_coor=lattice(g,1);
        z_coor=lattice(g,2);
        r1=(x_coor-donor_x)*(x_coor-donor_x)+(y_coor-donor_y)*
                        (y_coor-donor_y)+(z_coor-donor_z)*(z_coor-donor_z);
	r1=sqrt(r1);
	A=exp(-alpha*r1);
	AA=A*A;
	B=exp(-beta*r1);
	BB=B*B;
        for(int i=0;i<numAtoms;i++){
                x_coor_prime=lattice(i,0);
                y_coor_prime=lattice(i,1);
                z_coor_prime=lattice(i,2);
                r12=(x_coor-x_coor_prime)*(x_coor-x_coor_prime)+(y_coor-y_coor_prime)*
                        (y_coor-y_coor_prime)+(z_coor-z_coor_prime)*(z_coor-z_coor_prime);
                r12=sqrt(r12);
		r2=(x_coor-donor_x)*(x_coor-donor_x)+(y_coor-donor_y)*
                        (y_coor-donor_y)+(z_coor-donor_z)*(z_coor-donor_z);
		r2=sqrt(r2);
                nr[g]+=(AA*exp(-2*beta*r2)+BB*exp(-2*alpha*r2)+2*A*B*exp(-(alpha+beta)*r2))*
			(1+2*gamma*r12+gamma2*r12);
		wf[g]+=(A*exp(-beta*r2)+B*exp(-alpha*r2))*(1+gamma*r12);

        }
}
}

void Coulomb::getSCF(double *wf, double *HF, double *XC){

double x_coor=0.0;
double y_coor=0.0;
double z_coor=0.0;
double R=0.0;
double x_coor_prime=0.0;
double y_coor_prime=0.0;
double z_coor_prime=0.0;

double r0 = 4.28; //a.u.
double q = 1.10;
double bohr_a0=0.0529177208;        //nm
double r0_nm=r0*bohr_a0; //nm
//cout<<"Cut off distnace r0 = "<<r0_nm<<endl;
double dielectric_f=11.9;
double q_r0=r0*q;
double esqr =-1.4399766;
double R_atomic=0.0;
double argument =0.0;

for(int g=g_min;g<g_max;g++){
	x_coor=lattice(g,0);
	y_coor=lattice(g,1);
	z_coor=lattice(g,2);
	for(int i=0;i<numAtoms;i++){
		x_coor_prime=lattice(i,0);
		y_coor_prime=lattice(i,1);
		z_coor_prime=lattice(i,2);
		R=(x_coor-x_coor_prime)*(x_coor-x_coor_prime)+(y_coor-y_coor_prime)*
                        (y_coor-y_coor_prime)+(z_coor-z_coor_prime)*(z_coor-z_coor_prime);
                R=sqrt(R);
		if(R>r0_nm)
			dielectric_f=11.9;
		else{
			R_atomic=R/bohr_a0;
			argument=q_r0-q*R_atomic;
			dielectric_f=11.9*q_r0/((q*R_atomic) + sinh(argument));
		}	
		if(i==g)
			R=0.543095*sqrt(3.0)/8.0;
	 	HF[g]+=esqr*wf[i]/(R*dielectric_f);
		//XC[g]+=pow(wf[i], 1.0/3.0);
	}
}
}

void Coulomb::getSCF_new(double *wf_coul, double *wf_xc_r, double *wf_xc_i, double *HF, double *XC_r, double *XC_i){

double x_coor=0.0;
double y_coor=0.0;
double z_coor=0.0;
double R=0.0;
double x_coor_prime=0.0;
double y_coor_prime=0.0;
double z_coor_prime=0.0;

double r0 = 4.28; //a.u.
double q = 1.10;
double bohr_a0=0.0529177208;        //nm
double r0_nm=r0*bohr_a0; //nm
//cout<<"Cut off distnace r0 = "<<r0_nm<<endl;
double dielectric_f=11.9;
double q_r0=r0*q;
double esqr =-1.4399766;
double R_atomic=0.0;
double argument =0.0;

for(int g=g_min;g<g_max;g++){
        x_coor=lattice(g,0);
        y_coor=lattice(g,1);
        z_coor=lattice(g,2);
        for(int i=0;i<numAtoms;i++){
                x_coor_prime=lattice(i,0);
                y_coor_prime=lattice(i,1);
                z_coor_prime=lattice(i,2);
                R=(x_coor-x_coor_prime)*(x_coor-x_coor_prime)+(y_coor-y_coor_prime)*
                        (y_coor-y_coor_prime)+(z_coor-z_coor_prime)*(z_coor-z_coor_prime);
                R=sqrt(R);
                if(R>r0_nm)
                        dielectric_f=11.9;
                else{
                        R_atomic=R/bohr_a0;
                        argument=q_r0-q*R_atomic;
                        dielectric_f=11.9*q_r0/((q*R_atomic) + sinh(argument));
                }
                if(i==g)
                        R=0.543095*sqrt(3.0)/8.0;

                HF[g]+=esqr*wf_coul[i]/(R*dielectric_f);
                XC_r[g]+=esqr*wf_xc_r[i]/(R*dielectric_f);
		XC_i[g]+=esqr*wf_xc_i[i]/(R*dielectric_f);
        }
}
}


void Coulomb::getSCF_CI(double *n_r, double *HF){

double x_coor=0.0;
double y_coor=0.0;
double z_coor=0.0;
double R=0.0;
double x_coor_prime=0.0;
double y_coor_prime=0.0;
double z_coor_prime=0.0;

double r0 = 4.28; //a.u.
double q = 1.10;
double bohr_a0=0.0529177208;        //nm
double r0_nm=r0*bohr_a0; //nm
double dielectric_f=11.9;
double q_r0=r0*q;
double esqr =-1.4399766;
double R_atomic=0.0;
double argument =0.0;

for(int g=g_min;g<g_max;g++){
        x_coor=lattice(g,0);
        y_coor=lattice(g,1);
        z_coor=lattice(g,2);
        for(int i=0;i<numAtoms;i++){
                x_coor_prime=lattice(i,0);
                y_coor_prime=lattice(i,1);
                z_coor_prime=lattice(i,2);
                R=(x_coor-x_coor_prime)*(x_coor-x_coor_prime)+(y_coor-y_coor_prime)*
                        (y_coor-y_coor_prime)+(z_coor-z_coor_prime)*(z_coor-z_coor_prime);
                R=sqrt(R);

               /* if(R>r0_nm)
                        dielectric_f=11.9;
                else{
                        R_atomic=R/bohr_a0;
                        argument=q_r0-q*R_atomic;
                        dielectric_f=11.9*q_r0/((q*R_atomic) + sinh(argument));
                }
		*/
                if(i==g)
                        R=0.543095*sqrt(3.0)/8.0;

                HF[g]+=esqr*n_r[i]/(R*dielectric_f);
        }
}
}

void Coulomb::getSCF_image(double *wf, double *HF, double *XC){

double x_coor=0.0;
double y_coor=0.0;
double z_coor=0.0;
double R=0.0;
double x_coor_prime=0.0;
double y_coor_prime=0.0;
double z_coor_prime=0.0;
//partial metallic screening
//add -1 in front as charges are opposite in sign
double Q=0.5;

double r0 = 4.28; //a.u.
double q = 1.10;
double bohr_a0=0.0529177208;        //nm
double r0_nm=r0*bohr_a0; //nm
//cout<<"Cut off distnace r0 = "<<r0_nm<<endl;
double dielectric_f=11.9;
double q_r0=r0*q;
double esqr =-1.4399766;
double R_atomic=0.0;
double argument =0.0;

for(int g=g_min;g<g_max;g++){
        x_coor=lattice(g,0);
        y_coor=lattice(g,1);
        z_coor=lattice(g,2);
        for(int i=0;i<numAtoms;i++){
                x_coor_prime=lattice(i,0);
                y_coor_prime=lattice(i,1);
                z_coor_prime=lattice(i,2);
                R=(x_coor-x_coor_prime)*(x_coor-x_coor_prime)+(y_coor+y_coor_prime)*
                        (y_coor+y_coor_prime)+(z_coor-z_coor_prime)*(z_coor-z_coor_prime);
                R=sqrt(R);
                if(R>r0_nm)
                        dielectric_f=11.9;
                else{
                        R_atomic=R/bohr_a0;
                        argument=q_r0-q*R_atomic;
                        dielectric_f=11.9*q_r0/((q*R_atomic) + sinh(argument));
                }
                if(i==g)
                        R=0.543095*sqrt(3.0)/8.0;

                HF[g]+=-esqr*Q*wf[i]/(R*dielectric_f);
                //XC[g]+=pow(wf[i], 1.0/3.0);
        }
}
}

void Coulomb::AHF_wf(double *wf, double *A_aniso){
double x_coor=0.0;
double y_coor=0.0;
double z_coor=0.0;
double R_sqr=0.0;
double R =0.0;
double R5=0.0;

//calculate anisotropic tensor
for(int g=g_min;g<g_max;g++){
for(int i=0;i<numAtoms;i++){
        if(i!=g){
        	x_coor=lattice(i,0)-lattice(g,0);
        	y_coor=lattice(i,1)-lattice(g,1);
        	z_coor=lattice(i,2)-lattice(g,2);
        	R_sqr=x_coor*x_coor+y_coor*y_coor+z_coor*z_coor;
        	R=sqrt(R_sqr);
        	R5=R_sqr*R_sqr*R;
        	A_aniso[g*9]+=wf[i]*(3*x_coor*x_coor-R_sqr)/R5;
        	A_aniso[g*9+1]+=wf[i]*(3*x_coor*y_coor)/R5;
        	A_aniso[g*9+2]+=wf[i]*(3*x_coor*z_coor)/R5;
        	A_aniso[g*9+3]+=wf[i]*(3*y_coor*x_coor)/R5;
        	A_aniso[g*9+4]+=wf[i]*(3*y_coor*y_coor-R_sqr)/R5;
        	A_aniso[g*9+5]+=wf[i]*(3*y_coor*z_coor)/R5;
        	A_aniso[g*9+6]+=wf[i]*(3*z_coor*x_coor)/R5;
        	A_aniso[g*9+7]+=wf[i]*(3*z_coor*y_coor)/R5;
        	A_aniso[g*9+8]+=wf[i]*(3*z_coor*z_coor-R_sqr)/R5;
        }
}
}//end of g
}

void Coulomb::AHF(double *P_coor, double *Si0, double *A_aniso, double *A_iso){
double x_coor=0.0;
double y_coor=0.0;
double z_coor=0.0;
double R_sqr=0.0;
double R =0.0;
double R5=0.0;
for(int i=0;i<9;i++)
	A_aniso[i]=0.0;

int index=-1;
double previous_R_sqr=1000000;
for(int i=0;i<numAtoms;i++){
        x_coor=lattice(i,0)-Si0[0];
        y_coor=lattice(i,1)-Si0[1];
        z_coor=lattice(i,2)-Si0[2];
        R_sqr=x_coor*x_coor+y_coor*y_coor+z_coor*z_coor;
        if(R_sqr<=previous_R_sqr){
                index=i;
                previous_R_sqr=R_sqr;
        }
}
R_sqr=0;
cout<<"Index of the 29 Si atom : "<<index<<endl;
cout<<"x_coor : "<<lattice(index,0)<<endl;
cout<<"y_coor : "<<lattice(index,1)<<endl;
cout<<"z_coor : "<<lattice(index,2)<<endl;

//calculate anisotropic tensor
for(int i=0;i<numAtoms;i++){
	if(i!=index){
	x_coor=lattice(i,0)-Si0[0];
	y_coor=lattice(i,1)-Si0[1];
	z_coor=lattice(i,2)-Si0[2];
	R_sqr=x_coor*x_coor+y_coor*y_coor+z_coor*z_coor;
	R=sqrt(R_sqr);
	R5=R_sqr*R_sqr*R;
	Complex_t tmp=0.0;
	for(int s=0;s<2;s++)
	for(int o=0;o<o_max;o++)
		tmp+=conj(cf_e(0,i,s,o))*cf_e(0,i,s,o);
	A_aniso[0]+=tmp.real()*(3*x_coor*x_coor-R_sqr)/R5;
	A_aniso[1]+=tmp.real()*(3*x_coor*y_coor)/R5;
	A_aniso[2]+=tmp.real()*(3*x_coor*z_coor)/R5;
	A_aniso[3]+=tmp.real()*(3*y_coor*x_coor)/R5;
	A_aniso[4]+=tmp.real()*(3*y_coor*y_coor-R_sqr)/R5;
	A_aniso[5]+=tmp.real()*(3*y_coor*z_coor)/R5;
	A_aniso[6]+=tmp.real()*(3*z_coor*x_coor)/R5;
	A_aniso[7]+=tmp.real()*(3*z_coor*y_coor)/R5;
	A_aniso[8]+=tmp.real()*(3*z_coor*z_coor-R_sqr)/R5;
	}
}

//calculate isotropic contribution
A_iso[0]=0.0;
Complex_t tmp1=0.0;
if(index!=-1){
 for(int s=0;s<2;s++)
        for(int o=0;o<o_max;o++)
                tmp1+=conj(cf_e(0,index,s,o))*cf_e(0,index,s,o);
}
A_iso[0]=tmp1.real();
cout<<"tmp 1 : "<<tmp1<<endl;
cout<<endl<<endl;

//find the donor nearest neigbor coordinates
index=-1;
previous_R_sqr=1000000;
for(int i=0;i<numAtoms;i++){
        x_coor=lattice(i,0)-P_coor[0];
        y_coor=lattice(i,1)-P_coor[1];
        z_coor=lattice(i,2)-P_coor[2];
        R_sqr=x_coor*x_coor+y_coor*y_coor+z_coor*z_coor;
        if(R_sqr<=previous_R_sqr){
                index=i;
                previous_R_sqr=R_sqr;
        }
}

if(index!=-1){
	int nbr=-1;
	for(int i=0;i<4;i++){
		nbr=nn(index, i);
		cout<<"Nearest Neighbour "<<i<<"   :  ("<<lattice(nbr,0)<<", "<<lattice(nbr,1)<<", "<<lattice(nbr,2)<<")"<<endl;
	}
}

cout<<endl<<endl;
}

void Coulomb::HF(int c1, int *P_index, double *A_aniso, double *A_iso){
cout<<"Number of donors : "<<c1<<endl;
 for(int h=0;h<c1;h++){
        for(int i=0;i<9;i++)
                A_aniso[9*h+i]=0.0;
        A_iso[h]=0.0;
        }

 for(int h=0;h<c1;h++){
        for(int i=0;i<9;i++)
                cout<<A_aniso[9*h+i]<<endl;
        cout<<A_iso[h]<<endl<<endl;
        }

for(int h=0;h<c1;h++){

	double x_coor=0.0;
	double y_coor=0.0;
	double z_coor=0.0;
	double R_sqr=0.0;
	double R =0.0;
	double R5=0.0;
	double P_x=lattice(P_index[h],0);
	double P_y=lattice(P_index[h],1);
	double P_z=lattice(P_index[h],2);
	cout<<"Index of the P atom : "<<P_index[h]<<endl;
	cout<<"x_coor : "<<P_x<<endl;
	cout<<"y_coor : "<<P_y<<endl;
	cout<<"z_coor : "<<P_z<<endl;

	//calculate anisotropic tensor
	for(int i=0;i<numAtoms;i++){
		if(i!=P_index[h]){
			x_coor=lattice(i,0)-P_x;
			y_coor=lattice(i,1)-P_y;
			z_coor=lattice(i,2)-P_z;
			R_sqr=x_coor*x_coor+y_coor*y_coor+z_coor*z_coor;
			R=sqrt(R_sqr);
			R5=R_sqr*R_sqr*R;
			Complex_t tmp=0.0;
			for(int s=0;s<2;s++)
				for(int o=0;o<10;o++)
					tmp+=conj(cf_e(0,i,s,o))*cf_e(0,i,s,o);
			A_aniso[9*h]+=tmp.real()*(3*x_coor*x_coor-R_sqr)/R5;
			A_aniso[9*h+1]+=tmp.real()*(3*x_coor*y_coor)/R5;
			A_aniso[9*h+2]+=tmp.real()*(3*x_coor*z_coor)/R5;
			A_aniso[9*h+3]+=tmp.real()*(3*y_coor*x_coor)/R5;
			A_aniso[9*h+4]+=tmp.real()*(3*y_coor*y_coor-R_sqr)/R5;
			A_aniso[9*h+5]+=tmp.real()*(3*y_coor*z_coor)/R5;
			A_aniso[9*h+6]+=tmp.real()*(3*z_coor*x_coor)/R5;
			A_aniso[9*h+7]+=tmp.real()*(3*z_coor*y_coor)/R5;
			A_aniso[9*h+8]+=tmp.real()*(3*z_coor*z_coor-R_sqr)/R5;
		}
	}	

	//calculate isotropic contribution
	Complex_t tmp1=0.0;
    		for(int s=0;s<2;s++)
        		for(int o=0;o<10;o++)
                		tmp1+=conj(cf_e(0,P_index[h],s,o))*cf_e(0,P_index[h],s,o);
	
	A_iso[h]=tmp1.real();
}//end of h
}

void Coulomb::AHF_shell(double *Si0, double *A_aniso, double *A_iso){
double x_coor=0.0;
double y_coor=0.0;
double z_coor=0.0;
double R_sqr=0.0;
double R =0.0;
double R5=0.0;
for(int i=0;i<9;i++)
	A_aniso[i]=0.0;
A_iso[0]=0.0;
int index=-1;
double previous_R_sqr=10000000;
for(int i=0;i<numAtoms;i++){
        x_coor=lattice(i,0)-Si0[0];
        y_coor=lattice(i,1)-Si0[1];
        z_coor=lattice(i,2)-Si0[2];
        R_sqr=x_coor*x_coor+y_coor*y_coor+z_coor*z_coor;
        if(R_sqr<=previous_R_sqr){
                index=i;
                previous_R_sqr=R_sqr;
        }
}

R_sqr=0;
cout<<"Index of the 29 Si atom : "<<index<<endl;
cout<<"x_coor : "<<lattice(index,0)<<endl;
cout<<"y_coor : "<<lattice(index,1)<<endl;
cout<<"z_coor : "<<lattice(index,2)<<endl;

//calculate anisotropic tensor
for(int i=0;i<numAtoms;i++){
	if(i!=index){
		x_coor=lattice(i,0)-Si0[0];
		y_coor=lattice(i,1)-Si0[1];
		z_coor=lattice(i,2)-Si0[2];
		R_sqr=x_coor*x_coor+y_coor*y_coor+z_coor*z_coor;
		R=sqrt(R_sqr);
		R5=R_sqr*R_sqr*R;
		Complex_t tmp=0.0;
		for(int s=0;s<2;s++)
		for(int o=0;o<10;o++)
			tmp+=conj(cf_e(0,i,s,o))*cf_e(0,i,s,o);
		A_aniso[0]+=tmp.real()*(3*x_coor*x_coor-R_sqr)/R5;
		A_aniso[1]+=tmp.real()*(3*x_coor*y_coor)/R5;
		A_aniso[2]+=tmp.real()*(3*x_coor*z_coor)/R5;
		A_aniso[3]+=tmp.real()*(3*y_coor*x_coor)/R5;
		A_aniso[4]+=tmp.real()*(3*y_coor*y_coor-R_sqr)/R5;
		A_aniso[5]+=tmp.real()*(3*y_coor*z_coor)/R5;
		A_aniso[6]+=tmp.real()*(3*z_coor*x_coor)/R5;
		A_aniso[7]+=tmp.real()*(3*z_coor*y_coor)/R5;
		A_aniso[8]+=tmp.real()*(3*z_coor*z_coor-R_sqr)/R5;
	}
}

//calculate isotropic contribution
A_iso[0]=0.0;
Complex_t tmp1=0.0;
if(index!=-1){
    for(int s=0;s<2;s++)
        for(int o=0;o<10;o++)
                tmp1+=conj(cf_e(0,index,s,o))*cf_e(0,index,s,o);
}
A_iso[0]=tmp1.real();
}


double Coulomb::w_exch(int g1, int g2, int o1, int o2)
{
   double final;
   	
   if(g1==g2) final=exch_on(o1,o2,_atom[g1]) ; //consider this line only for onsite exchange
   else       final=exch_off(o1,o2,pair(g1,g2)); //of-site exchange is 0 for our case
   return final;
}

double Coulomb::w_coul(double dd, int g1, int o1, int o2)
{
  double final; 
  if(dd==0.0) {
    final=coul_on(o1,o2,_atom[g1]);
  }
  else {
    double tmp = coul_on_avg(o1,o2);	  
    final = 1.0/sqrt(1.0/(tmp*tmp)+dd*0.48268);
  }    
  return final;
}

double Coulomb::distance_sq(int g1, int g2)
{
  double dx = lattice(g1,0)-lattice(g2,0);
  double dy = lattice(g1,1)-lattice(g2,1);
  double dz = lattice(g1,2)-lattice(g2,2);

  double dd = dx*dx+dy*dy+dz*dz;
  return dd;
}  

void Coulomb::dotProduct(int wf1, int wf2, Complex_t *sum){
//computes the dot product beween two wavefunctions
wf1--;
wf2--;
Complex_t result=0.0;
for(int x=0;x<numAtoms;x++)
	for(int s=0;s<2;s++)
		for(int o=0;o<o_max;o++)
			result+=conj(cf_e(wf1,x,s,o))*cf_e(wf2,x,s,o);
*sum=result;
}

void Coulomb::x2(Complex_t *stdev_x, Complex_t *stdev_y, Complex_t *stdev_z, Complex_t *stdev_r, 
Complex_t *stdev_lx, Complex_t *stdev_ly, Complex_t *stdev_lz, double xP, double yP, double zP){
//computes the dot product beween two wavefunctions
cout<<"O_max : "<<o_max<<endl;
double x=0.0;
double y=0.0;
double z=0.0;
double r=0.0;

Complex_t std_x=0.0;
Complex_t std_y=0.0;
Complex_t std_z=0.0;
Complex_t std_r=0.0;
Complex_t std_lx=0.0;
Complex_t std_ly=0.0;
Complex_t std_lz=0.0;

cout<<"Coordinates :  "<<xP<<"    "<<yP<<"    "<<zP<<endl;

for(int i=0;i<numAtoms;i++){
	Complex_t tmp=0.0;
	x=lattice(i,0)-xP;
	y=lattice(i,1)-yP;
	z=lattice(i,2)-zP;
	r=x*x+y*y+z*z;
        for(int s=0;s<2;s++)
                for(int o=0;o<o_max;o++)
                        tmp+=conj(cf_e(0,i,s,o))*cf_e(0,i,s,o);
	std_x+=x*x*tmp;
	std_y+=y*y*tmp;
	std_z+=z*z*tmp;
	std_r+=r*tmp;
	std_lx+=x*tmp;
        std_ly+=y*tmp;
        std_lz+=z*tmp;
}
*stdev_x=std_x;
*stdev_y=std_y;
*stdev_z=std_z;
*stdev_r=std_r;
*stdev_lx=std_lx;
*stdev_ly=std_ly;
*stdev_lz=std_lz;
}


void Coulomb::dipole_moment(Complex_t *moment_x, Complex_t *moment_y, Complex_t *moment_z){
//computes the dot product beween two wavefunctions
double x=0.0;
double y=0.0;
double z=0.0;

Complex_t dipole_x=0.0;
Complex_t dipole_y=0.0;
Complex_t dipole_z=0.0;

for(int i=0;i<numAtoms;i++){
        Complex_t tmp=0.0;
        x=lattice(i,0);
        y=lattice(i,1);
        z=lattice(i,2);
        for(int s=0;s<2;s++)
                for(int o=0;o<o_max;o++)
                        tmp+=conj(cf_e(0,i,s,o))*cf_e(1,i,s,o);
        dipole_x+=x*tmp;
        dipole_y+=y*tmp;
        dipole_z+=z*tmp;
}
*moment_x=dipole_x;
*moment_y=dipole_y;
*moment_z=dipole_z;
}


void Coulomb::moment(Complex_t *moment_x, Complex_t *moment_y, Complex_t *moment_z, double xP, double yP, double zP){
//computes the dot product beween two wavefunctions
cout<<"O_max : "<<o_max<<endl;
double x=0.0;
double y=0.0;
double z=0.0;

Complex_t dipole_x=0.0;
Complex_t dipole_y=0.0;
Complex_t dipole_z=0.0;

cout<<"Coordinates :  "<<xP<<"    "<<yP<<"    "<<zP<<endl;

for(int i=0;i<numAtoms;i++){
        Complex_t tmp=0.0;
        x=lattice(i,0)-xP;
        y=lattice(i,1)-yP;
        z=lattice(i,2)-zP;
        for(int s=0;s<2;s++)
                for(int o=0;o<o_max;o++)
                        tmp+=conj(cf_e(0,i,s,o))*cf_e(0,i,s,o);
        dipole_x+=x*tmp;
        dipole_y+=y*tmp;
        dipole_z+=z*tmp;
}
*moment_x=dipole_x;
*moment_y=dipole_y;
*moment_z=dipole_z;
}

void Coulomb::r2(Complex_t *r2_x, Complex_t *r2_y, Complex_t *r2_z, double xP, double yP, double zP){
//computes the dot product beween two wavefunctions
cout<<"O_max : "<<o_max<<endl;
double x=0.0;
double y=0.0;
double z=0.0;

Complex_t r2_x_loc=0.0;
Complex_t r2_y_loc=0.0;
Complex_t r2_z_loc=0.0;

cout<<"Coordinates :  "<<xP<<"    "<<yP<<"    "<<zP<<endl;

for(int i=0;i<numAtoms;i++){
        Complex_t tmp=0.0;
        x=lattice(i,0)-xP;
        y=lattice(i,1)-yP;
        z=lattice(i,2)-zP;
        for(int s=0;s<2;s++)
                for(int o=0;o<o_max;o++)
                        tmp+=conj(cf_e(0,i,s,o))*cf_e(0,i,s,o);
        r2_x_loc+=(y*y+z*z)*tmp;
        r2_y_loc+=(x*x+z*z)*tmp;
        r2_z_loc+=(x*x+y*y)*tmp;
}
*r2_x=r2_x_loc;
*r2_y=r2_y_loc;
*r2_z=r2_z_loc;
}


void Coulomb::getImpurityInfo(int wf){
int g1=-1;
wf--;
for(int i=0;i<numAtoms;i++)
	if(_atom[i]==1)
		g1=i;

if(g1!=-1){

cout<<"Atom Index of the impurity: "<<g1<<endl;
cout<<"Coordinates of the impurity:";
cout<<"("<<lattice(g1,0)<<","<<lattice(g1,1)<<","<<lattice(g1,2)<<")"<<endl<<endl;

cout<<"Orbitals at the impurity site for wavefunction "<<wf+1<<":"<<endl;
cout<<"s_up       :   "<<cf_e(wf,g1,0,0)<<endl;
cout<<"px_up      :   "<<cf_e(wf,g1,0,1)<<endl;
cout<<"py_up      :   "<<cf_e(wf,g1,0,2)<<endl;
cout<<"pz_up      :   "<<cf_e(wf,g1,0,3)<<endl;
cout<<"dxy_up     :   "<<cf_e(wf,g1,0,4)<<endl;
cout<<"dyz_up     :   "<<cf_e(wf,g1,0,5)<<endl;
cout<<"dzx_up     :   "<<cf_e(wf,g1,0,6)<<endl;
cout<<"dx2-y2_up  :   "<<cf_e(wf,g1,0,7)<<endl;
cout<<"dz2_up     :   "<<cf_e(wf,g1,0,8)<<endl;
cout<<"s*_up      :   "<<cf_e(wf,g1,0,9)<<endl;
cout<<"s_down     :   "<<cf_e(wf,g1,1,0)<<endl;
cout<<"px_down    :   "<<cf_e(wf,g1,1,1)<<endl;
cout<<"py_down    :   "<<cf_e(wf,g1,1,2)<<endl;
cout<<"pz_down    :   "<<cf_e(wf,g1,1,3)<<endl;
cout<<"dxy_down   :   "<<cf_e(wf,g1,1,4)<<endl;
cout<<"dyz_down   :   "<<cf_e(wf,g1,1,5)<<endl;
cout<<"dzx_down   :   "<<cf_e(wf,g1,1,6)<<endl;
cout<<"dx2-y2_down:   "<<cf_e(wf,g1,1,7)<<endl;
cout<<"dz2_down   :   "<<cf_e(wf,g1,1,8)<<endl;
cout<<"s*_down    :   "<<cf_e(wf,g1,1,9)<<endl<<endl;
cout<<"Nearest neighbour information:"<<endl;
int nn_P;
for(int adj=0;adj<4;adj++){
      nn_P=nn(g1,adj);

cout<<"Atom Index of nearest Neigbour "<<adj<<" : "<<nn_P<<endl;
cout<<"Coordinates: ";
cout<<"("<<lattice(nn_P,0)<<","<<lattice(nn_P,1)<<","<<lattice(nn_P,2)<<")"<<endl<<endl;

cout<<"Orbitals at the impurity site for wavefunction "<<wf+1<<":"<<endl;
cout<<"s_up       :   "<<cf_e(wf,nn_P,0,0)<<endl;
cout<<"px_up      :   "<<cf_e(wf,nn_P,0,1)<<endl;
cout<<"py_up      :   "<<cf_e(wf,nn_P,0,2)<<endl;
cout<<"pz_up      :   "<<cf_e(wf,nn_P,0,3)<<endl;
cout<<"dxy_up     :   "<<cf_e(wf,nn_P,0,4)<<endl;
cout<<"dyz_up     :   "<<cf_e(wf,nn_P,0,5)<<endl;
cout<<"dzx_up     :   "<<cf_e(wf,nn_P,0,6)<<endl;
cout<<"dx2-y2_up  :   "<<cf_e(wf,nn_P,0,7)<<endl;
cout<<"dz2_up     :   "<<cf_e(wf,nn_P,0,8)<<endl;
cout<<"s*_up      :   "<<cf_e(wf,nn_P,0,9)<<endl;
cout<<"s_down     :   "<<cf_e(wf,nn_P,1,0)<<endl;
cout<<"px_down    :   "<<cf_e(wf,nn_P,1,1)<<endl;
cout<<"py_down    :   "<<cf_e(wf,nn_P,1,2)<<endl;
cout<<"pz_down    :   "<<cf_e(wf,nn_P,1,3)<<endl;
cout<<"dxy_down   :   "<<cf_e(wf,nn_P,1,4)<<endl;
cout<<"dyz_down   :   "<<cf_e(wf,nn_P,1,5)<<endl;
cout<<"dzx_down   :   "<<cf_e(wf,nn_P,1,6)<<endl;
cout<<"dx2-y2_down:   "<<cf_e(wf,nn_P,1,7)<<endl;
cout<<"dz2_down   :   "<<cf_e(wf,nn_P,1,8)<<endl;
cout<<"s*_down    :   "<<cf_e(wf,nn_P,1,9)<<endl<<endl;


}// end of for

cout<<"Coulomb table"<<endl;
cout<<"Onsite coulomb for Si"<<endl;

for(int o1=0; o1< o_max; o1++)
	for(int o2=0; o2<o_max; o2++)
 		cout<<"o1="<<o1<<"  o2="<<o2<<" coul: "<<coul_on(o1,o2,_atom[1])<<endl;
	


cout<<"Onsite Coulomb for P"<<endl;

for(int o1=0;o1<o_max;o1++)
	for(int o2=0;o2<o_max;o2++)
 		cout<<"o1="<<o1<<"  o2="<<o2<<" coul: "<<coul_on(o1,o2,_atom[g1])<<endl;
	


cout<<"Average Onsite Coulomb:"<<endl;

for(int o1=0;o1<o_max;o1++)
	for(int o2=0;o2<o_max;o2++)
 		cout<<"o1="<<o1<<"  o2="<<o2<<" coul: "<<coul_on_avg(o1,o2)<<endl;
	


}//end of if
else{
  cout<<"Impurity not found"<<endl;
}

}//end of function



void Coulomb::compute_coulomb_exchange(void)
{
/* for(int ne=0;ne<e_max;ne++){
   for(int nh=0;nh<h_max;nh++){
     H_coul(ne,nh)=0.0; 
     H_exch(ne,nh)=0.0;
     Complex_t t1=0.0,t2=0.0,t3=0.0,t4=0.0, t5=0.0, t6=0.0;
     for(int s1=0;s1<2;s1++){
       int s1f; if(s1==0) s1f = 1; else s1f = 0;		  
       for(int s2=0;s2<2;s2++){
         int s2f; if(s2==0) s2f = 1; else s2f = 0;       
         for(int x1=g_min;x1<g_max;x1++){
          for(int o1=0;o1<o_max;o1++){
           Complex_t ce1_s1=cf_e(ne,x1,s1,o1);
	   Complex_t ce1_s2f=cf_e(ne,x1,s2f,o1);
           Complex_t ch1_s1f=cf_h(nh,x1,s1f,o1);
           Complex_t ch1_s2=cf_h(nh,x1,s2,o1);
           for(int x2=0;x2<numAtoms;x2++){
              double dd = distance_sq(x1,x2);		   
              for(int o2=0;o2<o_max;o2++){
                Complex_t ce2_s2f=cf_e(ne,x2,s2f,o2);
                Complex_t ch2_s2=cf_h(nh,x2,s2,o2);
                double coul_tmp = w_coul(dd,x1,o1,o2);
                t1+=conj(ce1_s1)*ce1_s1*conj(ch2_s2)*ch2_s2*coul_tmp;
                t4+=conj(ce1_s1)*ch1_s1f*conj(ch2_s2)*ce2_s2f*coul_tmp;
              } //o2
           } // x2
           //for(int adj=0;adj<5;adj++) {
              //int x2=x1;		    
             // if(adj==4) x2=x1; else x2=nn(x1,adj);
             // if(x2!=-1) { 
                for(int o2=0;o2<o_max;o2++){
           	  if(o1!=o2) {
           	    double exch_tmp=w_exch(x1,x1,o1,o2);
                    if(exch_tmp!=0.0) {
                      Complex_t ce2_s1=cf_e(ne,x1,s1,o2);
                      Complex_t ce2_s2f=cf_e(ne,x1,s2f,o2);
           	      Complex_t ch2_s2=cf_h(nh,x1,s2,o2);
		      Complex_t ch2_s1f=cf_h(nh,x1,s1f,o2);
                      t2+=conj(ce1_s1)*ce2_s1*conj(ch1_s2)*ch2_s2*exch_tmp;
             	      t3+=conj(ce1_s1)*ce2_s1*conj(ch2_s2)*ch1_s2*exch_tmp;
                      t5+=conj(ce1_s1)*ch2_s1f*conj(ch2_s2)*ce1_s2f*exch_tmp;
             	      t6+=conj(ce1_s1)*ch2_s1f*conj(ch1_s2)*ce2_s2f*exch_tmp;
                    } / if(exch_tmp!=0.0)
                  } // if((x1!=x2)||(o1!=o2))
		} // o2 
              //} // if(x2!=-1)
           //} // adj
      }} // x1 and o1
    }} // s1 and s2   
   cout<<"T1: "<<t1<<"   T2: "<<t2<<"   T3: "<<t3<<endl;
   cout<<"T4: "<<t4<<"   T5: "<<t5<<"   T6: "<<t6<<endl;
   H_coul(ne,nh)-=(t1+t2+t3)*inv_dielect;
   H_exch(ne,nh)+=(t4+t5+t6)*inv_dielect;
//   cout<<"ne, nh, Coulomb, Exchange: "<<ne<<" "<<nh<<" "<<H_coul(ne,nh)<<" "<<H_exch(ne,nh)<<endl;
 }} // ne and nh
*/
}

void Coulomb::compute_exchange(int m, int n, int p, int q, int id)
{

//dump to file

#ifdef VERBAL 
	string b1;
	stringstream a1;
	a1 << id;
	b1 = a1.str();
	string filename="status_";
	filename+=b1;
	//cout<<filename<<endl;
	//const char *sr=filename.c_str();
	int interval=skip;
	FILE* fp=fopen(filename.c_str(),"w");
	if(fp==NULL){
		cout<<"Error opening file"<<endl;
	}
#endif
m--;
n--;
p--;
q--;

//<mn||pq>=[mp|nq]-[mq|np]
//<mn||mn>=<mn|mn>-<mn|nm>

H_coul()=0.0;
H_exch()=0.0;
   
int count_coul_ptr_o1=0;
int count_coul_ptr_o2=0;
int o_max_new=o_max-4;

Complex_t t1=0.0,t2=0.0,t3=0.0,t4=0.0, t5=0.0, t6=0.0;
     double tmp=0.0;
     double dd=0.0;
     //double dd_cutoff=9.0;	
	int counter1=0;
	complex<double>* dot1_s2_0=new complex<double>[numAtoms*o_max_new];
        complex<double>* dot2_s2_0=new complex<double>[numAtoms*o_max_new];
	complex<double>* dot1_s2_1=new complex<double>[numAtoms*o_max_new];
        complex<double>* dot2_s2_1=new complex<double>[numAtoms*o_max_new];

	complex<double>* mul_1_s1_0=new complex<double>[numAtoms*o_max_new];
        complex<double>* mul_2_s1_0=new complex<double>[numAtoms*o_max_new];
	complex<double>* mul_1_s1_1=new complex<double>[numAtoms*o_max_new];
        complex<double>* mul_2_s1_1=new complex<double>[numAtoms*o_max_new];

        complex<double>* dot3_s2_0=new complex<double>[numAtoms];
        complex<double>* dot4_s2_0=new complex<double>[numAtoms];
	complex<double>* dot3_s2_1=new complex<double>[numAtoms];
        complex<double>* dot4_s2_1=new complex<double>[numAtoms];


	   for(int i=0;i<numAtoms;i++){
		for(int j=0;j<o_max;j++){
			if(j<5){
				dot1_s2_0[counter1]=conj(cf_e(n,i,0,j))*cf_e(q,i,0,j);
				dot2_s2_0[counter1]=conj(cf_e(n,i, 0,j))*cf_e(p, i, 1,j);
				mul_1_s1_0[counter1]=conj(cf_e(m,i,0,j))*cf_e(p,i,0,j);
                          	mul_2_s1_0[counter1]=conj(cf_e(m,i,0,j))*cf_e(q,i,1,j);
				dot1_s2_1[counter1]=conj(cf_e(n,i,1,j))*cf_e(q,i,1,j);
                                dot2_s2_1[counter1]=conj(cf_e(n,i, 1,j))*cf_e(p, i, 0,j);
                                mul_1_s1_1[counter1]=conj(cf_e(m,i,1,j))*cf_e(p,i,1,j);
                                mul_2_s1_1[counter1]=conj(cf_e(m,i,1,j))*cf_e(q,i,0,j);
				counter1++;
			}
			else{
				dot1_s2_0[counter1]+=conj(cf_e(n,i, 0,j))*cf_e(q,i,0,j);
				dot2_s2_0[counter1]+=conj(cf_e(n,i, 0,j))*cf_e(p, i, 1,j);
				mul_1_s1_0[counter1]+=conj(cf_e(m,i,0,j))*cf_e(p,i,0,j);
                          	mul_2_s1_0[counter1]+=conj(cf_e(m,i,0,j))*cf_e(q,i,1,j);
				dot1_s2_1[counter1]+=conj(cf_e(n,i, 1,j))*cf_e(q,i,1,j);
                                dot2_s2_1[counter1]+=conj(cf_e(n,i, 1,j))*cf_e(p, i, 0,j);
                                mul_1_s1_1[counter1]+=conj(cf_e(m,i,1,j))*cf_e(p,i,1,j);
                                mul_2_s1_1[counter1]+=conj(cf_e(m,i,1,j))*cf_e(q,i,0,j);

			}//end of if
		}//end of j
		for(int k=0;k<o_max_new;k++){
			dot3_s2_0[i]+=dot1_s2_0[i*o_max_new+k];
			dot4_s2_0[i]+=dot2_s2_0[i*o_max_new+k];
			dot3_s2_1[i]+=dot1_s2_1[i*o_max_new+k];
                        dot4_s2_1[i]+=dot2_s2_1[i*o_max_new+k];

		}
		counter1++;
	   }//end of atom index i	

        //cout<<"Entering x1 block"<<endl;
	int counter_x1=g_min*o_max_new;
	int counter=0;	
	for(int x1=g_min;x1<g_max;x1++){
	  	count_coul_ptr_o1=0;	
          for(int o1=0;o1<o_max_new;o1++){
    
		Complex_t t1_temp_s2_0=0.0;
		Complex_t t4_temp_s2_0=0.0;
		Complex_t t1_temp_s2_1=0.0;
                Complex_t t4_temp_s2_1=0.0;

		counter=0;
           for(int x2=0;x2<numAtoms;x2++){
		if(x1!=x2) 
   	        	dd = distance_sq(x1,x2);
		else
			dd=0.0;
	//	if(dd<dd_cutoff){	
		for(int o2=0;o2<o_max_new;o2++){
			if(x1!=x2){
				tmp = _coul_on_rec_sq[count_coul_ptr_o1+o2];
				tmp = 1.0/sqrt(tmp+dd*0.48268);
			}
			else{
				tmp =_coul_on_avg[count_coul_ptr_o1+o2];
			}
		t1_temp_s2_0+=dot1_s2_0[counter+o2]*tmp;
		t4_temp_s2_0+=dot2_s2_0[counter+o2]*tmp;
		t1_temp_s2_1+=dot1_s2_1[counter+o2]*tmp;
                t4_temp_s2_1+=dot2_s2_1[counter+o2]*tmp;

	      } //o2
	 /* }else{
		tmp=1.439363044/sqrt(dd);
		t1_temp_s2_0+=dot3_s2_0[x2]*tmp;
		t4_temp_s2_0+=dot4_s2_0[x2]*tmp;
		t1_temp_s2_1+=dot3_s2_1[x2]*tmp;
                t4_temp_s2_1+=dot4_s2_1[x2]*tmp;
	  }*/
	
	   counter+=o_max_new;	
           } // x2

        t1+=(mul_1_s1_0[counter_x1]+mul_1_s1_1[counter_x1])*(t1_temp_s2_0+t1_temp_s2_1);
        t4+=(mul_2_s1_0[counter_x1]+mul_2_s1_1[counter_x1])*(t4_temp_s2_0+t4_temp_s2_1);

	count_coul_ptr_o1+=o_max;

        counter_x1++;
        } //o1
 

for(int s1=0;s1<2;s1++){
  int s1f; if(s1==0) s1f = 1; else s1f = 0;
       for(int s2=0;s2<2;s2++){
         int s2f; if(s2==0) s2f = 1; else s2f = 0;



	for(int o1=0;o1<o_max;o1++){
		Complex_t ce1_m_s1_o1=cf_e(m,x1,s1,o1); 

                Complex_t mul3=conj(ce1_m_s1_o1)*conj(cf_e(n,x1,s2,o1));
		Complex_t mul4=conj(ce1_m_s1_o1)*cf_e(q,x1,s2,o1);
  		Complex_t mul5=conj(ce1_m_s1_o1)*cf_e(p,x1,s2f,o1);	
		Complex_t t2_temp=0.0;
		Complex_t t3_temp=0.0;
		Complex_t t5_temp=0.0;
		Complex_t t6_temp=0.0;

                for(int o2=0;o2<o_max;o2++){
                  if(o1!=o2) {
                    double exch_tmp=exch_on(o1,o2,_atom[x1]);   //w_exch(x1,x2,o1,o2);
			if(exch_tmp!=0.0) {
                     	
                      Complex_t ce2_p_s2f_o2=cf_e(p,x1,s2f,o2); //t5 
		      Complex_t ce1_p_s1_o2=cf_e(p,x1,s1,o2); //t2
                      
		      Complex_t ce2_n_s2_o2=cf_e(n,x1,s2,o2); //t3, t6
               
	              Complex_t ce2_q_s2_o2=cf_e(q,x1,s2,o2); //t2
		      Complex_t ce1_q_s1f_o2=cf_e(q,x1,s1f,o2); //t5
			
	              //coulomb integrals

		      //n1=n3 & n2=n4	
                      //t2+=conj(ce1_m_s1_o1)*ce1_p_s1_o2*conj(ce2_n_s2_o1)*ce2_q_s2_o2*exch_tmp;
		  	t2_temp+=ce1_p_s1_o2*ce2_q_s2_o2*exch_tmp;		
	
		      //n1=n4 & n2=n3 	
                      t3_temp+=ce1_p_s1_o2*conj(ce2_n_s2_o2)*exch_tmp;

		      //exchange integrals

		      //n1=n3 & n2=n4		
                      t5_temp+=ce1_q_s1f_o2*ce2_p_s2f_o2*exch_tmp;

		      //n1=n4 & n2=n3	
                      t6_temp+=ce1_q_s1f_o2*conj(ce2_n_s2_o2)*exch_tmp;

                    } // if(exch_tmp!=0.0)
                  } // if((x1!=x2)||(o1!=o2))   
             } // o2

		t2+=t2_temp*mul3;
		t3+=t3_temp*mul4;
		t5+=t5_temp*mul3;
		t6+=t6_temp*mul5;

       	} //o1
}}//s1 and s2
#ifdef VERBAL
if((fp!=NULL) && (x1%interval==0)){
        fprintf(fp,"Finished %i atoms from %i to %i in processor %i\n", x1, g_min, g_max, id);
	fflush(fp);
}
#endif

} // x1 


delete[] dot1_s2_0;
delete[] dot2_s2_0;
delete[] dot1_s2_1;
delete[] dot2_s2_1;

delete[] mul_1_s1_0;
delete[] mul_2_s1_0;
delete[] mul_1_s1_1;
delete[] mul_2_s1_1;

delete[] dot3_s2_0;
delete[] dot4_s2_0;
delete[] dot3_s2_1;
delete[] dot4_s2_1;

   cout<<"T1: "<<t1<<"   T2: "<<t2<<"   T3: "<<t3<<endl;
   cout<<"T4: "<<t4<<"   T5: "<<t5<<"   T6: "<<t6<<endl;
   	
   H_coul()=(t1+t2+t3)*inv_dielect;
   H_exch()=(t4+t5+t6)*inv_dielect;
#ifdef VERBAL
if(fp!=NULL)
	fclose(fp);
#endif
}


void Coulomb::compute_exchange_cutoff(int m, int n, int p, int q, int id, double cutoff)
{

//dump to file

#ifdef VERBAL 
	string b1;
	stringstream a1;
	a1 << id;
	b1 = a1.str();
	string filename="status_";
	filename+=b1;
	//cout<<filename<<endl;
	//const char *sr=filename.c_str();
	int interval=skip;
	FILE* fp=fopen(filename.c_str(),"w");
	if(fp==NULL){
		cout<<"Error opening file"<<endl;
	}
#endif

m--;
n--;
p--;
q--;

//<mn||pq>=[mp|nq]-[mq|np]
//<mn||mn>=<mn|mn>-<mn|nm>

H_coul()=0.0;
H_exch()=0.0;
   
int count_coul_ptr_o1=0;
int count_coul_ptr_o2=0;
int o_max_new=o_max-4;

Complex_t t1=0.0,t2=0.0,t3=0.0,t4=0.0, t5=0.0, t6=0.0;
     double tmp=0.0;
     double dd=0.0;
     double dd_cutoff=cutoff; //9.0;	
	int counter1=0;
	complex<double>* dot1_s2_0=new complex<double>[numAtoms*o_max_new];
        complex<double>* dot2_s2_0=new complex<double>[numAtoms*o_max_new];
	complex<double>* dot1_s2_1=new complex<double>[numAtoms*o_max_new];
        complex<double>* dot2_s2_1=new complex<double>[numAtoms*o_max_new];

	complex<double>* mul_1_s1_0=new complex<double>[numAtoms*o_max_new];
        complex<double>* mul_2_s1_0=new complex<double>[numAtoms*o_max_new];
	complex<double>* mul_1_s1_1=new complex<double>[numAtoms*o_max_new];
        complex<double>* mul_2_s1_1=new complex<double>[numAtoms*o_max_new];

        complex<double>* dot3_s2_0=new complex<double>[numAtoms];
        complex<double>* dot4_s2_0=new complex<double>[numAtoms];
	complex<double>* dot3_s2_1=new complex<double>[numAtoms];
        complex<double>* dot4_s2_1=new complex<double>[numAtoms];


	   for(int i=0;i<numAtoms;i++){
		for(int j=0;j<o_max;j++){
			if(j<5){
				dot1_s2_0[counter1]=conj(cf_e(n,i,0,j))*cf_e(q,i,0,j);
				dot2_s2_0[counter1]=conj(cf_e(n,i, 0,j))*cf_e(p, i, 1,j);
				mul_1_s1_0[counter1]=conj(cf_e(m,i,0,j))*cf_e(p,i,0,j);
                          	mul_2_s1_0[counter1]=conj(cf_e(m,i,0,j))*cf_e(q,i,1,j);
				dot1_s2_1[counter1]=conj(cf_e(n,i,1,j))*cf_e(q,i,1,j);
                                dot2_s2_1[counter1]=conj(cf_e(n,i, 1,j))*cf_e(p, i, 0,j);
                                mul_1_s1_1[counter1]=conj(cf_e(m,i,1,j))*cf_e(p,i,1,j);
                                mul_2_s1_1[counter1]=conj(cf_e(m,i,1,j))*cf_e(q,i,0,j);
				counter1++;
			}
			else{
				dot1_s2_0[counter1]+=conj(cf_e(n,i, 0,j))*cf_e(q,i,0,j);
				dot2_s2_0[counter1]+=conj(cf_e(n,i, 0,j))*cf_e(p, i, 1,j);
				mul_1_s1_0[counter1]+=conj(cf_e(m,i,0,j))*cf_e(p,i,0,j);
                          	mul_2_s1_0[counter1]+=conj(cf_e(m,i,0,j))*cf_e(q,i,1,j);
				dot1_s2_1[counter1]+=conj(cf_e(n,i, 1,j))*cf_e(q,i,1,j);
                                dot2_s2_1[counter1]+=conj(cf_e(n,i, 1,j))*cf_e(p, i, 0,j);
                                mul_1_s1_1[counter1]+=conj(cf_e(m,i,1,j))*cf_e(p,i,1,j);
                                mul_2_s1_1[counter1]+=conj(cf_e(m,i,1,j))*cf_e(q,i,0,j);

			}//end of if
		}//end of j
		for(int k=0;k<o_max_new;k++){
			dot3_s2_0[i]+=dot1_s2_0[i*o_max_new+k];
			dot4_s2_0[i]+=dot2_s2_0[i*o_max_new+k];
			dot3_s2_1[i]+=dot1_s2_1[i*o_max_new+k];
                        dot4_s2_1[i]+=dot2_s2_1[i*o_max_new+k];

		}
		counter1++;
	   }//end of atom index i	

        //cout<<"Entering x1 block"<<endl;
	int counter_x1=g_min*o_max_new;
	int counter=0;	
	for(int x1=g_min;x1<g_max;x1++){
	  	count_coul_ptr_o1=0;	
          for(int o1=0;o1<o_max_new;o1++){
    
		Complex_t t1_temp_s2_0=0.0;
		Complex_t t4_temp_s2_0=0.0;
		Complex_t t1_temp_s2_1=0.0;
                Complex_t t4_temp_s2_1=0.0;

		counter=0;
           for(int x2=0;x2<numAtoms;x2++){
		if(x1!=x2) 
   	        	dd = distance_sq(x1,x2);
		else
			dd=0.0;
		if(dd<dd_cutoff){	
		for(int o2=0;o2<o_max_new;o2++){
			if(x1!=x2){
				tmp = _coul_on_rec_sq[count_coul_ptr_o1+o2];
				tmp = 1.0/sqrt(tmp+dd*0.48268);
			}
			else{
				tmp =_coul_on_avg[count_coul_ptr_o1+o2];
			}
		t1_temp_s2_0+=dot1_s2_0[counter+o2]*tmp;
		t4_temp_s2_0+=dot2_s2_0[counter+o2]*tmp;
		t1_temp_s2_1+=dot1_s2_1[counter+o2]*tmp;
                t4_temp_s2_1+=dot2_s2_1[counter+o2]*tmp;

	      } //o2
	  }else{
		tmp=1.439363044/sqrt(dd);
		t1_temp_s2_0+=dot3_s2_0[x2]*tmp;
		t4_temp_s2_0+=dot4_s2_0[x2]*tmp;
		t1_temp_s2_1+=dot3_s2_1[x2]*tmp;
                t4_temp_s2_1+=dot4_s2_1[x2]*tmp;
	  }
	
	   counter+=o_max_new;	
           } // x2

        t1+=(mul_1_s1_0[counter_x1]+mul_1_s1_1[counter_x1])*(t1_temp_s2_0+t1_temp_s2_1);
        t4+=(mul_2_s1_0[counter_x1]+mul_2_s1_1[counter_x1])*(t4_temp_s2_0+t4_temp_s2_1);

	count_coul_ptr_o1+=o_max;

        counter_x1++;
        } //o1
 

for(int s1=0;s1<2;s1++){
  int s1f; if(s1==0) s1f = 1; else s1f = 0;
       for(int s2=0;s2<2;s2++){
         int s2f; if(s2==0) s2f = 1; else s2f = 0;



	for(int o1=0;o1<o_max;o1++){
		Complex_t ce1_m_s1_o1=cf_e(m,x1,s1,o1); 

                Complex_t mul3=conj(ce1_m_s1_o1)*conj(cf_e(n,x1,s2,o1));
		Complex_t mul4=conj(ce1_m_s1_o1)*cf_e(q,x1,s2,o1);
  		Complex_t mul5=conj(ce1_m_s1_o1)*cf_e(p,x1,s2f,o1);	
		Complex_t t2_temp=0.0;
		Complex_t t3_temp=0.0;
		Complex_t t5_temp=0.0;
		Complex_t t6_temp=0.0;

                for(int o2=0;o2<o_max;o2++){
                  if(o1!=o2) {
                    double exch_tmp=exch_on(o1,o2,_atom[x1]);   //w_exch(x1,x2,o1,o2);
			if(exch_tmp!=0.0) {
                     	
                      Complex_t ce2_p_s2f_o2=cf_e(p,x1,s2f,o2); //t5 
		      Complex_t ce1_p_s1_o2=cf_e(p,x1,s1,o2); //t2
                      
		      Complex_t ce2_n_s2_o2=cf_e(n,x1,s2,o2); //t3, t6
               
	              Complex_t ce2_q_s2_o2=cf_e(q,x1,s2,o2); //t2
		      Complex_t ce1_q_s1f_o2=cf_e(q,x1,s1f,o2); //t5
			
	              //coulomb integrals

		      //n1=n3 & n2=n4	
                      //t2+=conj(ce1_m_s1_o1)*ce1_p_s1_o2*conj(ce2_n_s2_o1)*ce2_q_s2_o2*exch_tmp;
		  	t2_temp+=ce1_p_s1_o2*ce2_q_s2_o2*exch_tmp;		
	
		      //n1=n4 & n2=n3 	
                      t3_temp+=ce1_p_s1_o2*conj(ce2_n_s2_o2)*exch_tmp;

		      //exchange integrals

		      //n1=n3 & n2=n4		
                      t5_temp+=ce1_q_s1f_o2*ce2_p_s2f_o2*exch_tmp;

		      //n1=n4 & n2=n3	
                      t6_temp+=ce1_q_s1f_o2*conj(ce2_n_s2_o2)*exch_tmp;

                    } // if(exch_tmp!=0.0)
                  } // if((x1!=x2)||(o1!=o2))   
             } // o2

		t2+=t2_temp*mul3;
		t3+=t3_temp*mul4;
		t5+=t5_temp*mul3;
		t6+=t6_temp*mul5;

       	} //o1
}}//s1 and s2
#ifdef VERBAL
if((fp!=NULL) && (x1%interval==0)){
        fprintf(fp,"Finished %i atoms from %i to %i in processor %i\n", x1, g_min, g_max, id);
	fflush(fp);
}
#endif

} // x1 


delete[] dot1_s2_0;
delete[] dot2_s2_0;
delete[] dot1_s2_1;
delete[] dot2_s2_1;

delete[] mul_1_s1_0;
delete[] mul_2_s1_0;
delete[] mul_1_s1_1;
delete[] mul_2_s1_1;

delete[] dot3_s2_0;
delete[] dot4_s2_0;
delete[] dot3_s2_1;
delete[] dot4_s2_1;

   cout<<"T1: "<<t1<<"   T2: "<<t2<<"   T3: "<<t3<<endl;
   cout<<"T4: "<<t4<<"   T5: "<<t5<<"   T6: "<<t6<<endl;
   	
   H_coul()=(t1+t2+t3)*inv_dielect;
   H_exch()=(t4+t5+t6)*inv_dielect;
#ifdef VERBAL
if(fp!=NULL)
	fclose(fp);
#endif
}


void Coulomb::compute_HL_exchange(int m, int n, int p, int q)
{

/*m--;
n--;
p--;
q--;

//<mn||pq>=[mp|nq]-[mq|np]
//<mn||mn>=<mn|mn>-<mn|nm>

H_coul()=0.0;
H_exch()=0.0;

int count_coul_ptr_o1=0;
int count_coul_ptr_o2=0;
dot1=new complex<double>[numAtoms*o_max];
dot2=new complex<double>[numAtoms*o_max];


Complex_t t1=0.0,t2=0.0,t3=0.0,t4=0.0, t5=0.0, t6=0.0;
     double tmp=0.0;
     double dd=0.0;
     for(int s1=0;s1<2;s1++){
       int s1f; if(s1==0) s1f = 1; else s1f = 0;
       for(int s2=0;s2<2;s2++){
         int s2f; if(s2==0) s2f = 1; else s2f = 0;
        int counter1=0;
           for(int i=0;i<numAtoms;i++){
                for(int j=0;j<o_max;j++){
                        dot1[counter1]=conj(cf_e(n,i, s2,j))*cf_e(q,i,s2,j);
                        dot2[counter1]=conj(cf_e(n,i, s2,j))*cf_e(p, i, s2f,j);
                        counter1++;
                }
           }
         for(int x1=g_min;x1<g_max;x1++){
          count_coul_ptr_o1=0;
          for(int o1=0;o1<o_max;o1++){

                Complex_t ce1_m_s1_o1=cf_e(m,x1,s1,o1);  //t1, t4, t3, t4, t5, t6
                Complex_t mul1=conj(ce1_m_s1_o1)*cf_e(p,x1,s1,o1);
                Complex_t mul2=conj(ce1_m_s1_o1)*cf_e(q,x1,s1f,o1);
                Complex_t t1_temp=0.0;
                Complex_t t4_temp=0.0;
                int counter=0;

           for(int x2=0;x2<numAtoms;x2++){
                if(x1!=x2)
                        dd = distance_sq(x1,x2);
                count_coul_ptr_o2=0;
                for(int o2=0;o2<o_max;o2++){
                        if(x1!=x2){
                                tmp = _coul_on_rec_sq[count_coul_ptr_o1+count_coul_ptr_o2];
                                tmp = 1.0/sqrt(tmp+dd*0.48268);
                        }
                        else
                                tmp =_coul_on_avg[count_coul_ptr_o1+count_coul_ptr_o2];


                // Coulomb integral
                // n1=n2 & n3=n4
                t1_temp+=dot1[counter]*tmp;

                // Exchange integral
                //n1=n2 & n3=n4
                t4_temp+=dot2[counter]*tmp;
    		counter++;
                count_coul_ptr_o2++;
              } //o2
           } // x2
        t1+=mul1*t1_temp;
        t4+=mul2*t4_temp;

                Complex_t mul3=conj(ce1_m_s1_o1)*conj(cf_e(n,x1,s2,o1));
                Complex_t mul4=conj(ce1_m_s1_o1)*cf_e(q,x1,s2,o1);
                Complex_t mul5=conj(ce1_m_s1_o1)*cf_e(p,x1,s2f,o1);
                Complex_t t2_temp=0.0;
                Complex_t t3_temp=0.0;
                Complex_t t5_temp=0.0;
                Complex_t t6_temp=0.0;

                for(int o2=0;o2<o_max;o2++){
                  if(o1!=o2) {
                    double exch_tmp=exch_on(o1,o2,_atom[x1]);   //w_exch(x1,x2,o1,o2);
                        if(exch_tmp!=0.0) {

                      Complex_t ce2_p_s2f_o2=cf_e(p,x1,s2f,o2); //t5
                      Complex_t ce1_p_s1_o2=cf_e(p,x1,s1,o2); //t2

                      Complex_t ce2_n_s2_o2=cf_e(n,x1,s2,o2); //t3, t6

                      Complex_t ce2_q_s2_o2=cf_e(q,x1,s2,o2); //t2
                      Complex_t ce1_q_s1f_o2=cf_e(q,x1,s1f,o2); //t5

                      //coulomb integrals

                      //n1=n3 & n2=n4
                      //t2+=conj(ce1_m_s1_o1)*ce1_p_s1_o2*conj(ce2_n_s2_o1)*ce2_q_s2_o2*exch_tmp;
                        t2_temp+=ce1_p_s1_o2*ce2_q_s2_o2*exch_tmp;

                      //n1=n4 & n2=n3
                     t3_temp+=ce1_p_s1_o2*conj(ce2_n_s2_o2)*exch_tmp;

                      //exchange integrals

                      //n1=n3 & n2=n4
                      t5_temp+=ce1_q_s1f_o2*ce2_p_s2f_o2*exch_tmp;

                      //n1=n4 & n2=n3
                      t6_temp+=ce1_q_s1f_o2*conj(ce2_n_s2_o2)*exch_tmp;

                    } // if(exch_tmp!=0.0)
                  } // if((x1!=x2)||(o1!=o2))
             } // o2

                t2+=t2_temp*mul3;
                t3+=t3_temp*mul4;
                t5+=t5_temp*mul3;
                t6+=t6_temp*mul5;

        count_coul_ptr_o1+=o_max;
        } //o1

} // x1
    }} // s1 and s2
   cout<<"T1: "<<t1<<"   T2: "<<t2<<"   T3: "<<t3<<endl;
   cout<<"T4: "<<t4<<"   T5: "<<t5<<"   T6: "<<t6<<endl;

   H_coul()=(t1+t2+t3)*inv_dielect;
   H_exch()=(t4+t5+t6)*inv_dielect;

*/
}


