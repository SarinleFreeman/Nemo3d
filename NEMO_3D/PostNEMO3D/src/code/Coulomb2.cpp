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
  
/*
  cout<<"electron wavefunction"<<endl;
   for(int g=g_min; g<g_max; g++) 
     for(int s=0; s<2; s++) 
         for(int o=0; o<o_max; o++)
            cout<<cf_e(0,g,s,o)<<endl; 
  // cout<<endl<<"hole wavefunction"<<endl;
 //  for(int g=g_min; g<g_max; g++) 
 //    for(int s=0; s<2; s++) 
 //        for(int o=0; o<o_max; o++)
  //          cout<<cf_h(0,g,s,o)<<endl;
   cout<<endl<<"atomic id"<<endl;
   for(int g=g_min; g<g_max; g++) 
      cout<<_atom[g]<<endl;
   cout<<endl<<"lattice point"<<endl;
   for(int g=g_min; g<g_max; g++)
      cout<<lattice(g,0)<<" "<<lattice(g,1)<<" "<<lattice(g,2)<<endl;
   cout<<endl<<"nearest neighbors"<<endl;
   for(int g=g_min; g<g_max; g++)
      cout<<nn(g,0)<<" "<<nn(g,1)<<" "<<nn(g,2)<<" "<<nn(g,3)<<endl;
*/


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
    double tmp = _coul_on_rec_sq[o1*o_max+o2];
    final = 1.0/sqrt(tmp+dd*0.48226876);
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

  if (g_min1 == g_min2) {

    for(int x=0;x<g_max1-g_min1;x++)
        for(int s=0;s<2;s++)
	    for(int o=0;o<o_max;o++)
		result+=conj(cf_e(wf1,x,s,o))*cf_e(wf2,x,s,o);
  } 
  *sum=result;
}

void Coulomb::dotProduct_Slater(Complex_t *m_p, Complex_t *n_q, Complex_t *n_p, Complex_t *m_q){
//computes the dot product beween two wavefunctions

  Complex_t mp=0.0;
  Complex_t nq=0.0;
  Complex_t np=0.0;
  Complex_t mq=0.0;
  int m=0;
  int n=1;
  int p=2;
  int q=3;
  if (g_min1 == g_min2) {
    for(int x=0;x<g_max1-g_min1;x++)
        for(int s=0;s<2;s++)
            for(int o=0;o<o_max;o++){
                mp+=conj(cf_e(m,x,s,o))*cf_e(p,x,s,o);
		nq+=conj(cf_e(n,x,s,o))*cf_e(q,x,s,o);
		np+=conj(cf_e(n,x,s,o))*cf_e(p,x,s,o);
		mq+=conj(cf_e(m,x,s,o))*cf_e(q,x,s,o);
	    }//end of orbital loop
  }// end of if
*m_p=mp; *n_q=nq; *n_p=np; *m_q=mq;
}

void Coulomb::coulomb_Slater(int p1, int p2, Complex_t *mVp, Complex_t *nVq, Complex_t *nVp, Complex_t *mVq){
//computes the dot product beween two wavefunctions

  Complex_t mp=0.0;
  Complex_t nq=0.0;
  Complex_t np=0.0;
  Complex_t mq=0.0;
  int m=0;
  int n=1;
  int p=2;
  int q=3;
  double dd1=0.0;
  double dd2=0.0;
  double V_r1=0.0;
  double V_r2=0.0;
  double V_constant=1.4399766/11.9;
  double U0=3.782; // eV
  double V_r=0.0;

  if (g_min1 == g_min2) {
    for(int x=0;x<g_max1-g_min1;x++){
	dd1=distance_sq(x,p2);
	dd2=distance_sq(x,p1);
	if(x!=p2)
		V_r1=-V_constant/sqrt(dd1);
	else
		V_r1=-U0;
	 if(x!=p1)
                V_r2=-V_constant/sqrt(dd2);
        else
                V_r2=-U0;
	V_r=V_r1+V_r2;	
        for(int s=0;s<2;s++)
            for(int o=0;o<o_max;o++){
                mp+=conj(cf_e(m,x,s,o))*cf_e(p,x,s,o)*V_r;
                nq+=conj(cf_e(n,x,s,o))*cf_e(q,x,s,o)*V_r;
                np+=conj(cf_e(n,x,s,o))*cf_e(p,x,s,o)*V_r;
                mq+=conj(cf_e(m,x,s,o))*cf_e(q,x,s,o)*V_r;
	 }//end of orbital loop	
	}//end of atom loop
  }//end of if
	*mVp=mp; *nVq=nq; *nVp=np; *mVq=mq;
}

void Coulomb::spinPolarization(int wf1, Complex_t *sum_up, Complex_t *sum_down){
//computes the up and down spin contribution in a wavefunction

  wf1--;
  Complex_t result_up=0.0;
  Complex_t result_down=0.0;	

  if (g_min1 == g_min2) {
    for(int x=0;x<g_max1-g_min1;x++)
            for(int o=0;o<o_max;o++){
                result_up+=conj(cf_e(wf1,x,0,o))*cf_e(wf1,x,0,o);
		result_down+=conj(cf_e(wf1,x,1,o))*cf_e(wf1,x,1,o);
	    }	
  }
  *sum_up=result_up;
  *sum_down=result_down;	
}


void Coulomb::getImpurityInfo(int wf){
/*int g1=-1;
wf--;
for(int i=g_min1;i<g_max1;i++)
	if(_atom[i]==1)
		g1=i;

if(g1!=-1){

cout<<"Atom Index of the impurity: "<<g1<<endl;
cout<<"Coordinates of the impurity:";
cout<<"("<<lattice(g1,0)<<","<<lattice(g1,1)<<","<<lattice(g1,2)<<")"<<endl<<endl;

cout<<"Orbitals at the impurity site for wavefunction "<<wf+1<<":"<<endl;
cout<<"s_up       :   "<<cf_e(wf,g1-g_min1,0,0)<<endl;
cout<<"px_up      :   "<<cf_e(wf,g1-g_min1,0,1)<<endl;
cout<<"py_up      :   "<<cf_e(wf,g1-g_min1,0,2)<<endl;
cout<<"pz_up      :   "<<cf_e(wf,g1-g_min1,0,3)<<endl;
cout<<"dxy_up     :   "<<cf_e(wf,g1-g_min1,0,4)<<endl;
cout<<"dyz_up     :   "<<cf_e(wf,g1-g_min1,0,5)<<endl;
cout<<"dzx_up     :   "<<cf_e(wf,g1-g_min1,0,6)<<endl;
cout<<"dx2-y2_up  :   "<<cf_e(wf,g1-g_min1,0,7)<<endl;
cout<<"dz2_up     :   "<<cf_e(wf,g1-g_min1,0,8)<<endl;
cout<<"s*_up      :   "<<cf_e(wf,g1-g_min1,0,9)<<endl;
cout<<"s_down     :   "<<cf_e(wf,g1-g_min1,1,0)<<endl;
cout<<"px_down    :   "<<cf_e(wf,g1-g_min1,1,1)<<endl;
cout<<"py_down    :   "<<cf_e(wf,g1-g_min1,1,2)<<endl;
cout<<"pz_down    :   "<<cf_e(wf,g1-g_min1,1,3)<<endl;
cout<<"dxy_down   :   "<<cf_e(wf,g1-g_min1,1,4)<<endl;
cout<<"dyz_down   :   "<<cf_e(wf,g1-g_min1,1,5)<<endl;
cout<<"dzx_down   :   "<<cf_e(wf,g1-g_min1,1,6)<<endl;
cout<<"dx2-y2_down:   "<<cf_e(wf,g1-g_min1,1,7)<<endl;
cout<<"dz2_down   :   "<<cf_e(wf,g1-g_min1,1,8)<<endl;
cout<<"s*_down    :   "<<cf_e(wf,g1-g_min1,1,9)<<endl<<endl;
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
*/
}//end of function


void Coulomb::compute_coulomb_exchange2(double dd_cutoff_sq, const char* base, const int* eList, const int NeList) {

 H_coul() = 0.0;
 H_exch() = 0.0;

 complex<double>*cf_e_second = new complex<double>[20*NeList*(g_max2-g_min2)];

 FILE *fp; int count = 0;
 char filename[200]; char c;
 
 // read additional data.

 for (int f=0; f<NeList; f++) {

    sprintf(filename, "%s.nd_evec_%d", base, eList[f]);
    fp = fopen(filename,"rb");
    while(1) {
    c=fgetc(fp);
    if(c=='\n') count++;
    else if(c=='>' && count == 6) break;
    }

    fseek(fp,g_min2*20*2*sizeof(double),SEEK_CUR);
    fread(cf_e_second + f*(g_max2 - g_min2)*20, 2*sizeof(double), (g_max2 - g_min2)*20, fp);
    fclose(fp);count=0;
 
 }
 
 // start of calclulation.
 Complex_t t1=0.0,t2=0.0,t3=0.0,t4=0.0, t5=0.0, t6=0.0;
 Complex_t c0101=0.0, c0111=0.0, c1201=0.0, c1211=0.0, cR101=0.0, cR111=0.0, cR201=0.0, cR211=0.0;
 Complex_t c1202=0.0, c1212=0.0, c0102=0.0, c0112=0.0, cR202=0.0, cR212=0.0, cR102=0.0, cR112=0.0;
 double t2_Real=0.0, t3_Real=0.0, t5_Real=0.0, t6_Real=0.0;
 double coul_dd=0.0, coul_pphetero=0.0, coul_pphomo=0.0, coul_pd=0.0, coul_sss=0.0, coul_ssp =0.0, coul_ssd=0.0, coul_sp=0.0, coul_sd=0.0, coul_ss=0.0, coul_ssss=0.0;
 int o_max_new = 4;

 double *t1_o1_Real = new double[(g_max1-g_min1)*o_max_new]; double *t4_o1_Real = new double[(g_max1-g_min1)*o_max_new];         
 double *t4_o1_Imag = new double[(g_max1-g_min1)*o_max_new]; double *t1_o2_Real = new double[(g_max2-g_min2)*o_max_new];
 double *t4_o2_Real = new double[(g_max2-g_min2)*o_max_new]; double *t4_o2_Imag = new double[(g_max2-g_min2)*o_max_new];
 
// allocation of pre-computed values.

 count = 0;
 for(int x1=g_min1;x1<g_max1;x1++) {

     for(int o1=0;o1<2;o1++) {

         c0101 = _cf_e[(x1-g_min1)*2*o_max+o1];                   c0111 = _cf_e[((x1-g_min1)*2+1)*o_max+o1];
         cR101 = _cf_e[((g_max1 - g_min1)+x1-g_min1)*2*o_max+o1]; cR111 = _cf_e[(((g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o1]; 
         
         t1_o1_Real[count]   = c0101._Val[0]*c0101._Val[0] + c0101._Val[1]*c0101._Val[1] + c0111._Val[0]*c0111._Val[0] + c0111._Val[1]*c0111._Val[1];
         t4_o1_Real[count]   = c0101._Val[0]*cR101._Val[0] + c0101._Val[1]*cR101._Val[1] + c0111._Val[0]*cR111._Val[0] + c0111._Val[1]*cR111._Val[1];
         t4_o1_Imag[count++] = c0101._Val[0]*cR101._Val[1] - c0101._Val[1]*cR101._Val[0] + c0111._Val[0]*cR111._Val[1] - c0111._Val[1]*cR111._Val[0];

     }

     for(int o1=2;o1<5;o1++) {

         c0101 = _cf_e[(x1-g_min1)*2*o_max+o1];                   c0111 = _cf_e[((x1-g_min1)*2+1)*o_max+o1];
         cR101 = _cf_e[((g_max1 - g_min1)+x1-g_min1)*2*o_max+o1]; cR111 = _cf_e[(((g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o1]; 
 
         t1_o1_Real[count] += c0101._Val[0]*c0101._Val[0] + c0101._Val[1]*c0101._Val[1] + c0111._Val[0]*c0111._Val[0] + c0111._Val[1]*c0111._Val[1];
         t4_o1_Real[count] += c0101._Val[0]*cR101._Val[0] + c0101._Val[1]*cR101._Val[1] + c0111._Val[0]*cR111._Val[0] + c0111._Val[1]*cR111._Val[1];
         t4_o1_Imag[count] += c0101._Val[0]*cR101._Val[1] - c0101._Val[1]*cR101._Val[0] + c0111._Val[0]*cR111._Val[1] - c0111._Val[1]*cR111._Val[0];

     }
     
     count++;
     for(int o1=5;o1<o_max;o1++) {

         c0101 = _cf_e[(x1-g_min1)*2*o_max+o1];                   c0111 = _cf_e[((x1-g_min1)*2+1)*o_max+o1];
         cR101 = _cf_e[((g_max1 - g_min1)+x1-g_min1)*2*o_max+o1]; cR111 = _cf_e[(((g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o1]; 
 
         t1_o1_Real[count] += c0101._Val[0]*c0101._Val[0] + c0101._Val[1]*c0101._Val[1] + c0111._Val[0]*c0111._Val[0] + c0111._Val[1]*c0111._Val[1];
         t4_o1_Real[count] += c0101._Val[0]*cR101._Val[0] + c0101._Val[1]*cR101._Val[1] + c0111._Val[0]*cR111._Val[0] + c0111._Val[1]*cR111._Val[1];
         t4_o1_Imag[count] += c0101._Val[0]*cR101._Val[1] - c0101._Val[1]*cR101._Val[0] + c0111._Val[0]*cR111._Val[1] - c0111._Val[1]*cR111._Val[0];

     }
     
     count++;

 }

 count=0; 
 for(int x2=g_min2;x2<g_max2;x2++) {

     for(int o1=0;o1<2;o1++) {
 
         c1201 = cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o1]; c1211 = cf_e_second[(((g_max2 - g_min2)+x2-g_min2)*2+1)*o_max+o1];
         cR201 = cf_e_second[(x2-g_min2)*2*o_max+o1];                 cR211 = cf_e_second[((x2-g_min2)*2+1)*o_max+o1];

         t1_o2_Real[count]   = c1201._Val[0]*c1201._Val[0] + c1201._Val[1]*c1201._Val[1] + c1211._Val[0]*c1211._Val[0] + c1211._Val[1]*c1211._Val[1];
         t4_o2_Real[count]   = c1201._Val[0]*cR201._Val[0] + c1201._Val[1]*cR201._Val[1] + c1211._Val[0]*cR211._Val[0] + c1211._Val[1]*cR211._Val[1];
         t4_o2_Imag[count++] = c1201._Val[0]*cR201._Val[1] - c1201._Val[1]*cR201._Val[0] + c1211._Val[0]*cR211._Val[1] - c1211._Val[1]*cR211._Val[0];

     }

    for(int o1=2;o1<5;o1++) {

         c1201 = cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o1]; c1211 = cf_e_second[(((g_max2 - g_min2)+x2-g_min2)*2+1)*o_max+o1];
         cR201 = cf_e_second[(x2-g_min2)*2*o_max+o1];                 cR211 = cf_e_second[((x2-g_min2)*2+1)*o_max+o1];

         t1_o2_Real[count] += c1201._Val[0]*c1201._Val[0] + c1201._Val[1]*c1201._Val[1] + c1211._Val[0]*c1211._Val[0] + c1211._Val[1]*c1211._Val[1];
         t4_o2_Real[count] += c1201._Val[0]*cR201._Val[0] + c1201._Val[1]*cR201._Val[1] + c1211._Val[0]*cR211._Val[0] + c1211._Val[1]*cR211._Val[1];
         t4_o2_Imag[count] += c1201._Val[0]*cR201._Val[1] - c1201._Val[1]*cR201._Val[0] + c1211._Val[0]*cR211._Val[1] - c1211._Val[1]*cR211._Val[0];

     }

     count++;
     for(int o1=5;o1<o_max;o1++) {

         c1201 = cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o1]; c1211 = cf_e_second[(((g_max2 - g_min2)+x2-g_min2)*2+1)*o_max+o1];
         cR201 = cf_e_second[(x2-g_min2)*2*o_max+o1];                 cR211 = cf_e_second[((x2-g_min2)*2+1)*o_max+o1];

         t1_o2_Real[count] += c1201._Val[0]*c1201._Val[0] + c1201._Val[1]*c1201._Val[1] + c1211._Val[0]*c1211._Val[0] + c1211._Val[1]*c1211._Val[1];
         t4_o2_Real[count] += c1201._Val[0]*cR201._Val[0] + c1201._Val[1]*cR201._Val[1] + c1211._Val[0]*cR211._Val[0] + c1211._Val[1]*cR211._Val[1];
         t4_o2_Imag[count] += c1201._Val[0]*cR201._Val[1] - c1201._Val[1]*cR201._Val[0] + c1211._Val[0]*cR211._Val[1] - c1211._Val[1]*cR211._Val[0];

     }

     count++;
 }

 // start actual calculations 
 count = 0; 
 for(int x1=g_min1;x1<g_max1;x1++){
         
     for(int x2=g_min2;x2<g_max2;x2++){ 
        
         double dd, coul_tmp;
         
         if(x1 == x2) dd=0.0;
         else  dd = distance_sq(x1,x2);
         
         if(dd < dd_cutoff_sq) {
   
            coul_dd = w_coul(dd,x1,5,5); coul_pphetero = w_coul(dd,x1,2,3); coul_pphomo = w_coul(dd,x1,2,2); coul_pd = w_coul(dd,x1,2,5);
            coul_sss = w_coul(dd,x1,0,1); coul_ssp = w_coul(dd,x1,0,2); coul_ssd = w_coul(dd,x1,0,5); coul_sp = w_coul(dd,x1,1,2);
            coul_sd = w_coul(dd,x1,1,5); coul_ss = w_coul(dd,x1,1,1); coul_ssss = w_coul(dd,x1,0,0);

            for(int o1=0;o1<o_max_new;o1++){

                int count1 = 4*(x1-g_min1)+o1;
                int count2 = 4*(x2-g_min2)+o1;
                //Complex_t t1_o1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(0,x1-g_min1,0,o1) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(0,x1-g_min1,1,o1);
                //Complex_t t4_o1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e_second_x1[((x1-g_min1)*2+0)*o_max+o1] + conj(cf_e(0,x1-g_min1,1,o1))*cf_e_second_x1[((x1-g_min1)*2+1)*o_max+o1];
                //Complex_t t1_o2_Reverse = conj(cf_e(1,x2-g_min2,0,o1))*cf_e(1,x2-g_min2,0,o1) + conj(cf_e(1,x2-g_min2,1,o1))*cf_e(1,x2-g_min2,1,o1);
                //Complex_t t4_o2_Reverse = conj(cf_e(1,x2-g_min2,0,o1))*cf_e_first_x2[((x2-g_min2)*2+0)*o_max+o1] + conj(cf_e(1,x2-g_min2,1,o1))*cf_e_first_x2[((x2-g_min2)*2+1)*o_max+o1];
 
                for(int o2=o1+1;o2<o_max_new;o2++){

                    int count3 = 4*(x2-g_min2)+o2;
                    int count4 = 4*(x1-g_min1)+o2;

                    coul_tmp = (2-o2)*(3-o2)/2*coul_sss + (3-o2)*(o2-1)*((1-o1)*coul_ssp + o1*coul_sp) + (o2-2)*(o2-1)/2*((1-o1)*(2-o1)/2*coul_ssd + o1*(2-o1)*coul_sd + o1*(o1-1)/2*coul_pd);

                    t1._Val[0] += (t1_o1_Real[count1]*t1_o2_Real[count3] + t1_o1_Real[count4]*t1_o2_Real[count2])*coul_tmp;
                    t4._Val[0] += (t4_o1_Real[count1]*t4_o2_Real[count3] - t4_o1_Imag[count1]*t4_o2_Imag[count3] + t4_o1_Real[count4]*t4_o2_Real[count2] - t4_o1_Imag[count4]*t4_o2_Imag[count2])*coul_tmp;
                    t4._Val[1] += (t4_o1_Imag[count1]*t4_o2_Real[count3] + t4_o1_Real[count1]*t4_o2_Imag[count3] + t4_o1_Imag[count4]*t4_o2_Real[count2] + t4_o1_Real[count4]*t4_o2_Imag[count2])*coul_tmp;

                    //Complex_t t1_o2 = conj(cf_e(1,x2-g_min2,0,o2))*cf_e(1,x2-g_min2,0,o2) + conj(cf_e(1,x2-g_min2,1,o2))*cf_e(1,x2-g_min2,1,o2);
                    //Complex_t t4_o2 = conj(cf_e(1,x2-g_min2,0,o2))*cf_e_first_x2[((x2-g_min2)*2+0)*o_max+o2] + conj(cf_e(1,x2-g_min2,1,o2))*cf_e_first_x2[((x2-g_min2)*2+1)*o_max+o2];
                    //Complex_t t1_o1_Reverse = conj(cf_e(0,x1-g_min1,0,o2))*cf_e(0,x1-g_min1,0,o2) + conj(cf_e(0,x1-g_min1,1,o2))*cf_e(0,x1-g_min1,1,o2);
                    //Complex_t t4_o1_Reverse = conj(cf_e(0,x1-g_min1,0,o2))*cf_e_second_x1[((x1-g_min1)*2+0)*o_max+o2] + conj(cf_e(0,x1-g_min1,1,o2))*cf_e_second_x1[((x1-g_min1)*2+1)*o_max+o2];

                    //t1+=(t1_o1*t1_o2 + t1_o1_Reverse*t1_o2*Reverse)*coul_tmp;
                    //t4+=(t4_o1*t4_o2 + t4_o1_Reverse*t4_o2_Reverse)*coul_tmp;
              
                } //o2

                coul_tmp = (o1-2)*(o1-1)*o1/6*coul_dd + (3-o1)*(o1-1)*o1/2*(coul_pphomo - coul_pphetero) + (3-o1)*(2-o1)*o1/2*coul_ss + (3-o1)*(2-o1)*(1-o1)/6*coul_ssss;
           
                if(o1==2) {

                   double temp_t1_o1_Real=0.0, temp_t1_o2_Real_Reverse=0.0, temp_t4_o1_Real=0.0, temp_t4_o1_Imag=0.0, temp_t4_o2_Real_Reverse=0.0, temp_t4_o2_Imag_Reverse=0.0;

                   for (int oindx = 2; oindx < 5; oindx++) {
 
                        c0101 = _cf_e[(x1-g_min1)*2*o_max+oindx];                       c0111 = _cf_e[((x1-g_min1)*2+1)*o_max+oindx];
                        c1201 = cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+oindx]; c1211 = cf_e_second[(((g_max2 - g_min2)+x2-g_min2)*2+1)*o_max+oindx];
                        cR101 = _cf_e[((g_max1-g_min1)+x1-g_min1)*2*o_max+oindx];       cR111 = _cf_e[(((g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+oindx];
                        cR201 = cf_e_second[(x2-g_min2)*2*o_max+oindx];                 cR211 = cf_e_second[((x2-g_min2)*2+1)*o_max+oindx];

                        temp_t1_o1_Real         = c0101._Val[0]*c0101._Val[0] + c0101._Val[1]*c0101._Val[1] + c0111._Val[0]*c0111._Val[0] + c0111._Val[1]*c0111._Val[1];
                        temp_t1_o2_Real_Reverse = c1201._Val[0]*c1201._Val[0] + c1201._Val[1]*c1201._Val[1] + c1211._Val[0]*c1211._Val[0] + c1211._Val[1]*c1211._Val[1];
                        temp_t4_o1_Real         = c0101._Val[0]*cR101._Val[0] + c0101._Val[1]*cR101._Val[1] + c0111._Val[0]*cR111._Val[0] + c0111._Val[1]*cR111._Val[1];
                        temp_t4_o1_Imag         = c0101._Val[0]*cR101._Val[1] - c0101._Val[1]*cR101._Val[0] + c0111._Val[0]*cR111._Val[1] - c0111._Val[1]*cR111._Val[0];
                        temp_t4_o2_Real_Reverse = c1201._Val[0]*cR201._Val[0] + c1201._Val[1]*cR201._Val[1] + c1211._Val[0]*cR211._Val[0] + c1211._Val[1]*cR211._Val[1];
                        temp_t4_o2_Imag_Reverse = c1201._Val[0]*cR201._Val[1] - c1201._Val[1]*cR201._Val[0] + c1211._Val[0]*cR211._Val[1] - c1211._Val[1]*cR211._Val[0];

                        t1._Val[0] += (temp_t1_o1_Real*temp_t1_o2_Real_Reverse)*coul_tmp;
                        t4._Val[0] += (temp_t4_o1_Real*temp_t4_o2_Real_Reverse - temp_t4_o1_Imag*temp_t4_o2_Imag_Reverse)*coul_tmp;
                        t4._Val[1] += (temp_t4_o1_Imag*temp_t4_o2_Real_Reverse + temp_t4_o1_Real*temp_t4_o2_Imag_Reverse)*coul_tmp;

                   }

                   coul_tmp = coul_pphetero;

                } 

                t1._Val[0] += (t1_o1_Real[count1]*t1_o2_Real[count2])*coul_tmp;
                t4._Val[0] += (t4_o1_Real[count1]*t4_o2_Real[count2] - t4_o1_Imag[count1]*t4_o2_Imag[count2])*coul_tmp;
                t4._Val[1] += (t4_o1_Imag[count1]*t4_o2_Real[count2] + t4_o1_Real[count1]*t4_o2_Imag[count2])*coul_tmp; 

                //t1 += t1_o1*t1_o2_Reverse*coul_tmp;
                //t4 += t4_o1*t4_o2_Reverse*coul_tmp;

            } // o1 
         } // if(dd<dd_cutoff_sq)

         else { // if (dd<dd_cutoff_sq)

              coul_tmp = 1.439363044/sqrt(dd);
              double t1_o1_Real_Cutoff=0.0, t1_o2_Real_Cutoff=0.0, t4_o1_Real_Cutoff=0.0, t4_o2_Real_Cutoff=0.0, t4_o1_Imag_Cutoff=0.0, t4_o2_Imag_Cutoff=0.0;

              for (int o1=0; o1<o_max_new; o1++) {

                   int count1 = 4*(x1-g_min1) + o1; int count2 = 4*(x2-g_min2) + o1;
                   
                   t1_o1_Real_Cutoff += t1_o1_Real[count1]; t1_o2_Real_Cutoff += t1_o2_Real[count2];
                   t4_o1_Real_Cutoff += t4_o1_Real[count1]; t4_o2_Real_Cutoff += t4_o2_Real[count2];
                   t4_o1_Imag_Cutoff += t4_o1_Imag[count1]; t4_o2_Imag_Cutoff += t4_o2_Imag[count2];
               
              }

              t1._Val[0] += (t1_o1_Real_Cutoff*t1_o2_Real_Cutoff)*coul_tmp;
              t4._Val[0] += (t4_o1_Real_Cutoff*t4_o2_Real_Cutoff - t4_o1_Imag_Cutoff*t4_o2_Imag_Cutoff)*coul_tmp;
              t4._Val[1] += (t4_o1_Imag_Cutoff*t4_o2_Real_Cutoff + t4_o1_Real_Cutoff*t4_o2_Imag_Cutoff)*coul_tmp;

         } // else to if (dd<dd_cutoff_sq)

     } //x2
                
     if(g_min1 == g_min2) {
   
        Complex_t c1101=0.0, c1102=0.0, c1111=0.0, c1112=0.0;
        Complex_t t2_1=0.0, t2_2=0.0, t5_1=0.0, t6_2=0.0;

        for(int o1=0;o1<o_max;o1++) {
    
            c0101 = _cf_e[(x1-g_min1)*2*o_max+o1];                   c0111 = _cf_e[((x1-g_min1)*2+1)*o_max+o1]; 
            c1101 = _cf_e[(((g_max1-g_min1)+x1-g_min1)*2)*o_max+o1]; c1111 = _cf_e[(((g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o1];

            for(int o2=o1+1;o2<o_max;o2++){
                                
                double exch_tmp=w_exch(x1,x1,o1,o2);
                 
                if(exch_tmp!=0.0) {              

                   c0102 = _cf_e[(x1-g_min1)*2*o_max+o2];                   c0112 = _cf_e[((x1-g_min1)*2+1)*o_max+o2];
                   c1102 = _cf_e[(((g_max1-g_min1)+x1-g_min1)*2)*o_max+o2]; c1112 = _cf_e[(((g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o2];
                   
                   t2_1._Val[0] = c0101._Val[0]*c0102._Val[0] + c0101._Val[1]*c0102._Val[1] + c0111._Val[0]*c0112._Val[0] + c0111._Val[1]*c0112._Val[1];
                   t2_1._Val[1] = c0101._Val[0]*c0102._Val[1] - c0101._Val[1]*c0102._Val[0] + c0111._Val[0]*c0112._Val[1] - c0111._Val[1]*c0112._Val[0];

                   t2_2._Val[0] = c1101._Val[0]*c1102._Val[0] + c1101._Val[1]*c1102._Val[1] + c1111._Val[0]*c1112._Val[0] + c1111._Val[1]*c1112._Val[1];
                   t2_2._Val[1] = c1101._Val[0]*c1102._Val[1] - c1101._Val[1]*c1102._Val[0] + c1111._Val[0]*c1112._Val[1] - c1111._Val[1]*c1112._Val[0];

                   t5_1._Val[0] = c0101._Val[0]*c1102._Val[0] + c0101._Val[1]*c1102._Val[1] + c0111._Val[0]*c1112._Val[0] + c0111._Val[1]*c1112._Val[1];
                   t5_1._Val[1] = c0101._Val[0]*c1102._Val[1] - c0101._Val[1]*c1102._Val[0] + c0111._Val[0]*c1112._Val[1] - c0111._Val[1]*c1112._Val[0];

                   t6_2._Val[0] = c1101._Val[0]*c0102._Val[0] + c1101._Val[1]*c0102._Val[1] + c1111._Val[0]*c0112._Val[0] + c1111._Val[1]*c0112._Val[1];
                   t6_2._Val[1] = c1101._Val[0]*c0102._Val[1] - c1101._Val[1]*c0102._Val[0] + c1111._Val[0]*c0112._Val[1] - c1111._Val[1]*c0112._Val[0];

                   t2._Val[0] += (t2_1._Val[0]*t2_2._Val[0] - t2_1._Val[1]*t2_2._Val[1])*exch_tmp*2;
                   t3._Val[0] += (t2_1._Val[0]*t2_2._Val[0] + t2_1._Val[1]*t2_2._Val[1])*exch_tmp*2;
                   t5._Val[0] += (t5_1._Val[0]*t5_1._Val[0] + t5_1._Val[1]*t5_1._Val[1] + t6_2._Val[0]*t6_2._Val[0] + t6_2._Val[1]*t6_2._Val[1])*exch_tmp;
                   t6._Val[0] += (t5_1._Val[0]*t6_2._Val[0] - t5_1._Val[1]*t6_2._Val[1])*exch_tmp*2;
                   
                
                   //Complex_t t2_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(0,x1-g_min1,0,o2) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(0,x1-g_min1,1,o2);
                   //Complex_t t2_2 = conj(cf_e(1,x1-g_min1,0,o1))*cf_e(1,x1-g_min1,0,o2) + conj(cf_e(1,x1-g_min1,1,o1))*cf_e(1,x1-g_min1,1,o2);
                   //Complex_t t5_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(1,x1-g_min1,0,o2) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(1,x1-g_min1,1,o2);
                   //Complex_t t6_2 = conj(cf_e(1,x1-g_min1,0,o1))*cf_e(0,x1-g_min1,0,o2) + conj(cf_e(1,x1-g_min1,1,o1))*cf_e(0,x1-g_min1,1,o2);

                   //t2 += (t2_1*t2_2 + conj(t2_1)*conj(t2_2))*exch_tmp;
                   //t3 += (t2_1*conj(t2_2) + conj(t2_1)*t2_2)*exch_tmp;
                   //t5 += (t5_1*conj(t5_1) + conj(t6_2)*t6_2)*exch_tmp;
                   //t6 += (t5_1*t6_2 + conj(t6_2)*conj(t5_1))*exch_tmp;
 
                } // if(exch_tmp!=0.0)
            } // o2 
        } // o1
     } // if(g_min1 == g_min2)                  
 } // x1
 
 delete[] cf_e_second;   
 delete[] t1_o1_Real;         delete[] t4_o1_Real;         
 delete[] t4_o1_Imag;         delete[] t1_o2_Real;         
 delete[] t4_o2_Real;         delete[] t4_o2_Imag; 

 //cout<<"T1: "<<t1<<"   T2:"<<t2<<"   T3:"<<t3<<endl;
 //cout<<"T4: "<<t4<<"   T5:"<<t5<<"   T6:"<<t6<<endl;
 
 H_coul()+=(t1+t2+t3)*inv_dielect;
 H_exch()-=(t4+t5+t6)*inv_dielect;

}

void Coulomb::compute_coulomb_exchange4(double dd_cutoff_sq, const char* base, const int* eList, const int NeList) {
    
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
double dielectric_f_exch=argument=11.9*q_r0/(sinh(q_r0));
double R=0.0; 


  H_coul() = 0.0;
  H_exch() = 0.0;
         
  FILE *fp; int count = 0;
  char filename[200]; char c;
  Complex_t t1=0.0,t2=0.0,t3=0.0,t4=0.0, t5=0.0, t6=0.0;
  Complex_t ce1=0.0, ce2=0.0, ce3=0.0, ce4=0.0, ce5=0.0, ce6=0.0;
  double coul_dd=0.0, coul_pphetero=0.0, coul_pphomo=0.0, coul_pd=0.0, coul_sss=0.0, coul_ssp =0.0, coul_ssd=0.0, coul_sp=0.0, coul_sd=0.0, coul_ss=0.0, coul_ssss=0.0;
  double coul_tmp=0.0; int o_max_eff = 4;

  complex<double>*cf_e_second = new complex<double>[2*o_max*NeList*(g_max2-g_min2)]; // 2*o_max for spds* with spin.

  complex<double>*t1_1 = new complex<double>[o_max_eff*(g_max1-g_min1)];
  complex<double>*t4_1 = new complex<double>[o_max_eff*(g_max1-g_min1)];
  complex<double>*t1_2 = new complex<double>[o_max_eff*(g_max2-g_min2)];
  complex<double>*t4_2 = new complex<double>[o_max_eff*(g_max2-g_min2)];

  // read additional data.
  for (int f=0; f<NeList; f++) {

     sprintf(filename, "%s.nd_evec_%d", base, eList[f]);
     fp = fopen(filename,"rb");
     while(1) {
       c=fgetc(fp);
       if(c=='\n') count++;
       else if(c=='>' && count == 6) break;
     }

     fseek(fp,g_min2*20*2*sizeof(double),SEEK_CUR);
     fread(cf_e_second + f*(g_max2 - g_min2)*20, 2*sizeof(double), (g_max2 - g_min2)*20, fp);
     fclose(fp);count=0;

  }

  // pre-allocation of values.

  count = 0;
  for(int x1=g_min1;x1<g_max1;x1++){

     for(int o=0;o<2;o++) {

        ce1._Val[0] = _cf_e[(x1-g_min1)*2*o_max+o]._Val[0];         ce1._Val[1] = (-1)*_cf_e[(x1-g_min1)*2*o_max+o]._Val[1]; 
        ce3._Val[0] = _cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[0];     ce3._Val[1] = (-1)*_cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[1];
        ce2 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce4 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];       
        ce5 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce6 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];

        t1_1[count] = ce1*ce2 + ce3*ce4;        
        t4_1[count++] = ce1*ce5 + ce3*ce6; 
      
     }

     for(int o=2;o<5;o++) {

        ce1._Val[0] = _cf_e[(x1-g_min1)*2*o_max+o]._Val[0];         ce1._Val[1] = (-1)*_cf_e[(x1-g_min1)*2*o_max+o]._Val[1];  
        ce3._Val[0] = _cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[0];     ce3._Val[1] = (-1)*_cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[1];
        ce2 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce4 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];        
        ce5 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce6 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];

        t1_1[count] += ce1*ce2 + ce3*ce4;
        t4_1[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

     for(int o=5;o<o_max;o++) {

        ce1._Val[0] = _cf_e[(x1-g_min1)*2*o_max+o]._Val[0];         ce1._Val[1] = (-1)*_cf_e[(x1-g_min1)*2*o_max+o]._Val[1];  
        ce3._Val[0] = _cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[0];     ce3._Val[1] = (-1)*_cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[1];
        ce2 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce4 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];        
        ce5 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce6 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];
  
        t1_1[count] += ce1*ce2 + ce3*ce4;
        t4_1[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

  }

  count = 0;
  for(int x2=g_min2;x2<g_max2;x2++){

     for(int o=0;o<2;o++) {

        ce1 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o]); ce2 = cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce3 = conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);
        ce4 = cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; ce5 = cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

        t1_2[count] = ce1*ce2 + ce3*ce4;  
        t4_2[count++] = ce1*ce5 + ce3*ce6;

     }

     for(int o=2;o<5;o++) {
  
        ce1 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o]); ce2 = cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce3 = conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);
        ce4 = cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; ce5 = cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

        t1_2[count] += ce1*ce2 + ce3*ce4;
        t4_2[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

     for(int o=5;o<o_max;o++) {
  
        ce1 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o]); ce2 = cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce3 = conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);
        ce4 = cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; ce5 = cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

        t1_2[count] += ce1*ce2 + ce3*ce4;
        t4_2[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

  }

  for(int x1=g_min1;x1<g_max1;x1++){
     for(int x2=g_min2;x2<g_max2;x2++){
              
         double dd = distance_sq(x1,x2);
	
	 R=sqrt(dd);
         if(R>r0_nm)
                dielectric_f=11.9;
         else{
		
                R_atomic=R/bohr_a0;
                argument=q_r0-q*R_atomic;
                dielectric_f=11.9*q_r0/((q*R_atomic) + sinh(argument));
		//dielectric_f=11.9;
         }


         if (dd < dd_cutoff_sq) {

            coul_dd = w_coul(dd,x1,5,5); coul_pphetero = w_coul(dd,x1,2,3); coul_pphomo = w_coul(dd,x1,2,2); coul_pd = w_coul(dd,x1,2,5);
            coul_sss = w_coul(dd,x1,0,1); coul_ssp = w_coul(dd,x1,0,2); coul_ssd = w_coul(dd,x1,0,5); coul_sp = w_coul(dd,x1,1,2);
            coul_sd = w_coul(dd,x1,1,5); coul_ss = w_coul(dd,x1,1,1); coul_ssss = w_coul(dd,x1,0,0);
              
            for(int o1=0;o1<o_max_eff;o1++){
    
               int count1 = 4*(x1-g_min1)+o1;
               int count2 = 4*(x2-g_min2)+o1;

               for(int o2=o1+1;o2<o_max_eff;o2++){

                  int count3 = 4*(x2-g_min2)+o2;
                  int count4 = 4*(x1-g_min1)+o2;

                  coul_tmp = (2-o2)*(3-o2)/2*coul_sss + (3-o2)*(o2-1)*((1-o1)*coul_ssp + o1*coul_sp) + (o2-2)*(o2-1)/2*((1-o1)*(2-o1)/2*coul_ssd + o1*(2-o1)*coul_sd + o1*(o1-1)/2*coul_pd);

                  //Complex_t t1_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(2,x1-g_min1,0,o1) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(2,x1-g_min1,1,o1);
                  //Complex_t t4_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(3,x1-g_min1,0,o1) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(3,x1-g_min1,1,o1);

                  //Complex_t t1_2 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o2])*cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o2]+conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2])*cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2];
                  //Complex_t t4_2 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o2])*cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o2]+conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2])*cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2];

                  t1+=(t1_1[count1]*t1_2[count3] + t1_1[count4]*t1_2[count2])*coul_tmp/dielectric_f;
                  t4+=(t4_1[count1]*t4_2[count3] + t4_1[count4]*t4_2[count2])*coul_tmp/dielectric_f;

               } //o2

               coul_tmp = (o1-2)*(o1-1)*o1/6*coul_dd + (3-o1)*(o1-1)*o1/2*(coul_pphomo - coul_pphetero) + (3-o1)*(2-o1)*o1/2*coul_ss + (3-o1)*(2-o1)*(1-o1)/6*coul_ssss;

               if(o1==2) {

                  Complex_t temp_t1_1=0.0, temp_t1_2=0.0, temp_t4_1=0.0, temp_t4_2=0.0;

                  for(int o = 2; o < 5; o++) {

                      ce1 = conj(cf_e(0,x1-g_min1,0,o)); ce2 = cf_e(2,x1-g_min1,0,o); ce3 = conj(cf_e(0,x1-g_min1,1,o));
                      ce4 = cf_e(2,x1-g_min1,1,o);       ce5 = cf_e(3,x1-g_min1,0,o); ce6 = cf_e(3,x1-g_min1,1,o);
                    
                      temp_t1_1 = ce1*ce2 + ce3*ce4;
                      temp_t4_1 = ce1*ce5 + ce3*ce6;

                      ce1 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o]); ce2 = cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce3 = conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);
                      ce4 = cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; ce5 = cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

                      temp_t1_2 = ce1*ce2 + ce3*ce4;
                      temp_t4_2 = ce1*ce5 + ce3*ce6;

                      t1 += temp_t1_1*temp_t1_2*coul_tmp/dielectric_f;
                      t4 += temp_t4_1*temp_t4_2*coul_tmp/dielectric_f;

	          }

                  coul_tmp = coul_pphetero;

               }

               t1+= t1_1[count1]*t1_2[count2]*coul_tmp/dielectric_f;
               t4+= t4_1[count1]*t4_2[count2]*coul_tmp/dielectric_f;

            } // o1

         } 

         else { // if (dd<dd_cutoff_sq)

              coul_tmp = 1.439363044/sqrt(dd);
              Complex_t t1_1_Cutoff=0.0, t1_2_Cutoff=0.0, t4_1_Cutoff=0.0, t4_2_Cutoff=0.0;

              for (int o1=0; o1<o_max_eff; o1++) {

                   int count1 = 4*(x1-g_min1) + o1; int count2 = 4*(x2-g_min2) + o1;

                   t1_1_Cutoff += t1_1[count1]; t1_2_Cutoff += t1_2[count2];
                   t4_1_Cutoff += t4_1[count1]; t4_2_Cutoff += t4_2[count2];

              }

              t1 += t1_1_Cutoff*t1_2_Cutoff*coul_tmp/dielectric_f;
              t4 += t4_1_Cutoff*t4_2_Cutoff*coul_tmp/dielectric_f;

         } // else to if (dd<dd_cutoff_sq)


     } // x2

     if (g_min1 == g_min2) {
        for(int o1=0;o1<o_max;o1++){
            for(int o2=0;o2<o_max;o2++){
                       
               double exch_tmp=w_exch(x1,x1,o1,o2);
                       
               if (exch_tmp!=0.0) {
                          
                   Complex_t t23_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(2,x1-g_min1,0,o2) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(2,x1-g_min1,1,o2);
                   Complex_t t2_2 = conj(cf_e(1,x1-g_min1,0,o1))*cf_e(3,x1-g_min1,0,o2) + conj(cf_e(1,x1-g_min1,1,o1))*cf_e(3,x1-g_min1,1,o2);
                   Complex_t t3_2 = conj(cf_e(1,x1-g_min1,0,o2))*cf_e(3,x1-g_min1,0,o1) + conj(cf_e(1,x1-g_min1,1,o2))*cf_e(3,x1-g_min1,1,o1);
                   Complex_t t56_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(3,x1-g_min1,0,o2) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(3,x1-g_min1,1,o2);
                   Complex_t t5_2 = conj(cf_e(1,x1-g_min1,0,o2))*cf_e(2,x1-g_min1,0,o1) + conj(cf_e(1,x1-g_min1,1,o2))*cf_e(2,x1-g_min1,1,o1);
                   Complex_t t6_2 = conj(cf_e(1,x1-g_min1,0,o1))*cf_e(2,x1-g_min1,0,o2) + conj(cf_e(1,x1-g_min1,1,o1))*cf_e(2,x1-g_min1,1,o2);
             
                   t2+=t23_1*t2_2*exch_tmp/dielectric_f_exch;
                   t3+=t23_1*t3_2*exch_tmp/dielectric_f_exch;
                   t5+=t56_1*t5_2*exch_tmp/dielectric_f_exch;
                   t6+=t56_1*t6_2*exch_tmp/dielectric_f_exch;

               } // if(exch_tmp!=0.0)
            } // o2 
         } // o1
      } // g_min1 == g_min2
   } // x1

   //cout<<"T1: "<<t1<<"   T2: "<<t2<<"   T3: "<<t3<<endl;
   //cout<<"T4: "<<t4<<"   T5: "<<t5<<"   T6: "<<t6<<endl;

   //H_coul()+=(t1+t2+t3)*inv_dielect;
   //H_exch()-=(t4+t5+t6)*inv_dielect;
   H_coul()+=(t1+t2+t3);
   H_exch()-=(t4+t5+t6);

   
   delete[] cf_e_second;
   delete[] t1_1; delete[] t1_2; delete[] t4_1; delete[] t4_2;

//   cout<<"ne, nh, Coulomb, Exchange: "<<ne<<" "<<nh<<" "<<H_coul(ne,nh)<<" "<<H_exch(ne,nh)<<endl;

}


void Coulomb::compute_coulomb_exchange4_no_dielectric_correction(double dd_cutoff_sq, const char* base, const int* eList, const int NeList) {
    
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
//double dielectric_f_exch=argument=11.9*q_r0/(sinh(q_r0));
double dielectric_f_exch=argument=11.9;
double R=0.0; 


  H_coul() = 0.0;
  H_exch() = 0.0;
         
  FILE *fp; int count = 0;
  char filename[200]; char c;
  Complex_t t1=0.0,t2=0.0,t3=0.0,t4=0.0, t5=0.0, t6=0.0;
  Complex_t ce1=0.0, ce2=0.0, ce3=0.0, ce4=0.0, ce5=0.0, ce6=0.0;
  double coul_dd=0.0, coul_pphetero=0.0, coul_pphomo=0.0, coul_pd=0.0, coul_sss=0.0, coul_ssp =0.0, coul_ssd=0.0, coul_sp=0.0, coul_sd=0.0, coul_ss=0.0, coul_ssss=0.0;
  double coul_tmp=0.0; int o_max_eff = 4;

  complex<double>*cf_e_second = new complex<double>[2*o_max*NeList*(g_max2-g_min2)]; // 2*o_max for spds* with spin.

  complex<double>*t1_1 = new complex<double>[o_max_eff*(g_max1-g_min1)];
  complex<double>*t4_1 = new complex<double>[o_max_eff*(g_max1-g_min1)];
  complex<double>*t1_2 = new complex<double>[o_max_eff*(g_max2-g_min2)];
  complex<double>*t4_2 = new complex<double>[o_max_eff*(g_max2-g_min2)];

  // read additional data.
  for (int f=0; f<NeList; f++) {

     sprintf(filename, "%s.nd_evec_%d", base, eList[f]);
     fp = fopen(filename,"rb");
     while(1) {
       c=fgetc(fp);
       if(c=='\n') count++;
       else if(c=='>' && count == 6) break;
     }

     fseek(fp,g_min2*20*2*sizeof(double),SEEK_CUR);
     fread(cf_e_second + f*(g_max2 - g_min2)*20, 2*sizeof(double), (g_max2 - g_min2)*20, fp);
     fclose(fp);count=0;

  }

  // pre-allocation of values.

  count = 0;
  for(int x1=g_min1;x1<g_max1;x1++){

     for(int o=0;o<2;o++) {

        ce1._Val[0] = _cf_e[(x1-g_min1)*2*o_max+o]._Val[0];         ce1._Val[1] = (-1)*_cf_e[(x1-g_min1)*2*o_max+o]._Val[1]; 
        ce3._Val[0] = _cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[0];     ce3._Val[1] = (-1)*_cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[1];
        ce2 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce4 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];       
        ce5 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce6 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];

        t1_1[count] = ce1*ce2 + ce3*ce4;        
        t4_1[count++] = ce1*ce5 + ce3*ce6; 
      
     }

     for(int o=2;o<5;o++) {

        ce1._Val[0] = _cf_e[(x1-g_min1)*2*o_max+o]._Val[0];         ce1._Val[1] = (-1)*_cf_e[(x1-g_min1)*2*o_max+o]._Val[1];  
        ce3._Val[0] = _cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[0];     ce3._Val[1] = (-1)*_cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[1];
        ce2 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce4 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];        
        ce5 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce6 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];

        t1_1[count] += ce1*ce2 + ce3*ce4;
        t4_1[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

     for(int o=5;o<o_max;o++) {

        ce1._Val[0] = _cf_e[(x1-g_min1)*2*o_max+o]._Val[0];         ce1._Val[1] = (-1)*_cf_e[(x1-g_min1)*2*o_max+o]._Val[1];  
        ce3._Val[0] = _cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[0];     ce3._Val[1] = (-1)*_cf_e[((x1-g_min1)*2+1)*o_max+o]._Val[1];
        ce2 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce4 = _cf_e[((2*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];        
        ce5 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce6 = _cf_e[((3*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];
  
        t1_1[count] += ce1*ce2 + ce3*ce4;
        t4_1[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

  }

  count = 0;
  for(int x2=g_min2;x2<g_max2;x2++){

     for(int o=0;o<2;o++) {

        ce1 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o]); ce2 = cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce3 = conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);
        ce4 = cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; ce5 = cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

        t1_2[count] = ce1*ce2 + ce3*ce4;  
        t4_2[count++] = ce1*ce5 + ce3*ce6;

     }

     for(int o=2;o<5;o++) {
  
        ce1 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o]); ce2 = cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce3 = conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);
        ce4 = cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; ce5 = cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

        t1_2[count] += ce1*ce2 + ce3*ce4;
        t4_2[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

     for(int o=5;o<o_max;o++) {
  
        ce1 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o]); ce2 = cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce3 = conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);
        ce4 = cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; ce5 = cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

        t1_2[count] += ce1*ce2 + ce3*ce4;
        t4_2[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

  }

  for(int x1=g_min1;x1<g_max1;x1++){
     for(int x2=g_min2;x2<g_max2;x2++){
              
         double dd = distance_sq(x1,x2);
	
	 R=sqrt(dd);
/*         if(R>r0_nm)
                dielectric_f=11.9;
         else{
		
                R_atomic=R/bohr_a0;
                argument=q_r0-q*R_atomic;
                dielectric_f=11.9*q_r0/((q*R_atomic) + sinh(argument));
		//dielectric_f=11.9;
         }
*/
         if (dd < dd_cutoff_sq) {

            coul_dd = w_coul(dd,x1,5,5); coul_pphetero = w_coul(dd,x1,2,3); coul_pphomo = w_coul(dd,x1,2,2); coul_pd = w_coul(dd,x1,2,5);
            coul_sss = w_coul(dd,x1,0,1); coul_ssp = w_coul(dd,x1,0,2); coul_ssd = w_coul(dd,x1,0,5); coul_sp = w_coul(dd,x1,1,2);
            coul_sd = w_coul(dd,x1,1,5); coul_ss = w_coul(dd,x1,1,1); coul_ssss = w_coul(dd,x1,0,0);
              
            for(int o1=0;o1<o_max_eff;o1++){
    
               int count1 = 4*(x1-g_min1)+o1;
               int count2 = 4*(x2-g_min2)+o1;

               for(int o2=o1+1;o2<o_max_eff;o2++){

                  int count3 = 4*(x2-g_min2)+o2;
                  int count4 = 4*(x1-g_min1)+o2;

                  coul_tmp = (2-o2)*(3-o2)/2*coul_sss + (3-o2)*(o2-1)*((1-o1)*coul_ssp + o1*coul_sp) + (o2-2)*(o2-1)/2*((1-o1)*(2-o1)/2*coul_ssd + o1*(2-o1)*coul_sd + o1*(o1-1)/2*coul_pd);

                  //Complex_t t1_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(2,x1-g_min1,0,o1) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(2,x1-g_min1,1,o1);
                  //Complex_t t4_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(3,x1-g_min1,0,o1) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(3,x1-g_min1,1,o1);

                  //Complex_t t1_2 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o2])*cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o2]+conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2])*cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2];
                  //Complex_t t4_2 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o2])*cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o2]+conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2])*cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2];

                  t1+=(t1_1[count1]*t1_2[count3] + t1_1[count4]*t1_2[count2])*coul_tmp;
                  t4+=(t4_1[count1]*t4_2[count3] + t4_1[count4]*t4_2[count2])*coul_tmp;

               } //o2

               coul_tmp = (o1-2)*(o1-1)*o1/6*coul_dd + (3-o1)*(o1-1)*o1/2*(coul_pphomo - coul_pphetero) + (3-o1)*(2-o1)*o1/2*coul_ss + (3-o1)*(2-o1)*(1-o1)/6*coul_ssss;

               if(o1==2) {

                  Complex_t temp_t1_1=0.0, temp_t1_2=0.0, temp_t4_1=0.0, temp_t4_2=0.0;

                  for(int o = 2; o < 5; o++) {

                      ce1 = conj(cf_e(0,x1-g_min1,0,o)); ce2 = cf_e(2,x1-g_min1,0,o); ce3 = conj(cf_e(0,x1-g_min1,1,o));
                      ce4 = cf_e(2,x1-g_min1,1,o);       ce5 = cf_e(3,x1-g_min1,0,o); ce6 = cf_e(3,x1-g_min1,1,o);
                    
                      temp_t1_1 = ce1*ce2 + ce3*ce4;
                      temp_t4_1 = ce1*ce5 + ce3*ce6;

                      ce1 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o]); ce2 = cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce3 = conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);
                      ce4 = cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; ce5 = cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

                      temp_t1_2 = ce1*ce2 + ce3*ce4;
                      temp_t4_2 = ce1*ce5 + ce3*ce6;

                      t1 += temp_t1_1*temp_t1_2*coul_tmp;
                      t4 += temp_t4_1*temp_t4_2*coul_tmp;

	          }

                  coul_tmp = coul_pphetero;

               }

               t1+= t1_1[count1]*t1_2[count2]*coul_tmp;
               t4+= t4_1[count1]*t4_2[count2]*coul_tmp;

            } // o1

         } 

         else { // if (dd<dd_cutoff_sq)

              coul_tmp = 1.439363044/sqrt(dd);
              Complex_t t1_1_Cutoff=0.0, t1_2_Cutoff=0.0, t4_1_Cutoff=0.0, t4_2_Cutoff=0.0;

              for (int o1=0; o1<o_max_eff; o1++) {

                   int count1 = 4*(x1-g_min1) + o1; int count2 = 4*(x2-g_min2) + o1;

                   t1_1_Cutoff += t1_1[count1]; t1_2_Cutoff += t1_2[count2];
                   t4_1_Cutoff += t4_1[count1]; t4_2_Cutoff += t4_2[count2];

              }

              t1 += t1_1_Cutoff*t1_2_Cutoff*coul_tmp;
              t4 += t4_1_Cutoff*t4_2_Cutoff*coul_tmp;

         } // else to if (dd<dd_cutoff_sq)


     } // x2

     if (g_min1 == g_min2) {
        for(int o1=0;o1<o_max;o1++){
            for(int o2=0;o2<o_max;o2++){
                       
               double exch_tmp=w_exch(x1,x1,o1,o2);
                       
               if (exch_tmp!=0.0) {
                          
                   Complex_t t23_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(2,x1-g_min1,0,o2) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(2,x1-g_min1,1,o2);
                   Complex_t t2_2 = conj(cf_e(1,x1-g_min1,0,o1))*cf_e(3,x1-g_min1,0,o2) + conj(cf_e(1,x1-g_min1,1,o1))*cf_e(3,x1-g_min1,1,o2);
                   Complex_t t3_2 = conj(cf_e(1,x1-g_min1,0,o2))*cf_e(3,x1-g_min1,0,o1) + conj(cf_e(1,x1-g_min1,1,o2))*cf_e(3,x1-g_min1,1,o1);
                   Complex_t t56_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(3,x1-g_min1,0,o2) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(3,x1-g_min1,1,o2);
                   Complex_t t5_2 = conj(cf_e(1,x1-g_min1,0,o2))*cf_e(2,x1-g_min1,0,o1) + conj(cf_e(1,x1-g_min1,1,o2))*cf_e(2,x1-g_min1,1,o1);
                   Complex_t t6_2 = conj(cf_e(1,x1-g_min1,0,o1))*cf_e(2,x1-g_min1,0,o2) + conj(cf_e(1,x1-g_min1,1,o1))*cf_e(2,x1-g_min1,1,o2);
             
                   t2+=t23_1*t2_2*exch_tmp;
                   t3+=t23_1*t3_2*exch_tmp;
                   t5+=t56_1*t5_2*exch_tmp;
                   t6+=t56_1*t6_2*exch_tmp;

               } // if(exch_tmp!=0.0)
            } // o2 
         } // o1
      } // g_min1 == g_min2
   } // x1

   //cout<<"T1: "<<t1<<"   T2: "<<t2<<"   T3: "<<t3<<endl;
   //cout<<"T4: "<<t4<<"   T5: "<<t5<<"   T6: "<<t6<<endl;

   H_coul()+=(t1+t2+t3)*inv_dielect;
   H_exch()-=(t4+t5+t6)*inv_dielect;
   //H_coul()+=(t1+t2+t3);
   //H_exch()-=(t4+t5+t6);

   
   delete[] cf_e_second;
   delete[] t1_1; delete[] t1_2; delete[] t4_1; delete[] t4_2;

//   cout<<"ne, nh, Coulomb, Exchange: "<<ne<<" "<<nh<<" "<<H_coul(ne,nh)<<" "<<H_exch(ne,nh)<<endl;

}


void Coulomb::hartree_matrix(double dd_cutoff_sq, const char* base, const int* eList, int NeList) {

         
  FILE *fp; int count = 0;
  char filename[200]; char c;
  int o_max_eff = 4;

  complex<double>*cf_e_second = new complex<double>[2*o_max*NeList*(g_max2-g_min2)]; // 2*o_max for spds* with spin.

  complex<double>*t1_1 = new complex<double>[o_max_eff*(g_max1-g_min1)];
  complex<double>*t4_1 = new complex<double>[o_max_eff*(g_max1-g_min1)];
  complex<double>*t1_2 = new complex<double>[o_max_eff*(g_max2-g_min2)];
  complex<double>*t4_2 = new complex<double>[o_max_eff*(g_max2-g_min2)];

  // read additional data.
  for (int f=0; f<NeList; f++) {

     sprintf(filename, "%s.nd_evec_%d", base, eList[f]);
     fp = fopen(filename,"rb");
     while(1) {
       c=fgetc(fp);
       if(c=='\n') count++;
       else if(c=='>' && count == 6) break;
     }

     fseek(fp,g_min2*20*2*sizeof(double),SEEK_CUR);
     fread(cf_e_second + f*(g_max2 - g_min2)*20, 2*sizeof(double), (g_max2 - g_min2)*20, fp);
     fclose(fp);count=0;

  }

   count = 0;

   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 1, 0, 1,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 1, 0, 2,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 1, 0, 3,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 1, 1, 2,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 1, 1, 3,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 1, 2, 3,count);count++;

   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 2, 0, 2,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 2, 0, 3,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 2, 1, 2,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 2, 1, 3,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 2, 2, 3,count);count++;

   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 3, 0, 3,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 3, 1, 2,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 3, 1, 3,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 0, 3, 2, 3,count);count++;

   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 1, 2, 1, 2,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 1, 2, 1, 3,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 1, 2, 2, 3,count);count++;

   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 1, 3, 1, 3,count);count++;
   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 1, 3, 2, 3,count);count++;

   compute_one_component(dd_cutoff_sq, cf_e_second, t1_1, t4_1, t1_2, t4_2, 2, 3, 2, 3,count);count++;


   delete[] cf_e_second;
   delete[] t1_1; delete[] t1_2; delete[] t4_1; delete[] t4_2;


}

void Coulomb::compute_one_component(double dd_cutoff_sq, complex<double>*cf_e_second, complex<double>*t1_1, complex<double>*t4_1, complex<double>*t1_2, complex<double>*t4_2, int N1, int N2, int N3, int N4, int index) {

  int o_max_eff = 4;
  coul_Matrix[index] = 0.0; exch_Matrix[index] = 0.0;

  // pre-allocation of values.
  int count = 0;
  Complex_t t1=0.0,t2=0.0,t3=0.0,t4=0.0, t5=0.0, t6=0.0;
  Complex_t ce1=0.0, ce2=0.0, ce3=0.0, ce4=0.0, ce5=0.0, ce6=0.0;
  double coul_dd=0.0, coul_pphetero=0.0, coul_pphomo=0.0, coul_pd=0.0, coul_sss=0.0, coul_ssp =0.0, coul_ssd=0.0, coul_sp=0.0, coul_sd=0.0, coul_ss=0.0, coul_ssss=0.0;
  double coul_tmp=0.0;
 
  count = 0;
  for(int x1=g_min1;x1<g_max1;x1++){

     for(int o=0;o<2;o++) {

        t1_1[count] = 0.0;  t4_1[count] = 0.0;
        ce1._Val[0] = _cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o]._Val[0];       ce1._Val[1] = (-1)*_cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o]._Val[1]; 
        ce3._Val[0] = _cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o]._Val[0];     ce3._Val[1] = (-1)*_cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o]._Val[1];
        ce2 = _cf_e[((N3*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce4 = _cf_e[((N3*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];       
        ce5 = _cf_e[((N4*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];    ce6 = _cf_e[((N4*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];

        t1_1[count] = ce1*ce2 + ce3*ce4;        
        t4_1[count++] = ce1*ce5 + ce3*ce6; 
      
     }

     t1_1[count] = 0.0;  t4_1[count] = 0.0;
     for(int o=2;o<5;o++) {

        ce1._Val[0] = _cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o]._Val[0];       ce1._Val[1] = (-1)*_cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o]._Val[1];
        ce3._Val[0] = _cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o]._Val[0];     ce3._Val[1] = (-1)*_cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o]._Val[1];        
        ce2 = _cf_e[((N3*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce4 = _cf_e[((N3*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];        
        ce5 = _cf_e[((N4*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];    ce6 = _cf_e[((N4*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];

        t1_1[count] += ce1*ce2 + ce3*ce4;
        t4_1[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

     t1_1[count] = 0.0;  t4_1[count] = 0.0;
     for(int o=5;o<o_max;o++) {

        ce1._Val[0] = _cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o]._Val[0];       ce1._Val[1] = (-1)*_cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o]._Val[1];
        ce3._Val[0] = _cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o]._Val[0];     ce3._Val[1] = (-1)*_cf_e[((N1*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o]._Val[1];        
        ce2 = _cf_e[((N3*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];   ce4 = _cf_e[((N3*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];        
        ce5 = _cf_e[((N4*(g_max1 - g_min1)+x1-g_min1)*2)*o_max+o];    ce6 = _cf_e[((N4*(g_max1 - g_min1)+x1-g_min1)*2+1)*o_max+o];

        t1_1[count] += ce1*ce2 + ce3*ce4;
        t4_1[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

  }

  count = 0;
  for(int x2=g_min2;x2<g_max2;x2++){

     
     for(int o=0;o<2;o++) {

        t1_2[count] = 0.0;  t4_2[count] = 0.0;
        ce1 = conj(cf_e_second[(N2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]);      ce2 = cf_e_second[(N4*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; 
        ce3 = conj(cf_e_second[((N2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);  ce4 = cf_e_second[((N4*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; 
        ce5 = cf_e_second[(N3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((N3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

        t1_2[count] = ce1*ce2 + ce3*ce4;  
        t4_2[count++] = ce1*ce5 + ce3*ce6;

     }

     t1_2[count] = 0.0;  t4_2[count] = 0.0;
     for(int o=2;o<5;o++) {
  
        ce1 = conj(cf_e_second[(N2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]);       ce2 = cf_e_second[(N4*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; 
        ce3 = conj(cf_e_second[((N2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);   ce4 = cf_e_second[((N4*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; 
        ce5 = cf_e_second[(N3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((N3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

        t1_2[count] += ce1*ce2 + ce3*ce4;
        t4_2[count] += ce1*ce5 + ce3*ce6;

     }
     count++;
     t1_2[count] = 0.0;  t4_2[count] = 0.0;
     for(int o=5;o<o_max;o++) {
  
        ce1 = conj(cf_e_second[(N2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]);       ce2 = cf_e_second[(N4*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; 
        ce3 = conj(cf_e_second[((N2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);   ce4 = cf_e_second[((N4*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; 
        ce5 = cf_e_second[(N3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; ce6 = cf_e_second[((N3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

        t1_2[count] += ce1*ce2 + ce3*ce4;
        t4_2[count] += ce1*ce5 + ce3*ce6;

     }
     count++;

  }

  for(int x1=g_min1;x1<g_max1;x1++){
     for(int x2=g_min2;x2<g_max2;x2++){
              
         double dd = distance_sq(x1,x2);

         if (dd < dd_cutoff_sq) {

            coul_dd = w_coul(dd,x1,5,5); coul_pphetero = w_coul(dd,x1,2,3); coul_pphomo = w_coul(dd,x1,2,2); coul_pd = w_coul(dd,x1,2,5);
            coul_sss = w_coul(dd,x1,0,1); coul_ssp = w_coul(dd,x1,0,2); coul_ssd = w_coul(dd,x1,0,5); coul_sp = w_coul(dd,x1,1,2);
            coul_sd = w_coul(dd,x1,1,5); coul_ss = w_coul(dd,x1,1,1); coul_ssss = w_coul(dd,x1,0,0);
              
            for(int o1=0;o1<o_max_eff;o1++){
    
               int count1 = 4*(x1-g_min1)+o1;
               int count2 = 4*(x2-g_min2)+o1;

               for(int o2=o1+1;o2<o_max_eff;o2++){

                  int count3 = 4*(x2-g_min2)+o2;
                  int count4 = 4*(x1-g_min1)+o2;

                  coul_tmp = (2-o2)*(3-o2)/2*coul_sss + (3-o2)*(o2-1)*((1-o1)*coul_ssp + o1*coul_sp) + (o2-2)*(o2-1)/2*((1-o1)*(2-o1)/2*coul_ssd + o1*(2-o1)*coul_sd + o1*(o1-1)/2*coul_pd);

                  //Complex_t t1_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(2,x1-g_min1,0,o1) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(2,x1-g_min1,1,o1);
                  //Complex_t t4_1 = conj(cf_e(0,x1-g_min1,0,o1))*cf_e(3,x1-g_min1,0,o1) + conj(cf_e(0,x1-g_min1,1,o1))*cf_e(3,x1-g_min1,1,o1);

                  //Complex_t t1_2 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o2])*cf_e_second[(3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o2]+conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2])*cf_e_second[((3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2];
                  //Complex_t t4_2 = conj(cf_e_second[((g_max2-g_min2)+x2-g_min2)*2*o_max+o2])*cf_e_second[(2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o2]+conj(cf_e_second[(((g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2])*cf_e_second[((2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o2];

                  t1+=(t1_1[count1]*t1_2[count3] + t1_1[count4]*t1_2[count2])*coul_tmp;
                  t4+=(t4_1[count1]*t4_2[count3] + t4_1[count4]*t4_2[count2])*coul_tmp;

               } //o2

               coul_tmp = (o1-2)*(o1-1)*o1/6*coul_dd + (3-o1)*(o1-1)*o1/2*(coul_pphomo - coul_pphetero) + (3-o1)*(2-o1)*o1/2*coul_ss + (3-o1)*(2-o1)*(1-o1)/6*coul_ssss;

               if(o1==2) {

                  Complex_t temp_t1_1=0.0, temp_t1_2=0.0, temp_t4_1=0.0, temp_t4_2=0.0;

                  for(int o = 2; o < 5; o++) {

                      ce1 = conj(cf_e(N1,x1-g_min1,0,o)); ce2 = cf_e(N3,x1-g_min1,0,o); ce3 = conj(cf_e(N1,x1-g_min1,1,o));
                      ce4 = cf_e(N3,x1-g_min1,1,o);       ce5 = cf_e(N4,x1-g_min1,0,o); ce6 = cf_e(N4,x1-g_min1,1,o);
                    
                      temp_t1_1 = ce1*ce2 + ce3*ce4;
                      temp_t4_1 = ce1*ce5 + ce3*ce6;

                      ce1 = conj(cf_e_second[(N2*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]);       ce2 = cf_e_second[(N4*(g_max2-g_min2)+x2-g_min2)*2*o_max+o]; 
                      ce3 = conj(cf_e_second[((N2*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]);   ce4 = cf_e_second[((N4*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o]; 
                      ce5 = cf_e_second[(N3*(g_max2-g_min2)+x2-g_min2)*2*o_max+o];             ce6 = cf_e_second[((N3*(g_max2-g_min2)+x2-g_min2)*2+1)*o_max+o];

                      temp_t1_2 = ce1*ce2 + ce3*ce4;
                      temp_t4_2 = ce1*ce5 + ce3*ce6;

                      t1 += temp_t1_1*temp_t1_2*coul_tmp;
                      t4 += temp_t4_1*temp_t4_2*coul_tmp;

	          }

                  coul_tmp = coul_pphetero;

               }

               t1+= t1_1[count1]*t1_2[count2]*coul_tmp;
               t4+= t4_1[count1]*t4_2[count2]*coul_tmp;

            } // o1

         } 

         else { // if (dd<dd_cutoff_sq)

              coul_tmp = 1.439363044/sqrt(dd);
              Complex_t t1_1_Cutoff=0.0, t1_2_Cutoff=0.0, t4_1_Cutoff=0.0, t4_2_Cutoff=0.0;

              for (int o1=0; o1<o_max_eff; o1++) {

                   int count1 = 4*(x1-g_min1) + o1; int count2 = 4*(x2-g_min2) + o1;

                   t1_1_Cutoff += t1_1[count1]; t1_2_Cutoff += t1_2[count2];
                   t4_1_Cutoff += t4_1[count1]; t4_2_Cutoff += t4_2[count2];

              }

              t1 += t1_1_Cutoff*t1_2_Cutoff*coul_tmp;
              t4 += t4_1_Cutoff*t4_2_Cutoff*coul_tmp;

         } // else to if (dd<dd_cutoff_sq)


     } // x2

     if (g_min1 == g_min2) {
        for(int o1=0;o1<o_max;o1++){
            for(int o2=0;o2<o_max;o2++){
                       
               double exch_tmp=w_exch(x1,x1,o1,o2);
                       
               if (exch_tmp!=0.0) {
                          
                   Complex_t t23_1 = conj(cf_e(N1,x1-g_min1,0,o1))*cf_e(N3,x1-g_min1,0,o2) + conj(cf_e(N1,x1-g_min1,1,o1))*cf_e(N3,x1-g_min1,1,o2);
                   Complex_t t2_2 = conj(cf_e(N2,x1-g_min1,0,o1))*cf_e(N4,x1-g_min1,0,o2) + conj(cf_e(N2,x1-g_min1,1,o1))*cf_e(N4,x1-g_min1,1,o2);
                   Complex_t t3_2 = conj(cf_e(N2,x1-g_min1,0,o2))*cf_e(N4,x1-g_min1,0,o1) + conj(cf_e(N2,x1-g_min1,1,o2))*cf_e(N4,x1-g_min1,1,o1);
                   Complex_t t56_1 = conj(cf_e(N1,x1-g_min1,0,o1))*cf_e(N4,x1-g_min1,0,o2) + conj(cf_e(N1,x1-g_min1,1,o1))*cf_e(N4,x1-g_min1,1,o2);
                   Complex_t t5_2 = conj(cf_e(N2,x1-g_min1,0,o2))*cf_e(N3,x1-g_min1,0,o1) + conj(cf_e(N2,x1-g_min1,1,o2))*cf_e(N3,x1-g_min1,1,o1);
                   Complex_t t6_2 = conj(cf_e(N2,x1-g_min1,0,o1))*cf_e(N3,x1-g_min1,0,o2) + conj(cf_e(N2,x1-g_min1,1,o1))*cf_e(N3,x1-g_min1,1,o2);
             
                   t2+=t23_1*t2_2*exch_tmp;
                   t3+=t23_1*t3_2*exch_tmp;
                   t5+=t56_1*t5_2*exch_tmp;
                   t6+=t56_1*t6_2*exch_tmp;

               } // if(exch_tmp!=0.0)
            } // o2 
         } // o1
      } // g_min1 == g_min2
   } // x1

   //cout<<"T1: "<<t1<<"   T2: "<<t2<<"   T3: "<<t3<<endl;
   //cout<<"T4: "<<t4<<"   T5: "<<t5<<"   T6: "<<t6<<endl;

   coul_Matrix[index] += (t1+t2+t3)*inv_dielect;
   exch_Matrix[index] -= (t4+t5+t6)*inv_dielect;

}

void Coulomb::compute_exchange(int m, int n, int p, int q, int id)
{
/*
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
	 // }else{
	//	tmp=1.439363044/sqrt(dd);
	//	t1_temp_s2_0+=dot3_s2_0[x2]*tmp;
	//	t4_temp_s2_0+=dot4_s2_0[x2]*tmp;
	//	t1_temp_s2_1+=dot3_s2_1[x2]*tmp;
         //       t4_temp_s2_1+=dot4_s2_1[x2]*tmp;
	//  }
	
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
*/
}


void Coulomb::compute_exchange_cutoff(int m, int n, int p, int q, int id, double cutoff)
{
/*
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
*/
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


