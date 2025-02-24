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

  for(int m=0;m<o_total;m++) {
    for(int iatom=0; iatom<num_atom_type; iatom++) 
      _coul_on_avg[m]+=_coul_on[m*num_atom_type+iatom];
    _coul_on_avg[m]/=num_atom_type;
  }

}

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


void Coulomb::get_data(Complex_t* coef_e, Complex_t* coef_h, 
		       int* atomid, int* nn, double* lat) 
{
   _cf_e = coef_e;
   _cf_h = coef_h;
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
   cout<<endl<<"hole wavefunction"<<endl;
   for(int g=g_min; g<g_max; g++) 
     for(int s=0; s<2; s++) 
         for(int o=0; o<o_max; o++)
            cout<<cf_h(0,g,s,o)<<endl; 
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


double Coulomb::w_exch(int g1, int g2, int o1, int o2)
{
   double final;
   if(g1==g2) final=exch_on(o1,o2,_atom[g1]) ; 
   else       final=exch_off(o1,o2,pair(g1,g2)); 
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

void Coulomb::compute_coulomb_exchange(void)
{
 for(int ne=0;ne<e_max;ne++){
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
           for(int adj=0;adj<5;adj++) {
              int x2;		    
              if(adj==4) x2=x1; else x2=nn(x1,adj);
              if(x2!=-1) { 
                for(int o2=0;o2<o_max;o2++){
           	  if((x1!=x2)||(o1!=o2)) {
           	    double exch_tmp=w_exch(x1,x2,o1,o2);
                    if(exch_tmp!=0.0) {
                      Complex_t ce2_s1=cf_e(ne,x2,s1,o2);
                      Complex_t ce2_s2f=cf_e(ne,x2,s2f,o2);
           	      Complex_t ch2_s2=cf_h(nh,x2,s2,o2);
		      Complex_t ch2_s1f=cf_h(nh,x2,s1f,o2);
                      t2+=conj(ce1_s1)*ce2_s1*conj(ch1_s2)*ch2_s2*exch_tmp;
             	      t3+=conj(ce1_s1)*ce2_s1*conj(ch2_s2)*ch1_s2*exch_tmp;
                      t5+=conj(ce1_s1)*ch2_s1f*conj(ch2_s2)*ce1_s2f*exch_tmp;
             	      t6+=conj(ce1_s1)*ch2_s1f*conj(ch1_s2)*ce2_s2f*exch_tmp;
                    } // if(exch_tmp!=0.0)
                  } // if((x1!=x2)||(o1!=o2))
		} // o2 
              } // if(x2!=-1)
           } // adj
       }} // x1 and o1
    }} // s1 and s2   

   H_coul(ne,nh)-=(t1+t2+t3)*inv_dielect;
   H_exch(ne,nh)+=(t4+t5+t6)*inv_dielect;
//   cout<<"ne, nh, Coulomb, Exchange: "<<ne<<" "<<nh<<" "<<H_coul(ne,nh)<<" "<<H_exch(ne,nh)<<endl;
 }} // ne and nh
}


