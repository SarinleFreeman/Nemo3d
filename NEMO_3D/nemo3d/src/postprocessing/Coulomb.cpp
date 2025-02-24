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
$Header: /repo/nemo3d/src/postprocessing/Coulomb.cpp,v 1.1 2004/08/26 21:28:07 swlee Exp $
*****************************************************************************/

#include "Coulomb.h"
#include "SimpleParser.h"

const double
  Coulomb::CONSTANT1 = 0.482085578, // the factor to convert nm^2 to 1/eV^2.
  Coulomb::CONSTANT2 = 1.440250154; // the factor to convert 1/nm to 1/eV.

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


void Coulomb::get_data(complex* coef_e, complex* coef_h, 
		       int* atomid, int* nn, double* lat) 
{
   _cf_e = coef_e;
   _cf_h = coef_h;
   _atom = atomid;
   _nn   = nn;  
   _lattice = lat; 
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
    final = 1.0/sqrt(1.0/(tmp*tmp)+dd*CONSTANT1); 
  }   

  return final;
}

double Coulomb::distance_sq(int g1, int g2)
{
  double dx = lattice(g1,0)-lattice(g2,0); 
  double dy = lattice(g1,1)-lattice(g2,1); 
  double dz = lattice(g1,2)-lattice(g2,2); 
  double dd = dx*dx + dy*dy + dz*dz; 
  return dd;
}  


double Coulomb::coulomb_cell(int cell1, int cell2)
{
  int g1=cell1*d->geo.AtomsPerCellMax();
  int g2=cell2*d->geo.AtomsPerCellMax();
  double dx = lattice(g1,0)-lattice(g2,0); 
  double dy = lattice(g1,1)-lattice(g2,1); 
  double dz = lattice(g1,2)-lattice(g2,2); 
  double dd = dx*dx + dy*dy + dz*dz; 
  double answer = CONSTANT2/sqrt(dd); 
  return answer;
}


void Coulomb::compute_coulomb_exchange(void)
{
/* not yet implemented in this embedded version.
   please use the self-standing version for this function */
}

void Coulomb::print_phase()
{
   ofstream fout("phase.dat");
   int cell_start = 0;
   int cell_end = d->geo.N_Cell;
   for(int cell=cell_start; cell<cell_end; cell++){
/*
     int cell_x = d->geo.cell__ijk[cell][0]; 
     int cell_y = d->geo.cell__ijk[cell][1]; 
     int cell_z = d->geo.cell__ijk[cell][2];
*/
       for(int atom=0; atom<d->geo.AtomsPerCellMax(); atom++){
	 int x=cell*d->geo.AtomsPerCellMax()+atom;
         complex wf1=cf_e(0,x,1,2); //spin down and p_x orbital
         complex wf2=cf_e(0,x,1,3); //spin down and p_y orbital
         complex wf3=cf_e(0,x,1,4); //spin down and p_z orbital
/*
         complex wf4=cf_e(0,x,0,2); //spin up   and p_x orbital
         complex wf5=cf_e(0,x,0,3); //spin up   and p_y orbital
         complex wf6=cf_e(0,x,0,4); //spin up   and p_z orbital
         double ratio1 = (wf1.r*wf1.r+wf1.i*wf1.i)/(wf4.r*wf4.r+wf4.i*wf4.i);
         double ratio2 = (wf2.r*wf2.r+wf2.i*wf2.i)/(wf5.r*wf5.r+wf5.i*wf5.i);
         double ratio3 = (wf3.r*wf3.r+wf3.i*wf3.i)/(wf6.r*wf6.r+wf6.i*wf6.i);
*/
         complex wf4=cf_e(0,x,1,5); //spin down and d_xy orbital
         complex wf5=cf_e(0,x,1,6); //spin down and d_yz orbital
         complex wf6=cf_e(0,x,1,7); //spin down and d_zx orbital
         complex wf7=cf_e(0,x,1,8); //spin down and d_xx_yy orbital
         complex wf8=cf_e(0,x,1,9); //spin down and d_zz-3rr orbital

         double cos1=wf1.r/sqrt((wf1.r*wf1.r)+(wf1.i*wf1.i));
         double cos2=wf2.r/sqrt((wf2.r*wf2.r)+(wf2.i*wf2.i));
         double cos3=wf3.r/sqrt((wf3.r*wf3.r)+(wf3.i*wf3.i));
         double cos4=wf4.r/sqrt((wf4.r*wf4.r)+(wf4.i*wf4.i));
         double cos5=wf5.r/sqrt((wf5.r*wf5.r)+(wf5.i*wf5.i));
         double cos6=wf6.r/sqrt((wf6.r*wf6.r)+(wf6.i*wf6.i));
         double cos7=wf7.r/sqrt((wf7.r*wf7.r)+(wf7.i*wf7.i));
         double cos8=wf8.r/sqrt((wf8.r*wf8.r)+(wf8.i*wf7.i));

         double theta1 = acos(cos1); if(wf1.i<0.0) theta1=-theta1;
         double theta2 = acos(cos2); if(wf2.i<0.0) theta2=-theta2;
         double theta3 = acos(cos3); if(wf3.i<0.0) theta3=-theta3;
         double theta4 = acos(cos4); if(wf4.i<0.0) theta4=-theta4;
         double theta5 = acos(cos5); if(wf5.i<0.0) theta5=-theta5;
         double theta6 = acos(cos6); if(wf6.i<0.0) theta6=-theta6;
         double theta7 = acos(cos7); if(wf7.i<0.0) theta7=-theta7;
         double theta8 = acos(cos8); if(wf8.i<0.0) theta8=-theta8;

         fout<<lattice(x,0)<<" "<<lattice(x,1)<<" "<<lattice(x,2)
             <<" "<<theta1<<" "<<theta2<<" "<<theta3
             <<" "<<theta4<<" "<<theta5<<" "<<theta6
             <<" "<<theta7<<" "<<theta8<<endl;
       } 
   }
}

void Coulomb::compute_overlap_vs_dot_separation(int separation, int direction)
{
  complex Overlap;
  Overlap.r=0.0e0; Overlap.i=0.0e0;
  int cell_start = d->cell_s[my_id];
  int cell_end   = cell_start+d->cell_ln[my_id];
  int count_cell=0;
  for(int cell1=cell_start; cell1<cell_end; cell1++)
  {
    int cell1_overlap = find_cell_overlap(cell1,separation, direction);
    if(cell1_overlap>=0)
    {
      for(int atom1=0; atom1<d->geo.AtomsPerCellMax(); atom1++)
      {
        int xL=cell1*d->geo.AtomsPerCellMax()+atom1;
        int xR=cell1_overlap*d->geo.AtomsPerCellMax()+atom1;
        for(int orbital1=0;orbital1<o_max;orbital1++)
        {
          complex wf_L1=cf_e(0,xL,1,orbital1); // spin down
          complex wf_R1=cf_e(1,xR,1,orbital1); // spin down
//          complex wf_L1=cf_e(0,xL,0,orbital1);  // spin up
//          complex wf_R1=cf_e(0,xR,0,orbital1);  // spin up
          complex wf_tmp1=multiply(wf_L1,wf_R1);
          Overlap.r+=wf_tmp1.r;
          Overlap.i+=wf_tmp1.i;
        }
      }
      count_cell++;
    }
  }

  complex T_Overlap;
  MPI_Allreduce(&Overlap.r, &T_Overlap.r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Overlap.i, &T_Overlap.i, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double overlap_sq = T_Overlap.r*T_Overlap.r + T_Overlap.i*T_Overlap.i; 
  int T_count;
  MPI_Allreduce(&count_cell, &T_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  cout_master<<separation<<" "<<overlap_sq<<" "<<T_count<<endl;

}


void Coulomb::compute_exchange_vs_dot_separation(int separation, int direction)
//For now, the separation is given by the number of unit cells 
//over which the two wave functions overlap. 
{
   complex H_exch,t1,t2,t3;
   H_exch.r=0.0e0; H_exch.i=0.0e0;
   t1.r=0.0e0, t1.i=0.0e0, t2.r=0.0e0, t2.i=0.0e0, t3.r=0.0e0, t3.i=0.0e0;

   int cell_start = d->cell_s[my_id];
   int cell_end   = cell_start+d->cell_ln[my_id];

   for(int cell1=cell_start; cell1<cell_end; cell1++){
     int cell1_overlap = find_cell_overlap(cell1,separation, direction);
     if(cell1_overlap==-1) continue;
     for(int atom1=0; atom1<d->geo.AtomsPerCellMax(); atom1++){
          int xL=cell1*d->geo.AtomsPerCellMax()+atom1;
          int xR=cell1_overlap*d->geo.AtomsPerCellMax()+atom1; 
          for(int orbital1=0;orbital1<o_max;orbital1++){
            complex wf_L1=cf_e(0,xL,1,orbital1);
            complex wf_R1=cf_e(1,xR,1,orbital1);
            complex wf_tmp1=multiply(wf_R1,wf_L1);
           // conjugate of this will be mutiplied to another wf_tmp;

            for(int orbital2=0;orbital2<o_max;orbital2++){
              complex wf_L2=cf_e(0,xL,1,orbital2);
              complex wf_R2=cf_e(1,xR,1,orbital2);
              double  coul_tmp=w_coul(0.0,xL,orbital1, orbital2);
              double  exch_tmp=w_exch(xL,xL,orbital1, orbital2);
              complex wf_tmp2=multiply(wf_R2,wf_L2);
              complex wf_tmp3=multiply(wf_R2,wf_L1);
              complex wf_tmp4=multiply(wf_R1,wf_L2);
              complex wf_tmp5=multiply(wf_R2,wf_L1);
              complex wf_tmp6=multiply(wf_tmp1, wf_tmp2);
              complex wf_tmp7=multiply(wf_tmp3, wf_tmp4);
              complex wf_tmp8=multiply(wf_tmp3, wf_tmp5);
              
              t1.r+=wf_tmp6.r*coul_tmp;
//              t1.i+=wf_tmp6.i*coul_tmp;
              t2.r+=wf_tmp7.r*exch_tmp;
//              t2.i+=wf_tmp7.i*exch_tmp;
              t3.r+=wf_tmp8.r*exch_tmp;
//              t3.i+=wf_tmp8.i*exch_tmp;
           } // orbital2 for loop exchange interaction between orbitals from the same atom 

           for(int atom2=0; atom2<d->geo.AtomsPerCellMax(); atom2++){
             if(atom2==atom1) continue; 
             int x2L=cell1*d->geo.AtomsPerCellMax()+atom2;
             int x2R=cell1_overlap*d->geo.AtomsPerCellMax()+atom2; 
             double dd = distance_sq(xL,x2L);  
             for(int orbital2=0;orbital2<o_max;orbital2++){
                 complex wf_L2 = cf_e(0,x2L,1,orbital2);
                 complex wf_R2 = cf_e(1,x2R,1,orbital2);
                 double coul_tmp=w_coul(dd,x2L,orbital1, orbital2);
                 complex wf_tmp2=multiply(wf_R2, wf_L2);
                 complex wf_tmp6=multiply(wf_tmp1, wf_tmp2);
                 t1.r+=wf_tmp6.r*coul_tmp;
//               t1.i+=wf_tmp6.i*coul_tmp;
             } //orbital2 loop
           } // atom2 for loop exchange interaction between atoms from the same cell

           for(int cell2=0; cell2<d->geo.N_Cell; cell2++){
             bool close_cell = AreCellsClose(cell1, cell2);           
             if(!close_cell) continue; 
             int cell2_overlap = find_cell_overlap(cell2, separation, direction);
             if(cell2_overlap==-1) continue;
             for(int atom2=0; atom2<d->geo.AtomsPerCellMax(); atom2++){
                   int x2L=cell2*d->geo.AtomsPerCellMax()+atom2;
                   int x2R=cell2_overlap*d->geo.AtomsPerCellMax()+atom2;
                   double dd = distance_sq(xL,x2L);
                   for(int orbital2=0;orbital2<o_max;orbital2++){
                     complex wf_L2 = cf_e(0,x2L,1,orbital2);
                     complex wf_R2 = cf_e(1,x2R,1,orbital2);
                     double coul_tmp=w_coul(dd,x2L,orbital1, orbital2);
                     complex wf_tmp2=multiply(wf_R2, wf_L2);
                     complex wf_tmp6=multiply(wf_tmp1, wf_tmp2);
                     t1.r+=wf_tmp6.r*coul_tmp;
//                     t1.i+=wf_tmp6.i*coul_tmp;
                   } //orbital2 loop
             } // atom2 for loop
           } //cell2 loop exchange interaction between atoms from different cells.


         } // orbital1 for loop
     } // atom1 for loop
   } // cell for loop

   complex t4 = compute_distant_exchange_interaction(separation,direction);
   
   H_exch.r+=(t1.r+t2.r+t3.r+t4.r)*inv_dielect;
//   H_exch.i+=(t1.i+t2.i+t3.i+t4.i)*inv_dielect;
   
   double Total_exch_real;
   double Total_exch_imag;
    
   double T1[2], T2[2], T3[2], T4[2];
   MPI_Allreduce(&H_exch.r, &Total_exch_real, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&H_exch.i, &Total_exch_imag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&t1.r, &T1, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&t2.r, &T2, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&t3.r, &T3, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&t4.r, &T4, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
   cout_master<<"overlap, exch : "<<separation<<" "<<Total_exch_real <<" "<<Total_exch_imag<<endl;
   cout_master<<"overlap, t1   : "<<separation<<" "<<T1[0]*inv_dielect<<" "<<T1[1]*inv_dielect<<endl;
   cout_master<<"overlap, t2   : "<<separation<<" "<<T2[0]*inv_dielect<<" "<<T2[1]*inv_dielect<<endl;
   cout_master<<"overlap, t3   : "<<separation<<" "<<T3[0]*inv_dielect<<" "<<T3[1]*inv_dielect<<endl;
   cout_master<<"overlap, t4   : "<<separation<<" "<<T4[0]*inv_dielect<<" "<<T4[1]*inv_dielect<<endl<<endl;
}

complex Coulomb::compute_distant_exchange_interaction(int separation, int direction)
{
  complex t4; t4.r=0.0e0; t4.i=0.0e0;

  int cell_start = d->cell_s[my_id];
  int cell_end   = cell_start+d->cell_ln[my_id];

  for(int cell1=cell_start; cell1<cell_end; ++cell1){
     int cell1_overlap = find_cell_overlap(cell1,separation, direction);
     if(cell1_overlap==-1) continue;

     complex wf1; wf1.r=0.0e0; wf1.i=0.0e0;
     for(int atom1=0; atom1<d->geo.AtomsPerCellMax(); ++atom1){
         int x1L=cell1*d->geo.AtomsPerCellMax()+atom1;
         int x1R=cell1_overlap*d->geo.AtomsPerCellMax()+atom1;
         for(int orbital1=0;orbital1<o_max;++orbital1){
           complex wf_L1=cf_e(0,x1L,1,orbital1);
           complex wf_R1=cf_e(1,x1R,1,orbital1);
           complex wf1_tmp = multiply(wf_R1, wf_L1);
           wf1.r+=wf1_tmp.r; wf1.i+=wf1_tmp.i;
         }
     }// sum up all the atom and orbital to wave function for a given cell1.

     for(int cell2=0; cell2<d->geo.N_Cell; ++cell2){
       bool close_cell = AreCellsClose(cell1, cell2);
       if(close_cell) continue; 
       int cell2_overlap = find_cell_overlap(cell2,separation, direction);
       if( cell2_overlap==-1) continue;

       complex wf2; wf2.r=0.0e0; wf2.i=0.0e0;
       for(int atom2=0; atom2<d->geo.AtomsPerCellMax(); atom2++){
           int x2L=cell1*d->geo.AtomsPerCellMax()+atom2;
           int x2R=cell1_overlap*d->geo.AtomsPerCellMax()+atom2;
           for(int orbital2=0;orbital2<o_max;orbital2++){
             complex wf_L2=cf_e(0,x2L,1,orbital2);
             complex wf_R2=cf_e(1,x2R,1,orbital2);
             complex wf2_tmp = multiply(wf_R2, wf_L2);
             wf2.r+=wf2_tmp.r; wf2.i+=wf2_tmp.i;
           }
       }
//       complex wf = multiply(wf1, wf2); 
       double real_wf = wf1.r*wf2.r+wf1.i*wf2.i; 
       double coul_tmp = coulomb_cell(cell1, cell2);
       t4.r+=real_wf*coul_tmp;
//       t4.r+=wf.r*coul_tmp; 
//       t4.i+=wf.i*coul_tmp;
    }
  }
  return t4;
}

bool Coulomb::AreCellsClose(int cell1, int cell2)
{
  int cell1_x = d->geo.cell__ijk[cell1][0];
  int cell1_y = d->geo.cell__ijk[cell1][1];
  int cell1_z = d->geo.cell__ijk[cell1][2];
  int cell2_x = d->geo.cell__ijk[cell2][0];
  int cell2_y = d->geo.cell__ijk[cell2][1];
  int cell2_z = d->geo.cell__ijk[cell2][2];
  bool answer;
  if(abs(cell1_x-cell2_x)<=1 && abs(cell1_y-cell2_y)<=1 && abs(cell1_z-cell2_z)<=1) answer=true;
  else answer=false;
  return answer;
}

int Coulomb::find_cell_overlap(int cell, int separation, int direction)
{
/*
     cout_master<<"cell xmin, xmax, ymin, ymax, zmin, zmax "
         <<d->geo.cell_xmin<<" "<<d->geo.cell_xmax<<" "
         <<d->geo.cell_ymin<<" "<<d->geo.cell_ymax<<" "
         <<d->geo.cell_zmin<<" "<<d->geo.cell_zmax<<endl;
*/
     int cell_x = d->geo.cell__ijk[cell][0];
     int cell_y = d->geo.cell__ijk[cell][1];
     int cell_z = d->geo.cell__ijk[cell][2];
     int cell_overlap_x, cell_overlap_y, cell_overlap_z, cell_overlap_index;
     int cell_overlap_min, cell_overlap_max;
     if(direction==1) {
       cell_overlap_x = cell_x+separation-d->geo.cell_xmax+(d->geo.cell_xmin-1);
       cell_overlap_index = cell_overlap_x;
       cell_overlap_min = d->geo.cell_xmin;
       cell_overlap_max = d->geo.cell_xmax;
     }
     else if(direction==2) {
       cell_overlap_y = cell_y+separation-d->geo.cell_ymax+(d->geo.cell_ymin-1);
       cell_overlap_index = cell_overlap_y;
       cell_overlap_min = d->geo.cell_ymin;
       cell_overlap_max = d->geo.cell_ymax;
     }
     else if(direction==3) {
       cell_overlap_z = cell_z+separation-d->geo.cell_zmax+(d->geo.cell_zmin-1);
       cell_overlap_index = cell_overlap_z;
       cell_overlap_min = d->geo.cell_zmin;
       cell_overlap_max = d->geo.cell_zmax;
     }

     int cell_overlap;
     if(cell_overlap_index>=cell_overlap_min && cell_overlap_index<=cell_overlap_max) {
       //there is an overlap between the two wave functions
       if(direction==1)
         cell_overlap = d->geo.ijk__cell[cell_overlap_x][cell_y][cell_z];
       else if(direction==2)
         cell_overlap = d->geo.ijk__cell[cell_x][cell_overlap_y][cell_z];
       else if(direction==3)
         cell_overlap = d->geo.ijk__cell[cell_x][cell_y][cell_overlap_z];
     }
     else cell_overlap = -1;
/* 
     if(cell_overlap!=-1) cout_master<<"cell_y, cell_overlap_y, num_overlap_cell: "
                                     <<cell_y<<" "<<cell_overlap_y<<" "<<separation<<endl;
*/
     return cell_overlap;
}
