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
$Header: /repo/nemo3d/src/postprocessing/Decoherence.cpp,v 1.3 2007/02/23 19:00:36 hoonryu Exp $
*****************************************************************************/

#include "Decoherence.h"
#define magnetic_moment_As 1.86 
#define magnetic_moment_Ga 2.89 
#define magnetic_moment_In 6.13 
#define nuclear_spin_As 1.5 
#define nuclear_spin_Ga 1.5 
#define nuclear_spin_In 4.5 
#define nuclear_spin_square_As 3.75
#define nuclear_spin_square_Ga 3.75
#define nuclear_spin_square_In 24.75 
#define electron_density_As 9.8 // in unit of 10e25 /cm3
#define electron_density_Ga 5.8 // in unit of 10e25 /cm3
#define electron_density_In 9.4 // in unit of 10e25 /cm3
#define hbar 1.055e-27 //in unit of erg*sec
#define hf_coupling_coefficient 7.849e-18 //in unit of erg with electron density 10e25 /cm3
                                          // 16*Pi/3*bohr_magneton*nuclear_magneton*10e25/cm3 
#define hbar_over_hf_coupling_coefficient 1.4095e-10 //in unit of second 
#define convertor_erg_to_eV 6.241450383e11 //  1 erg = 6.241450383e11 eV
#define electron_density_As_s 18.4 // in unit of 10e25 /cm3
#define electron_density_Ga_s 5.2  // in unit of 10e25 /cm3
#define electron_density_In_s 7.9  // in unit of 10e25 /cm3
#define ratio_sstar_to_s_orbital_As 0.3 
#define ratio_sstar_to_s_orbital_Ga 0.44
#define ratio_sstar_to_s_orbital_In 0.53
#define bohr_magneton 5.788381749e-5 // in unit of eV/T
#define bohr_magneton_erg_over_gauss 9.2741e-21 // in unit of erg/gauss
#define nuclear_magneton 5.051e-24 // in unit of erg/gauss     
#define hbar_in_eV_sec 6.582183e-16 // in unit of eV sec  

void Decoherence::compute_decoherence_time (void)
{
  double A[3], I[3], II[3], Rs[3];
  A[0] = magnetic_moment_As/nuclear_spin_As*electron_density_As_s;  
  A[1] = magnetic_moment_Ga/nuclear_spin_Ga*electron_density_Ga_s;  
  A[2] = magnetic_moment_In/nuclear_spin_In*electron_density_In_s; 
  I[0] = nuclear_spin_As;
  I[1] = nuclear_spin_Ga;
  I[2] = nuclear_spin_In;
  II[0] = nuclear_spin_square_As; 
  II[1] = nuclear_spin_square_Ga; 
  II[2] = nuclear_spin_square_In;
  Rs[0] = ratio_sstar_to_s_orbital_As; 
  Rs[1] = ratio_sstar_to_s_orbital_Ga; 
  Rs[2] = ratio_sstar_to_s_orbital_In; 
  double AI_interface_Variance = 0.5*0.5*(I[2]-I[1])*(I[2]-I[1]); 
  double AI_alloy_Variance = Composition*(1.0-Composition)*(A[2]*I[2]-A[1]*I[1])*(A[2]*I[2]-A[1]*I[1]); 
 
  double* HF_coupling = new double[g_max-g_min]; 

  char *filename = NULL;
  char label[100];
  sprintf(label, "hf_coupling_%d", my_id);
  filename = n3d_strdup_n(d->inputfile);
  n3d_FileTypeSet(&filename, label, true);
  ofstream fout(filename);

  str_free(filename);

  int n_max = e_max+h_max;
  double A_sum[n_max], AI_sum[n_max], AA_sum[n_max], AAII_sum[n_max];
  double AAII_alloy_sum[n_max], AAII_interface_sum[n_max];
  int count_interface_atom[n_max];
  int count_dot_atom[n_max];
  for(int n=0;n<n_max;n++){
    A_sum[n]= 0.0; AI_sum[n]=0.0; AA_sum[n]= 0.0; AAII_sum[n] = 0.0; 
    AAII_alloy_sum[n]=0.0; AAII_interface_sum[n]=0.0;
    count_interface_atom[n]=0;
    count_dot_atom[n]=0;
    double A_max=0.0;
    for(int g=g_min;g<g_max;g++){
      double amplitude=0.0;
      if(n<e_max) {
        int m=n;
        complex coef_up, coef_down;
        coef_up.r = cf_e(m,g,0,0).r*Rs[_atom[g]]+ cf_e(m,g,0,1).r; 
        coef_up.i = cf_e(m,g,0,0).i*Rs[_atom[g]]+ cf_e(m,g,0,1).i; 
        coef_down.r = cf_e(m,g,1,0).r*Rs[_atom[g]]+ cf_e(m,g,1,1).r; 
        coef_down.i = cf_e(m,g,1,0).i*Rs[_atom[g]]+ cf_e(m,g,1,1).i; 
        amplitude += coef_up.r*coef_up.r + coef_up.i*coef_up.i; 
        amplitude += coef_down.r*coef_down.r + coef_down.i*coef_down.i;
      }
      else {
        int m=n-e_max;
        complex coef_up, coef_down;
        coef_up.r = cf_h(m,g,0,0).r*Rs[_atom[g]]+ cf_h(m,g,0,1).r; 
        coef_up.i = cf_h(m,g,0,0).i*Rs[_atom[g]]+ cf_h(m,g,0,1).i; 
        coef_down.r = cf_h(m,g,1,0).r*Rs[_atom[g]]+ cf_h(m,g,1,1).r; 
        coef_down.i = cf_h(m,g,1,0).i*Rs[_atom[g]]+ cf_h(m,g,1,1).i; 
        amplitude += coef_up.r*coef_up.r + coef_up.i*coef_up.i; 
        amplitude += coef_down.r*coef_down.r + coef_down.i*coef_down.i;
      }
      HF_coupling[g-g_min] = amplitude*A[_atom[g]];
      if(HF_coupling[g-g_min]> A_max) A_max=HF_coupling[g-g_min];
      double hf_hf = HF_coupling[g-g_min]*HF_coupling[g-g_min];
      A_sum[n]  += HF_coupling[g-g_min];
      AI_sum[n] += HF_coupling[g-g_min]*I[_atom[g]];
      AA_sum[n] += hf_hf; 
      AAII_sum[n] += hf_hf*II[_atom[g]];

      if(CationNearInterface(g)) {
        AAII_interface_sum[n] += amplitude*amplitude*AI_interface_Variance; 
        count_interface_atom[n]++;
      }

      if(CationInsideDot(g)) { 
        AAII_alloy_sum[n] += amplitude*amplitude*AI_alloy_Variance; 
        count_dot_atom[n]++;
      } 
    } 

    double tmp_convertor = hf_coupling_coefficient*convertor_erg_to_eV;
    double Total_A_max;
    MPI_Allreduce(&A_max, &Total_A_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if(my_id==0) cout<<"A_max is "<<Total_A_max*tmp_convertor<<" eV"<<endl;
    double A_min=0.1*Total_A_max;

    int count_nuclei=0;
    int total_nuclei=0;
    for(int g=g_min; g<g_max; g++)
      if(HF_coupling[g-g_min]> A_min) { 
          count_nuclei++;
          fout<<I[_atom[g]]<<"\t"<<HF_coupling[g-g_min]*tmp_convertor<<endl;
      }    
    MPI_Allreduce(&count_nuclei, &total_nuclei, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(my_id==0) cout<<"Number of Nuclei interacting with electron "<<total_nuclei<<endl;

/*
    char *filename1 = NULL;
    char label[100];
    sprintf(label, "A_along_z_%d", my_id);
    filename1 = n3d_strdup_n(d->inputfile);
    n3d_FileTypeSet(&filename1, label, true);
    ofstream HFzout(filename1);
    HFzout.setf(ios::scientific);
    for(int g=g_min;g<g_max;g++){
      if ( lattice(g,0) > 10.2 && lattice(g,0)<10.8 && lattice(g,1)>10.2 && lattice(g,1)<10.8 ) 
        HFzout<< lattice(g,2) << "\t" << HF_coupling[g-g_min]*tmp_convertor<<"\t"<<_atom[g]<< "\t"<< HF_coupling[g-g_min]/A[_atom[g]]<<endl;
    }

    char *filename2 = NULL;
    sprintf(label, "A_along_x_%d", my_id);
    filename2 = n3d_strdup_n(d->inputfile);
    n3d_FileTypeSet(&filename2, label, true);
    ofstream HFxout(filename2);
    HFxout.setf(ios::scientific);
    for(int g=g_min;g<g_max;g++){
      if ( lattice(g,2) > 5 && lattice(g,2)<5.4 && lattice(g,1)>10.3 && lattice(g,1)<10.7 ) 
        HFxout<< lattice(g,0) << "\t" << HF_coupling[g-g_min]*tmp_convertor<<"\t"<<_atom[g]<<"\t"<< HF_coupling[g-g_min]/A[_atom[g]]<<endl;
    }
*/

    
  } 

  double total_A_sum[n_max];
  double total_AI_sum[n_max];
  double total_AA_sum[n_max];
  double total_AAII_sum[n_max];
  double total_AAII_alloy_sum[n_max];
  double total_AAII_interface_sum[n_max];
  int total_interface_atom[n_max];
  int total_dot_atom[n_max];
  MPI_Allreduce(&A_sum, &total_A_sum, n_max, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&AI_sum, &total_AI_sum, n_max, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&AA_sum, &total_AA_sum, n_max, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&AAII_sum, &total_AAII_sum, n_max, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&AAII_alloy_sum, &total_AAII_alloy_sum, n_max, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&AAII_interface_sum, &total_AAII_interface_sum, n_max, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_interface_atom, &total_interface_atom, n_max, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_dot_atom, &total_dot_atom, n_max, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(my_id==0) {
    cout.setf(ios::scientific);
    cout<<setprecision(20);
    double tmp_convertor = hf_coupling_coefficient*convertor_erg_to_eV;
    cout<<"A: hyperfine coupling without electron density "<<endl;
    cout<<"As nucleus: "<<A[0]*tmp_convertor<<" eV"<<endl;
    cout<<"Ga nucleus: "<<A[1]*tmp_convertor<<" eV"<<endl;
    cout<<"In nucleus: "<<A[2]*tmp_convertor<<" eV"<<endl;
    for(int n=0; n<n_max; n++) {
      dephasing_time[n] = 1.0/sqrt(total_AAII_sum[n]); 
      dephasing_time[n] *= hbar_over_hf_coupling_coefficient; 
      alloy_dephasing_time[n] = 1.0/sqrt(total_AAII_alloy_sum[n]);
      alloy_dephasing_time[n] *= hbar_over_hf_coupling_coefficient; 
      interface_dephasing_time[n] = 1.0/sqrt(total_AAII_interface_sum[n]);
      interface_dephasing_time[n] *= hbar_over_hf_coupling_coefficient; 

      double denominator = (total_AA_sum[n]/numAtoms) - (total_A_sum[n]*total_A_sum[n])
                           /numAtoms/numAtoms;
      decoherence_time[n] = 2.0* hbar_over_hf_coupling_coefficient/sqrt(denominator); 
      if(n<e_max)      cout<<"Electron Level: "<<n<<endl;
      else             cout<<"Hole Level: "<<n-e_max<<endl;
      cout<<"Number of Interface Atoms: "<<total_interface_atom[n]<<endl;
      cout<<"Number of Dot Atoms: "<<total_dot_atom[n]<<endl;
      cout<<"A_sum (eV): "<< total_A_sum[n]*hf_coupling_coefficient*convertor_erg_to_eV<<endl;
      cout<<"AI_sum (eV) with polarized nuclear spins: "
          << total_AI_sum[n]*hf_coupling_coefficient*convertor_erg_to_eV<<endl;
      cout<<"AA_sum (eV_square): "<< 
        total_AA_sum[n]*hf_coupling_coefficient*convertor_erg_to_eV
        *hf_coupling_coefficient*convertor_erg_to_eV <<endl;
      cout<<"Dephasing_time due to nuclear spin configuration (second):  "<<dephasing_time[n]<<endl;
      cout<<"Dephasing_time due to alloy disorder (second):  "<<alloy_dephasing_time[n]<<endl;
      cout<<"Dephasing_time due to interface disorder (second):  "<<interface_dephasing_time[n]<<endl;
      cout<<"Decoherence_time (second): "<<decoherence_time[n]<<endl<<endl;
    } 
  }
}

void Decoherence::get_data(complex* coef_e, complex* coef_h, int* atom_id, int* shape_rank, 
                           int* nn, double* lat)
{
   _cf_e = coef_e;
   _cf_h = coef_h;
   _atom = atom_id;
   _shape = shape_rank;
   _nn = nn;
   _lattice = lat;
}

void Decoherence::get_geometry(double x0, double y0, double z0, 
         double xL, double yL, double zL, double composition)
{
   Xorigin = x0;
   Yorigin = y0;
   Zorigin = z0;
   Xdim = xL;
   Ydim = yL;
   Zdim = zL;
   Xradius2 = xL*xL*0.25;
   Yradius2 = yL*yL*0.25;
   Zradius2 = zL*zL;
   Xcenter = Xorigin + Xdim*0.5;
   Ycenter = Yorigin + Ydim*0.5;
   Zcenter = Zorigin;   
   Composition = composition;
}

bool Decoherence::CationInsideDot(int g)
{
  if(_atom[g]==0) return false; //anion atom
  else { //cation atom
    double x=lattice(g,0), y=lattice(g,1), z=lattice(g,2);
    if(z < Zorigin || z > Zorigin+Zdim ||
       x < Xorigin || x > Xorigin+Xdim || 
       y < Yorigin || y > Yorigin+Ydim )
       return false; // not inside a dot
    else { 
      double sum = 0.0;
      sum+=(x-Xcenter)*(x-Xcenter)/Xradius2;
      sum+=(y-Ycenter)*(y-Ycenter)/Yradius2;
      sum+=(z-Zcenter)*(z-Zcenter)/Zradius2;
      if(sum>1.0)  return false; // not inside a dot
      else return true; // inside a dot
    }
  }
}

bool Decoherence::CationNearInterface(int g)
{
  double f1 = 0.1;
  double f2 = f1+1.0;
  double f3 = (1.0-f1)*(1.0-f1);
  double f4 = (1.0+f1)*(1.0+f1);
  if(_atom[g]==0) return false; //anion atom
  else if(lattice(g,2)<(Zorigin-Zdim*f1)) return false; // cation atom at the bottom of the dot
  else if(lattice(g,2)<(Zorigin+Zdim*f1) &&
          lattice(g,0)>(Xorigin-Xdim*f1) &&
          lattice(g,0)<(Xorigin+Xdim*f2) &&
          lattice(g,1)>(Yorigin-Ydim*f1) &&
          lattice(g,1)>(Yorigin-Ydim*f2)) return true; // cation atom near the substrate 
  else {
    double x=lattice(g,0), y=lattice(g,1), z=lattice(g,2);
    double sum = 0.0;
    sum+=(x-Xcenter)*(x-Xcenter)/Xradius2;
    sum+=(y-Ycenter)*(y-Ycenter)/Yradius2;
    sum+=(z-Zcenter)*(z-Zcenter)/Zradius2;
    if(sum>f3 || sum<f4)  return false; // cation atom away from the interface 
    else return true; // cation atom near the interface 
  }
}

void Decoherence::compute_spin_correlator(const char& carrier, int n) 
{
  double A[3], I[3], II[3], Rs[3];
  A[0] = magnetic_moment_As/nuclear_spin_As*electron_density_As_s;  
  A[1] = magnetic_moment_Ga/nuclear_spin_Ga*electron_density_Ga_s;  
  A[2] = magnetic_moment_In/nuclear_spin_In*electron_density_In_s; 
  I[0] = nuclear_spin_As;
  I[1] = nuclear_spin_Ga;
  I[2] = nuclear_spin_In;
  II[0] = nuclear_spin_square_As; 
  II[1] = nuclear_spin_square_Ga; 
  II[2] = nuclear_spin_square_In;
  Rs[0] = ratio_sstar_to_s_orbital_As; 
  Rs[1] = ratio_sstar_to_s_orbital_Ga; 
  Rs[2] = ratio_sstar_to_s_orbital_In; 
  
  double VV[g_num]; // |<f|V|i>^2 where V = sum of A(SxIx+SyIy)
  double omega[g_num]; // (Ei - Ef)/hbar where Ei and Ef are the energies 
                            // (sum of ASzIz) of initial and final states
  double AI_sum=0.0;
  RandRandom myRandom(my_id, 1, true);
  for(int g=g_min;g<g_max;g++){
      double amplitude=0.0;
      if(carrier == 'e') {
        complex coef_up, coef_down;
        coef_up.r = cf_e(n,g,0,0).r*Rs[_atom[g]]+ cf_e(n,g,0,1).r; 
        coef_up.i = cf_e(n,g,0,0).i*Rs[_atom[g]]+ cf_e(n,g,0,1).i; 
        coef_down.r = cf_e(n,g,1,0).r*Rs[_atom[g]]+ cf_e(n,g,1,1).r; 
        coef_down.i = cf_e(n,g,1,0).i*Rs[_atom[g]]+ cf_e(n,g,1,1).i; 
        amplitude += coef_up.r*coef_up.r + coef_up.i*coef_up.i; 
        amplitude += coef_down.r*coef_down.r + coef_down.i*coef_down.i;
      }
      else {
        complex coef_up, coef_down;
        coef_up.r = cf_h(n,g,0,0).r*Rs[_atom[g]]+ cf_h(n,g,0,1).r; 
        coef_up.i = cf_h(n,g,0,0).i*Rs[_atom[g]]+ cf_h(n,g,0,1).i; 
        coef_down.r = cf_h(n,g,1,0).r*Rs[_atom[g]]+ cf_h(n,g,1,1).r; 
        coef_down.i = cf_h(n,g,1,0).i*Rs[_atom[g]]+ cf_h(n,g,1,1).i; 
        amplitude += coef_up.r*coef_up.r + coef_up.i*coef_up.i; 
        amplitude += coef_down.r*coef_down.r + coef_down.i*coef_down.i;
      }
      double hyperfine_coupling = amplitude*A[_atom[g]]*hf_coupling_coefficient; //in unit of erg
      omega[g-g_min] = hyperfine_coupling*0.5/hbar; //in unit of 1/second 
      double random_number = myRandom.getRandom();
      double Iz;  

      if(RandomOrientation){ //nuclear spins are randomly oriented
        if(_atom[g]==2) { // for In atom case
          if(random_number<0.1)      Iz = -4.5;
          else if(random_number<0.2) Iz = -3.5;
          else if(random_number<0.3) Iz = -2.5;
          else if(random_number<0.4) Iz = -1.5;
          else if(random_number<0.5) Iz = -0.5;
          else if(random_number<0.6) Iz =  0.5;
          else if(random_number<0.7) Iz =  1.5;
          else if(random_number<0.8) Iz =  2.5;
          else if(random_number<0.9) Iz =  3.5;
          else                       Iz =  4.5;
        }
        else { // for Ga or As atom case
          if(random_number<0.25)      Iz = -1.5; 
          else if(random_number<0.5)  Iz = -0.5; 
          else if(random_number<0.75) Iz =  0.5; 
          else                        Iz =  1.5; 
        }
     }
     else { //nuclear spins are polarized 
        if(_atom[g]==2) Iz=-4.5; else Iz=-1.5;
     }

      VV[g-g_min] = hyperfine_coupling*hyperfine_coupling*0.25; // in unit of erg*erg
      VV[g-g_min] *= II[_atom[g]]-Iz*(Iz+1); // in unit of erg*erg 
      AI_sum += hyperfine_coupling*Iz; // in unit of erg
//      cout<<"g, Iz, AIz "<<g<<" "<<Iz<<" "<<hyperfine_coupling*Iz<<endl;
  } 

  double total_AI_sum;
  MPI_Allreduce(&AI_sum, &total_AI_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double total_AI_sum_over_hbar = total_AI_sum/hbar; //in unit of 1/second
  double two_external_magnetic_field_omega =  2*bohr_magneton*MagneticField/hbar_in_eV_sec;
                                              //in unit of 1/second,  g*bohr_magneton*B/hbar
  double constant_omega = total_AI_sum_over_hbar+two_external_magnetic_field_omega;
  for(int g=0; g<g_num; g++){
    omega[g] += constant_omega; // in unit of 1/second
    VV[g] /=(hbar*hbar); // V*V/(hbar*hbar) in unit of 1/(second*second)  
  }
  double V_over_hbar_omega = 0.0;
  for(int g=0; g<g_num; g++){
    V_over_hbar_omega += VV[g]/(omega[g]*omega[g]);
  }
  double total_V_over_hbar_omega;
  MPI_Allreduce(&V_over_hbar_omega, &total_V_over_hbar_omega, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(my_id==0) {
    cout.setf(ios::scientific);
    cout<<setprecision(10);
    cout<<"Sum of AI: "<<total_AI_sum*convertor_erg_to_eV<<" eV"<<endl;
    cout<<"Zeeman Energy: "
        <<2*bohr_magneton*MagneticField+total_AI_sum*convertor_erg_to_eV<<" eV"<<endl;
    cout<<"Electron precession time: "<<2*M_PI/fabs(constant_omega)<<" second"<<endl;
    cout<<"The sum of the amplitude of cosine function : "<<total_V_over_hbar_omega<<endl;
  }


  double time_step = 0.2*M_PI/fabs(constant_omega); //in unit of second;
  double time_length = 1.0e-5; //in unit of second;
  int long num_time_step = long(time_length/time_step);
  ofstream* fout;
  if(my_id==0) {
    if(RandomOrientation) fout = new ofstream("Sz.unpolarized.coarse.dat");
    else                  fout = new ofstream("Sz.polarized.coarse.dat");
    *fout<<"#Total simulation time: "<<time_length<<endl;
    *fout<<"#Simulation time step: "<<time_step<<endl;
    *fout<<"#Number of time steps: "<<num_time_step<<endl;
    *fout<<"#Time (in second),  Sz (unitless)"<<endl;
  }
  for(int long it=0; it<num_time_step; it++) {
    double Sz=0.0;  // expectation value of electron spin z component, unitless.
    double time = it*time_step; // in unit of second;
    for(int g=0; g<g_num; g++){ 
//      Sz += VV[g] / (omega[g]*omega[g]) *  ( cos(omega[g]*time) - 1.0 ); // unitless
      Sz +=  cos(omega[g]*time) - 1.0 ; // unitless
    }
    double total_Sz;
    MPI_Allreduce(&Sz, &total_Sz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(my_id==0) *fout<<time<<" "<<total_Sz<<endl;
  }
  if(my_id==0) delete fout;

/*
  time_step = 0.2/fabs(constant_omega);
  num_time_step = long(time_length/time_step);
  ofstream* gout;
  if(my_id==0) {
    if(RandomOrientation) gout = new ofstream("Sz.unpolarized.fine.dat");
    else                  gout = new ofstream("Sz.polarized.fine.dat");
    *gout<<"#Total simulation time: "<<time_length<<endl;
    *gout<<"#Simulation time step: "<<time_step<<endl;
    *gout<<"#Number of time steps: "<<num_time_step<<endl;
    *gout<<"#Time (in second),  Sz (unitless)"<<endl;
  }
  for(int long it=0; it<num_time_step; it++) {
    double Sz=0.0;  // expectation value of electron spin z component, unitless.
    double time = it*time_step; // in unit of second;
    for(int g=0; g<g_num; g++){
      Sz += VV[g] / (omega[g]*omega[g]) *  ( cos(omega[g]*time) - 1.0 ); // unitless
    }
    double total_Sz;
    MPI_Allreduce(&Sz, &total_Sz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(my_id==0) *gout<<time<<" "<<total_Sz<<endl;
  }
  if(my_id==0) delete gout;
*/




}

