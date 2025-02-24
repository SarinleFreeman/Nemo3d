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

#include "FakeMPI.h"
#if (defined MPI3d && !defined FAKE_MPI)
#include <mpi.h>
#endif
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include "nml_dcvector.h"
#include "Coulomb.h"
#include "Dipole.h"
#include "Projection.h"
#include "ReadData.h"
#include "SimpleParser.h"
//#include "DispersionUnfold.h"
#include "FCCproj.h"
#include "CubicProj.h"
#include "FCCguess.h"
#include "CubicGuess.h"
using namespace std;

int main(int argc, char** argv) {

  int rank=0, size=1;
#if (defined MPI3d && !defined FAKE_MPI)
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  char data[123];
  int NumElectrons = 0;
  int NumHoles = 0;
  int NumAtoms = 0;
  int NumOrbitals = 0;
  int* eList;
  int* hList;

//Read Options =================================================
  vector<string> options;
  for(int i=0; i<argc; i++) options.push_back(argv[i]);
  unsigned int iopt = 0;
  bool ComputeCoulomb = false;
  bool ComputeCoulomb1=false;
  bool ComputeDipole = false;
  bool ComputeProjection = false;
//  bool ComputeUnfold = false;
  bool ComputeFCCproj = false;
  bool ComputeFCCguess = false;
  bool ComputeCubicProj = false;
  bool ComputeCubicGuess = false;

  while(iopt < options.size()) {
    if(options[iopt] == "--data") {
      strcpy(data,options[++iopt].c_str());
    }      
    else if(options[iopt] == "--electron") {
      NumElectrons = atoi(options[++iopt].c_str());
      eList = new int[NumElectrons];
      for(int ne=0; ne<NumElectrons; ne++){ 
         eList[ne] = atoi(options[++iopt].c_str());
      }
    }
    else if(options[iopt] == "--hole") {
      NumHoles = atoi(options[++iopt].c_str());
      hList = new int[NumHoles];
      for(int nh=0; nh<NumHoles; nh++) 
         hList[nh] = atoi(options[++iopt].c_str());
    }
    else if(options[iopt] == "--coulomb") {
      ComputeCoulomb1 = true;
    }
    else if(options[iopt] == "--hl") {
      ComputeCoulomb = true;
    }
    else if(options[iopt] == "--dipole") {
      ComputeDipole = true;
    }
    else if(options[iopt] == "--projection") {
      ComputeProjection = true;
    }
//    else if(options[iopt] == "--unfold") {
//      ComputeUnfold = true;
//    }

    else if(options[iopt] == "--fccproj") {
      ComputeFCCproj = true;
    }
    else if(options[iopt] == "--fccguess") {
      ComputeFCCguess = true;
    }

    else if(options[iopt] == "--cubicproj") {
      ComputeCubicProj = true;
    }
    else if(options[iopt] == "--cubicguess") {
      ComputeCubicGuess = true;
    }
    iopt++;
  }   

//  if (ComputeUnfold){
//    unfold_dispersion();
//    exit(-1);
//  }

   if (ComputeFCCproj){
    fccproj();
    exit(-1);
  }

  if (ComputeFCCguess){
    fccguess();
    exit(-1);
  }
  
   if (ComputeCubicProj){
    cubicproj();
    exit(-1);
  }

  if (ComputeCubicGuess){
    cubicguess();
    exit(-1);
  }


  cout.setf(ios::scientific);
  cout<<"NumElectrons, NumHoles "<<NumElectrons<<" " <<NumHoles<<endl;
  cout<<"eList: ";
  for(int ne=0; ne<NumElectrons; ne++) 
    cout<<eList[ne]<<" ";
  cout<<endl;
  cout<<"hList: ";
  for(int nh=0; nh<NumHoles; nh++) 
    cout<<hList[nh]<<" ";
  cout<<endl<<endl;
//===================================================================

//Read NEMO3D Output Files ==========================================  
  complex<double>* wf_e;
  complex<double>* wf_h;
  int* atomid;
  int* neighbor;
  double* lattice;
  double* Ee = new double[NumElectrons];
  double* Eh = new double[NumHoles];

  read_natoms_norbitals(data, &NumAtoms, &NumOrbitals);

  cout<<"Number of atoms: "<<NumAtoms<<endl;
  cout<<"Number of orbitals: "<<NumOrbitals<<endl;
//cout<<"Reached here"<<endl;

  read_data_for_postprocessing(data, NumElectrons, eList, NumHoles, hList,
                               &wf_e, &wf_h, &atomid, &neighbor, &lattice, &Ee, &Eh);
//cout<<"Assigned nml"<<endl;
 
  for(int iatom=0; iatom<NumAtoms; iatom++) {
    if(atomid[iatom]==7) atomid[iatom]=0; //Si atom
    else if(atomid[iatom]==8) atomid[iatom]=0; //P atom
    else cout<<"atomid does not match with Si, P atom"<<endl;
  }

/*
cout<<"Reached here"<<endl;
int dim_wf=NumElectrons*NumAtoms*20;
nml_dcvector *cf=nml_dcv_new(dim_wf);
for(int j=0;j<dim_wf;j++){
  nml_dscalar *r2_n_j = &(cf[j].r);
  *r2_n_j=wf_e[j].real();
  r2_n_j++; 	
  *r2_n_j=wf_e[j].imag();	
}
cout<<"Assigned nml"<<endl;
*/

//===================================================================

  int bat= NumAtoms/size;
  int residue = NumAtoms%size;
  int myAtoms_begin, myAtoms_end;
  if(rank<residue){
     myAtoms_begin = rank*(bat+1);
     myAtoms_end   = myAtoms_begin + bat + 1;
  }
  else {
     myAtoms_begin = residue*(bat+1) + (rank - residue)*bat;
     myAtoms_end   = myAtoms_begin + bat;
  }

  if(ComputeCoulomb) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb = 
	  Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    //ifstream coulomb_input("coulomb.table");
    //data_coulomb.assign_table(coulomb_input);
       data_coulomb.assign_table1();
	
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);

	FILE* fp=fopen("exchange_parameters","r");
	float dd_cutoff=12.0;
	double dd_cutoff1=dd_cutoff;
	int use_dd_cutoff=1;
	if(fp!=NULL){
		fscanf(fp, "%i\n", &use_dd_cutoff);
		fscanf(fp, "%f\n", &dd_cutoff);
		fclose(fp);
		dd_cutoff1=dd_cutoff;
	}else{
		if(rank==0)
			cout<<"File not found: Using default parameters. dd_cutoff=12 nm"<<endl;
	}
	if(rank==0){
		cout<<"Use dd_cutoff: "<<use_dd_cutoff<<endl;
		cout<<"dd_cutoff: "<<dd_cutoff1<<"  nm"<<endl;
	}


     complex<double> S=0.0;
     complex<double> S_sqr=0.0;

     data_coulomb.dotProduct(1, 2, &S);
     S_sqr=conj(S)*S;

    double elapsed_time=0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time=-MPI_Wtime();	
    if(use_dd_cutoff==0)	
    	data_coulomb.compute_exchange(1,2,1,2,rank);
    else
	data_coulomb.compute_exchange_cutoff(1,2,1,2,rank, dd_cutoff1);	
   elapsed_time+=MPI_Wtime();

    complex<double>* local_coul = data_coulomb.get_H_coul();
    complex<double>* local_exch = data_coulomb.get_H_exch();
    int csize = 1;
    complex<double>* total_coul;
    total_coul = new complex<double>[csize]; 
    complex<double>* total_exch;
    total_exch = new complex<double>[csize]; 
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(local_coul, total_coul, csize*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce(local_exch, total_exch, csize*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
#else 
    for(int cs=0; cs<csize; cs++) {
       total_coul[cs] = local_coul[cs];
       total_exch[cs] = local_exch[cs];
    }
#endif
    if (rank==0) {
      cout<<endl<<"====================================================="<<endl;
          cout<<"Electron and Hole Pair: "<<endl;
          cout<<"Coulomb:  "<<setw(10)<<total_coul[0]<<endl;
          cout<<"Exchange: "<<setw(10)<<total_exch[0]<<endl;
	  cout<<"Overlap_square:  "<<S_sqr<<endl;
          cout<<endl;
	  cout<<"Time to compute a matrix element:"<<elapsed_time<<endl;
   	  complex<double> exchange=(2/(1-S_sqr.real()*S_sqr.real()))*(S_sqr*total_coul[0]-total_exch[0].real());
	  cout<<"Exchange Energy="<<exchange<<endl;
	 }
 }

 
 if(ComputeCoulomb1) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    ifstream coulomb_input("coulomb.table");
    data_coulomb.assign_table(coulomb_input);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double elapsed_time=0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time=-MPI_Wtime();
    data_coulomb.compute_exchange(1,2,1,2, rank);
    elapsed_time+=MPI_Wtime();

    complex<double>* local_coul = data_coulomb.get_H_coul();
    complex<double>* local_exch = data_coulomb.get_H_exch();
    int csize = 1;
    complex<double>* total_coul;
    total_coul = new complex<double>[csize];
    complex<double>* total_exch;
    total_exch = new complex<double>[csize];
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(local_coul, total_coul, csize*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_exch, total_exch, csize*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for(int cs=0; cs<csize; cs++) {
       total_coul[cs] = local_coul[cs];
       total_exch[cs] = local_exch[cs];
    }
#endif
    if (rank==0) {
      cout<<endl<<"====================================================="<<endl;
          cout<<"Electron and Hole Pair: "<<endl;
          cout<<"Coulomb:  "<<setw(10)<<total_coul[0]<<endl;
          cout<<"Exchange: "<<setw(10)<<total_exch[0]<<endl;
          cout<<endl;
          cout<<"Time to compute a matrix element:"<<elapsed_time<<endl;
         }
 }

 
  if(ComputeDipole) {
    //Computing Dipole moments =========================================
    Dipole data_dipole = 
          Dipole(NumElectrons, NumHoles, NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    ifstream dipole_input("dipole.table");
    ofstream dipole_output("dipole.output");
    dipole_output.setf(ios::scientific);
    data_dipole.assign_table(dipole_input);
    data_dipole.get_data(wf_e, wf_h, atomid, neighbor, lattice);
    data_dipole.compute_dipole();
    complex<double>* local_dipole_x = data_dipole.get_H_dipole_X();
    complex<double>* local_dipole_y = data_dipole.get_H_dipole_Y();
    complex<double>* local_dipole_z = data_dipole.get_H_dipole_Z();
    int dsize = NumElectrons*NumHoles;
    complex<double>* total_dipole_x;
    complex<double>* total_dipole_y;
    complex<double>* total_dipole_z;
    total_dipole_x = new complex<double>[dsize];
    total_dipole_y = new complex<double>[dsize];
    total_dipole_z = new complex<double>[dsize];
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(local_dipole_x, total_dipole_x, dsize*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce(local_dipole_y, total_dipole_y, dsize*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce(local_dipole_z, total_dipole_z, dsize*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
#else 
    for(int ds=0; ds<dsize; ds++){
      total_dipole_x[ds] = local_dipole_x[ds];
      total_dipole_y[ds] = local_dipole_y[ds];
      total_dipole_z[ds] = local_dipole_z[ds];
    }
#endif
    if (rank==0) {
      cout<<endl<<"====================================================="<<endl;
      for(int ne=0; ne<NumElectrons; ne++){
        for(int nh=0; nh<NumHoles; nh++){ 
          complex<double> dipole_moment_x = total_dipole_x[NumHoles*ne+nh];
          complex<double> dipole_moment_y = total_dipole_y[NumHoles*ne+nh];
          complex<double> dipole_moment_z = total_dipole_z[NumHoles*ne+nh];
          double oscillator_strength_x = (dipole_moment_x* conj(dipole_moment_x)).real()
                                     * (Ee[ne]-Eh[nh]);
          double oscillator_strength_y = (dipole_moment_y* conj(dipole_moment_y)).real()
                                     * (Ee[ne]-Eh[nh]);
          double oscillator_strength_z = (dipole_moment_z* conj(dipole_moment_z)).real()
                                     * (Ee[ne]-Eh[nh]);
          cout<<"Electron and Hole Pair: "<<setw(5)<<ne+1<<setw(5)<<nh+1<<endl;
          cout<<"Dipole moment in x direction:   "<<setw(10)<<dipole_moment_x<<endl;
          cout<<"Dipole moment in y direction:   "<<setw(10)<<dipole_moment_y<<endl;
          cout<<"Dipole moment in z direction:   "<<setw(10)<<dipole_moment_z<<endl;
          cout<<"Oscillator Strength for x-direction linearly polarized light: "
            <<setw(10)<<oscillator_strength_x<<endl;
          cout<<"Oscillator Strength for y-direction linearly polarized light: "
            <<setw(10)<<oscillator_strength_y<<endl;
          cout<<"Oscillator Strength for z-direction linearly polarized light: "
            <<setw(10)<<oscillator_strength_z<<endl;
          cout<<endl;
          dipole_output<<setw(5)<<ne+1<<setw(5)<<nh+1
                       <<setw(15)<<Ee[ne]-Eh[nh]
                       <<setw(15)<<oscillator_strength_x
                       <<setw(15)<<oscillator_strength_y
                       <<setw(15)<<oscillator_strength_z
                       <<endl;
        }
      }  
    }
  }

  if(ComputeProjection) {
    //Projecting wavefunctions onto bulk-band Bloch functions ===========
    Projection data_projection = 
      Projection(NumElectrons, NumHoles, NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_projection.get_data(wf_e, wf_h, atomid, lattice);
    data_projection.compute_projection();
    double* local_projection = data_projection.get_projection();
    int psize = (NumElectrons+NumHoles)*8;
    double  total_projection[psize];
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(local_projection, total_projection, psize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
#else 
    for(int ps=0; ps<psize; ps++){
      total_projection[ps] = local_projection[ps];
    }
#endif
    if (rank==0) {
      cout<<endl<<"====================================================="<<endl;
      for(int n=0; n<NumElectrons+NumHoles; n++){
        if(n<NumElectrons) cout<<"Electron Level: "<<setw(5)<<n+1<<", Energy: "<<setw(10)<<Ee[n]<<endl;
        else cout<<"Hole Level: "<<setw(5)<<n-NumElectrons+1<<", Energy: "<<setw(10)<<Eh[n-NumElectrons]<<endl;
        for(int p=0; p<8; p++) {
          if(p==0) cout<<"Projection onto s* band:          "<<setw(10)<<total_projection[n*8]<<endl;
          if(p==1) cout<<"Projection onto s  band:          "<<setw(10)<<total_projection[n*8+1]<<endl;
          if(p==2) cout<<"Projection onto heavy-hole band:  "<<setw(10)<<total_projection[n*8+2]<<endl;
          if(p==3) cout<<"Projection onto light-hole band:  "<<setw(10)<<total_projection[n*8+3]<<endl;
          if(p==4) cout<<"Projection onto split-off  band:  "<<setw(10)<<total_projection[n*8+4]<<endl;
          if(p==5) cout<<"Projection onto d  band:          "<<setw(10)<<total_projection[n*8+5]<<endl;
          if(p==6) cout<<"Projection onto spin up:          "<<setw(10)<<total_projection[n*8+6]<<endl;
          if(p==7) cout<<"Projection onto spin down:        "<<setw(10)<<total_projection[n*8+7]<<endl;
        }
        cout<<endl;
      }  
    }
  }
 
  delete [] eList; delete [] hList; delete [] Ee; delete [] Eh;
#if (defined MPI3d && !defined FAKE_MPI)
  MPI_Finalize();
#endif
  return 0;
}	
