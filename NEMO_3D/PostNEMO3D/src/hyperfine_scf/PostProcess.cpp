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
#include <math.h>
#include <cstring>
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
  bool ComputeAHF=false;
  bool ComputeSlice=false;
  bool ComputeAHF_wf = false;
  bool Compute_SCF = false;
  bool Compute_SCF_D2 = false;
  bool Compute_SCF_D3 = false;
  bool Compute_SCF_D2_single = false;
  bool Compute_SCF_D3_single = false;
  bool Compute_SCF_image = false;
  bool Compute_SCF_CI=false;
  bool Compute_SCF_avg=false;
  bool Compute_avg = false;
  bool Compute_SCF_Elect=false;
  bool Compute_proj =false;
  bool Compute_add = false;
  bool Compute_AHF_shell=false; 
  bool Compute_HF_HF=false;
  bool save_electron_density=false;
  bool Compute_Dot = false;
  bool Compute_x2 = false;
  bool Compute_moment=false;
  bool Compute_dipole_moment=false;
  bool Compute_r2=false;
  bool Compute_moment1=false;
  bool Compute_moment2=false;
  bool Compute_moment3=false;
  bool Compute_cut=false;
  bool Compute_LS=false;
  bool Compute_CTAP=false;
  bool Compute_SiGe=false;
  bool charge=false;
  bool Compute_SCF_new=false;
  bool GenerateSpin_wfs=false;
  bool Compute_Chandra=false;
  bool Compute_wf_diff=false; 
  bool Compute_E_field=false;
  bool skip=false;	
  bool normalize=false; 
  int num_Electrons=0;

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
   else if(options[iopt] == "--ahf") {
      ComputeAHF = true;
    }
   else if(options[iopt] == "--normalize") {
      normalize = true;
    }	
   else if(options[iopt] == "--ahf_wf") {
      ComputeAHF_wf = true;
    }
   else if(options[iopt] == "--ahf_shell") {
      Compute_AHF_shell = true;
    }
   else if(options[iopt] == "--hf_hf") {
      Compute_HF_HF = true;
    }
   else if(options[iopt] == "--skip") {
      skip = true;
    }
   else if(options[iopt] == "--wf_diff") {
      Compute_wf_diff = true;
    }	
   else if(options[iopt] == "--spin_wfs") {
      GenerateSpin_wfs = true;
    }
    else if(options[iopt] == "--charge") {
      charge = true;
    }
   else if(options[iopt] == "--slice") {
      ComputeSlice = true;
    }
   else if(options[iopt] == "--scf") {
      Compute_SCF = true;
    }
   else if(options[iopt] == "--scf_d2") {
      Compute_SCF_D2 = true;
    }
   else if(options[iopt] == "--scf_d3") {
      Compute_SCF_D3 = true;
    }
   else if(options[iopt] == "--scf_elect") {
      Compute_SCF_Elect = true;
	cout<<"2017 Hartree calculation set"<<endl;
    }
   else if(options[iopt] == "--elect") {
      num_Electrons = atoi(options[++iopt].c_str());
	cout<<"Number of electrons: "<<num_Electrons<<endl;
    }
    else if(options[iopt] == "--scf_new") {
      Compute_SCF_new = true;
    }
    else if(options[iopt] == "--chandra") {
      Compute_Chandra = true;
    }		
   else if(options[iopt] == "--scf_d2_single") {
      Compute_SCF_D2_single = true;
    }
   else if(options[iopt] == "--scf_d3_single") {
      Compute_SCF_D3_single = true;
    }
   else if(options[iopt] == "--scf_image") {
      Compute_SCF_image = true;
    }
   else if(options[iopt] == "--dot") {
      Compute_Dot = true;
    }
   else if(options[iopt] == "--scf_ci") {
      Compute_SCF_CI = true;
    }
   else if(options[iopt] == "--x2") {
      Compute_x2 = true;
    }
   else if(options[iopt] == "--moment") {
      Compute_moment = true;
    }
   else if(options[iopt] == "--dipole_moment") {
      Compute_dipole_moment = true;
    }
   else if(options[iopt] == "--moment1") {
      Compute_moment1 = true;
    }
    else if(options[iopt] == "--moment2") {
      Compute_moment2 = true;
    }
    else if(options[iopt] == "--moment3") {
      Compute_moment3 = true;
    }	
   else if(options[iopt] == "--r2") {
      Compute_r2 = true;
    }		
   else if(options[iopt] == "--cut") {
      Compute_cut = true;
    } 
  else if(options[iopt] == "--nr") {
      save_electron_density = true;
   }
   else if(options[iopt] == "--scf_avg") {
      Compute_SCF_avg = true;
    }
   else if(options[iopt] == "--avg") {
      Compute_avg = true;
    }
      else if(options[iopt] == "--add") {
      Compute_add = true;
    }	
    else if(options[iopt] == "--dipole") {
      ComputeDipole = true;
    }
    else if(options[iopt] == "--LS") {
      Compute_LS = true;
    }
    else if(options[iopt] == "--ctap") {
      Compute_CTAP = true;
    }
   else if(options[iopt] == "--sige") {
      Compute_SiGe = true;
    }
   else if(options[iopt] == "--e_field") {
      Compute_E_field = true;
    }
    else if(options[iopt] == "--proj") {
      Compute_proj = true;
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
 int counter=0;
  for(int iatom=0; iatom<NumAtoms; iatom++) {
    if(atomid[iatom]==7) atomid[iatom]=0; //Si atom
    else if(atomid[iatom]==8){ 
		counter++;
		atomid[iatom]=1;
    }
    else if(atomid[iatom]==11) atomid[iatom]=2; //Ge atom
    else cout<<"atomid does not match with Si, P atom"<<endl;
  }
  cout<<"Counted number of P atoms : "<<counter<<endl;
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

  if(Compute_proj) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    //ifstream coulomb_input("coulomb.table");
    //data_coulomb.assign_table(coulomb_input);
      data_coulomb.assign_table1();

     data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
     complex<double> S=0.0;
     complex<double> S_sqr=0.0;
     data_coulomb.dotProduct(1, 2, &S);
     S_sqr=conj(S)*S;
     cout<<"Dot product between wfs "<<eList[0]<<"    "<<eList[1]<<endl;	
     cout<<S<<endl;
     cout<<S_sqr<<endl;		
}	

if(Compute_E_field) {

  cout<<"Number of atoms : "<<NumAtoms<<endl; 
  for(int i=0; i<6;i++){
	if(i<5){
	cout<<"Atom id of "<<i<<"   :  "<<atomid[i]<<endl;
	cout<<"Lattice coordinates : ("<<lattice[i*3]<<", "<<lattice[i*3+1]<<", "<<lattice[i*3+2]<<")"<<endl;
	cout<<"Nearest neighbors : ("<<neighbor[i*4]<<", "<<neighbor[i*4+1]<<", "<<neighbor[i*4+2]<<", "<<neighbor[i*4+3]<<")"<<endl;
}else{
int j=NumAtoms-1;

cout<<"Atom id of "<<j<<"  :  "<<atomid[j]<<endl;        

cout<<"Lattice coordinates : ("<<lattice[j*3]<<", "<<lattice[j*3+1]<<", "<<lattice[j*3+2]<<")"<<endl;        

cout<<"Nearest neighbors : ("<<neighbor[j*4]<<", "<<neighbor[j*4+1]<<", "<<neighbor[j*4+2]<<", "<<neighbor[j*4+3]<<")"<<endl;

}
}
 	double* V_t=new double[NumAtoms];
	cout<<"Reading potential : "<<endl;

        FILE* fp=fopen("V_total","r");
        double pot=0.0;
        if(fp!=NULL){
		for(int i=0;i<NumAtoms;i++){
                	fscanf(fp, "%lf\n", &pot);
			V_t[i]=pot;
			pot=0.0;
		}
		cout<<"Printing potential"<<endl;
		for(int i=0;i<6;i++)
			if(i<5)
				cout<<V_t[i]<<endl;
			else
				cout<<V_t[NumAtoms-1]<<endl;	
        }else{
               cout<<"Potential File not found"<<endl;
        }
	fclose(fp);

//begin computing E-field
int nn1=0;
int nn2=0;
int nn3=0;
int nn4=0;
double dx1=0.0;
double dx2=0.0;
double dx3=0.0;
double dx4=0.0;
double dy1=0.0;
double dy2=0.0;
double dy3=0.0;
double dy4=0.0;
double dz1=0.0;
double dz2=0.0;
double dz3=0.0;
double dz4=0.0;
double dVx1=0.0;
double dVy1=0.0;
double dVz1=0.0;
double dVx2=0.0;
double dVy2=0.0;
double dVz2=0.0;
double dVx3=0.0;
double dVy3=0.0;
double dVz3=0.0;
double dVx4=0.0;
double dVy4=0.0;
double dVz4=0.0;

double x0=0.0;
double y0=0.0;
double z0=0.0;
double x1=0.0;
double y1=0.0;
double z1=0.0;
double x2=0.0;
double y2=0.0;
double z2=0.0;
double x3=0.0;
double y3=0.0;
double z3=0.0;
double x4=0.0;
double y4=0.0;
double z4=0.0;
double V0=0.0;
double V1=0.0;
double V2=0.0;
double V3=0.0;
double V4=0.0;
int counter=0;

double Fx1=0.0;
double Fy1=0.0;
double Fz1=0.0;
double Fx2=0.0;
double Fy2=0.0;
double Fz2=0.0;
double Fx3=0.0;
double Fy3=0.0;
double Fz3=0.0;
double Fx4=0.0;
double Fy4=0.0;
double Fz4=0.0;

double F_x=0.0;
double F_y=0.0;
double F_z=0.0;
double F_r=0.0;

 FILE *fp1;
 fp1=fopen("Fx", "w");
  if (fp1==NULL)
    {
      printf("Can't open input file 1\n");
      exit(-1);
    }

 FILE *fp2;
 fp2=fopen("Fy", "w");
  if (fp2==NULL)
    {
      printf("Can't open input file 1\n");
      exit(-1);
    } 

 FILE *fp3;
 fp3=fopen("Fz", "w");
  if (fp3==NULL)
    {
      printf("Can't open input file 1\n");
      exit(-1);
    }

FILE *fp4;
 fp4=fopen("Fr", "w");
  if (fp4==NULL)
    {
      printf("Can't open input file 1\n");
      exit(-1);
    }
 

for(int i=0;i<NumAtoms;i++){
	x0=lattice[3*i];
	y0=lattice[3*i+1];
	z0=lattice[3*i+2];
	V0=V_t[i];

	nn1=neighbor[i*4];
	nn2=neighbor[i*4+1];
	nn3=neighbor[i*4+2];
	nn4=neighbor[i*4+3];

	dVx1=0.0;
	dVy1=0.0;
	dVz1=0.0;
	dVx2=0.0;
	dVy2=0.0;
	dVz2=0.0;
	dVx3=0.0;
	dVy3=0.0;
	dVz3=0.0;
	dVx4=0.0;
	dVy4=0.0;
	dVz4=0.0;
	counter=0;
	Fx1=0.0;
	Fy1=0.0;
	Fz1=0.0;
	Fx2=0.0;
        Fy2=0.0;
        Fz2=0.0;
	Fx3=0.0;
        Fy3=0.0;
        Fz3=0.0;
	Fx4=0.0;
        Fy4=0.0;
        Fz4=0.0;
	F_x=0.0;
	F_y=0.0;
	F_z=0.0;
	F_r=0.0;

	if(nn1>=0){
		x1=lattice[3*nn1];
		y1=lattice[3*nn1+1];
		z1=lattice[3*nn1+2];
		V1=V_t[nn1];
		if(x0>x1){
		      dx1=x0-x1;
		      dVx1=V0-V1;	
		}else{
		      dx1=x1-x0;
		      dVx1=V1-V0; 		
		}
		Fx1=-dVx1/dx1;
		if(y0>y1){
                      dy1=y0-y1;
                      dVy1=V0-V1;
                }else{
                      dy1=y1-y0;
                      dVy1=V1-V0;
                }
		Fy1=-dVy1/dy1;
		if(z0>z1){
                      dz1=z0-z1;
                      dVz1=V0-V1;
                }else{
                      dz1=z1-z0;
                      dVz1=V1-V0;
                }
		Fz1=-dVz1/dz1;
		counter++;
	}
	if(nn2>=0){
                x2=lattice[3*nn2];
                y2=lattice[3*nn2+1];
                z2=lattice[3*nn2+2];
                V2=V_t[nn2];
                if(x0>x2){
                      dx2=x0-x2;
                      dVx2=V0-V2;
                }else{
                      dx2=x2-x0;
                      dVx2=V2-V0;
                }
		Fx2=-dVx2/dx2;
                if(y0>y2){
                      dy2=y0-y2;
                      dVy2=V0-V2;
                }else{
                      dy2=y2-y0;
                      dVy2=V2-V0;
                }
		Fy2=-dVy2/dy2;
                if(z0>z2){
                      dz2=z0-z2;
                      dVz2=V0-V2;
                }else{
                      dz2=z2-z0;
                      dVz2=V2-V0;
                }
		Fz2=-dVz2/dz2;
                counter++;
        }
	if(nn3>=0){
                x3=lattice[3*nn3];
                y3=lattice[3*nn3+1];
                z3=lattice[3*nn3+2];
                V3=V_t[nn3];
                if(x0>x3){
                      dx3=x0-x3;
                      dVx3=V0-V3;
                }else{
                      dx3=x3-x0;
                      dVx3=V3-V0;
                }
		Fx3=-dVx3/dx3;
                if(y0>y3){
                      dy3=y0-y3;
                      dVy3=V0-V3;
                }else{
                      dy3=y3-y0;
                      dVy3=V3-V0;
                }
		Fy3=-dVy3/dy3;
                if(z0>z3){
                      dz3=z0-z3;
                      dVz3=V0-V3;
                }else{
                      dz3=z3-z0;
                      dVz3=V3-V0;
                }
		Fz3=-dVz3/dz3;
                counter++;
        }
 
	if(nn4>=0){
                x4=lattice[3*nn4];
                y4=lattice[3*nn4+1];
                z4=lattice[3*nn4+2];
                V4=V_t[nn4];
                if(x0>x4){
                      dx4=x0-x4;
                      dVx4=V0-V4;
                }else{
                      dx4=x4-x0;
                      dVx4=V4-V0;
                }
		Fx4=-dVx4/dx4;
                if(y0>y4){
                      dy4=y0-y4;
                      dVy4=V0-V4;
                }else{
                      dy4=y4-y0;
                      dVy4=V4-V0;
                }
      		Fy4=-dVy4/dy4; 
                if(z0>z4){
                      dz4=z0-z4;
                      dVz4=V0-V4;
                }else{
                      dz4=z4-z0;
                      dVz4=V4-V0;
                }
		Fz4=-dVz4/dz4;
                counter++;
        }


/*cout<<"Number of nearest neighbors : "<<counter<<endl;
cout<<"Atom index : "<<i<<endl;
cout<<"Nearest neighbors : "<<nn1<<"  "<<nn2<<"  "<<nn3<<"  "<<nn4<<endl;
cout<<"Potential : "<<V0<<"  "<<V1<<"  "<<V2<<"  "<<V3<<"  "<<V4<<endl;
cout<<"dVx : "<<dVx1<<"  "<<dVx2<<"  "<<dVx3<<"  "<<dVx4<<endl;
cout<<"dx : "<<dx1<<"  "<<dx2<<"  "<<dx3<<"  "<<dx4<<endl;
cout<<"Fx : "<<Fx1<<"  "<<Fx2<<"  "<<Fx3<<"  "<<Fx4<<endl;
cout<<"dVy : "<<dVy1<<"  "<<dVy2<<"  "<<dVy3<<"  "<<dVy4<<endl;
cout<<"dy : "<<dy1<<"  "<<dy2<<"  "<<dy3<<"  "<<dy4<<endl;
cout<<"Fy : "<<Fy1<<"  "<<Fy2<<"  "<<Fy3<<"  "<<Fy4<<endl;
cout<<"dVz : "<<dVz1<<"  "<<dVz2<<"  "<<dVz3<<"  "<<dVz4<<endl;
cout<<"dz : "<<dz1<<"  "<<dz2<<"  "<<dz3<<"  "<<dz4<<endl;
cout<<"Fz : "<<Fz1<<"  "<<Fz2<<"  "<<Fz3<<"  "<<Fz4<<endl;
*/
F_x=(Fx1+Fx2+Fx3+Fx4)/counter;
F_y=(Fy1+Fy2+Fy3+Fy4)/counter;
F_z=(Fz1+Fz2+Fz3+Fz4)/counter;
F_r=sqrt(F_x*F_x+F_y*F_y+F_z*F_z);
fprintf(fp1,"%12.10le\n",F_x);
fprintf(fp2,"%12.10le\n",F_y);
fprintf(fp3,"%12.10le\n",F_z);
fprintf(fp4,"%12.10le\n",F_r);
//cout<<"Average field  Fx : "<<F_x<<"   Fy : "<<F_y<<" Fz : "<<F_z<<endl;
}
fclose(fp1);
fclose(fp2);
fclose(fp3);
fclose(fp4);
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

if(ComputeSlice) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb = Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    //ifstream coulomb_input("coulomb.table");
    //data_coulomb.assign_table(coulomb_input);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *r_init=new double[3];
    double a0=0.543095;
    r_init[0]=35;
    r_init[1]=35;
    r_init[2]=0;
    double *dr=new double[3];
    dr[0]=0.1;
    dr[1]=0.1;
    dr[2]=0.1;
    data_coulomb.oneDSlice(0,0,1,1,r_init,dr);
}

if(normalize){
	Coulomb data_coulomb = Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
        data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
	double sum=0.0;
	int n=0;
	int next_index=20*NumAtoms;
	for(int k=0; k<NumElectrons;k++){
		sum=0.0;
		n=k*next_index;	
        	for(int x=0;x<NumAtoms;x++){
                	for(int s=0;s<2;s++)
                	for(int o=0;o<10;o++){
                        	sum+=(wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag());   
                        	n++;
                	}
        	}
	cout<<"Normalization factor of wf "<<k<<" : "<<sum<<endl;
	}
}


if(Compute_wf_diff){
FILE *fp1=fopen("wf_diff","w");
if(fp1==NULL){
	cout<<"Error opening file for output"<<endl;
	exit(-1);
}else{
double* wf_1=new double[NumAtoms];
double* wf_2=new double[NumAtoms];
double* wf_12=new double[NumAtoms];
read_wfs_prob(NumAtoms, data, 1 , wf_1);
read_wfs_prob(NumAtoms, data, 2 , wf_2);
double sum1=0.0;
double sum2=0.0;
for(int i=0;i<NumAtoms;i++){
	sum1+=wf_1[i];
	sum2+=wf_2[i];
}
cout<<"sum1 = "<<sum1<<endl;
cout<<"sum2 = "<<sum2<<endl;
for(int i=0;i<NumAtoms;i++){
        wf_1[i]=wf_1[i]/sum1;
        wf_2[i]=wf_2[i]/sum2;
	wf_12[i]=wf_2[i]-wf_1[i];
}

}
}

if(Compute_avg) {

FILE* fp1=fopen("HF0","r");
FILE* fp2=fopen("HF1","r");
FILE* fp3=fopen("HF","w");

double a1=0.0;
double a2=0.0;
double a3=0.0;
if((fp1!=NULL)&&(fp2!=NULL)&&(fp3!=NULL)){
	for(int i=0;i<NumAtoms;i++){
		a1=a2=a3=0.0;
		fscanf(fp1,"%lf\n",&a1);
		fscanf(fp2,"%lf\n",&a2);
	 	a3=(a1+a2)/2.0;
	 	fprintf(fp3,"%12.10le\n", a3);	
		if(i<5)
			cout<<"a1 = "<<a1<<"    a2 = "<<a2<<"    a3 = "<<a3<<endl;
				
	}
fclose(fp1);
fclose(fp2);
fclose(fp3);
}else{
   if(rank==0)
     	cout<<"Error opening file."<<endl;
}
}

if(Compute_add) {

FILE* fp1=fopen("HF","r");
FILE* fp2=fopen("TCAD_atomic_ref.data","r");
FILE* fp3=fopen("potential","w");

double a1=0.0;
double a2=0.0;
double a3=0.0;
if((fp1!=NULL)&&(fp2!=NULL)&&(fp3!=NULL)){
        for(int i=0;i<NumAtoms;i++){
                a1=a2=a3=0.0;
                fscanf(fp1,"%lf\n",&a1);
                fscanf(fp2,"%lf\n",&a2);
                a3=a1+a2;
                fprintf(fp3,"%12.10le\n", a3);
                if(i<5)
                        cout<<"a1 = "<<a1<<"    a2 = "<<a2<<"    a3 = "<<a3<<endl;

        }
fclose(fp1);
fclose(fp2);
fclose(fp3);
}else{
   if(rank==0)
        cout<<"Error opening file."<<endl;
}
}

if(Compute_Dot){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    complex<double> dot_result=0.0;
    data_coulomb.dotProduct(1,2,&dot_result);	
    cout<<dot_result.real()<<endl;
    cout<<dot_result.imag()<<endl;
    double mag=dot_result.real()*dot_result.real()+dot_result.imag()*dot_result.imag();
    cout<<mag<<endl;				   	
}
if(Compute_Chandra){
 Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *nr = new double[NumAtoms];
    double *nr_total = new double[NumAtoms];
    double *wf_total = new double[NumAtoms];
    double *wf = new double[NumAtoms];

    data_coulomb.density_Chandra(nr, wf);

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(nr, nr_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(wf, wf_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                nr_total[i]=nr[i];
		wf_total[i]=wf[i];	
        }
#endif

if(rank==0){
    int n=0;
    double max1=-1;
    double max2=-1;
        
for(int x=0;x<NumAtoms;x++){
        if(wf_total[x]>max2)
                max2=wf_total[x];
        if(nr_total[x]>max1)
                max1=nr_total[x];
}
for(int x=0;x<NumAtoms;x++){
        nr_total[x]=nr_total[x]/max1;
	wf_total[x]=wf_total[x]/max2;
}
	cout<<"Finished calculating Chandrashekhar wavefunction density"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("nr_chandra", "w");
        for(int i=0;i<NumAtoms;i++){
                fprintf(fp1,"%12.10le\n", nr_total[i]);
        }
        fclose(fp1);
	FILE* fp2=fopen("wf_chandra", "w");
        for(int i=0;i<NumAtoms;i++){
                fprintf(fp2,"%12.10le\n", wf_total[i]);
        }
        fclose(fp2);
}//end of rank
}

if(Compute_moment1){

int basis=9;
double xP=0.0;
double yP=0.0;
double zP=0.0;
double x=0.0;
double y=0.0;
double z=0.0;

FILE* fp=fopen("coordinates","r");
if(fp!=NULL){
        fscanf(fp,"%lf %lf %lf",&xP, &yP, &zP);
        cout<<"Coordinates of the impurity : "<<xP<<"    "<<yP<<"     "<<zP<<endl;
}else{
        cout<<"Error opening file. "<<endl;
        exit(-1);
}
fclose(fp);

for(int i=0;i<basis;i++){	
    	complex<double> *wf1 = new complex<double>[20*NumAtoms]; 
    	read_wfs1(NumAtoms, data, i+1 , wf1);

	complex<double> dipole_x=0.0;
	complex<double> dipole_y=0.0;
	complex<double> dipole_z=0.0;

	for(int i=0;i<NumAtoms;i++){
        	complex<double> tmp=0.0;
        	x=lattice[3*i]-xP;
        	y=lattice[3*i+1]-yP;
        	z=lattice[3*i+2]-zP;
        	for(int s=0;s<2;s++)
                	for(int o=0;o<10;o++)
                        	tmp+=conj(wf1[i*20+s*10+o])*wf1[i*20+s*10+o];
        	dipole_x+=x*tmp;
        	dipole_y+=y*tmp;
        	dipole_z+=z*tmp;
	}
    	cout<<"X   :   "<<dipole_x<<endl;
    	cout<<"Y   :   "<<dipole_y<<endl;
    	cout<<"Z   :   "<<dipole_z<<endl;

	FILE* fp1=fopen("moment_x","a");
	FILE* fp2=fopen("moment_y","a");
	FILE* fp3=fopen("moment_z","a");
	if((fp1!=NULL)||(fp2!=NULL)||(fp3!=NULL)){
        	fprintf(fp1,"%12.10le\n", dipole_x.real());
        	fprintf(fp2,"%12.10le\n", dipole_y.real());
        	fprintf(fp3,"%12.10le\n", dipole_z.real());
	}else{
        	cout<<"Error opening file. "<<endl;
        	exit(-1);
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	delete[] wf1;
}//end of basis
}

if(Compute_moment3){

int basis=1;
double xP=0.0;
double yP=0.0;
double zP=0.0;
double x=0.0;
double y=0.0;
double z=0.0;

FILE* fp=fopen("coordinates","r");
if(fp!=NULL){
        fscanf(fp,"%lf %lf %lf",&xP, &yP, &zP);
        cout<<"Coordinates of the impurity : "<<xP<<"    "<<yP<<"     "<<zP<<endl;
}else{
        cout<<"Error opening file. "<<endl;
        exit(-1);
}
fclose(fp);

for(int i=0;i<basis;i++){
        double* wf1 = new double[NumAtoms];
        read_wfs_prob(NumAtoms, data, i+1 , wf1);

        double dipole_x=0.0;
        double dipole_y=0.0;
        double dipole_z=0.0;
	double rx2=0.0;
        double ry2=0.0;
        double rz2=0.0;

        for(int i=0;i<NumAtoms;i++){
                double tmp=0.0;
                x=lattice[3*i]-xP;
                y=lattice[3*i+1]-yP;
                z=lattice[3*i+2]-zP;
                dipole_x+=x*wf1[i];
                dipole_y+=y*wf1[i];
                dipole_z+=z*wf1[i];
		rx2+=sqrt(y*y+z*z)*wf1[i];
                ry2+=sqrt(x*x+z*z)*wf1[i];
                rz2+=sqrt(x*x+y*y)*wf1[i];
        }
        cout<<"X   :   "<<dipole_x<<endl;
        cout<<"Y   :   "<<dipole_y<<endl;
        cout<<"Z   :   "<<dipole_z<<endl;
	cout<<"r_x   :   "<<rx2<<endl;
        cout<<"r_y   :   "<<ry2<<endl;
        cout<<"r_z   :   "<<rz2<<endl;

        FILE* fp1=fopen("moment","a");
        if(fp1!=NULL){
                fprintf(fp1,"%12.10le\n", dipole_x);
                fprintf(fp1,"%12.10le\n", dipole_y);
                fprintf(fp1,"%12.10le\n", dipole_z);
		fprintf(fp1,"%12.10le\n", rx2);
                fprintf(fp1,"%12.10le\n", ry2);
                fprintf(fp1,"%12.10le\n", rz2);
        }else{
                cout<<"Error opening file. "<<endl;
                exit(-1);
        }
        fclose(fp1);
        delete[] wf1;
}//end of basis
}

if(Compute_cut){

int basis=9;
double xP=0.0;
double yP=0.0;
double zP=0.0;
double x=0.0;
double y=0.0;
double z=0.0;

double x0=0.0;
double z0=0.0;
double dr=0.1;

FILE* fp2=fopen("slope","r");
if(fp2!=NULL){
        fscanf(fp2,"%lf %lf",&x0, &z0);
        cout<<"Coordinates of the origin : "<<x0<<"    "<<z0<<endl;
}else{
        cout<<"Error opening file. "<<endl;
        exit(-1);
}
fclose(fp2);
FILE* WF=fopen("WF","w");
if(WF==NULL){
	cout<<"Error opening file"<<endl;
	exit(-1);
}
for(int i=0;i<basis;i++){
        double* wf1 = new double[NumAtoms];
        read_wfs_prob(NumAtoms, data, i+1 , wf1);
	cout<<"WF "<<i+1<<endl;
	fprintf(WF,"\n%i\n",i);
	for(int j=0;j<NumAtoms;j++){
		x=lattice[3*j];
		y=lattice[3*j+1];
		z=lattice[3*j+2];
		if((x<=x0+dr)&&(x>=x0-dr)&&(z<=z0+dr)&&(z>=z0-dr)){
			cout<<x<<"  "<<y<<"  "<<z<<"  "<<wf1[j]<<endl;
			fprintf(WF,"%12.10le %12.10le\n", y, wf1[j]);		
		}	
	}
	delete[] wf1;
}//end of basis
fclose(WF);
}

if(Compute_moment2){

int basis=9;
double xP=0.0;
double yP=0.0;
double zP=0.0;
double x=0.0;
double y=0.0;
double z=0.0;

FILE* fp=fopen("coordinates","r");
if(fp!=NULL){
        fscanf(fp,"%lf %lf %lf",&xP, &yP, &zP);
        cout<<"Coordinates of the impurity : "<<xP<<"    "<<yP<<"     "<<zP<<endl;
}else{
        cout<<"Error opening file. "<<endl;
        exit(-1);
}
fclose(fp);

for(int i=0;i<basis;i++){
        complex<double> *wf1 = new complex<double>[20*NumAtoms];
        read_wfs1(NumAtoms, data, i+1 , wf1);

        complex<double> rx2=0.0;
        complex<double> ry2=0.0;
        complex<double> rz2=0.0;

        for(int i=0;i<NumAtoms;i++){
                complex<double> tmp=0.0;
                x=lattice[3*i]-xP;
                y=lattice[3*i+1]-yP;
                z=lattice[3*i+2]-zP;
                for(int s=0;s<2;s++)
                        for(int o=0;o<10;o++)
                                tmp+=conj(wf1[i*20+s*10+o])*wf1[i*20+s*10+o];
                rx2+=sqrt(y*y+z*z)*tmp;
                ry2+=sqrt(x*x+z*z)*tmp;
                rz2+=sqrt(x*x+y*y)*tmp;
        }
        cout<<"r_x   :   "<<rx2<<endl;
        cout<<"r_y   :   "<<ry2<<endl;
        cout<<"r_z   :   "<<rz2<<endl;

        FILE* fp1=fopen("rx","a");
        FILE* fp2=fopen("ry","a");
        FILE* fp3=fopen("rz","a");
        if((fp1!=NULL)||(fp2!=NULL)||(fp3!=NULL)){
                fprintf(fp1,"%12.10le\n", rx2.real());
                fprintf(fp2,"%12.10le\n", ry2.real());
                fprintf(fp3,"%12.10le\n", rz2.real());
        }else{
                cout<<"Error opening file. "<<endl;
                exit(-1);
        }
        fclose(fp1);
        fclose(fp2);
        fclose(fp3);
        delete[] wf1;
}//end of basis
}




if(Compute_moment){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);

double xP=0.0;
double yP=0.0;
double zP=0.0;
	
FILE* fp=fopen("coordinates","r");
if(fp!=NULL){
	fscanf(fp,"%lf %lf %lf",&xP, &yP, &zP); 
	cout<<"Coordinates of the impurity : "<<xP<<"    "<<yP<<"     "<<zP<<endl;
}else{
	cout<<"Error opening file. "<<endl;
	exit(-1);
}
fclose(fp);

    complex<double> moment_x=0.0;
    complex<double> moment_y=0.0;
    complex<double> moment_z=0.0;	
	
    data_coulomb.moment(&moment_x, &moment_y, &moment_z, xP, yP, zP);
    cout<<"X   :   "<<moment_x<<endl;
    cout<<"Y   :   "<<moment_y<<endl;	
    cout<<"Z   :   "<<moment_z<<endl;

FILE* fp1=fopen("moment_x","a");
FILE* fp2=fopen("moment_y","a");
FILE* fp3=fopen("moment_z","a");
if((fp1!=NULL)||(fp2!=NULL)||(fp3!=NULL)){
        fprintf(fp1,"%12.10le\n", moment_x.real());
	fprintf(fp2,"%12.10le\n", moment_y.real());
	fprintf(fp3,"%12.10le\n", moment_z.real());
}else{
        cout<<"Error opening file. "<<endl;
        exit(-1);
}
fclose(fp1);
fclose(fp2);
fclose(fp3);
}


//Dipole moment between evec_1 and 2
if(Compute_dipole_moment){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);

    complex<double> moment_x=0.0;
    complex<double> moment_y=0.0;
    complex<double> moment_z=0.0;

    data_coulomb.dipole_moment(&moment_x, &moment_y, &moment_z);
    cout<<"Dipole X   :   "<<moment_x<<endl;
    cout<<"Dipole Y   :   "<<moment_y<<endl;
    cout<<"Dipole Z   :   "<<moment_z<<endl;

FILE* fp1=fopen("moment","w");
if(fp1!=NULL){
        fprintf(fp1,"%12.10le\n", moment_x.real());
        fprintf(fp1,"%12.10le\n", moment_y.real());
        fprintf(fp1,"%12.10le\n", moment_z.real());
}else{
        cout<<"Error opening file. "<<endl;
        exit(-1);
}
fclose(fp1);
}

if(Compute_r2){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);

double xP=0.0;
double yP=0.0;
double zP=0.0;

FILE* fp=fopen("coordinates","r");
if(fp!=NULL){
        fscanf(fp,"%lf %lf %lf",&xP, &yP, &zP);
        cout<<"Coordinates of the impurity : "<<xP<<"    "<<yP<<"     "<<zP<<endl;
}else{
        cout<<"Error opening file. "<<endl;
        exit(-1);
}
fclose(fp);

    complex<double> r2_x=0.0;
    complex<double> r2_y=0.0;
    complex<double> r2_z=0.0;

    data_coulomb.r2(&r2_x, &r2_y, &r2_z, xP, yP, zP);
    cout<<"r2_x   :   "<<r2_x<<endl;
    cout<<"r2_y   :   "<<r2_y<<endl;
    cout<<"r2_z   :   "<<r2_z<<endl;

FILE* fp1=fopen("r2","w");
if(fp1!=NULL){
        fprintf(fp1,"%12.10le\n", r2_x.real());
        fprintf(fp1,"%12.10le\n", r2_y.real());
        fprintf(fp1,"%12.10le\n", r2_z.real());
}else{
        cout<<"Error opening file. "<<endl;
        exit(-1);
}
fclose(fp1);
}

/*if(Compute_LS){

int dim=20*NumAtoms;
nml_dcvector *wf=nml_dcv_new(2*dim);
int m=0;
for(int iat=0; iat<NumAtoms; iat++){
	for(int iorb=0; iorb<10; iorb++){
		wf[m].r=wf_e[m].real();
		wf[m].i=wf_e[m].imag();
		wf[m+10].r=wf_e[m+10].real();
                wf[m+10].i=wf_e[m+10].imag();
	
        	wf[dim+m+10].r=wf_e[m].real();
                wf[dim+m+10].i=-wf_e[m].imag();
                wf[dim+m].r=-wf_e[m+10].real();
                wf[dim+m].i=wf_e[m+10].imag();
                m++;
         }//end of orbital
         m+=10;
}//end of iat

double dot_product1=0.0;
double dot_product2=0.0;
double dp3_r=0.0;
double dp3_i=0.0;

m=0;
for(int iat=0; iat<NumAtoms; iat++){
	for(int s=0;s<2;s++)
        for(int iorb=0; iorb<10; iorb++){
		dot_product1+=wf[m].r*wf[m].r+wf[m].i*wf[m].i;
		dot_product2+=wf[dim+m].r*wf[dim+m].r+wf[dim+m].i*wf[dim+m].i;
		dp3_r+=wf[m].r*wf[dim+m].r+wf[m].i*wf[dim+m].i;
		dp3_i+=wf[m].r*wf[dim+m].i-wf[m].i*wf[dim+m].r;
		m++;
	}                                                                                                                                           
}
cout<<"dot product: wf1*wf1  "<<dot_product1<<endl;
cout<<"dot product: wf2*wf2  "<<dot_product2<<endl;
cout<<"dot product: wf1*wf2  "<<dp3_r<<"   "<<dp3_i<<endl;

nml_dcvector *Sx=nml_dcv_new(4);
nml_dcvector *Sy=nml_dcv_new(4);
nml_dcvector *Sz=nml_dcv_new(4);

for(int i=0;i<2;i++){
	for(int j=0;j<2;j++){

		for(int at=0; at<NumAtoms; at++){
      			for(int iorb=0; iorb<10; iorb++){
          			nml_dcscalar c1_up, c1_down, c2_up, c2_down;
          			c1_up=wf[i*dim+at*20+iorb];
          			c2_up=wf[j*dim+at*20+iorb];
          			c1_down=wf[i*dim+at*20+iorb+10];
          			c2_down=wf[j*dim+at*20+iorb+10];


	  			Sx[i*2+j].r +=c1_up.r*c2_down.r+c1_up.i*c2_down.i
                   			+c1_down.r*c2_up.r+c1_down.i*c2_up.i;
          			Sx[i*2+j].i +=c1_up.r*c2_down.i-c1_up.i*c2_down.r
                   			+c1_down.r*c2_up.i-c1_down.i*c2_up.r;
          			Sy[i*2+j].r +=c1_up.r*c2_down.i-c1_up.i*c2_down.r
                   			+c1_down.i*c2_up.r-c1_down.r*c2_up.i;
          			Sy[i*2+j].i +=c1_down.r*c2_up.r+c1_down.i*c2_up.i
                   			-c1_up.r*c2_down.r-c1_up.i*c2_down.i;
          			Sz[i*2+j].r +=c1_up.r*c2_up.r+c1_up.i*c2_up.i
                   			-c1_down.r*c2_down.r-c1_down.i*c2_down.i;
          			Sz[i*2+j].i +=c1_up.r*c2_up.i-c1_up.i*c2_up.r
                    			-c1_down.r*c2_down.i+c1_down.i*c2_down.r;
        	}
	}
}}

FILE* fp_x=fopen("Sx","w");
FILE* fp_y=fopen("Sy","w");
FILE* fp_z=fopen("Sz","w");

for(int i=0;i<2;i++)
        for(int j=0;j<2;j++){
		Sx[i*2+j].r *=0.5; Sx[i*2+j].i *=0.5;
		Sy[i*2+j].r *=0.5; Sy[i*2+j].i *=0.5;
		Sz[i*2+j].r *=0.5; Sz[i*2+j].i *=0.5;
	
		fprintf(fp_x, "%12.10le %12.10le\n", Sx[i*2+j].r, Sx[i*2+j].i);
		fprintf(fp_y, "%12.10le %12.10le\n", Sy[i*2+j].r, Sy[i*2+j].i);
		fprintf(fp_z, "%12.10le %12.10le\n", Sz[i*2+j].r, Sz[i*2+j].i);
	}
fclose(fp_x);
fclose(fp_y);
fclose(fp_z);
cout<<"Done with spin calculations"<<endl;

nml_dcvector *X_WF=nml_dcv_new(NumAtoms*20);
nml_dcvector *Y_WF=nml_dcv_new(NumAtoms*20);
nml_dcvector *Z_WF=nml_dcv_new(NumAtoms*20);
nml_dcvector *H_X_WF=nml_dcv_new(NumAtoms*20);
nml_dcvector *H_Y_WF=nml_dcv_new(NumAtoms*20);
nml_dcvector *H_Z_WF=nml_dcv_new(NumAtoms*20);

// read H_X_WF
// read H_Y_WF
// read H_Z_WF

FILE *hx=fopen("H_X_WF1","w");
if(hx!=NULL)
	fread(H_X_WF, 2*sizeof(double), 20*Natoms, hx);
else
	cout<<"Error reading file"<<endl;
fclose(hx);

FILE *hy=fopen("H_Y_WF1","w");
if(hy!=NULL)
        fread(H_Y_WF, 2*sizeof(double), 20*Natoms, hy);
else
        cout<<"Error reading file"<<endl;
fclose(hy);

FILE *hz=fopen("H_Z_WF1","w");
if(hz!=NULL)
        fread(H_Z_WF, 2*sizeof(double), 20*Natoms, hz);
else
        cout<<"Error reading file"<<endl;
fclose(hz);

double X0=15.20666;
double Y0=15.20666;
double Z0=15.20666;
double X=0;
double Y=0;
double Z=0;

for(int iat=0; iat<NumAtoms; iat++){
      X=lattice[3*iat]-X0;
      Y=lattice[3*iat+1]-Y0;
      Z=lattice[3*iat+2]-Z0;
      for(int spin=0; spin<2; spin++){
        for(int iorb=0; iorb<10; iorb++){
		cf = cf_e(n1,iat,spin,iorb);
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
}
*/

if(Compute_SCF){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *HF = new double[NumAtoms];
    double *XC = new double[NumAtoms];
    double *HF_total = new double[NumAtoms];
    double *wf = new double[NumAtoms];

    double inv_eps = 1.0 / 11.9;
    double esqr = 1.4399766; // eV-nm
    double coulomb_constant=-1*inv_eps*esqr;

    int n=0;
    int index=-1; 	
    double max=-1;

for(int x=0;x<NumAtoms;x++){
        for(int s=0;s<2;s++)
        for(int o=0;o<10;o++){
		wf[x]+=wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag();
                n++;
        }
        /*if(wf[x]>max){
                max=wf[x];
		index=x;
	}*/	
}	
/*
   cout<<"Index of the impurity : "<<index<<endl;
   cout<<"Coordinates :  ("<<lattice[index*3]<<", "<<lattice[index*3+1]<<", "<<lattice[index*3+2]<<")"<<endl;		
*/  
	
    data_coulomb.getSCF(wf, HF, XC);

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(HF, HF_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                HF_total[i]=HF[i];
	}
#endif
    if (rank==0) {
        cout<<"Finished calculating SCF potential"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("HF", "w");
        for(int i=0;i<NumAtoms;i++){
                fprintf(fp1,"%12.10le\n", HF_total[i]);
        }
        fclose(fp1);
    }
}

if(GenerateSpin_wfs){

read_write_wf_10(data, NumElectrons, eList, NumAtoms);

}

if(Compute_SCF_new){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *HF = new double[NumAtoms];
    double *XC_r = new double[NumAtoms];
    double *XC_i = new double[NumAtoms];
    double *HF_total = new double[NumAtoms];
    double *XC_total_r = new double[NumAtoms];
    double *XC_total_i = new double[NumAtoms];
    double *XC_total= new double[NumAtoms];	
    double *wf_coul = new double[NumAtoms];
    double *wf_xc_r = new double [NumAtoms];
    double *wf_xc_i = new double [NumAtoms];		

    int n=0;
    int index=-1;
    double max=-1;
    int skip=NumAtoms*20;

for(int x=0;x<NumAtoms;x++){
        for(int s=0;s<2;s++)
        for(int o=0;o<10;o++){
                wf_coul[x]+=wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag();
		wf_xc_r[x]+=wf_e[n].real()*wf_e[skip+n].real()+wf_e[n].imag()*wf_e[skip+n].imag();
		wf_xc_i[x]+=wf_e[n].real()*wf_e[skip+n].imag()-wf_e[n].imag()*wf_e[skip+n].real();
                n++;
        }
}

data_coulomb.getSCF_new(wf_coul, wf_xc_r, wf_xc_i, HF, XC_r, XC_i);
double V=0.0;

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(HF, HF_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(XC_r, XC_total_r, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(XC_i, XC_total_i, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                HF_total[i]=HF[i];
		XC_total_r[i]=XC_r[i];
		XC_total_i[i]=XC_i[i];		
        }
#endif
    if (rank==0) {
        cout<<"Finished calculating SCF potential"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("HF", "w");
	FILE* fp2=fopen("XC_r", "w");
	FILE* fp3=fopen("XC_i", "w");
	FILE* fp4=fopen("XC_total", "w");
	FILE* fp5=fopen("V_total", "w");
        for(int i=0;i<NumAtoms;i++){
		fprintf(fp1,"%12.10le\n", HF_total[i]);
                fprintf(fp2,"%12.10le\n", XC_total_r[i]);
		fprintf(fp3,"%12.10le\n", XC_total_i[i]);
		V=sqrt(XC_total_r[i]*XC_total_r[i]+XC_total_i[i]*XC_total_i[i]);
		fprintf(fp4,"%12.10le\n", V);
		V=HF_total[i]-V;
		fprintf(fp5,"%12.10le\n", V);
        }
        fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);
    }
}

if(Compute_SCF_Elect){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *HF = new double[NumAtoms];
    double *XC = new double[NumAtoms];
    double *HF_total = new double[NumAtoms];
    double *XC_total = new double[NumAtoms];
    double *wf = new double[NumAtoms];
//Should again normaliz the evecs

if(!charge){
cout<<"Consider spin filling"<<endl;
if (num_Electrons==1){
	int n=0;
	cout<<"Calculating 1e occupation density"<<endl;
	for(int x=0;x<NumAtoms;x++){
        	for(int s=0;s<2;s++)
        	for(int o=0;o<10;o++){
			wf[x]+=(wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag());	
        		n++;
		}
	}	
}else if (num_Electrons==2){
	int n=0;
        cout<<"Calculating 2e occupation density"<<endl;
        for(int x=0;x<NumAtoms;x++){
                for(int s=0;s<2;s++)
                for(int o=0;o<10;o++){
                        wf[x]+=2.0*(wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag());    
                        n++;
                }
        }
} 
else if (num_Electrons==3){
        int n=0;
	int next_index=20*NumAtoms;
        cout<<"Calculating 3e occupation density"<<endl;
        for(int x=0;x<NumAtoms;x++){
                for(int s=0;s<2;s++)
                for(int o=0;o<10;o++){
                        wf[x]+=2.0*(wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag());
			wf[x]+=(wf_e[next_index+n].real()*wf_e[next_index+n].real()+wf_e[next_index+n].imag()*wf_e[next_index+n].imag());
                        n++;
                }
        }
}else if (num_Electrons==4){
        int n=0;
        int next_index=20*NumAtoms;
        cout<<"Calculating 4e occupation density"<<endl;
        for(int x=0;x<NumAtoms;x++){
                for(int s=0;s<2;s++)
                for(int o=0;o<10;o++){
                        wf[x]+=2.0*(wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag());
                        wf[x]+=2.0*(wf_e[next_index+n].real()*wf_e[next_index+n].real()+wf_e[next_index+n].imag()*wf_e[next_index+n].imag());
                        n++;
                }
        }
}else if (num_Electrons==5){
        int n=0;
        int next_index=20*NumAtoms;
        cout<<"Calculating 5e occupation density"<<endl;
        for(int x=0;x<NumAtoms;x++){
                for(int s=0;s<2;s++)
                for(int o=0;o<10;o++){
                        wf[x]+=2.0*(wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag());
                        wf[x]+=2.0*(wf_e[next_index+n].real()*wf_e[next_index+n].real()+wf_e[next_index+n].imag()*wf_e[next_index+n].imag());
			wf[x]+=(wf_e[2*next_index+n].real()*wf_e[2*next_index+n].real()+wf_e[2*next_index+n].imag()*wf_e[2*next_index+n].imag());
			n++;
                }
        }
}else if (num_Electrons==6){
        int n=0;
        int next_index=20*NumAtoms;
        cout<<"Calculating 6e occupation density"<<endl;
        for(int x=0;x<NumAtoms;x++){
                for(int s=0;s<2;s++)
                for(int o=0;o<10;o++){
                        wf[x]+=2.0*(wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag());
                        wf[x]+=2.0*(wf_e[next_index+n].real()*wf_e[next_index+n].real()+wf_e[next_index+n].imag()*wf_e[next_index+n].imag());
                        wf[x]+=2.0*(wf_e[2*next_index+n].real()*wf_e[2*next_index+n].real()+wf_e[2*next_index+n].imag()*wf_e[2*next_index+n].imag());
                        n++;
                }       
        }
}else
 cout<<"Cannot do this case: electron number"<<num_Electrons<<endl;	
} else{
   cout<<"Just use charge degrees of freedom"<<endl;	
	int n=0;
	for(int i=0;i<num_Electrons;i++){
        for(int x=0;x<NumAtoms;x++){
                for(int s=0;s<2;s++)
                for(int o=0;o<10;o++){
                        wf[x]+=wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag();
                        n++;
                }
        }
	}
}
 	
    data_coulomb.getSCF(wf, HF, XC);

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(HF, HF_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //MPI_Reduce(XC, XC_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                HF_total[i]=HF[i];
	}
#endif
    if (rank==0) {
        cout<<"Finished calculating SCF potential"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("HF", "w");
        for(int i=0;i<NumAtoms;i++){
                fprintf(fp1,"%12.10le\n", HF_total[i]);
        }
        fclose(fp1);
  }
}

if(Compute_SCF_D2){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *HF = new double[NumAtoms];
    double *XC = new double[NumAtoms];
    double *HF_total = new double[NumAtoms];
    double *XC_total = new double[NumAtoms];
    double *wf = new double[NumAtoms];

    int n=0;
for(int x=0;x<NumAtoms;x++){
        for(int s=0;s<2;s++)
        for(int o=0;o<10;o++){
		wf[x]+=2.0*(wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag());
                n++;
        }
    }	
  	
    data_coulomb.getSCF(wf, HF, XC);

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(HF, HF_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //MPI_Reduce(XC, XC_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                HF_total[i]=HF[i];
	}
#endif
    if (rank==0) {
        cout<<"Finished calculating SCF potential"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("HF", "w");
        for(int i=0;i<NumAtoms;i++){
                fprintf(fp1,"%12.10le\n", HF_total[i]);
        }
        fclose(fp1);
  }
}


if(Compute_SCF_D2_single){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *HF = new double[NumAtoms];
    double *XC = new double[NumAtoms];
    double *HF_total = new double[NumAtoms];
    double *wf = new double[NumAtoms];

    double inv_eps = 1.0 / 11.9;
    double esqr = 1.4399766; // eV-nm
    double coulomb_constant=-1*inv_eps*esqr;

    int n=0;
    int index=-1; 	
    double max=-1;

for(int i=0;i<NumElectrons;i++){
	cout<<"evec number : "<<i<<endl;	
	int cur_index=i*NumAtoms*20;
	n=0;
for(int x=0;x<NumAtoms;x++){
        for(int s=0;s<2;s++)
        for(int o=0;o<10;o++){
                wf[x]+=wf_e[cur_index+n].real()*wf_e[cur_index+n].real()+wf_e[cur_index+n].imag()*wf_e[cur_index+n].imag();
                n++;
        }
    }	
}//end of i
    data_coulomb.getSCF(wf, HF, XC);

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(HF, HF_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                HF_total[i]=HF[i];
	}
#endif
    if (rank==0) {
        cout<<"Finished calculating SCF potential"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("HF", "w");
        for(int i=0;i<NumAtoms;i++){
                fprintf(fp1,"%12.10le\n", HF_total[i]);
        }
        fclose(fp1);
  }
}


if(Compute_SCF_D3_single){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *HF = new double[NumAtoms];
    double *XC = new double[NumAtoms];
    double *HF_total = new double[NumAtoms];
    double *wf = new double[NumAtoms];

    double inv_eps = 1.0 / 11.9;
    double esqr = 1.4399766; // eV-nm
    double coulomb_constant=-1*inv_eps*esqr;

    int n=0;
    int index=-1; 	
    double max=-1;

for(int i=0;i<NumElectrons;i++){
	cout<<"evec number : "<<i<<endl;	
	int cur_index=i*NumAtoms*20;
	n=0;
for(int x=0;x<NumAtoms;x++){
        for(int s=0;s<2;s++)
        for(int o=0;o<10;o++){
                wf[x]+=wf_e[cur_index+n].real()*wf_e[cur_index+n].real()+wf_e[cur_index+n].imag()*wf_e[cur_index+n].imag();
                n++;
        }
    }	
}//end of i
    data_coulomb.getSCF(wf, HF, XC);

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(HF, HF_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                HF_total[i]=HF[i];
	}
#endif
    if (rank==0) {
        cout<<"Finished calculating SCF potential"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("HF", "w");
        for(int i=0;i<NumAtoms;i++){
                fprintf(fp1,"%12.10le\n", HF_total[i]);
        }
        fclose(fp1);
  }
}


if(Compute_SCF_D3){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *HF = new double[NumAtoms];
    double *XC = new double[NumAtoms];
    double *HF_total = new double[NumAtoms];
    double *wf = new double[NumAtoms];

    double inv_eps = 1.0 / 11.9;
    double esqr = 1.4399766; // eV-nm
    double coulomb_constant=-1*inv_eps*esqr;

    int n=0;
    int index=-1; 	

for(int i=0;i<NumElectrons;i++){
	cout<<"evec number : "<<i<<endl;	
	int cur_index=i*NumAtoms*20;
	n=0;
for(int x=0;x<NumAtoms;x++){
        for(int s=0;s<2;s++)
        for(int o=0;o<10;o++){
		if(i==0)
                	wf[x]+=2.0*(wf_e[cur_index+n].real()*wf_e[cur_index+n].real()+wf_e[cur_index+n].imag()*wf_e[cur_index+n].imag());
		else
			wf[x]+=wf_e[cur_index+n].real()*wf_e[cur_index+n].real()+wf_e[cur_index+n].imag()*wf_e[cur_index+n].imag();
                n++;
        }
    }	
}//end of i

    data_coulomb.getSCF(wf, HF, XC);

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(HF, HF_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                HF_total[i]=HF[i];
	}
#endif
    if (rank==0) {
        cout<<"Finished calculating SCF potential"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("HF", "w");
        for(int i=0;i<NumAtoms;i++){
                fprintf(fp1,"%12.10le\n", HF_total[i]);
        }
        fclose(fp1);
  }
}

if(Compute_SiGe){
    //Coulomb data_coulomb =
     //     Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    //data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
	double a0=0.543095;
    	double max_z=0.0;
	double z_sige=0.0;
	for(int i=0;i<NumAtoms;i++){
		if(atomid[i]==2){
			z_sige=lattice[3*i+2];
			if((z_sige>max_z)&&(z_sige<20*a0)){
				max_z=z_sige;
				//break;
			}	
		}
	}
	
	cout<<"Maximum z coordinate of Ge: "<<max_z<<endl;

     	double* wf1 = new double[NumAtoms];
     	read_wfs_prob(NumAtoms, data, 1 , wf1);
     	//double a0=0.543095;
     	double Lz=max_z;

     double sum=0.0;
     double norm=0.0;
     double x_coor=0.0;
     double y_coor=0.0;
     double z_coor=0.0;
     double ge_density=0.0;

     for(int i=0;i<NumAtoms;i++){
        x_coor=lattice[3*i];
        y_coor=lattice[3*i+1];
        z_coor=lattice[3*i+2];
        norm+=wf1[i];
        if(z_coor<=Lz){
                sum+=wf1[i];
        }
	if(atomid[i]==2)
		ge_density+=wf1[i];	
     }
  cout<<"Overlap : "<<sum<<endl;
  cout<<"Norm: "<<norm<<endl;
  cout<<"SiGe barrier Leakage: "<<sum/norm<<endl;
  cout<<"Total Ge density:  "<<ge_density/norm<<endl;
}

if(Compute_CTAP){
    //Coulomb data_coulomb =
     //     Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    //data_coulomb.get_data(wf_e, atomid, neighbor, lattice);

     double* wf1 = new double[NumAtoms];
     read_wfs_prob(NumAtoms, data, 2 , wf1);
     double a0=0.543095;
     double L=10;	
     double Lx1=(56-L)*a0;
     double Lx2=(56+L)*a0;
     double Ly1=(28-L)*a0;
     double Ly2=(28+L)*a0;
     double Lz1=(28-L)*a0;
     double Lz2=(28+L)*a0;
	
     double sum=0.0;
     double norm=0.0;	
     double x_coor=0.0;	
     double y_coor=0.0;
     double z_coor=0.0;

     for(int i=0;i<NumAtoms;i++){	 
	x_coor=lattice[3*i];
	y_coor=lattice[3*i+1];
	z_coor=lattice[3*i+2];
	norm+=wf1[i];
	if((x_coor<=Lx2) && (x_coor>=Lx1) && (y_coor<=Ly2) && (y_coor>=Ly1) && (z_coor<=Lz2) && (z_coor>=Lz1)){
		sum+=wf1[i];	
	}
     }	
  cout<<"Sum : "<<sum<<endl;
  cout<<"Norm: "<<norm<<endl; 	
  cout<<"Middle donor Leakage: "<<sum/norm<<endl;	
}

if(Compute_SCF_CI){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *HF = new double[NumAtoms];
    double *HF_total=new double[NumAtoms];
    double *n_r = new double[NumAtoms];

    double inv_eps = 1.0 / 11.9;
    double esqr = 1.4399766; // eV-nm
    double coulomb_constant=-1*inv_eps*esqr;

 //Read CI ground state
  int num_Basis=0; 	
  double *coef;
  FILE* fp_gs=fopen("CI_GS","r");
  if(!fp_gs){
	die("CI_GS file not found\n.");
	exit(-1);
   }else{
	fscanf(fp_gs,"%d\n",&num_Basis);
	cout<<"Length of the eigenvector: "<<num_Basis<<endl;
        coef = new double[num_Basis];
	for(int i=0;i<num_Basis;i++)
		fscanf(fp_gs, "%lf\n", &coef[i]);	
   }	
	cout<<"CI ground state"<<endl;
  	for(int i=0;i<num_Basis;i++){
		cout<<coef[i]<<endl;
	}

FILE* fp_states=fopen("States","r");
int *first_state=new int[num_Basis];
int *second_state=new int[num_Basis];
if(fp_states){
cout<<"Reading States: "<<endl;
for(int i=0;i<num_Basis;i++){
	fscanf(fp_states, "%d %d\n", &first_state[i], &second_state[i]);
	cout<<first_state[i]<<"   "<<second_state[i]<<endl;
}
}else{
	cout<<"File not found : States "<<endl;
	exit(-1);
}

for(int i=0;i<num_Basis;i++){
  if(coef[i]>0.01){
	complex<double> *wf1 = new complex<double>[20*NumAtoms];
    	complex<double> *wf2 = new complex<double>[20*NumAtoms];
	read_wfs(NumAtoms, data, first_state[i], second_state[i], wf1, wf2);

        cout<<"Reading wf1 "<<endl;
        for(int k=0;k<5;k++){
	for(int j=0;j<20;j++){	
        if(k<4)

                cout<<k<<" :  "<<wf1[k*20+j]<<endl;
        else
                cout<<NumAtoms-1<<"   :  "<<wf1[20*(NumAtoms-1)+j]<<endl;
   	}
	}
	cout<<"Reading wf2 "<<endl;
   	for(int k=0;k<5;k++){          
	for(int j=0;j<20;j++){
        if(k<4)                   
                cout<<k<<" :  "<<wf2[k*20+j]<<endl;
        else                      
                cout<<NumAtoms-1<<"   :  "<<wf2[20*(NumAtoms-1)+j]<<endl;
   	}
	}

	for(int j=0;j<NumAtoms;j++){
		double tmp1=0.0;
		double tmp2=0.0;
		for(int l=0;l<20;l++){
			tmp1+=wf1[j*20+l].real()*wf1[j*20+l].real()+wf1[j*20+l].imag()*wf1[j*20+l].imag();
			tmp2+=wf2[j*20+l].real()*wf2[j*20+l].real()+wf2[j*20+l].imag()*wf2[j*20+l].imag();		
		}	
		n_r[j]+=coef[i]*0.5*(tmp1+tmp2);
	}//end of j
	delete[] wf1;
	delete[] wf2;		
  }else{
	cout<<"Skipping coefficient as it is too small:  "<<i<<endl;	
  }
}//end of i

    data_coulomb.getSCF_CI(n_r, HF);

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(HF, HF_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                HF_total[i]=HF[i];
        }
#endif
    if (rank==0) {
        cout<<"Finished calculating SCF potential"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("HF", "w");
        for(int i=0;i<NumAtoms;i++){
                fprintf(fp1,"%12.10le\n", HF_total[i]);
        }
        fclose(fp1);
}
if(rank==0){
	if(save_electron_density){
		FILE* nr_density=fopen("nr_GS","w");
		if(nr_density){
        		for(int m=0;m<NumAtoms;m++)
                		fprintf(nr_density,"%12.10le\n",n_r[m]);
		}else{
        		cout<<"Could not open file for printing: nr_GS"<<endl;
        		exit(-1);
		}
		fclose(nr_density);
	}//end of save_electron_density
}//end of rank

}

if(save_electron_density){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *n_r = new double[NumAtoms];

 //Read CI ground state
  int num_Basis=0;
  double *coef;
  FILE* fp_gs=fopen("CI_GS","r");
  if(!fp_gs){
        die("CI_GS file not found\n.");
        exit(-1);
   }else{
        fscanf(fp_gs,"%d\n",&num_Basis);
        cout<<"Length of the eigenvector: "<<num_Basis<<endl;
        coef = new double[num_Basis];
        for(int i=0;i<num_Basis;i++)
                fscanf(fp_gs, "%lf\n", &coef[i]);
   }
        cout<<"CI ground state"<<endl;
        for(int i=0;i<num_Basis;i++){
                cout<<coef[i]<<endl;
        }

FILE* fp_states=fopen("States","r");
int *first_state=new int[num_Basis];
int *second_state=new int[num_Basis];
if(fp_states){
cout<<"Reading States: "<<endl;
for(int i=0;i<num_Basis;i++){
    fscanf(fp_states, "%d %d\n", &first_state[i], &second_state[i]);
        cout<<first_state[i]<<"   "<<second_state[i]<<endl;
}
}else{
        cout<<"File not found : States "<<endl;
        exit(-1);
}

for(int i=0;i<num_Basis;i++){
  if(coef[i]>0.01){
        complex<double> *wf1 = new complex<double>[20*NumAtoms];
        complex<double> *wf2 = new complex<double>[20*NumAtoms];
        read_wfs(NumAtoms, data, first_state[i], second_state[i], wf1, wf2);

        cout<<"Reading wf1 "<<endl;
        for(int k=0;k<5;k++){
        for(int j=0;j<20;j++){
        if(k<4)

                cout<<k<<" :  "<<wf1[k*20+j]<<endl;
        else
                cout<<NumAtoms-1<<"   :  "<<wf1[20*(NumAtoms-1)+j]<<endl;
        }
        }
        cout<<"Reading wf2 "<<endl;
        for(int k=0;k<5;k++){
        for(int j=0;j<20;j++){
        if(k<4)
                cout<<k<<" :  "<<wf2[k*20+j]<<endl;
        else
                cout<<NumAtoms-1<<"   :  "<<wf2[20*(NumAtoms-1)+j]<<endl;
        }
        }

        for(int j=0;j<NumAtoms;j++){
                double tmp1=0.0;
                double tmp2=0.0;
                for(int l=0;l<20;l++){
                        tmp1+=wf1[j*20+l].real()*wf1[j*20+l].real()+wf1[j*20+l].imag()*wf1[j*20+l].imag();
                        tmp2+=wf2[j*20+l].real()*wf2[j*20+l].real()+wf2[j*20+l].imag()*wf2[j*20+l].imag();
                }
                n_r[j]+=coef[i]*0.5*(tmp1+tmp2);
        }//end of j
       delete[] wf1;
        delete[] wf2;
  }else{
        cout<<"Skipping coefficient as it is too small:  "<<i<<endl;
  }
}//end of i

//save electron desnity of each many-body state in file
        if(save_electron_density){
                FILE* nr_density=fopen("nr_GS","w");
                if(nr_density){
                        for(int m=0;m<NumAtoms;m++)
                                fprintf(nr_density,"%12.10le\n",n_r[m]);
                }else{
                        cout<<"Could not open file for printing: nr_GS"<<endl;
                        exit(-1);
                }
                fclose(nr_density);
        }//end of save_electron_density

}





if(Compute_SCF_avg){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *HF = new double[NumAtoms];
    double *XC = new double[NumAtoms];
    double *HF_total = new double[NumAtoms];
    double *XC_total = new double[NumAtoms];
    double *wf = new double[NumAtoms];

    double inv_eps = 1.0 / 11.9;
    double esqr = 1.4399766; // eV-nm
    double coulomb_constant=-1*inv_eps*esqr;

    double C=-1;	
    double correlation_constant=coulomb_constant*C;	

    int n=0;
    int index=-1; 	
    double max=-1;
	for(int x=0;x<NumAtoms;x++){
		double tmp1=0.0;
		double tmp2=0.0;
        	for(int s=0;s<2;s++){
        	for(int o=0;o<10;o++){
                	tmp1+=wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag();
			tmp2+=wf_e[NumAtoms*20+n].real()*wf_e[NumAtoms*20+n].real()+wf_e[NumAtoms*20+n].imag()*wf_e[NumAtoms*20+n].imag();
                	n++;
        	}
		}
		wf[x]=0.5*(tmp1+tmp2);
        	if(wf[x]>max){
                	max=wf[x];
			index=x;
		}	
    	}
    cout<<"Index of the impurity : "<<index<<endl;
    cout<<"Coordinates :  ("<<lattice[index*3]<<", "<<lattice[index*3+1]<<", "<<lattice[index*3+2]<<")"<<endl;		
    data_coulomb.getSCF(wf, HF, XC);

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(HF, HF_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //MPI_Reduce(XC, XC_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                HF_total[i]=HF[i];
		//XC_total[i]=XC[i];
	}
#endif
    if (rank==0) {
        cout<<"Finished calculating SCF potential"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("HF", "w");
        //FILE* fp2=fopen("XC", "w");
        //FILE* fp3=fopen("V_total", "w");
        for(int i=0;i<NumAtoms;i++){
		//HF_total[i]=coulomb_constant*HF_total[i];
		//HF_total[i]=HF_total[i];
		//XC_total[i]=correlation_constant*XC_total[i];
                fprintf(fp1,"%12.10le\n", HF_total[i]);
                //fprintf(fp2,"%12.10le\n", XC_total[i]);
                //fprintf(fp3,"%12.10le\n", HF_total[i]+XC_total[i]);
        }
        fclose(fp1);
        //fclose(fp2);
        //fclose(fp3);
}
}

if(Compute_SCF_image){
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *HF = new double[NumAtoms];
    double *XC = new double[NumAtoms];
    double *HF_total = new double[NumAtoms];
    double *XC_total = new double[NumAtoms];
    double *wf = new double[NumAtoms];

    double inv_eps = 1.0 / 11.9;
    double esqr = 1.4399766; // eV-nm
    double coulomb_constant=-1*inv_eps*esqr;

    double C=-1;
    double correlation_constant=coulomb_constant*C;
/*
    double r0 = 4.28; //a.u.
    double q = 1.10;
    double bohr_a0=0.0529177208;        //nm
    double r0_nm=r0*bohr_a0; //nm
    cout<<"Cut off distnace r0 = "<<r0_nm<<endl;
    double q_r0=q*r0;
    double R=0.27;
    double dielectric_f=11.9;
    double R_atomic=0.0;
    double argument=0.0;

        if(R>r0_nm)
                dielectric_f=11.9;
        else{
                R_atomic=R/bohr_a0;
                argument=q_r0-q*R_atomic;
                dielectric_f=11.9*q_r0/((q*R_atomic) + sinh(argument));
        }
        cout<<"Dielectric value at R = "<<R<<"   is  :  "<<dielectric_f<<endl;
*/

    int n=0;
    int index=-1;
    double max=-1;
        for(int x=0;x<NumAtoms;x++){
        for(int s=0;s<2;s++)
        for(int o=0;o<10;o++){
                wf[x]+=wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag();
                n++;
        }
        if(wf[x]>max){
                max=wf[x];
                index=x;
        }
    }
    cout<<"Index of the impurity : "<<index<<endl;
    cout<<"Coordinates :  ("<<lattice[index*3]<<", "<<lattice[index*3+1]<<", "<<lattice[index*3+2]<<")"<<endl;
    data_coulomb.getSCF_image(wf, HF, XC);

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(HF, HF_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //MPI_Reduce(XC, XC_total, NumAtoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for(int i=0;i<NumAtoms;i++){
                HF_total[i]=HF[i];
                //XC_total[i]=XC[i];
        }
#endif
    if (rank==0) {
        cout<<"Finished calculating SCF potential"<<endl;
        cout<<"Writing to file"<<endl;
        FILE* fp1=fopen("HF", "w");
        //FILE* fp2=fopen("XC", "w");
        //FILE* fp3=fopen("V_total", "w");
        for(int i=0;i<NumAtoms;i++){
                //HF_total[i]=coulomb_constant*HF_total[i];
                //HF_total[i]=HF_total[i];
                //XC_total[i]=correlation_constant*XC_total[i];
                fprintf(fp1,"%12.10le\n", HF_total[i]);
                //fprintf(fp2,"%12.10le\n", XC_total[i]);
                //fprintf(fp3,"%12.10le\n", HF_total[i]+XC_total[i]);
        }
        fclose(fp1);
        //fclose(fp2);
        //fclose(fp3);
}
}

if(Compute_HF_HF) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *A_aniso= new double[9*counter];
    double *A_iso = new double[counter]; 	
    int *P_index=new int[counter]; 	
    cout<<"Number of P atoms : "<<counter<<endl;
    cout<<"Atomic coordinates : "<<endl;	
    int count=0;	
	
    for(int i=0; i<NumAtoms;i++){
		if(atomid[i]==1){
			P_index[count]=i;
			count++;
			cout<<i<<endl;
			cout<<"("<<lattice[i*3]<<", "<<lattice[i*3+1]<<", "<<lattice[i*3+2]<<")"<<endl;
		}
	}

    data_coulomb.HF(counter, P_index, A_aniso, A_iso);
	
FILE* fp1=fopen("A_ij", "w");
for(int j=0;j<counter;j++){
for(int i=0;i<9;i++)
	fprintf(fp1,"%12.10le\n", A_aniso[9*j+i]);

	fprintf(fp1, "\n");
}	
fclose(fp1);
cout<<"1P density is 1.7741717755e-02 and it corresponds to 117 MHZ"<<endl;

FILE* fp2=fopen("A_0","w");
for(int i=0;i<counter;i++)
	fprintf(fp2,"%12.10le\n", (A_iso[i]/1.7741717755e-02)*117);
fclose(fp2);

FILE* fp3=fopen("P_coordinates","w");
for(int i=0;i<counter;i++){
fprintf(fp3,"%12.10le\n", lattice[3*P_index[i]]);
fprintf(fp3,"%12.10le\n", lattice[3*P_index[i]+1]);
fprintf(fp3,"%12.10le\n\n", lattice[3*P_index[i]+2]);
}
fclose(fp3);

}


if(ComputeAHF_wf) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
    double *A_aniso= new double[9*NumAtoms];
    double *wf=new double[NumAtoms];

    if(rank==0)	
    	cout<<"Printing wf "<<endl;	
    int n=0;
	double max=-1;
    for(int x=0;x<NumAtoms;x++){
        for(int s=0;s<2;s++)
        for(int o=0;o<10;o++){
        	wf[x]+=wf_e[n].real()*wf_e[n].real()+wf_e[n].imag()*wf_e[n].imag();
		n++;
	}
	if(wf[x]>max)
		max=wf[x];	
    }		
	
//end of file wf generation

/*	
    	//generate Gaussian wf
	double x_coor=0.0;
	double y_coor=0.0;
	double z_coor=0.0;
	double R=0.0;
	double a_bohr=3;	//nm	
	double xP=15.20666;
	double yP=15.20666;
	double zP=15.20666;

	//need to check Gaussian formula
	for(int i=0;i<NumAtoms;i++){ 	
		x_coor=lattice[i*3]-xP;
		y_coor=lattice[i*3+1]-yP;
		z_coor=lattice[i*3+2]-zP;
		R=sqrt(x_coor*x_coor+y_coor*y_coor+z_coor*z_coor);
		
		//s-type
		//wf[i]=exp(-R/a_bohr);
		
		//px-type
                //wf[i]=x_coor*exp(-R/a_bohr);
		
		//py-type
                //wf[i]=y_coor*exp(-R/a_bohr);
		
		//pz-type
                //wf[i]=z_coor*exp(-R/a_bohr);
	}
	//Normalize wf
	double sum=0.0;
	double norm=0.0;
	for(int i=0;i<NumAtoms;i++)
		sum+=wf[i]*wf[i];
	norm=sqrt(sum);	
	cout<<"Normalization constant: "<<norm<<endl;
	for(int i=0;i<NumAtoms;i++)
		wf[i]=wf[i]*wf[i]/sum;	

	if(rank==0){	
		//save wf
		FILE*fp=fopen("wf_4","w");
		for(int i=0;i<NumAtoms;i++)
             		fprintf(fp,"%12.10le\n", wf[i]);
		fclose(fp);
	}
*/

    data_coulomb.AHF_wf(wf,A_aniso);
    double *A_aniso_final= new double[9*NumAtoms];	

MPI_Barrier(MPI_COMM_WORLD);
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(A_aniso, A_aniso_final, NumAtoms*9, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	for(int i=0;i<NumAtoms*9;i++)
		A_aniso_final[i]=A_aniso[i];	
#endif
    if (rank==0) {
        cout<<"Finished calculating Anisotropic hyperfine wavefunction"<<endl;
	cout<<"Writing to file"<<endl; 
	FILE* fp1=fopen("B_xx", "w");
	FILE* fp2=fopen("B_xy", "w");
	FILE* fp3=fopen("B_xz", "w");
	//FILE* fp4=fopen("B_yx", "w");
	FILE* fp5=fopen("B_yy", "w");
	FILE* fp6=fopen("B_yz", "w");
	//FILE* fp7=fopen("B_zx", "w");
	//FILE* fp8=fopen("B_zy", "w");
	FILE* fp9=fopen("B_zz", "w");
	FILE* fp10=fopen("A_xx", "w");
	FILE* fp11=fopen("A_yy", "w");
	FILE* fp12=fopen("A_yz", "w");

  	for(int i=0;i<NumAtoms;i++){
		fprintf(fp1,"%12.10le\n", A_aniso_final[i*9]);
		fprintf(fp2,"%12.10le\n", A_aniso_final[i*9+1]);
		fprintf(fp3,"%12.10le\n", A_aniso_final[i*9+2]);
		//fprintf(fp4,"%12.10le\n", A_aniso_final[i*9+3]);
		fprintf(fp5,"%12.10le\n", A_aniso_final[i*9+4]);
		fprintf(fp6,"%12.10le\n", A_aniso_final[i*9+5]);
		//fprintf(fp7,"%12.10le\n", A_aniso_final[i*9+6]);
                //fprintf(fp8,"%12.10le\n", A_aniso_final[i*9+7]);
                fprintf(fp9,"%12.10le\n", A_aniso_final[i*9+8]);
		fprintf(fp10,"%12.10le\n", (A_aniso_final[i*9]+8.0*3.1415926535/3.0*wf[i]));
		fprintf(fp11,"%12.10le\n", (A_aniso_final[i*9+4]+8.0*3.1415926535/3.0*wf[i]));
		fprintf(fp12,"%12.10le\n", (A_aniso_final[i*9+8]+8.0*3.1415926535/3.0*wf[i]));
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
        //fclose(fp4);
	fclose(fp5);
        fclose(fp6);
	//fclose(fp7);
        //fclose(fp8);
        fclose(fp9);
	fclose(fp10);
	 fclose(fp11);
	 fclose(fp12);
	}
 
}

if(ComputeAHF) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
	double a0=0.543095;
	double* P_coor=new double[3];
	double* Si_coor=new double[3];
	double xP=28*a0;
	double yP=28*a0;
	double zP=28*a0;
	double xSi=30*a0;
	double ySi=28*a0;
	double zSi=28*a0;

/*	FILE *fp=fopen("coordinates","r");
	if(fp!=NULL){
		fscanf(fp, "%lf %lf %lf\n", &xP, &yP, &zP);
		//fscanf(fp, "%f %f %f\n", &xSi, &ySi, &zSi);                
		
	}else{
		cout<<"File not found."<<endl;
	}

	fclose(fp);
*/
        P_coor[0]=xP;
        P_coor[1]=yP;
        P_coor[2]=zP;
        Si_coor[0]=xSi;
        Si_coor[1]=ySi;
        Si_coor[2]=zSi;

	double *A_aniso= new double[9];
	double *A_iso=new double[1];	
	cout<<"P_coordinates : ("<<P_coor[0]<<", "<<P_coor[1]<<", "<<P_coor[2]<<")"<<endl;
	cout<<"Si_coordinates : ("<<Si_coor[0]<<", "<<Si_coor[1]<<", "<<Si_coor[2]<<")"<<endl;

	data_coulomb.AHF(P_coor, Si_coor, A_aniso, A_iso);

	cout<<"Anisotropic hyperfine components"<<endl;
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		cout<<"A_"<<i<<j<<"    :   "<<A_aniso[i*3+j]<<endl;
	cout<<"Isotropic hyperfine component: "<<A_iso[0]<<endl;

  	for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
                cout<<A_aniso[i*3+j]<<endl;
	cout<<A_iso[0]<<endl;
/*
	//Total hyperfine contribution;
	double *A_total=new double[3];
	//need to know the correct hyperfine constant gamma_s*gamma_n*h_bar^2;
	double h_bar=1.054571628e-34; 		//J-s
	double Gamma_S=1.76085977e11;    	//1/T-s
	double Gamma_I=-8.458e6;   		//1/T-s
	double mu0_over_4pi=1e-7; 		//N/A-A
	double e_charge=1.602176487e-19;	// J/eV
	double h=4.13566733e-15; 		//eV-s
	double m2nm=1e27;			//nm/m
	double hyp_constant=Gamma_S*Gamma_I*h_bar*h_bar*mu0_over_4pi*m2nm/(e_charge*h);
	cout<<"Multiplicative factor   :   "<<hyp_constant<<endl;
	double iso_constant = (8.0/3.0)*3.1415926535;
	double hyp_constant_simplified=-8.458*1.76085977*1000*1.054571628*1.054571628/(1.602176487*4.13566733);
	cout<<"Multiplicative factor (simplified)  :   "<<hyp_constant_simplified<<endl;
	hyp_constant=1;
	for(int i=0;i<3;i++)
        	for(int j=0;j<3;j++){
        	 if(i==j)
                	A_total[3*i+j]=hyp_constant*(iso_constant*A_iso[0]+A_aniso[i*3+j]);
        	else
                	A_total[3*i+j]=hyp_constant*A_aniso[i*3+j];
        	cout<<"Total Hyperine Contribution : A_total_"<<i<<j<<"   "<<A_total[3*i+j]<<endl;
}
  for(int i=0;i<3;i++)
                for(int j=0;j<3;j++){
                 if(i==j)
                        A_total[3*i+j]=hyp_constant*(iso_constant*A_iso[0]+A_aniso[i*3+j]);
                else
                        A_total[3*i+j]=hyp_constant*A_aniso[i*3+j];
                cout<<A_total[3*i+j]<<endl;
}
*/

}

if(Compute_AHF_shell) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb =
          Coulomb(NumAtoms, myAtoms_begin, myAtoms_end, NumOrbitals);
    data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
	
	int numPoints=0;
	double* Si_coor;
	double xSi, ySi, zSi = 0.0;

	FILE *fp=fopen("coordinates","r");
	if(fp!=NULL){
		fscanf(fp,"%d\n", &numPoints);
		cout<<"Number of Points in this shell: "<<numPoints<<endl;
		Si_coor=new double[3*numPoints];
		for(int i=0;i<numPoints;i++){
			fscanf(fp, "%lf %lf %lf\n", &xSi, &ySi, &zSi);
			cout<<"Point "<<i<<"   :   ("<<xSi<<", "<<ySi<<", "<<zSi<<")"<<endl;
			Si_coor[3*i]=xSi;
			Si_coor[3*i+1]=ySi;
			Si_coor[3*i+2]=zSi;                
		}
	}else{
		cout<<"File not found."<<endl;
	}
	fclose(fp);

for(int k=0; k<numPoints; k++){
	
	double *A_aniso= new double[9];
	double *A_iso=new double[1];	
	double *A_total=new double[9];
	double *Si_pos=new double[3];

	Si_pos[0] = Si_coor[3*k];
	Si_pos[1] = Si_coor[3*k+1];
	Si_pos[2] = Si_coor[3*k+2];

	cout<<k<<"  -   "<<"Si_coordinates : ("<<Si_pos[0]<<", "<<Si_pos[1]<<", "<<Si_pos[2]<<")"<<endl;

	data_coulomb.AHF_shell(Si_pos, A_aniso, A_iso);

	cout<<"Anisotropic hyperfine components : "<<endl;
	for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		cout<<"A_"<<i<<j<<"    :   "<<A_aniso[i*3+j]<<endl;
	cout<<"Isotropic hyperfine component: "<<A_iso[0]<<endl;
	FILE* B_tensor=fopen("B_tensor","a");
        if(!B_tensor) {
		cout<<"Could not open file: B_tensor"<<endl;
		exit(-1);
	}else{		
		fprintf(B_tensor, "%12.10le %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le\n", 
		A_aniso[0], A_aniso[1], A_aniso[2], A_aniso[3], A_aniso[4], A_aniso[5], A_aniso[6], A_aniso[7], A_aniso[8], A_iso[0]);		
	}
	fclose(B_tensor);

	FILE* B_xx=fopen("B_xx","a");
        if(!B_xx) {
                cout<<"Could not open file: B_xx"<<endl;
                exit(-1);
        }else{
                fprintf(B_xx, "%12.10le\n", A_aniso[0]);         
        }
        fclose(B_xx);

	FILE* B_xy=fopen("B_xy","a");
	if(!B_xy) {
                cout<<"Could not open file: B_xy"<<endl;
                exit(-1);
        }else{                fprintf(B_xy, "%12.10le\n", A_aniso[1]);                                                                             
        }
        fclose(B_xy);

	FILE* B_xz=fopen("B_xz","a");
	if(!B_xz) {
                cout<<"Could not open file: B_xz"<<endl;
                exit(-1);
        }else{               
		fprintf(B_xz, "%12.10le\n", A_aniso[2]);                                                                         
        }
        fclose(B_xz);

	FILE* B_yy=fopen("B_yy","a");
        if(!B_yy) {
                cout<<"Could not open file: B_yy"<<endl;
                exit(-1);
        }else{               
                fprintf(B_yy, "%12.10le\n", A_aniso[4]);
        }
        fclose(B_yy);
 	
	FILE* B_yz=fopen("B_yz","a");
	if(!B_yz) {
                cout<<"Could not open file: B_yz"<<endl;
                exit(-1);
        }else{
                fprintf(B_yz, "%12.10le\n", A_aniso[5]);
        }
        fclose(B_yz);

	FILE* B_zz=fopen("B_zz","a");
	if(!B_zz) {
                cout<<"Could not open file: B_zz"<<endl;
                exit(-1);
        }else{
                fprintf(B_zz, "%12.10le\n", A_aniso[8]);
        }
        fclose(B_zz);


	for(int i=0;i<9;i++)
		A_total[i]=A_aniso[i];	
	A_total[0]+=(8.0/3.0)*3.1415926535*A_iso[0];
	A_total[4]+=(8.0/3.0)*3.1415926535*A_iso[0];
	A_total[8]+=(8.0/3.0)*3.1415926535*A_iso[0];
	
	FILE* A_tensor=fopen("A_tensor","a");
        if(!A_tensor){
                cout<<"Could not open file: A_tensor"<<endl;
                exit(-1);
        }else{
                fprintf(A_tensor, "%12.10le %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le %12.10le\n", A_total[0], A_total[1],
                        A_total[2], A_total[3], A_total[4], A_total[5], A_total[6], A_total[7], A_total[8]);
        }
        fclose(A_tensor);

	FILE* A_isotropic=fopen("A_iso","a");
        if(!A_isotropic){
                cout<<"Could not open file: A_isotropic"<<endl;
                exit(-1);
        }else{
                fprintf(A_isotropic, "%12.10le\n", A_iso[0]);
        }
        fclose(A_isotropic);
	
	char filename[200];
	sprintf(filename, "AHF_B_%d", k);
	FILE* fp1=fopen(filename, "w");
	if(!fp1){
		cout<<"Error opening file: AHF_B "<<endl;
		exit(-1);
	}else{
		for(int i=0;i<9;i++)
			fprintf(fp1,"%12.10le\n", A_aniso[i]);
		fprintf(fp1, "%12.10le\n\n", A_iso[0]);				
	}
	fclose(fp1);
	
	char filename1[200];
        sprintf(filename1, "AHF_A_%d", k);
        FILE* fp2=fopen(filename1, "w");
        if(!fp2){
                cout<<"Error opening file: AHF_A "<<endl;
                exit(-1);
        }else{
                for(int i=0;i<9;i++)
                        fprintf(fp2,"%12.10le\n", A_total[i]);
        }
        fclose(fp2);

	delete[] A_aniso;
	delete[] A_iso;
	delete[] A_total;
	delete[] Si_pos;
} //end of k-loop	
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
