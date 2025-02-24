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

#include "nml_dcvector.h"
#include "Coulomb.h"
#include "ReadData.h"
#include "SimpleParser.h"
//#include "DispersionUnfold.h"
using namespace std;

#ifdef  __cplusplus
extern "C" {
#endif

extern void zheev_(char *JOBZ, char *UPLO, int *N, nml_dcvector *A, int *LDA, double *W, nml_dcvector *WORK, int *LWORK,
                double *RWORK, int *INFO);

extern void zhegv_(int *ITYPE, char *JOBZ, char *UPLO, int *N, nml_dcvector *A, int *LDA, nml_dcvector *B, int *LDB, double *W, 
		nml_dcvector *WORK, int *LWORK, double *RWORK, int *INFO);

extern void zheevr_(char *JOBZ, char *RANGE, char *UPLO, int *N, nml_dcvector *A, int *LDA, double *VL, double *VU, int *IL, 
		int *IU, double *ABSTOL, int *M, double *W, nml_dcvector *Z, int *LDZ, int *ISUPPZ, nml_dcvector *WORK, 
		int *LWORK, double *RWORK, int *LRWORK, int *IWORK, int *LIWORK, int *INFO);

#ifdef  __cplusplus
}
#endif


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
  int steps=0;
  int start_elem=0;
  int end_elem=0;
  int basis=0; 	
  int step_num=0;
  int selected=0; 

//Read Options =================================================
  vector<string> options;
  for(int i=0; i<argc; i++) options.push_back(argv[i]);
  unsigned int iopt = 0;
  bool ComputeCoulomb = false;
  bool ComputeCoulombHM = false;
  bool ComputeHartree = false;
  bool ComputeHL20 = false;
  bool ComputeSpin = false;
  bool Compute_CI = false;
  bool Compute_CI_Matrix = false;
  bool Compute_CI_Additional = false;
  bool file_dump = false;
  bool find_additional = false;
  bool ComputeCoulomb_new=false;
  bool Compute_hund_mulliken=false;	
  bool Compute_Kramer_vectors=false;
  bool Compute_CI_Matrix_Selected=false;
  bool GenerateSpin_wfs=false;
  bool dielectric=false;
  bool Compute_CI_fast=false; 

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
    else if(options[iopt] == "--basis") {
      basis = atoi(options[++iopt].c_str());	
    }		
    else if(options[iopt] == "--steps") {      
      steps = atoi(options[++iopt].c_str());
    }
    else if(options[iopt] == "--start") {
      start_elem = atoi(options[++iopt].c_str());
    }
    else if(options[iopt] == "--end") {
      end_elem = atoi(options[++iopt].c_str());
    }
    else if(options[iopt] == "--stepnumber") {
      step_num = atoi(options[++iopt].c_str());
    }	
    else if(options[iopt] == "--selected") {
      selected = atoi(options[++iopt].c_str());
    }
    if(options[iopt] == "--hl") {
      ComputeCoulomb = true;
    }
    if(options[iopt] == "--hl_new") {
      ComputeCoulomb_new = true;
    }
    else if(options[iopt] == "--HM") {
      ComputeCoulombHM = true;
    }
    else if(options[iopt] == "--hm") {
      ComputeHartree = true;
    }
    else if(options[iopt] == "--hl20") {
      ComputeHL20 = true;
    }
    else if(options[iopt] == "--spin") {
      ComputeSpin = true;
    }
    else if(options[iopt] == "--spin_wfs") {
      GenerateSpin_wfs = true;
    }	
    else if(options[iopt] == "--kramer") {
      Compute_Kramer_vectors = true;
    }
    else if(options[iopt] == "--ci_matrix_selected") {
      Compute_CI_Matrix_Selected = true;
    }
    else if(options[iopt] == "--file") {
      file_dump = true;
    }
    else if(options[iopt] == "--ci") {
      Compute_CI = true;
    }
    else if(options[iopt] == "--ci_fast") {
      Compute_CI_fast = true;
    }
    else if(options[iopt] == "--hund") {
      Compute_hund_mulliken = true;
    }
    else if(options[iopt] == "--ci_matrix") {
      Compute_CI_Matrix = true;
    }
    else if(options[iopt] == "--ci_additional") {
      Compute_CI_Additional = true;
    }
    else if(options[iopt] == "--ci_add") {
      find_additional = true;
    }
     else if(options[iopt] == "--dielectric") {
      dielectric = true;
      cout<<"Chose the dielectric option"<<endl;	
    }

    iopt++;
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
  //complex<double>* wf_h;
  int* atomid;
  int* neighbor;
  double* lattice;
  double* Ee = new double[NumElectrons];
  double* Eh = new double[NumHoles];

  read_natoms_norbitals(data, &NumAtoms, &NumOrbitals, eList);

  cout<<"Number of atoms: "<<NumAtoms<<endl;
  cout<<"Number of orbitals: "<<NumOrbitals<<endl;

//Parallelization of Data ============================================
  
  int truncate = sqrt((double)size); // need to use (integer)^2 processrors. 
  int bat= NumAtoms/truncate;
  int residue = NumAtoms%truncate;
  int *myAtoms_begin,*myAtoms_end;

  myAtoms_begin = new int [truncate];
  myAtoms_end = new int [truncate];

  for(int indx = 0; indx < truncate; indx++) {

     if(indx<residue) {
        myAtoms_begin[indx] = indx*(bat+1);
        myAtoms_end[indx] = myAtoms_begin[indx] + bat + 1;
     }

     else {
        myAtoms_begin[indx] = residue*(bat+1) + (indx - residue)*bat;
        myAtoms_end[indx]   = myAtoms_begin[indx] + bat;
     }

  }

  int x1_start = 0, x1_end = 0, x2_start = 0, x2_end = 0;
  int flag = 0;
 
  for (int indx = 0; indx < truncate; indx++) {
    for (int indx2 = 0; indx2 < truncate; indx2++) {
        if(rank == (truncate*indx + indx2))  {
          flag = 1;
          x1_start = myAtoms_begin[indx]; x1_end = myAtoms_end[indx];
          x2_start = myAtoms_begin[indx2]; x2_end = myAtoms_end[indx2];
          break;
        }
    }

    if(flag == 1) break;
  }

  if (rank >= truncate*truncate) {

	  x1_start = 0; x1_end = 0;
          x2_start = 0; x2_end = 0;
         
  }

 // read_data_for_postprocessing(data, NumElectrons, eList, NumHoles, hList,
 //                              &wf_e, &wf_h, &atomid, &neighbor, &lattice, &Ee, &Eh, x1_start, x1_end, x2_start, x2_end);
 
read_data_for_postprocessing_new(data, NumElectrons, eList, NumHoles, hList, &atomid, &neighbor, &lattice, &Ee, &Eh);


  for(int iatom=0; iatom<NumAtoms; iatom++) {
	
     if(atomid[iatom]==12) atomid[iatom]=0; //As atom
     else if(atomid[iatom]==10) atomid[iatom]=1; //Ga atom
     else if(atomid[iatom]==14) atomid[iatom]=2; // In atom
     else if(atomid[iatom]==6) atomid[iatom]=3; // Al atom

     if(atomid[iatom]==12) atomid[iatom]=0; //As atom
     else if(atomid[iatom]==10) atomid[iatom]=1; //Ga atom
     else if(atomid[iatom]==6) atomid[iatom]=2; // Al atom

     if(atomid[iatom] == 7) atomid[iatom]=0; // Si atom
     else if(atomid[iatom] == 8) atomid[iatom]=1; // P atom

     //else cout<<"atomid does not match with As, Ga, In, Al atom"<<endl;

  } 
//===================================================================

if(GenerateSpin_wfs){

read_write_wf_10(data, NumElectrons, eList, NumAtoms);

}


 if(Compute_CI) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb =
          Coulomb(NumAtoms, x1_start, x1_end, x2_start, x2_end, NumOrbitals);
    ifstream coulomb_input("coulomb.table");
    data_coulomb.assign_table(coulomb_input);
       //data_coulomb.assign_table1();
    //data_coulomb.get_data(wf_e, atomid, neighbor, lattice);

    int numStates=(basis*(basis-1))/2; 		
    int num_Matrix_Elements=(numStates*numStates-numStates)/2+numStates;	
     	
if(rank==0){
    cout<<"Single electron basis : "<<basis<<endl;
    cout<<"Number of steps : "<<steps<<endl;
    cout<<"Starting element : "<<start_elem<<endl;
    cout<<"Ending element : "<<end_elem<<endl;
    cout<<"Step Number : "<<step_num<<endl;
    cout<<"Number of Two electron States : "<<numStates<<endl;
    cout<<"Number of Matrix Elements : "<<num_Matrix_Elements<<endl;
}//end of rank
   
	int *first_state=new int[numStates];
	int *second_state=new int[numStates];
	int count=0;
	for(int i=1;i<basis;i++)
        for(int j=i+1;j<=basis;j++){
                first_state[count]=i;
                second_state[count]=j;
                count++;
        }
if(rank==0){	
	for(int i=0;i<numStates;i++)
	cout<<"CI state "<<i<<"    :   "<<first_state[i]<<"  "<<second_state[i]<<endl;

	if(file_dump==true){
		FILE* fp=fopen("States","w");
        	for(int i=0;i<numStates;i++)
        	fprintf(fp,"%i %i\n", first_state[i], second_state[i]);
		fclose(fp);
	}
} //end of rank

int *first_elem=new int[num_Matrix_Elements];
int *second_elem=new int[num_Matrix_Elements];
int *third_elem=new int[num_Matrix_Elements];
int *fourth_elem=new int[num_Matrix_Elements];
count=0;

for(int i=0;i<numStates;i++){
	for(int j=i;j<numStates;j++){
		first_elem[count]=first_state[i];
		second_elem[count]=second_state[i];
		third_elem[count]=first_state[j];
		fourth_elem[count]=second_state[j];
		count++;
	}
}

if(rank==0){
	for(int i=0;i<num_Matrix_Elements;i++)
		cout<<"Index "<<i<<"  :  "<<first_elem[i]<<"  "<<second_elem[i]<<"  "<<third_elem[i]<<"  "<<fourth_elem[i]<<endl;
	if(file_dump==true){
		FILE* fp1=fopen("Matrix_Elements","w");
        	for(int i=0;i<num_Matrix_Elements;i++)
                	fprintf(fp1, "%i %i %i %i\n", first_elem[i], second_elem[i], third_elem[i], fourth_elem[i]);
        	fclose(fp1);
}
}//end of rank

 	FILE* fp1=fopen("exchange_parameters","r");
        float dd_cutoff=12.0;
        double dd_cutoff1=dd_cutoff;
        int use_dd_cutoff=1;
        if(fp1!=NULL){
                fscanf(fp1, "%i\n", &use_dd_cutoff);
                fscanf(fp1, "%f\n", &dd_cutoff);
                fclose(fp1);
                dd_cutoff1=dd_cutoff;
        }else{
                if(rank==0)
                        cout<<"File not found: Using default parameters. dd_cutoff=12 nm"<<endl;
        }
        if(rank==0){
                cout<<"Use dd_cutoff: "<<use_dd_cutoff<<endl;
                cout<<"dd_cutoff: "<<dd_cutoff1<<"  nm"<<endl;
        }

    double elapsed_time=0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time=-MPI_Wtime();
    int *electron_states=new int[4];	
if(rank==0){
	char filename_out[200];
        sprintf(filename_out, "steps_%i",step_num);
	FILE* fp_out=fopen(filename_out,"w");
	if(fp_out!=NULL)
		fprintf(fp_out,"%i\n",end_elem-start_elem+1);
	else{
		cout<<"Error Opening file"<<endl;
		exit(-1);
	}
	fclose(fp_out);
}

cout<<"Entering matrix element loop"<<endl;

for(int i=start_elem;i<=end_elem;i++){

	electron_states[0]=first_elem[i];
	electron_states[1]=second_elem[i];
	electron_states[2]=third_elem[i];
	electron_states[3]=fourth_elem[i];

	cout<<"Calculating matrix element : "<<i<<endl;

	complex<double>* wf_e1;
	read_wf(data, 4, electron_states, &wf_e1, x1_start, x1_end, x2_start, x2_end);

	data_coulomb.get_data(wf_e1, atomid, neighbor, lattice);
	if(dielectric)
		data_coulomb.compute_coulomb_exchange4(dd_cutoff1,data,electron_states,4);
	else
		data_coulomb.compute_coulomb_exchange4_no_dielectric_correction(dd_cutoff1,data,electron_states,4);

	complex<double>* local_coul = data_coulomb.get_H_coul();
    	complex<double>* local_exch = data_coulomb.get_H_exch();

    	int csize = 1;
    	complex<double>* total_coul=new complex<double>[csize];
    	complex<double>* total_exch= new complex<double>[csize];

#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
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
       cout<<"Matrix Element : "<<i<<endl;
       cout<<"States  : "<<first_elem[i]<<"  "<<second_elem[i]<<"  "<<third_elem[i]<<"  "<<fourth_elem[i]<<endl;	
       cout<<"Coulomb:  "<<setw(10)<<total_coul[0]<<endl;
       cout<<"Exchange: "<<setw(10)<<total_exch[0]<<endl;
       cout<<endl;
       char filename_out[200];
       sprintf(filename_out, "steps_%i", step_num);
       FILE* fp_out=fopen(filename_out,"a");
       if(fp_out!=NULL){
       	fprintf(fp_out, "%i %i %i %i %12.10le %12.10le %12.10le %12.10le\n", first_elem[i], second_elem[i],
	third_elem[i], fourth_elem[i], total_coul[0].real(), total_coul[0].imag(), total_exch[0].real(), total_exch[0].imag());        
       }else{
	 cout<<"Error opening file for output."<<endl;
	 exit(-1);
       }
       fclose(fp_out);
   }//end of rank	
	delete(wf_e1);
}//end of for

elapsed_time+=MPI_Wtime();
if(rank==0)
cout<<"Time to compute the matrix elements:"<<elapsed_time<<endl;
}

if(find_additional){
if(rank==0){
	cout<<"steps  :  "<<steps<<endl;
	int numStates=(basis*(basis-1))/2;
    	int num_Matrix_Elements=(numStates*numStates-numStates)/2+numStates;
        int *first_state=new int[numStates];
        int *second_state=new int[numStates];
        int count=0;
	cout<<"Basis = "<<basis<<endl;
	cout<<"Number of States :  "<<numStates<<endl;
	cout<<"Number of matrix elements : "<<num_Matrix_Elements<<endl;
	cout<<"States :  "<<endl;
        for(int i=1;i<basis;i++)
        for(int j=i+1;j<=basis;j++){
                first_state[count]=i;
                second_state[count]=j;
		cout<<first_state[count]<<"  "<<second_state[count]<<endl;
     		count++;
	   }

	FILE* fp_states=fopen("New_States","w");
	if(fp_states==NULL)
		cout<<"Error Opening file : New States "<<endl;
	else{
		fprintf(fp_states, "%i\n",numStates);	
                for(int i=0;i<numStates;i++)
                	fprintf(fp_states,"%i %i\n", first_state[i], second_state[i]);
	}             
   	fclose(fp_states);


int *first_elem=new int[num_Matrix_Elements];
int *second_elem=new int[num_Matrix_Elements];
int *third_elem=new int[num_Matrix_Elements];
int *fourth_elem=new int[num_Matrix_Elements];
bool *mat_elem=new bool[num_Matrix_Elements];
for(int i=0;i<num_Matrix_Elements;i++)
	mat_elem[i]=false;

count=0;
cout<<endl;
cout<<"Matrix Elements  :  "<<endl;
for(int i=0;i<numStates;i++){
        for(int j=i;j<numStates;j++){
                first_elem[count]=first_state[i];
                second_elem[count]=second_state[i];
                third_elem[count]=first_state[j];
                fourth_elem[count]=second_state[j];
		cout<<first_elem[count]<<"  "<<second_elem[count]<<"  "<<third_elem[count]<<"  "<<fourth_elem[count]<<endl;
		count++;
        }
}
cout<<endl;
int count_old_states=0;

for(int i=0;i<steps;i++){
        char filename[200];
        sprintf(filename, "steps_%i", i);
        cout<<"Reading file : "<<filename<<endl;
        double coul_tmp_r=0.0;
        double coul_tmp_i=0.0;
        double exch_tmp_r=0.0;
        double exch_tmp_i=0.0;
        int first = 0;
        int second = 0;
        int third = 0;
        int fourth = 0;
        int count = 0;
        FILE* fp_tmp = fopen(filename, "rb");
        if(fp_tmp==NULL)
                cout<<"Error Opening File : "<<filename<<endl;
        else{
                fscanf(fp_tmp,"%i\n", &count);
                cout<<"Step : "<<i<<"      Elements :  "<<count<<endl<<endl;
                for(int j=0;j<count;j++){
                        first, second, third, fourth=0;
                        coul_tmp_r, coul_tmp_i, exch_tmp_r, exch_tmp_i = 0.0;
                        fscanf(fp_tmp,"%i %i %i %i %lf %lf %lf %lf\n",&first, &second,
                        &third, &fourth, &coul_tmp_r, &coul_tmp_i, &exch_tmp_r, &exch_tmp_i);
                        //cout<<first<<"  "<<second<<"   "<<third<<"   "<<fourth<<"   "<<endl;
                        //cout<<coul_tmp_r<<"   "<<coul_tmp_i<<"   "<<exch_tmp_r<<"    "<<exch_tmp_i<<endl<<endl;
                        for(int k=0;k<num_Matrix_Elements;k++){
  				if((first_elem[k]==first)&&(second_elem[k]==second)&&(third_elem[k]==third)&&(fourth_elem[k]==fourth)){
					mat_elem[k]=true;	
					count_old_states++;
					//cout<<"Found element at "<<k<<endl;
					break; 
				}                     
			}//end of k	
                }//end of count
        }//end of if
        fclose(fp_tmp);
}//end of steps

//int new_States=numStates-count_old_states;
//cout<<"Number of additional states"<<new_States<<endl;

count=0;
for(int i=0;i<num_Matrix_Elements;i++)
	if(mat_elem[i]==false)
		count++;
cout<<"Number of new matrix elements : "<<count<<endl;

FILE *new_matrix_elements=fopen("New_Matrix_Elements","w");
if(new_matrix_elements==NULL){
        cout<<"Error opening file for writing : New Matrix Elements"<<endl;
        exit(-1);
}else{
        fprintf(new_matrix_elements, "%i\n", count);
        for(int i=0;i<num_Matrix_Elements;i++){
                if(mat_elem[i]==false){
                        fprintf(new_matrix_elements, "%i %i %i %i\n", first_elem[i], second_elem[i], third_elem[i], fourth_elem[i]);
                }
        }
}//end of if
fclose(new_matrix_elements);

}//end of rank
}

 if(Compute_CI_fast) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb =
          Coulomb(NumAtoms, x1_start, x1_end, x2_start, x2_end, NumOrbitals);
    ifstream coulomb_input("coulomb.table");
    data_coulomb.assign_table(coulomb_input);
       //data_coulomb.assign_table1();
    //data_coulomb.get_data(wf_e, atomid, neighbor, lattice);

    int numStates=(basis*(basis-1))/2; 		
    int num_Matrix_Elements=(numStates*numStates-numStates)/2+numStates;	
     	
if(rank==0){
    cout<<"Single electron basis : "<<basis<<endl;
    cout<<"Number of steps : "<<steps<<endl;
    cout<<"Starting element : "<<start_elem<<endl;
    cout<<"Ending element : "<<end_elem<<endl;
    cout<<"Step Number : "<<step_num<<endl;
    cout<<"Number of Two electron States : "<<numStates<<endl;
    cout<<"Number of Matrix Elements : "<<num_Matrix_Elements<<endl;
}//end of rank
   
	int *first_state=new int[numStates];
	int *second_state=new int[numStates];
	int count=0;
	for(int i=1;i<basis;i++)
        for(int j=i+1;j<=basis;j++){
                first_state[count]=i;
                second_state[count]=j;
                count++;
        }
if(rank==0){	
	for(int i=0;i<numStates;i++)
	cout<<"CI state "<<i<<"    :   "<<first_state[i]<<"  "<<second_state[i]<<endl;

	if(file_dump==true){
		FILE* fp=fopen("States","w");
        	for(int i=0;i<numStates;i++)
        	fprintf(fp,"%i %i\n", first_state[i], second_state[i]);
		fclose(fp);
	}
} //end of rank

int *first_elem=new int[num_Matrix_Elements];
int *second_elem=new int[num_Matrix_Elements];
int *third_elem=new int[num_Matrix_Elements];
int *fourth_elem=new int[num_Matrix_Elements];
count=0;

for(int i=0;i<numStates;i++){
	for(int j=i;j<numStates;j++){
		first_elem[count]=first_state[i];
		second_elem[count]=second_state[i];
		third_elem[count]=first_state[j];
		fourth_elem[count]=second_state[j];
		count++;
	}
}

if(rank==0){
	for(int i=0;i<num_Matrix_Elements;i++)
		cout<<"Index "<<i<<"  :  "<<first_elem[i]<<"  "<<second_elem[i]<<"  "<<third_elem[i]<<"  "<<fourth_elem[i]<<endl;
	if(file_dump==true){
		FILE* fp1=fopen("Matrix_Elements","w");
        	for(int i=0;i<num_Matrix_Elements;i++)
                	fprintf(fp1, "%i %i %i %i\n", first_elem[i], second_elem[i], third_elem[i], fourth_elem[i]);
        	fclose(fp1);
}
}//end of rank

 	FILE* fp1=fopen("exchange_parameters","r");
        float dd_cutoff=12.0;
        double dd_cutoff1=dd_cutoff;
        int use_dd_cutoff=1;
        if(fp1!=NULL){
                fscanf(fp1, "%i\n", &use_dd_cutoff);
                fscanf(fp1, "%f\n", &dd_cutoff);
                fclose(fp1);
                dd_cutoff1=dd_cutoff;
        }else{
                if(rank==0)
                        cout<<"File not found: Using default parameters. dd_cutoff=12 nm"<<endl;
        }
        if(rank==0){
                cout<<"Use dd_cutoff: "<<use_dd_cutoff<<endl;
                cout<<"dd_cutoff: "<<dd_cutoff1<<"  nm"<<endl;
        }

    double elapsed_time=0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time=-MPI_Wtime();
    int *electron_states=new int[4];	
if(rank==0){
	char filename_out[200];
        sprintf(filename_out, "steps_%i",step_num);
	FILE* fp_out=fopen(filename_out,"w");
	if(fp_out!=NULL)
		fprintf(fp_out,"%i\n",end_elem-start_elem+1);
	else{
		cout<<"Error Opening file"<<endl;
		exit(-1);
	}
	fclose(fp_out);
}

cout<<"Entering matrix element loop"<<endl;

for(int i=start_elem;i<=end_elem;i++){

	electron_states[0]=first_elem[i];
	electron_states[1]=second_elem[i];
	electron_states[2]=third_elem[i];
	electron_states[3]=fourth_elem[i];
	
	cout<<"Calculating matrix element : "<<i<<endl;
	bool skip=false;
	int skip_test1=first_elem[i]%2 + second_elem[i]%2 + third_elem[i]%2 + fourth_elem[i]%2;
	int group1=first_elem[i]%2 + second_elem[i]%2;
	int group2=third_elem[i]%2 + fourth_elem[i]%2;
 	//cout<<"test 1 : "<<skip_test1<<"   group 1 : "<<group1<<"    group 2 : "<<group2<<endl;	
	if((skip_test1==1) || (skip_test1==3) || (group1==0 && group2==2) || (group1==2 && group2==0))
		skip=true;

	//cout<<first_elem[i]<<"  "<<second_elem[i]<<"    "<<third_elem[i]<<"   "<<fourth_elem[i]<<endl;

if(!skip){
	//cout<<"Need to calculate this element"<<endl;
	complex<double>* wf_e1;
	read_wf(data, 4, electron_states, &wf_e1, x1_start, x1_end, x2_start, x2_end);

	data_coulomb.get_data(wf_e1, atomid, neighbor, lattice);
	if(dielectric)
		data_coulomb.compute_coulomb_exchange4(dd_cutoff1,data,electron_states,4);
	else
		data_coulomb.compute_coulomb_exchange4_no_dielectric_correction(dd_cutoff1,data,electron_states,4);

	complex<double>* local_coul = data_coulomb.get_H_coul();
    	complex<double>* local_exch = data_coulomb.get_H_exch();

    	int csize = 1;
    	complex<double>* total_coul=new complex<double>[csize];
    	complex<double>* total_exch= new complex<double>[csize];

#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
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
       cout<<"Matrix Element : "<<i<<endl;
       cout<<"States  : "<<first_elem[i]<<"  "<<second_elem[i]<<"  "<<third_elem[i]<<"  "<<fourth_elem[i]<<endl;	
       cout<<"Coulomb:  "<<setw(10)<<total_coul[0]<<endl;
       cout<<"Exchange: "<<setw(10)<<total_exch[0]<<endl;
       cout<<endl;
       char filename_out[200];
       sprintf(filename_out, "steps_%i", step_num);
       FILE* fp_out=fopen(filename_out,"a");
       if(fp_out!=NULL){
       	fprintf(fp_out, "%i %i %i %i %12.10le %12.10le %12.10le %12.10le\n", first_elem[i], second_elem[i],
	third_elem[i], fourth_elem[i], total_coul[0].real(), total_coul[0].imag(), total_exch[0].real(), total_exch[0].imag());        
       }else{
	 cout<<"Error opening file for output."<<endl;
	 exit(-1);
       }
       fclose(fp_out);
   }//end of rank	
	delete(wf_e1);
}else{
 cout<<"Matrix element zero by symmetry. Skipping calculation"<<endl;
if (rank==0) {
       cout<<endl<<"====================================================="<<endl;
       cout<<"Matrix Element : "<<i<<endl;
       cout<<"States  : "<<first_elem[i]<<"  "<<second_elem[i]<<"  "<<third_elem[i]<<"  "<<fourth_elem[i]<<endl;
	complex<double> Jt=0.0;
	complex<double> Kt=0.0;
       cout<<"Coulomb:  "<<setw(10)<<Jt<<endl;
       cout<<"Exchange: "<<setw(10)<<Kt<<endl;
       cout<<endl;
       char filename_out[200];
       sprintf(filename_out, "steps_%i", step_num);
       FILE* fp_out=fopen(filename_out,"a");
       if(fp_out!=NULL){
        fprintf(fp_out, "%i %i %i %i %12.10le %12.10le %12.10le %12.10le\n", first_elem[i], second_elem[i],
        third_elem[i], fourth_elem[i], Jt.real(), Jt.imag(), Kt.real(), Kt.imag());
       }else{
         cout<<"Error opening file for output."<<endl;
         exit(-1);
       }
       fclose(fp_out);
   }//end of rank
}//end of skip
}//end of for

elapsed_time+=MPI_Wtime();
if(rank==0)
cout<<"Time to compute the matrix elements:"<<elapsed_time<<endl;
}


//Compute additional elements related to basis expansion

if(Compute_CI_Additional) {
    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb =
          Coulomb(NumAtoms, x1_start, x1_end, x2_start, x2_end, NumOrbitals);
    ifstream coulomb_input("coulomb.table");
    data_coulomb.assign_table(coulomb_input);
       //data_coulomb.assign_table1();
    //data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
  int count=0;
  int first, second, third, fourth = 0;
  int *first_elem=new int[end_elem-start_elem+1];
  int *second_elem=new int[end_elem-start_elem+1];
  int *third_elem=new int[end_elem-start_elem+1];    
  int *fourth_elem=new int[end_elem-start_elem+1];
 			
  FILE *fp=fopen("New_Matrix_Elements","r");
  if(fp==NULL){
	cout<<"Error opening file : New_Matrix_Elements"<<endl;
	exit(-1);
  }else{
	fscanf(fp, "%i\n", &count);
	cout<<"Reading file : elements="<<count<<endl;
	for(int i=0;i<count;i++){
		first, second, third, fourth = 0;
		fscanf(fp,"%i %i %i %i\n",&first, &second, &third, &fourth);
		if((i>=start_elem)&&(i<=end_elem)){
			first_elem[i-start_elem]=first;
			second_elem[i-start_elem]=second;
			third_elem[i-start_elem]=third;
                        fourth_elem[i-start_elem]=fourth;
		}
	}
  }
  fclose(fp);

//test printing
cout<<"Matrix elements to be computed : "<<endl;
for(int i=0;i<end_elem-start_elem+1;i++)
	cout<<first_elem[i]<<"  "<<second_elem[i]<<"  "<<third_elem[i]<<"  "<<fourth_elem[i]<<endl;			
cout<<endl;

        FILE* fp1=fopen("exchange_parameters","r");
        float dd_cutoff=12.0;
        double dd_cutoff1=dd_cutoff;
        int use_dd_cutoff=1;
        if(fp1!=NULL){
                fscanf(fp1, "%i\n", &use_dd_cutoff);
                fscanf(fp1, "%f\n", &dd_cutoff);
                fclose(fp1);
                dd_cutoff1=dd_cutoff;
        }else{
                if(rank==0)
                        cout<<"File not found: Using default parameters. dd_cutoff=12 nm"<<endl;
        }
        if(rank==0){
                cout<<"Use dd_cutoff: "<<use_dd_cutoff<<endl;
                cout<<"dd_cutoff: "<<dd_cutoff1<<"  nm"<<endl;
        }

	if(rank==0){        
		char filename_out[200];        
		sprintf(filename_out, "steps_%i",step_num);        
		FILE* fp_out=fopen(filename_out,"w");        
		if(fp_out!=NULL)                
			fprintf(fp_out,"%i\n",end_elem-start_elem+1);        
		else{           
     			cout<<"Error Opening file"<<endl;
     			exit(-1);        
		}        
		fclose(fp_out);
	}

    double elapsed_time=0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time=-MPI_Wtime();

int *electron_states=new int[4];
cout<<"Entering matrix element loop"<<endl;

for(int i=0;i<end_elem-start_elem+1;i++){

        electron_states[0]=first_elem[i];
        electron_states[1]=second_elem[i];
        electron_states[2]=third_elem[i];
        electron_states[3]=fourth_elem[i];

        cout<<"Calculating matrix element : "<<i<<endl;

        complex<double>* wf_e1;
        read_wf(data, 4, electron_states, &wf_e1, x1_start, x1_end, x2_start, x2_end);

        data_coulomb.get_data(wf_e1, atomid, neighbor, lattice);
	if(dielectric)
        	data_coulomb.compute_coulomb_exchange4(dd_cutoff1,data,electron_states,4);
	else
		data_coulomb.compute_coulomb_exchange4_no_dielectric_correction(dd_cutoff1,data,electron_states,4);
	

        complex<double>* local_coul = data_coulomb.get_H_coul();
        complex<double>* local_exch = data_coulomb.get_H_exch();

        int csize = 1;
        complex<double>* total_coul=new complex<double>[csize];
        complex<double>* total_exch= new complex<double>[csize];

#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
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
	cout<<"Matrix Element : "<<i<<endl;
       cout<<"States  : "<<first_elem[i]<<"  "<<second_elem[i]<<"  "<<third_elem[i]<<"  "<<fourth_elem[i]<<endl;
       cout<<"Coulomb:  "<<setw(10)<<total_coul[0]<<endl;
       cout<<"Exchange: "<<setw(10)<<total_exch[0]<<endl;
       cout<<endl;
       char filename_out[200];
       sprintf(filename_out, "steps_%i", step_num);
       FILE* fp_out=fopen(filename_out,"a");
       if(fp_out!=NULL){
        fprintf(fp_out, "%i %i %i %i %12.10le %12.10le %12.10le %12.10le\n", first_elem[i], second_elem[i],
        third_elem[i], fourth_elem[i], total_coul[0].real(), total_coul[0].imag(), total_exch[0].real(), total_exch[0].imag());
       }else{
         cout<<"Error opening file for output."<<endl;
         exit(-1);
       }
       fclose(fp_out);
   }//end of rank
        delete(wf_e1);
}//end of for

elapsed_time+=MPI_Wtime();
if(rank==0)
cout<<"Time to compute the matrix elements:"<<elapsed_time<<endl;
}


if(Compute_CI_Matrix) {
if(rank==0){
	double E_ref=1.1317;
	//Compute Many body Eigen States  
	int numStates=(basis*(basis-1))/2;
    	int num_Matrix_Elements=(numStates*numStates-numStates)/2+numStates;
	nml_dcvector *coul_matrix=nml_dcv_new(numStates*numStates);
	nml_dcvector *exch_matrix=nml_dcv_new(numStates*numStates);
	double *Ek=new double[basis];
	cout<<"Basis : "<<basis<<endl;
	cout<<"Number of Matrix Elements : "<<num_Matrix_Elements<<endl;
	cout<<"Steps : "<<steps<<endl;

  	int *first_state=new int[numStates];
        int *second_state=new int[numStates];
        int count=0;
        for(int i=1;i<basis;i++)
        for(int j=i+1;j<=basis;j++){
                first_state[count]=i;
                second_state[count]=j;
                count++;
        }
	// read in eigenvalues
   	FILE* fp;
   	char filename[200];

   	sprintf(filename, "Ek");
   	fp=fopen(filename,"r");
   	cout << "Reading eigenvalues:  " << filename << endl;
   	if (!fp) die("Could not find file %s\n", filename);
   	for(int i=0;i<basis;i++)
        	fscanf(fp, "%lf\n", &Ek[i]);
   	fclose(fp);

	//print Eigenvalues
	cout<<"Single electron eigenvalues   "<<endl;
	for(int i=0;i<basis;i++){
        	Ek[i]=(Ek[i]-E_ref)*1000;
		cout<<Ek[i]<<endl;
	}

//Read Coulomb and Exchange Matrix Elements Here 
for(int i=0;i<steps;i++){
	char filename[200];
        sprintf(filename, "steps_%i", i);
	cout<<"Reading file : "<<filename<<endl;
	double coul_tmp_r=0.0;
	double coul_tmp_i=0.0;
	double exch_tmp_r=0.0;
	double exch_tmp_i=0.0;
	int first = 0;
	int second = 0;
	int third = 0;
	int fourth = 0;
	int count = 0;
	int row_index=0;
	int col_index=0;
        FILE* fp_tmp = fopen(filename, "rb");
    	if(fp_tmp==NULL)
		cout<<"Error Opening File : "<<filename<<endl;
    	else{
		fscanf(fp_tmp,"%i\n", &count);
		cout<<"Step : "<<i<<"      Elements :  "<<count<<endl<<endl;
		for(int j=0;j<count;j++){
			first, second, third, fourth=0;
			coul_tmp_r, coul_tmp_i, exch_tmp_r, exch_tmp_i = 0.0;	
			row_index, col_index=-1;
			fscanf(fp_tmp,"%i %i %i %i %lf %lf %lf %lf\n",&first, &second, 
			&third, &fourth, &coul_tmp_r, &coul_tmp_i, &exch_tmp_r, &exch_tmp_i);
			cout<<first<<"  "<<second<<"   "<<third<<"   "<<fourth<<"   "<<endl;
			cout<<coul_tmp_r<<"   "<<coul_tmp_i<<"   "<<exch_tmp_r<<"    "<<exch_tmp_i<<endl<<endl;
			for(int k=0;k<numStates;k++){
				if((first_state[k]==first)&&(second_state[k]==second))
					row_index=k;
				 if((first_state[k]==third)&&(second_state[k]==fourth))
                                        col_index=k;
			}
			cout<<"Row  = "<<row_index<<"   Col = "<<col_index<<endl;
			if((row_index!=-1)||(col_index!=-1)){
				coul_matrix[row_index*numStates+col_index].r=coul_tmp_r*1000;
				coul_matrix[row_index*numStates+col_index].i=coul_tmp_i*1000;
				exch_matrix[row_index*numStates+col_index].r=exch_tmp_r*1000;
				exch_matrix[row_index*numStates+col_index].i=exch_tmp_i*1000;
			}					
		}//end of count		
	}//end of if
	fclose(fp_tmp);
}//end of steps

//Generate the total Many-body Hamiltonian
//reverse sign of off-diagonal elements for FORTRAN
nml_dcvector *H=nml_dcv_new(numStates*numStates);
for(int i=0;i<numStates;i++)
        for(int j=i;j<numStates;j++){
                if(j>i){
                        H[i*numStates+j].r=+(coul_matrix[i*numStates+j].r+exch_matrix[i*numStates+j].r);
                        H[i*numStates+j].i=+(coul_matrix[i*numStates+j].i+exch_matrix[i*numStates+j].i);
                        H[j*numStates+i].r=+(coul_matrix[i*numStates+j].r+exch_matrix[i*numStates+j].r);
                        H[j*numStates+i].i=-(coul_matrix[i*numStates+j].i+exch_matrix[i*numStates+j].i);
                }else{
                        H[i*numStates+i].r=+(coul_matrix[i*numStates+i].r+exch_matrix[i*numStates+i].r);
                        H[i*numStates+i].i=0.0;
                        H[i*numStates+i].r+=Ek[first_state[i]-1]+Ek[second_state[i]-1];
                }
}
cout<<"Hartree Matrix"<<endl;
FILE* Hartree=fopen("Hartree_Matrix","w");
for(int i=0;i<numStates;i++)
        for(int j=0;j<numStates;j++){
                cout<<i<<"   "<<j<<"   :   "<<H[i*numStates+j].r<<"    "<<H[i*numStates+j].i<<endl;
		fprintf(Hartree,"%i %i %12.10le %12.10le\n", i, j, H[i*numStates+j].r, H[i*numStates+j].i);
}	
fclose(Hartree);
cout<<"Starting Eigenvalue Calculations"<<endl;
//test solver first
int num_Basis=numStates;

  char JOBZ='V';
  char UPLO='U';
  int N=num_Basis;
//  nml_dcvector *A=nml_dcv_new(N*N);


//test purposes
/*   A[0].r=1;
   A[0].i=0;
   A[1].r=0;
   A[1].i=0;
   A[2].r=0;
   A[2].i=0;

   A[3].r=0;
   A[3].i=0;
   A[4].r=0;
   A[4].i=0;
   A[5].r=1;
   A[5].i=0;

   A[6].r=0;
   A[6].i=0;
   A[7].r=1;
   A[7].i=0;
   A[8].r=1;
   A[8].i=0;
*/
  int LDA = N;
  double* W=new double[N];
  int LWORK=2*N;
  nml_dcvector *WORK=nml_dcv_new(LWORK);
  double* RWORK=new double[3*N];
  int INFO=0;
  zheev_( &JOBZ, &UPLO, &N, H, &LDA, W, WORK, &LWORK, RWORK, &INFO);

  int num_eigs=num_Basis;
  if(selected!=0)
	num_eigs=selected;
  cout<<"Eigenvalues : "<<endl;
  for(int i=0;i<num_eigs;i++)
        cout<<W[i]<<endl;

  cout<<"Eigenvalues : "<<endl;
  FILE *CI_eigs=fopen("CI_Eigenvalues","w");
  for(int i=0;i<num_eigs;i++){
        printf("%lf\n", W[i]);
	fprintf(CI_eigs, "%12.10le\n", W[i]);	
  }
  fclose(CI_eigs);

  FILE *CI_evecs=fopen("CI_Evecs","w");
  for(int i=0;i<num_eigs;i++){
        cout<<endl;
        cout<<"Eigenvector  "<<i<<endl;
	fprintf(CI_evecs,"Eigenvector %i\n",i);
        for(int j=0;j<num_Basis;j++){
                cout<<H[i*num_Basis+j].r<<" "<<H[i*num_Basis+j].i<<endl;
		fprintf(CI_evecs,"%12.10le %12.10le\n", H[i*num_Basis+j].r, H[i*num_Basis+j].i);
        }
	fprintf(CI_evecs,"\n");	
        }
   fclose(CI_evecs);
	
//Compute probability densities of each Slater determinant contribution for each state.

  double p_real=0.0;
  double sum=0.0;
  cout<<"Probability: "<<endl<<endl;
  FILE* CI_wfs=fopen("CI_wfs","w");	
  FILE* fp_gs=fopen("CI_GS","w");
  fprintf(fp_gs,"%d\n",num_Basis);	

  for(int i=0;i<num_eigs;i++){
        cout<<endl;
        cout<<"Eigenvector  "<<i<<endl;
	fprintf(CI_wfs, "Wavefunction %i\n",i);
        sum=0;
        for(int j=0;j<num_Basis;j++){
                p_real=H[i*num_Basis+j].r*H[i*num_Basis+j].r+H[i*num_Basis+j].i*H[i*num_Basis+j].i;
                sum+=p_real;
                cout<<p_real<<endl;
		fprintf(CI_wfs, "%12.10le\n", p_real);
		if(i==0)			
			fprintf(fp_gs, "%12.10le\n", p_real);
        }
	fprintf(CI_wfs,"\n");
        //cout<<"Sum : "<<sum<<endl;
  }

fclose(fp_gs);
fclose(CI_wfs);

  double singlet_triplet_diff=W[1]-W[0];
  cout<<"Exchange Energy (2-1 difference)  :   "<<singlet_triplet_diff<<"   (meV)"<<endl;
}// end of rank
}

if(Compute_CI_Matrix_Selected) {
if(rank==0){
	//Compute Many body Eigen States  
	int numStates=(basis*(basis-1))/2;
    	int num_Matrix_Elements=(numStates*numStates-numStates)/2+numStates;
	nml_dcvector *coul_matrix=nml_dcv_new(numStates*numStates);
	nml_dcvector *exch_matrix=nml_dcv_new(numStates*numStates);
	double *Ek=new double[basis];
	cout<<"Basis : "<<basis<<endl;
	cout<<"Number of Matrix Elements : "<<num_Matrix_Elements<<endl;
	cout<<"Steps : "<<steps<<endl;

  	int *first_state=new int[numStates];
        int *second_state=new int[numStates];
        int count=0;
        for(int i=1;i<basis;i++)
        for(int j=i+1;j<=basis;j++){
                first_state[count]=i;
                second_state[count]=j;
                count++;
        }
	// read in eigenvalues
   	FILE* fp;
   	char filename[200];

   	sprintf(filename, "Ek");
   	fp=fopen(filename,"r");
   	cout << "Reading eigenvalues:  " << filename << endl;
   	if (!fp) die("Could not find file %s\n", filename);
   	for(int i=0;i<basis;i++)
        	fscanf(fp, "%lf\n", &Ek[i]);
   	fclose(fp);

	//print Eigenvalues
	cout<<"Single electron eigenvalues   "<<endl;
	for(int i=0;i<basis;i++)
        	cout<<Ek[i]<<endl;

//Read Coulomb and Exchange Matrix Elements Here 
for(int i=0;i<steps;i++){
	char filename[200];
        sprintf(filename, "steps_%i", i);
	cout<<"Reading file : "<<filename<<endl;
	double coul_tmp_r=0.0;
	double coul_tmp_i=0.0;
	double exch_tmp_r=0.0;
	double exch_tmp_i=0.0;
	int first = 0;
	int second = 0;
	int third = 0;
	int fourth = 0;
	int count = 0;
	int row_index=0;
	int col_index=0;
        FILE* fp_tmp = fopen(filename, "rb");
    	if(fp_tmp==NULL)
		cout<<"Error Opening File : "<<filename<<endl;
    	else{
		fscanf(fp_tmp,"%i\n", &count);
		cout<<"Step : "<<i<<"      Elements :  "<<count<<endl<<endl;
		for(int j=0;j<count;j++){
			first, second, third, fourth=0;
			coul_tmp_r, coul_tmp_i, exch_tmp_r, exch_tmp_i = 0.0;	
			row_index, col_index=-1;
			fscanf(fp_tmp,"%i %i %i %i %lf %lf %lf %lf\n",&first, &second, 
			&third, &fourth, &coul_tmp_r, &coul_tmp_i, &exch_tmp_r, &exch_tmp_i);
			cout<<first<<"  "<<second<<"   "<<third<<"   "<<fourth<<"   "<<endl;
			cout<<coul_tmp_r<<"   "<<coul_tmp_i<<"   "<<exch_tmp_r<<"    "<<exch_tmp_i<<endl<<endl;
			for(int k=0;k<numStates;k++){
				if((first_state[k]==first)&&(second_state[k]==second))
					row_index=k;
				 if((first_state[k]==third)&&(second_state[k]==fourth))
                                        col_index=k;
			}
			cout<<"Row  = "<<row_index<<"   Col = "<<col_index<<endl;
			if((row_index!=-1)||(col_index!=-1)){
				coul_matrix[row_index*numStates+col_index].r=coul_tmp_r;
				coul_matrix[row_index*numStates+col_index].i=coul_tmp_i;
				exch_matrix[row_index*numStates+col_index].r=exch_tmp_r;
				exch_matrix[row_index*numStates+col_index].i=exch_tmp_i;
			}					
		}//end of count		
	}//end of if
	fclose(fp_tmp);
}//end of steps

//Generate the total Many-body Hamiltonian
//reverse sign of off-diagonal elements for FORTRAN
nml_dcvector *H=nml_dcv_new(numStates*numStates);
for(int i=0;i<numStates;i++)
        for(int j=i;j<numStates;j++){
                if(j>i){
                        H[i*numStates+j].r=+(coul_matrix[i*numStates+j].r+exch_matrix[i*numStates+j].r);
                        H[i*numStates+j].i=+(coul_matrix[i*numStates+j].i+exch_matrix[i*numStates+j].i);
                        H[j*numStates+i].r=+(coul_matrix[i*numStates+j].r+exch_matrix[i*numStates+j].r);
                        H[j*numStates+i].i=-(coul_matrix[i*numStates+j].i+exch_matrix[i*numStates+j].i);
                }else{
                        H[i*numStates+i].r=+(coul_matrix[i*numStates+i].r+exch_matrix[i*numStates+i].r);
                        H[i*numStates+i].i=0.0;
                        H[i*numStates+i].r+=Ek[first_state[i]-1]+Ek[second_state[i]-1];
                }
}
cout<<"Hartree Matrix"<<endl;
FILE* Hartree=fopen("Hartree_Matrix","w");
for(int i=0;i<numStates;i++)
        for(int j=0;j<numStates;j++){
                cout<<i<<"   "<<j<<"   :   "<<H[i*numStates+j].r<<"    "<<H[i*numStates+j].i<<endl;
		fprintf(Hartree,"%i %i %12.10le %12.10le\n", i, j, H[i*numStates+j].r, H[i*numStates+j].i);
}	
fclose(Hartree);
cout<<"Starting Eigenvalue Calculations"<<endl;
//test solver first
int num_Basis=numStates;

  char JOBZ='V';
  char UPLO='U';
  int N=num_Basis;
  char RANGE='I';

//  nml_dcvector *A=nml_dcv_new(N*N);


//test purposes
/*   A[0].r=1;
   A[0].i=0;
   A[1].r=0;
   A[1].i=0;
   A[2].r=0;
   A[2].i=0;

   A[3].r=0;
   A[3].i=0;
   A[4].r=0;
   A[4].i=0;
   A[5].r=1;
   A[5].i=0;

   A[6].r=0;
   A[6].i=0;
   A[7].r=1;
   A[7].i=0;
   A[8].r=1;
   A[8].i=0;
*/

  int LDA = N;
  double VL=0.0;
  double VU=3.0;
  int IL=1;
  int IU=selected;
  double ABSTOL=1e-8;
  int M=0;
  double* W=new double[N];
  int LDZ=N;
  nml_dcvector *Z=nml_dcv_new(LDZ*N); 
  int LWORK=2*N;
  nml_dcvector *WORK=nml_dcv_new(LWORK);
  int INFO=0;
  int *ISUPPZ=new int[2*N];  
  int LRWORK=24*N;
  double* RWORK=new double[LRWORK];	
  int LIWORK=10*N;
  int *IWORK=new int[LIWORK];
  zheevr_(&JOBZ, &RANGE, &UPLO, &N, H, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, 
		ISUPPZ, WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);
  cout<<"Number of Eigenvalues Found : "<<M<<endl;
  cout<<"Lapack Status : "<<INFO<<endl;
	
  int num_eigs=M;

  cout<<"Eigenvalues : "<<endl;
  for(int i=0;i<num_eigs;i++)
        cout<<W[i]<<endl;

  cout<<"Eigenvalues : "<<endl;
  FILE *CI_eigs=fopen("CI_Eigenvalues","w");
  for(int i=0;i<num_eigs;i++){
        printf("%lf\n", W[i]);
	fprintf(CI_eigs, "%12.10le\n", W[i]);	
  }
  fclose(CI_eigs);

  FILE *CI_evecs=fopen("CI_Evecs","w");
  for(int i=0;i<num_eigs;i++){
        cout<<endl;
        cout<<"Eigenvector  "<<i<<endl;
	fprintf(CI_evecs,"Eigenvector %i\n",i);
        for(int j=0;j<num_Basis;j++){
                cout<<Z[i*num_Basis+j].r<<" "<<Z[i*num_Basis+j].i<<endl;
		fprintf(CI_evecs,"%12.10le %12.10le\n", Z[i*num_Basis+j].r, Z[i*num_Basis+j].i);
        }
	fprintf(CI_evecs,"\n");	
        }
   fclose(CI_evecs);
	
//Compute probability densities of each Slater determinant contribution for each state.

  double p_real=0.0;
  double sum=0.0;
  cout<<"Probability: "<<endl<<endl;
  FILE* CI_wfs=fopen("CI_wfs","w");	
  for(int i=0;i<num_eigs;i++){
        cout<<endl;
        cout<<"Eigenvector  "<<i<<endl;
	fprintf(CI_wfs, "Wavefunction %i\n",i);
        sum=0;
        for(int j=0;j<num_Basis;j++){
                p_real=Z[i*num_Basis+j].r*Z[i*num_Basis+j].r+Z[i*num_Basis+j].i*Z[i*num_Basis+j].i;
                sum+=p_real;
                cout<<p_real<<endl;
		fprintf(CI_wfs, "%12.10le\n", p_real);
        }
	fprintf(CI_wfs,"\n");
        //cout<<"Sum : "<<sum<<endl;
  }
  double singlet_triplet_diff=(W[1]-W[0])*1000;
  cout<<"Exchange Energy (2-1 difference)  :   "<<singlet_triplet_diff<<"   (meV)"<<endl;
}
}

if(Compute_Kramer_vectors){
cout<<"Generating Kramer's degenerate spin states"<<endl;

}

 if(ComputeCoulomb) {

    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb =
          Coulomb(NumAtoms, x1_start, x1_end, x2_start, x2_end, NumOrbitals);
    ifstream coulomb_input("coulomb.table");
    data_coulomb.assign_table(coulomb_input);
       //data_coulomb.assign_table1();

    //data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
  
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

    //complex<double> S=0.0;
    //complex<double> S_sqr=0.0;
    //if(x1_end != 0) data_coulomb.dotProduct(1, 2, &S);
    //S_sqr=conj(S)*S;

    double elapsed_time=0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time=-MPI_Wtime();
    if(use_dd_cutoff==0)
        data_coulomb.compute_exchange(1,2,1,2,rank);
    else {

       if(x1_end != 0) {
          //data_coulomb.compute_exchange_cutoff(1,2,1,2,rank, dd_cutoff1);
          if(NumElectrons == 2) data_coulomb.compute_coulomb_exchange2(dd_cutoff1,data, eList, NumElectrons);
          else if(NumElectrons == 4) data_coulomb.compute_coulomb_exchange4(dd_cutoff1,data,eList,NumElectrons);
          else {
             cout <<"Choose 2 or 4 as NumElectrons"<<endl;
             exit(1);
          }
       }
    }

    elapsed_time+=MPI_Wtime();
    complex<double>* local_coul = data_coulomb.get_H_coul();
    complex<double>* local_exch = data_coulomb.get_H_exch();
    //complex<double>* local_S = &S;

    int csize = 1;
    complex<double>* total_coul;
    total_coul = new complex<double>[csize];
    complex<double>* total_exch;
    total_exch = new complex<double>[csize];
    //complex<double>* total_S;
    //total_S = new complex<double>[csize];

#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Reduce(local_coul, total_coul, csize*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_exch, total_exch, csize*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //MPI_Reduce(local_S, total_S, csize*2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for(int cs=0; cs<csize; cs++) {
       total_coul[cs] = local_coul[cs];
       total_exch[cs] = local_exch[cs];
       //total_S[cs] = local_S[cs];
    }

#endif

    if (rank==0) {

       //S_sqr = conj(total_S[0])*total_S[0];

       cout<<endl<<"====================================================="<<endl;
       cout<<"Electron and Hole Pair: "<<endl;
       cout<<"Coulomb:  "<<setw(10)<<total_coul[0]<<endl;
       cout<<"Exchange: "<<setw(10)<<total_exch[0]<<endl;
       //cout<<"Overlap_square:  "<<S_sqr<<endl;
       cout<<endl;
       cout<<"Time to compute a matrix element:"<<elapsed_time<<endl;
        
	char filename1[200];
	if(NumElectrons==2)
        	sprintf(filename1, "Coulomb_%d_%d_%d_%d", eList[0], eList[1], eList[0], eList[1]);
	else
        	sprintf(filename1, "Coulomb_%d_%d_%d_%d", eList[0], eList[1], eList[2], eList[3]);

        FILE* fp1=fopen(filename1,"w");
  	if(fp1!=NULL){
        	fprintf(fp1, "%12.10le %12.10le\n", total_coul[0].real(), total_coul[0].imag());
                fprintf(fp1, "%12.10le %12.10le\n", total_exch[0].real(), total_exch[0].imag());
        }
	fclose(fp1);
   }
}

if(ComputeCoulomb_new) {

    //Computing Coulomb-Exchange energies ===============================
    	Coulomb data_coulomb =
        	Coulomb(NumAtoms, x1_start, x1_end, x2_start, x2_end, NumOrbitals);
    	ifstream coulomb_input("coulomb.table");
    	data_coulomb.assign_table(coulomb_input);
       	//data_coulomb.assign_table1();

    	//data_coulomb.get_data(wf_e, atomid, neighbor, lattice);
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

    double elapsed_time=0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time=-MPI_Wtime();
    complex<double>* wf_e1;

     if(x1_end != 0) {
          //data_coulomb.compute_exchange_cutoff(1,2,1,2,rank, dd_cutoff1);
          if(NumElectrons == 2){
		read_wf(data, 2, eList, &wf_e1, x1_start, x1_end, x2_start, x2_end);
		data_coulomb.get_data(wf_e1, atomid, neighbor, lattice);
		data_coulomb.compute_coulomb_exchange2(dd_cutoff1,data, eList, NumElectrons);
	  }
          else if(NumElectrons == 4){
		read_wf(data, 4, eList, &wf_e1, x1_start, x1_end, x2_start, x2_end);
		data_coulomb.get_data(wf_e1, atomid, neighbor, lattice);
		data_coulomb.compute_coulomb_exchange4(dd_cutoff1,data,eList,NumElectrons);
	  }	
          else {
             cout <<"Choose 2 or 4 as NumElectrons"<<endl;
             exit(1);
          }
       }
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

       char filename1[200];
       if(NumElectrons==2)
                sprintf(filename1, "Coulomb_%d_%d_%d_%d", eList[0], eList[1], eList[0], eList[1]);
       else
                sprintf(filename1, "Coulomb_%d_%d_%d_%d", eList[0], eList[1], eList[2], eList[3]);
       FILE* fp1=fopen(filename1,"w");
       if(fp1!=NULL){
                fprintf(fp1, "%12.10le %12.10le\n", total_coul[0].real(), total_coul[0].imag());
                fprintf(fp1, "%12.10le %12.10le\n", total_exch[0].real(), total_exch[0].imag());
       }
       fclose(fp1);
   }
}

//Identify spin polarization

if(ComputeSpin){
Coulomb data_coulomb =
        Coulomb(NumAtoms, x1_start, x1_end, x2_start, x2_end, NumOrbitals);

int *electron_states=new int[1];
FILE *fp=fopen("Spin_polarization","w");
if(fp==NULL){
	cout<<"Error opening file "<<endl;
	exit(-1);	
}	

for(int i=1;i<=basis;i++){
    	complex<double>* wf_e1;
        electron_states[0]=i;	
        read_wf(data, 1, electron_states, &wf_e1, x1_start, x1_end, x2_start, x2_end);
        data_coulomb.get_data(wf_e1, atomid, neighbor, lattice);
	complex<double> sum_up= 0.0;
	complex<double> sum_down=0.0; 
	data_coulomb.spinPolarization(1, &sum_up, &sum_down);
        complex<double>* total_up_spin=new complex<double>[1];
        complex<double>* total_down_spin= new complex<double>[1];

#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&sum_up, total_up_spin, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_down, total_down_spin, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    total_up_spin[0] = sum_up[0];
    total_down_spin[0] = sum_down[0];
#endif
    if (rank==0) {
       	cout<<endl;
	cout<<"Spin decomposition of wf "<<i<<"  :  "<<endl;
	cout<<"Up spin    :   "<<total_up_spin[0]<<endl;
	cout<<"Down spin  :   "<<total_down_spin[0]<<endl;
       	if(fp!=NULL){
        	fprintf(fp, "%i\n%12.10le %12.10le\n%12.10le %12.10le\n", i, 
		total_up_spin[0].real(), total_up_spin[0].imag(), total_down_spin[0].real(), total_down_spin[0].imag());
	}
    }//end of rank
    delete(wf_e1);
}//end of basis
fclose(fp);
}	

//Heitler-London formalism with spin
//Uses non-orthogonal wave functions
//Make sure: evec_1 is L_up, evec_2 is L_down, evec_3 is R_up, evec_4 is R_down

if(ComputeHL20) {
Coulomb data_coulomb =
         Coulomb(NumAtoms, x1_start, x1_end, x2_start, x2_end, NumOrbitals);
    ifstream coulomb_input("coulomb.table");
    data_coulomb.assign_table(coulomb_input);

        FILE* fp=fopen("exchange_parameters","r");
        float dd_cutoff=12.0;
        double dd_cutoff1=dd_cutoff;
        if(fp!=NULL){
                fscanf(fp, "%f\n", &dd_cutoff);
                fclose(fp);
                dd_cutoff1=dd_cutoff;
        }else{
                if(rank==0)
                        cout<<"File not found: Using default parameters. dd_cutoff=12 nm"<<endl;
        }
        if(rank==0)
                cout<<"dd_cutoff: "<<dd_cutoff1<<"  nm"<<endl;

        double *coordinates1=new double[3];
        double *coordinates2=new double[3];

        FILE* fp_impurity=fopen("Impurity_Parameters","r");
        if(fp_impurity==NULL){
                if(rank==0)
                        cout<<"File not found : Impurity_Parameters"<<endl;
                exit(1);
        }else{
                fscanf(fp_impurity,"%lf %lf %lf\n", &coordinates1[0], &coordinates1[1], &coordinates1[2]);
                fscanf(fp_impurity,"%lf %lf %lf\n", &coordinates2[0], &coordinates2[1], &coordinates2[2]);

                cout<<"Impurity 1 is at : ("<<coordinates1[0]<<", "<<coordinates1[1]<<", "<<coordinates1[2]<<")"<<endl;
                cout<<"Impurity 2 is at : ("<<coordinates2[0]<<", "<<coordinates2[1]<<", "<<coordinates2[2]<<")"<<endl;
        }
        fclose(fp_impurity);

        int p1=-1;
	int p2=-1;
        double x1, y1, z1, x2, y2, z2 = 0.0;
        double R1, R2 = 0.0;
        double R1_max=9999999.0;
	double R2_max=9999999.0;
        for(int i=0;i<NumAtoms;i++){
                x1=lattice[3*i]-coordinates1[0];
                y1=lattice[3*i+1]-coordinates1[1];
                z1=lattice[3*i+2]-coordinates1[2];
                x2=lattice[3*i]-coordinates2[0];
                y2=lattice[3*i+1]-coordinates2[1];
                z2=lattice[3*i+2]-coordinates2[2];
                R1=x1*x1+y1*y1+z1*z1;
                R2=x2*x2+y2*y2+z2*z2;
                if(R1<R1_max){
                        R1_max=R1;
                        p1=i;
                }
                if(R2<R2_max){
                        R2_max=R2;
                        p2=i;
                }
        }
	 if((p1==-1)||(p2==-1)){
                cout<<"Impurity not found in the lattice"<<endl;
                exit(1);
        }else{

        	if(rank==0){
                	cout<<"Impurity 1  : Index="<<p1<<endl;
                	cout<<"Coordinates : "<<lattice[3*p1]<<"  "<<lattice[3*p1+1]<<"  "<<lattice[3*p1+2]<<endl;
           		cout<<"Impurity 2  : Index="<<p2<<endl;
                	cout<<"Coordinates : "<<lattice[3*p2]<<"  "<<lattice[3*p2+1]<<"  "<<lattice[3*p2+2]<<endl;
        	}
       } 

        double* Ek=new double[4];
        FILE* fp_Ek_L=fopen("Ek_L","r");
        FILE* fp_Ek_R=fopen("Ek_R","r");
        if((fp_Ek_L==NULL)||(fp_Ek_R==NULL)){
                cout<<"Error in file opening : Ek_L, Ek_R"<<endl;
                exit(1);
        }else{
                fscanf(fp_Ek_L, "%lf\n%lf\n", &Ek[0], &Ek[1]);
                fscanf(fp_Ek_R, "%lf\n%lf\n", &Ek[2], &Ek[3]);
        }
        if(rank==0){
                cout<<"Ek  :  "<<endl;
                for(int i=0;i<4;i++)
                        cout<<Ek[i]<<endl;
        }

    int basis = 4;
    int numStates=2;

    int *first_state =new int[numStates];
    int *second_state=new int[numStates];
    first_state[0]=1;
    first_state[1]=2;
    second_state[0]=4;
    second_state[1]=3;		

    complex<double>* Overlap_Matrix = new complex<double>[numStates*numStates];
    complex<double>* V_imp = new complex<double>[numStates*numStates];
    complex<double>* two_e = new complex<double>[numStates*numStates];
    complex<double>* one_e = new complex<double>[numStates*numStates];

    int* electron_states=new int[4];
    FILE* fp_overlap=fopen("Overlap_Matrix_HL","w");
    FILE* fp_V_imp=fopen("V_imp_HL","w");
    FILE* fp_two_e=fopen("Two_e_HL","w");
    FILE* fp_one_e=fopen("One_e_HL","w");
    if((fp_overlap==NULL)||(fp_V_imp==NULL)||(fp_two_e==NULL)||(fp_one_e==NULL)){
                cout<<"Error opening file"<<endl;
                exit(1);
    }
for(int i=0;i<numStates;i++){
        for(int j=i;j<numStates;j++){
                complex<double>* wf_e1;
                electron_states[0]=first_state[i];
                electron_states[1]=second_state[i];
                electron_states[2]=first_state[j];
                electron_states[3]=second_state[j];
                if(rank==0){
                        cout<<"Calculating matrix element between |"<<electron_states[0]<<" "<<electron_states[1]<<"> and |"<<
                        electron_states[2]<<" "<<electron_states[3]<<">"<<endl;
                }
                read_wf(data, 4, electron_states, &wf_e1, x1_start, x1_end, x2_start, x2_end);
                data_coulomb.get_data(wf_e1, atomid, neighbor, lattice);
                complex<double> mp=0.0;
                complex<double> np=0.0;
                complex<double> mq=0.0;
                complex<double> nq=0.0;
                complex<double> mp_total=0.0;
                complex<double> np_total=0.0;
                complex<double> mq_total=0.0;
                complex<double> nq_total=0.0;
                complex<double> mVp=0.0;
                complex<double> nVp=0.0;
                complex<double> mVq=0.0;
                complex<double> nVq=0.0;
                complex<double> mVp_total=0.0;
                complex<double> nVp_total=0.0;
                complex<double> mVq_total=0.0;
                complex<double> nVq_total=0.0;

                data_coulomb.dotProduct_Slater(&mp, &nq, &np, &mq);
                data_coulomb.coulomb_Slater(p1, p2, &mVp, &nVq, &nVp, &mVq);
                data_coulomb.compute_coulomb_exchange4(dd_cutoff1,data,electron_states,4);

                complex<double>* local_coul = data_coulomb.get_H_coul();
                complex<double>* local_exch = data_coulomb.get_H_exch();
                complex<double> total_coul=0.0;
                complex<double> total_exch=0.0;

#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&mp, &mp_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nq, &nq_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&np, &np_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mq, &mq_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mVp, &mVp_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nVq, &nVq_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nVp, &nVp_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mVq, &mVq_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_coul, &total_coul, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_exch, &total_exch, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    mp_total = mp;
    nq_total = nq;
    np_total = np;
    mq_total = mq;
    mVp_total = mVp;
    nVq_total = nVq;
    nVp_total = nVp;
    mVq_total = mVq;
    total_coul = local_coul[0];
    total_exch = local_exch[0];
#endif
if(rank==0){

//Overlap Matrix Element
Overlap_Matrix[i*numStates+j]=2.0*(mp_total*nq_total-mq_total*np_total);
fprintf(fp_overlap, "%i %i %i %i %12.10le %12.10le\n", electron_states[0], electron_states[1], electron_states[2], electron_states[3],
Overlap_Matrix[i*numStates+j].real(), Overlap_Matrix[i*numStates+j].imag());
//debug
cout<<"Overlap: S"<<electron_states[0]<<electron_states[1]<<" and S"<<electron_states[2]<<electron_states[3]<<"  :  "
<<2.0*(mp_total*nq_total-mq_total*np_total)<<endl;

//Impurity potential Expectation
V_imp[i*numStates+j]=mVp_total*nq_total-mVq_total*np_total-nVp_total*mq_total+nVq_total*mp_total;
fprintf(fp_V_imp, "%i %i %i %i %12.10le %12.10le\n", electron_states[0], electron_states[1], electron_states[2], electron_states[3],
V_imp[i*numStates+j].real(), V_imp[i*numStates+j].imag());
//debug
cout<<"V_imp: S"<<electron_states[0]<<electron_states[1]<<" and S"<<electron_states[2]<<electron_states[3]<<"  :  "
<<mVp_total*nq_total-mVq_total*np_total-nVp_total*mq_total+nVq_total*mp_total<<endl;

//Two electron values
two_e[i*numStates+j]=total_coul+total_exch;
fprintf(fp_two_e, "%i %i %i %i %12.10le %12.10le %12.10le %12.10le\n", electron_states[0], electron_states[1],
electron_states[2], electron_states[3], total_coul.real(), total_coul.imag(), total_exch.real(), total_exch.imag());
//debug
cout<<"Two_e: S"<<electron_states[0]<<electron_states[1]<<" and S"<<electron_states[2]<<electron_states[3]<<"  :  "
<<total_coul+total_exch<<endl;

//Single electron values
one_e[i*numStates+j]=2.0*(Ek[first_state[j]-1]+Ek[second_state[j]-1])*(mp_total*nq_total-np_total*mq_total);
fprintf(fp_one_e, "%i %i %i %i %12.10le %12.10le\n", electron_states[0], electron_states[1], electron_states[2], electron_states[3],
one_e[i*numStates+j].real(), one_e[i*numStates+j].imag());
//debug
cout<<"One_electron : h1 + h2"<<electron_states[0]<<electron_states[1]<<" and S"<<electron_states[2]<<electron_states[3]<<"  :  "
<<2.0*(Ek[first_state[j]-1]+Ek[second_state[j]-1])*(mp_total*nq_total-np_total*mq_total)<<endl;

}
delete(wf_e1);
}//end of for
}//end of for
fclose(fp_overlap);
fclose(fp_V_imp);
fclose(fp_two_e);
fclose(fp_one_e);

if(rank==0){
//2x2 HL matrix
nml_dcvector *H=nml_dcv_new(numStates*numStates);
nml_dcvector *O=nml_dcv_new(numStates*numStates);

FILE* Hartree=fopen("Hartree_Matrix_HL","w");
if(Hartree==NULL)
        cout<<"Error opening Hartree File"<<endl;
        for(int i=0;i<numStates;i++)
          for(int j=i;j<numStates;j++){
                        H[i*numStates+j].r=one_e[i*numStates+j].real()+V_imp[i*numStates+j].real()+two_e[i*numStates+j].real();
                        H[i*numStates+j].i=one_e[i*numStates+j].imag()+V_imp[i*numStates+j].imag()+two_e[i*numStates+j].imag();
			O[i*numStates+j].r=Overlap_Matrix[i*numStates+j].real();
			O[i*numStates+j].i=Overlap_Matrix[i*numStates+j].imag();
                if(i!=j){
                        H[j*numStates+i].r=one_e[i*numStates+j].real()+V_imp[i*numStates+j].real()+two_e[i*numStates+j].real();
                        H[j*numStates+i].i=-(one_e[i*numStates+j].imag()+V_imp[i*numStates+j].imag()+two_e[i*numStates+j].imag());
			O[j*numStates+i].r=Overlap_Matrix[i*numStates+j].real();
                        O[j*numStates+i].i=-Overlap_Matrix[i*numStates+j].imag();
                }
		/*if(i==j){
			H[i*numStates+j].i=0.0;
			O[i*numStates+j].i=0.0;
		}
		*/
                if(Hartree!=NULL)
                        fprintf(Hartree,"%i %i %12.10le %12.10le\n", i, j, H[i*numStates+j].r, H[i*numStates+j].i);
          }
fclose(Hartree);
cout<<"Starting Eigenvalue Calculations"<<endl;
int ITYPE=1; //Ax=EBx
char JOBZ='V';
char UPLO='U';
int N=numStates;
int LDA=N;
int LDB=N;
double *W=new double[N];
int LWORK=2*N-1;
int INFO=0;
nml_dcvector *WORK = nml_dcv_new(LWORK);
double *RWORK = new double[3*N-2];
nml_dcvector *A=nml_dcv_new(numStates*numStates);
nml_dcvector *B=nml_dcv_new(numStates*numStates);
double tolerance = 1.0e-12;
for(int i=0;i<numStates;i++)
	for(int j=0;j<numStates;j++){
		A[i*numStates+j].r=H[j*numStates+i].r;
		A[i*numStates+j].i=H[j*numStates+i].i;
		B[i*numStates+j].r=O[j*numStates+i].r;
		B[i*numStates+j].i=O[j*numStates+i].i;		
		if(abs(A[i*numStates+j].r) < tolerance)
			A[i*numStates+j].r=0.0;
		if(abs(A[i*numStates+j].i) < tolerance)
                        A[i*numStates+j].i=0.0;
		if(abs(B[i*numStates+j].r) < tolerance)
                        B[i*numStates+j].r=0.0;
		if(abs(B[i*numStates+j].i) < tolerance)
                        B[i*numStates+j].i=0.0;
	}

for(int i=0;i<numStates;i++){
        for(int j=0;j<numStates;j++){
		cout<<"row = "<<i<<"     col = "<<j<<endl;
		cout<<"A =   "<<A[i*numStates+j].r<<"   "<<A[i*numStates+j].i<<endl<<endl;
		cout<<"B =   "<<B[i*numStates+j].r<<"   "<<B[i*numStates+j].i<<endl<<endl;
        }
}

zhegv_(&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, RWORK, &INFO);
cout<<"Eigenvalues : "<<endl;
for(int i=0;i<numStates;i++)
	cout<<W[i]<<endl;
cout<<"Singlet-Triplet Splitting (meV) : "<<(W[1]-W[0])*1000<<endl;
for(int i=0;i<numStates;i++){
	cout<<"Eigenvector "<<i<<"  :  "<<endl;	
        for(int j=0;j<numStates;j++)
		cout<<A[j*numStates+i].r<<"  "<<A[j*numStates+i].i<<endl;
}
double p_real=0.0;
cout<<"Computing Probablities : "<<endl;
for(int i=0;i<numStates;i++){
        cout<<"WF "<<i<<"  :  "<<endl;
        for(int j=0;j<numStates;j++){
		p_real=A[j*numStates+i].r*A[j*numStates+i].r+A[j*numStates+i].i*A[j*numStates+i].i;
		cout<<p_real<<endl;
	}	
}
//Triplet state : S14+S23
//Singlet state : S14-S23
/*nml_dcvector *HL_matrix=nml_dcv_new(numStates*numStates);
nml_dcvector *OL_matrix=nml_dcv_new(numStates*numStates);
HL_matrix[0].r=H[0].r+H[3].r-H[1];

cout<<"HL_matrix :  "<<HL_matrix[0].r<<"  "<<HL_matrix[0].i<<endl;
*/

}//end of rank
}


if(Compute_hund_mulliken){
Coulomb data_coulomb =
         Coulomb(NumAtoms, x1_start, x1_end, x2_start, x2_end, NumOrbitals);
    ifstream coulomb_input("coulomb.table");
    data_coulomb.assign_table(coulomb_input);
 
 	FILE* fp=fopen("exchange_parameters","r");
        float dd_cutoff=12.0;
        double dd_cutoff1=dd_cutoff;
        if(fp!=NULL){
                fscanf(fp, "%f\n", &dd_cutoff);
                fclose(fp);
                dd_cutoff1=dd_cutoff;
        }else{
                if(rank==0)
                        cout<<"File not found: Using default parameters. dd_cutoff=12 nm"<<endl;
        }
        if(rank==0)
                cout<<"dd_cutoff: "<<dd_cutoff1<<"  nm"<<endl;
	
	double *coordinates1=new double[3];
	double *coordinates2=new double[3];

	FILE* fp_impurity=fopen("Impurity_Parameters","r");
	if(fp_impurity==NULL){
		if(rank==0)			
			cout<<"File not found : Impurity_Parameters"<<endl;
		exit(1);
	}else{
		fscanf(fp_impurity,"%lf %lf %lf\n", &coordinates1[0], &coordinates1[1], &coordinates1[2]);
		fscanf(fp_impurity,"%lf %lf %lf\n", &coordinates2[0], &coordinates2[1], &coordinates2[2]);

		cout<<"Impurity 1 is at : ("<<coordinates1[0]<<", "<<coordinates1[1]<<", "<<coordinates1[2]<<")"<<endl;
		cout<<"Impurity 2 is at : ("<<coordinates2[0]<<", "<<coordinates2[1]<<", "<<coordinates2[2]<<")"<<endl;
	}
	fclose(fp_impurity);


        int p1=-1;
        int p2=-1;
        double x1, y1, z1, x2, y2, z2 = 0.0;
        double R1, R2 = 0.0;
        double R1_max=9999999.0;
        double R2_max=9999999.0;
        for(int i=0;i<NumAtoms;i++){
                x1=lattice[3*i]-coordinates1[0];
                y1=lattice[3*i+1]-coordinates1[1];
                z1=lattice[3*i+2]-coordinates1[2];
                x2=lattice[3*i]-coordinates2[0];
                y2=lattice[3*i+1]-coordinates2[1];
                z2=lattice[3*i+2]-coordinates2[2];
                R1=x1*x1+y1*y1+z1*z1;
                R2=x2*x2+y2*y2+z2*z2;
                if(R1<R1_max){
                        R1_max=R1;
                        p1=i;
                }
                if(R2<R2_max){
                        R2_max=R2;
                        p2=i;
                }
        }
         if((p1==-1)||(p2==-1)){
                cout<<"Impurity not found in the lattice"<<endl;
                exit(1);
        }else{

                if(rank==0){
                        cout<<"Impurity 1  : Index="<<p1<<endl;
                        cout<<"Coordinates : "<<lattice[3*p1]<<"  "<<lattice[3*p1+1]<<"  "<<lattice[3*p1+2]<<endl;
                        cout<<"Impurity 2  : Index="<<p2<<endl;
                        cout<<"Coordinates : "<<lattice[3*p2]<<"  "<<lattice[3*p2+1]<<"  "<<lattice[3*p2+2]<<endl;
                }
       }
	
	double* Ek=new double[4];
	FILE* fp_Ek_L=fopen("Ek_L","r");
	FILE* fp_Ek_R=fopen("Ek_R","r");
	if((fp_Ek_L==NULL)||(fp_Ek_R==NULL)){
		cout<<"Error in file opening : Ek_L, Ek_R"<<endl;
		exit(1);
	}else{
		fscanf(fp_Ek_L, "%lf\n%lf\n", &Ek[0], &Ek[1]);
		fscanf(fp_Ek_R, "%lf\n%lf\n", &Ek[2], &Ek[3]);
	}
	if(rank==0){
		cout<<"Ek  :  "<<endl;
		for(int i=0;i<4;i++)
			cout<<Ek[i]<<endl;
	}

    int basis = 4; 	
    int numStates=(basis*(basis-1))/2;

    int *first_state=new int[numStates];
    int *second_state=new int[numStates];
    int count=0;
    for(int i=1;i<basis;i++)
    	for(int j=i+1;j<=basis;j++){
        	first_state[count]=i;
                second_state[count]=j;
                count++;
    }

    complex<double>* Overlap_Matrix = new complex<double>[numStates*numStates];
    complex<double>* V_imp = new complex<double>[numStates*numStates];
    complex<double>* two_e = new complex<double>[numStates*numStates];
    complex<double>* one_e = new complex<double>[numStates*numStates];
    	
    int* electron_states=new int[4];
    FILE* fp_overlap=fopen("Overlap_Matrix","w");
    FILE* fp_V_imp=fopen("V_imp","w");
    FILE* fp_two_e=fopen("Two_e","w");
    FILE* fp_one_e=fopen("One_e","w");
    if((fp_overlap==NULL)||(fp_V_imp==NULL)||(fp_two_e==NULL)||(fp_one_e==NULL)){			
		cout<<"Error opening file"<<endl;
		exit(1);
    } 		
for(int i=0;i<numStates;i++){
	for(int j=i;j<numStates;j++){
		complex<double>* wf_e1;
  		electron_states[0]=first_state[i];
        	electron_states[1]=second_state[i];
        	electron_states[2]=first_state[j];
        	electron_states[3]=second_state[j];
		if(rank==0){
			cout<<"Calculating matrix element between |"<<electron_states[0]<<" "<<electron_states[1]<<"> and |"<<
			electron_states[2]<<" "<<electron_states[3]<<">"<<endl;	
		}
        	read_wf(data, 4, electron_states, &wf_e1, x1_start, x1_end, x2_start, x2_end);
        	data_coulomb.get_data(wf_e1, atomid, neighbor, lattice);

		complex<double> mp=0.0;
		complex<double> np=0.0;
		complex<double> mq=0.0;
		complex<double> nq=0.0;
		complex<double> mp_total=0.0;
        	complex<double> np_total=0.0;
        	complex<double> mq_total=0.0;
        	complex<double> nq_total=0.0;
		complex<double> mVp=0.0;
        	complex<double> nVp=0.0;
        	complex<double> mVq=0.0;
        	complex<double> nVq=0.0;
        	complex<double> mVp_total=0.0;
        	complex<double> nVp_total=0.0;
        	complex<double> mVq_total=0.0;
        	complex<double> nVq_total=0.0;

		data_coulomb.dotProduct_Slater(&mp, &nq, &np, &mq);
		data_coulomb.coulomb_Slater(p1, p2, &mVp, &nVq, &nVp, &mVq);
        	data_coulomb.compute_coulomb_exchange4(dd_cutoff1,data,electron_states,4);

        	complex<double>* local_coul = data_coulomb.get_H_coul();
        	complex<double>* local_exch = data_coulomb.get_H_exch();
		complex<double> total_coul=0.0;
        	complex<double> total_exch=0.0;
	
#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&mp, &mp_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nq, &nq_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&np, &np_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mq, &mq_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mVp, &mVp_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nVq, &nVq_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nVp, &nVp_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&mVq, &mVq_total, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_coul, &total_coul, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_exch, &total_exch, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    mp_total = mp;
    nq_total = nq;
    np_total = np;
    mq_total = mq;
    mVp_total = mVp;
    nVq_total = nVq;
    nVp_total = nVp;
    mVq_total = mVq;
    total_coul = local_coul[0];
    total_exch = local_exch[0];
#endif
if(rank==0){

//Overlap Matrix Element
Overlap_Matrix[i*numStates+j]=2.0*(mp_total*nq_total-mq_total*np_total);
fprintf(fp_overlap, "%i %i %i %i %12.10le %12.10le\n", electron_states[0], electron_states[1], electron_states[2], electron_states[3], 
Overlap_Matrix[i*numStates+j].real(), Overlap_Matrix[i*numStates+j].imag());
//debug
cout<<"Overlap: S"<<electron_states[0]<<electron_states[1]<<" and S"<<electron_states[2]<<electron_states[3]<<"  :  "
<<2.0*(mp_total*nq_total-mq_total*np_total)<<endl;

//Impurity potential Expectation
V_imp[i*numStates+j]=mVp_total*nq_total-mVq_total*np_total-nVp_total*mq_total+nVq_total*mp_total;
fprintf(fp_V_imp, "%i %i %i %i %12.10le %12.10le\n", electron_states[0], electron_states[1], electron_states[2], electron_states[3], 
V_imp[i*numStates+j].real(), V_imp[i*numStates+j].imag());
//debug
cout<<"V_imp: S"<<electron_states[0]<<electron_states[1]<<" and S"<<electron_states[2]<<electron_states[3]<<"  :  "
<<mVp_total*nq_total-mVq_total*np_total-nVp_total*mq_total+nVq_total*mp_total<<endl;

//Two electron values
two_e[i*numStates+j]=total_coul+total_exch;
fprintf(fp_two_e, "%i %i %i %i %12.10le %12.10le %12.10le %12.10le\n", electron_states[0], electron_states[1], 
electron_states[2], electron_states[3], total_coul.real(), total_coul.imag(), total_exch.real(), total_exch.imag());
//debug
cout<<"Two_e: S"<<electron_states[0]<<electron_states[1]<<" and S"<<electron_states[2]<<electron_states[3]<<"  :  "
<<total_coul+total_exch<<endl;

//Single electron values
one_e[i*numStates+j]=2.0*(Ek[first_state[j]-1]+Ek[second_state[j]-1])*(mp_total*nq_total-np_total*mq_total);
fprintf(fp_one_e, "%i %i %i %i %12.10le %12.10le\n", electron_states[0], electron_states[1], electron_states[2], electron_states[3],
one_e[i*numStates+j].real(), one_e[i*numStates+j].imag());
//debug
cout<<"One electron: S"<<electron_states[0]<<electron_states[1]<<" and S"<<electron_states[2]<<electron_states[3]<<"  :  "
<<2.0*(Ek[first_state[j]-1]+Ek[second_state[j]-1])*(mp_total*nq_total-np_total*mq_total)<<endl;

}
delete(wf_e1);
}//end of for
}//end of for
fclose(fp_overlap);
fclose(fp_V_imp);
fclose(fp_two_e);
fclose(fp_one_e);

if(rank==0){
	//Diagonalize the 6 x 6 Hamiltonian
	nml_dcvector *H=nml_dcv_new(numStates*numStates);
	nml_dcvector *O=nml_dcv_new(numStates*numStates);
	FILE* Hartree=fopen("Hartree_Matrix_HM","w");
	if(Hartree==NULL)
		cout<<"Error opening Hartree File"<<endl;
	for(int i=0;i<numStates;i++)
          for(int j=i;j<numStates;j++){
                        H[i*numStates+j].r=one_e[i*numStates+j].real()+V_imp[i*numStates+j].real()+two_e[i*numStates+j].real();
			H[i*numStates+j].i=one_e[i*numStates+j].imag()+V_imp[i*numStates+j].imag()+two_e[i*numStates+j].imag();
			O[i*numStates+j].r=Overlap_Matrix[i*numStates+j].real();
			O[i*numStates+j].i=Overlap_Matrix[i*numStates+j].imag();		
		if(i!=j){
			H[j*numStates+i].r=one_e[i*numStates+j].real()+V_imp[i*numStates+j].real()+two_e[i*numStates+j].real();
                        H[j*numStates+i].i=-(one_e[i*numStates+j].imag()+V_imp[i*numStates+j].imag()+two_e[i*numStates+j].imag());
			O[j*numStates+i].r=Overlap_Matrix[j*numStates+i].real();
                        O[j*numStates+i].i=-Overlap_Matrix[j*numStates+i].imag();
		}
		if(i==j){
			H[i*numStates+j].i=0.0;
			O[i*numStates+j].i=0.0;
		}
		if(Hartree!=NULL)				
			fprintf(Hartree,"%i %i %12.10le %12.10le\n", i, j, H[i*numStates+j].r, H[i*numStates+j].i);
          }

fclose(Hartree);
cout<<"Starting Eigenvalue Calculations"<<endl;
int ITYPE=1; //Ax=EBx
char JOBZ='V';
char UPLO='U';
int N=numStates;
int LDA=N;
int LDB=N;
double *W=new double[N];
int LWORK=2*N-1;
int INFO=0;
nml_dcvector *WORK = nml_dcv_new(LWORK);
double *RWORK = new double[3*N-2];
nml_dcvector *A=nml_dcv_new(numStates*numStates);
nml_dcvector *B=nml_dcv_new(numStates*numStates);
for(int i=0;i<numStates;i++)
        for(int j=0;j<numStates;j++){
                A[i*numStates+j].r=H[j*numStates+i].r;
                A[i*numStates+j].i=H[j*numStates+i].i;
                B[i*numStates+j].r=O[j*numStates+i].r;
                B[i*numStates+j].i=O[j*numStates+i].i;
        }
zhegv_(&ITYPE, &JOBZ, &UPLO, &N, A, &LDA, B, &LDB, W, WORK, &LWORK, RWORK, &INFO);
cout<<"Eigenvalues : "<<endl;
for(int i=0;i<numStates;i++)
        cout<<W[i]<<endl;
cout<<"Singlet-Triplet Splitting (meV): "<<(W[1]-W[0])*1000<<endl;
for(int i=0;i<numStates;i++){
        cout<<"Eigenvector "<<i<<"  :  "<<endl;
        for(int j=0;j<numStates;j++)
                cout<<A[j*numStates+i].r<<"  "<<A[j*numStates+i].i<<endl;
}

double p_real=0.0;
cout<<"Computing Probablities : "<<endl;
for(int i=0;i<numStates;i++){
	cout<<endl;
        cout<<"WF "<<i<<"  :  "<<endl;
        for(int j=0;j<numStates;j++){
                p_real=A[j*numStates+i].r*A[j*numStates+i].r+A[j*numStates+i].i*A[j*numStates+i].i;
                cout<<p_real<<endl;
        }
}
cout<<"States  :  "<<endl;
for(int i=0;i<numStates;i++)
	cout<<i<<"    "<<first_state[i]<<"  "<<second_state[i]<<endl;
   
}//end of rank
}


if(ComputeCoulombHM) {

    //Computing Coulomb-Exchange energies ===============================
    Coulomb data_coulomb = 
	  Coulomb(NumAtoms, x1_start, x1_end, x2_start, x2_end, NumOrbitals);
    ifstream coulomb_input("coulomb.table");
    data_coulomb.assign_table(coulomb_input);
       //data_coulomb.assign_table1();
	
    //data_coulomb.get_data(wf_e, atomid, neighbor, lattice);

	FILE* fp=fopen("exchange_parameters","r");
	float dd_cutoff=12.0;
	double dd_cutoff1=dd_cutoff;
	if(fp!=NULL){
		fscanf(fp, "%f\n", &dd_cutoff);
		fclose(fp);
		dd_cutoff1=dd_cutoff;
	}else{
		if(rank==0)
			cout<<"File not found: Using default parameters. dd_cutoff=12 nm"<<endl;
	}
	if(rank==0){
		cout<<"dd_cutoff: "<<dd_cutoff1<<"  nm"<<endl;
	}

    double elapsed_time=0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time=-MPI_Wtime();	
    
    if(NumElectrons != 4) {
       cout << "Choose 4 as NumElectrons for Matrix Construction\n" << endl;
       exit(1);
    }

   if(x1_end != 0) data_coulomb.hartree_matrix(dd_cutoff1, data, eList, NumElectrons);

    elapsed_time+=MPI_Wtime();

    complex<double>* local_coul = data_coulomb.get_coul_Matrix();
    complex<double>* local_exch = data_coulomb.get_exch_Matrix();

    int csize = 55;
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

int *first_state=new int[6];
int *second_state=new int[6];
int count=0;
for(int i=1;i<NumElectrons;i++)
	for(int j=i+1;j<=NumElectrons;j++){
		first_state[count]=eList[i-1];
		second_state[count]=eList[j-1];
		count++;
	}
for(int i=0;i<6;i++)
cout<<"CI state "<<i<<"    :   "<<first_state[i]<<"  "<<second_state[i]<<endl;

    if (rank==0) {

	char filename1[200];
        char filename2[200];
        sprintf(filename1, "Coulomb_%d_%d_%d_%d", eList[0], eList[1], eList[2], eList[3]);
        sprintf(filename2, "Exchange_%d_%d_%d_%d", eList[0], eList[1], eList[2], eList[3]);
        //cout<<"Filename 1 : "<<filename1<<endl;
        //cout<<"Filename 2 : "<<filename2<<endl;


	FILE* fp1=fopen(filename1,"w");
	FILE* fp2=fopen(filename2,"w");
 
       int count = 0;
       cout<<endl<<"====================================================="<<endl;
       cout<<"Electron and Hole Pair: "<<endl;
       for (int rindx=1; rindx <= 6; rindx++){
          for (int cindx = rindx; cindx <= 6; cindx++){  
		cout<<" ["<<rindx<<","<<cindx<<"] "<<"CI:  "<<total_coul[count]<<"   "<<"EI:  "<<total_exch[count]<<endl;
 		 if(fp1!=NULL){
                	fprintf(fp1, "%i %i %i %i %e %e\n", first_state[rindx-1], second_state[rindx-1], 
				first_state[cindx-1], second_state[cindx-1], total_coul[count].real(), total_coul[count].imag());
			fprintf(fp2, "%i %i %i %i %e %e\n", first_state[rindx-1], second_state[rindx-1],
				first_state[cindx-1], second_state[cindx-1], total_exch[count].real(), total_exch[count].imag());
        	}
		count++;
	  }     
	}
       cout<<endl;
       cout<<"Time to compute the matrix :"<<elapsed_time<<endl;
       fclose(fp1);
       fclose(fp2);	 
  } //end of rank
} 

if(ComputeHartree) {
int numBasis=0;
FILE* readBasis=fopen("Basis","r");
if(readBasis!=NULL)
fscanf(readBasis, "%i\n", &numBasis);
cout<<"Number of Single Particle states  :  "<<numBasis<<endl; 
fclose(readBasis);
int numStates=numBasis*(numBasis-1)/2;
int *first_state=new int[numStates];
int *second_state=new int[numStates];
int count=0;
for(int i=1;i<numBasis;i++)
        for(int j=i+1;j<=numBasis;j++){
                first_state[count]=i;
                second_state[count]=j;
                count++;
        }
for(int i=0;i<numStates;i++)
cout<<"CI state "<<i<<"    :   "<<first_state[i]<<"  "<<second_state[i]<<endl;

nml_dcvector *coul_matrix=nml_dcv_new(numStates*numStates);
nml_dcvector *exch_matrix=nml_dcv_new(numStates*numStates);
double *Ek=new double[numBasis];

int num_distinct_groups=numBasis/2;
int num_Files=num_distinct_groups*(num_distinct_groups-1)/2;
cout<<"Number of Coulomb or exchange files to read  :   "<<num_Files<<endl;
int *first_group=new int[num_Files];
int *second_group=new int[num_Files];
count=0;
for(int i=1;i<num_distinct_groups;i++)
        for(int j=i+1;j<=num_distinct_groups;j++){
                first_group[count]=2*i-1;
                second_group[count]=2*j-1;
                count++;
        }
for(int i=0;i<num_Files;i++)
cout<<"Files: "<<i<<"    :   "<<first_group[i]<<"  "<<first_group[i]+1<<"  "<<second_group[i]<<"  "<<second_group[i]+1<<endl;

char filename1[200];
char filename2[200];
int numElements=21;
int tmp1, tmp2, tmp3, tmp4=0;
double coul_real, coul_imag=0.0;
int tmp5, tmp6, tmp7, tmp8=0;
double exch_real, exch_imag=0.0;
int row_index=0;
int col_index=0;

for(int i=0;i<num_Files;i++){
 	sprintf(filename1, "Coulomb_%d_%d_%d_%d", first_group[i], first_group[i]+1, second_group[i], second_group[i]+1);
        sprintf(filename2, "Exchange_%d_%d_%d_%d",first_group[i], first_group[i]+1, second_group[i], second_group[i]+1);
        cout<<"Filename 1 : "<<filename1<<endl;
        cout<<"Filename 2 : "<<filename2<<endl;
        FILE* fp1=fopen(filename1,"r");
        FILE* fp2=fopen(filename2,"r");
	if(fp1!=NULL & fp2!=NULL){
		for(int j=0;j<numElements;j++){				
			fscanf(fp1, "%i %i %i %i %lf %lf\n", &tmp1, &tmp2, &tmp3, &tmp4, &coul_real, &coul_imag); 
			//cout<<"coulomb: "<<tmp1<<" "<<tmp2<<" "<<tmp3<<" "<<tmp4<<"   "<<coul_real<<"   "<<coul_imag<<endl;
			fscanf(fp2, "%i %i %i %i %lf %lf\n", &tmp5, &tmp6, &tmp7, &tmp8, &exch_real, &exch_imag);
                        //cout<<"Exchange: "<<tmp5<<" "<<tmp6<<" "<<tmp7<<" "<<tmp8<<"   "<<exch_real<<"   "<<exch_imag<<endl;
			row_index=-1;
			col_index=-1;
			for(int k=0;k<numStates;k++){
				if((tmp1==first_state[k])&&(tmp2==second_state[k]))
					row_index=k;
				if((tmp3==first_state[k])&&(tmp4==second_state[k]))
                                        col_index=k;
			}
			cout<<"Row_index   :   "<<row_index<<endl;
			cout<<"Col_index   :   "<<col_index<<endl;
			coul_matrix[row_index*numStates+col_index].r=coul_real;
			coul_matrix[row_index*numStates+col_index].i=coul_imag;
			exch_matrix[row_index*numStates+col_index].r=exch_real;
                        exch_matrix[row_index*numStates+col_index].i=exch_imag;
		}
	}else{
		cout<<"Error Opening Files"<<endl;
	}

	fclose(fp1);
	fclose(fp2);
}


//print coulomb and exchange matrices
cout<<"Coulomb Matrix"<<endl;
for(int i=0;i<numStates;i++)
	for(int j=0;j<numStates;j++)
		cout<<i<<"   "<<j<<"   :   "<<coul_matrix[i*numStates+j].r<<"    "<<coul_matrix[i*numStates+j].i<<endl;
cout<<"Exchange Matrix"<<endl;
for(int i=0;i<numStates;i++)
        for(int j=0;j<numStates;j++)
                cout<<i<<"   "<<j<<"   :   "<<exch_matrix[i*numStates+j].r<<"    "<<exch_matrix[i*numStates+j].i<<endl;


// read in eigenvalues
   FILE* fp;
   char filename[200];   

   sprintf(filename, "Ek");
   fp=fopen(filename,"r");	
   cout << "Reading eigenvalues:  " << filename << endl;
   if (!fp) die("Could not find file %s\n", filename);
   for(int i=0;i<numBasis;i++)
	fscanf(fp, "%lf\n", &Ek[i]);	
   fclose(fp);

//print Eigenvalues
cout<<"Single electron eigenvalues   "<<endl;
for(int i=0;i<numBasis;i++)
	cout<<Ek[i]<<endl;


//Generate the total Many-body Hamiltonian
//reverse sign of off-diagonal elements for FORTRAN
nml_dcvector *H=nml_dcv_new(numStates*numStates);
for(int i=0;i<numStates;i++)
        for(int j=i;j<numStates;j++){
		if(j>i){
			H[i*numStates+j].r=+(coul_matrix[i*numStates+j].r+exch_matrix[i*numStates+j].r);
			H[i*numStates+j].i=+(coul_matrix[i*numStates+j].i+exch_matrix[i*numStates+j].i);
			H[j*numStates+i].r=+(coul_matrix[i*numStates+j].r+exch_matrix[i*numStates+j].r);
                        H[j*numStates+i].i=-(coul_matrix[i*numStates+j].i+exch_matrix[i*numStates+j].i);
		}else{
			H[i*numStates+i].r=+(coul_matrix[i*numStates+i].r+exch_matrix[i*numStates+i].r);
                        H[i*numStates+i].i=0.0;
			H[i*numStates+i].r+=Ek[first_state[i]-1]+Ek[second_state[i]-1];
		}
}
cout<<"Hartree Matrix"<<endl;
for(int i=0;i<numStates;i++)
        for(int j=0;j<numStates;j++)
                cout<<i<<"   "<<j<<"   :   "<<H[i*numStates+j].r<<"    "<<H[i*numStates+j].i<<endl;


cout<<"Starting Eiegnvalue Calculations"<<endl;

//test solver first
int num_Basis=numStates;

  char JOBZ='V';
  char UPLO='U';
  int N=num_Basis;
//  nml_dcvector *A=nml_dcv_new(N*N);


//test purposes
/*   A[0].r=1;
   A[0].i=0;
   A[1].r=0;
   A[1].i=0;
   A[2].r=0;
   A[2].i=0;

   A[3].r=0;
   A[3].i=0;
   A[4].r=0;
   A[4].i=0;
   A[5].r=1;
   A[5].i=0;

   A[6].r=0;
   A[6].i=0;
   A[7].r=1;
   A[7].i=0;
   A[8].r=1;
   A[8].i=0;
*/
  int LDA = N;
  double* W=new double[N];
  int LWORK=2*N;
  nml_dcvector *WORK=nml_dcv_new(LWORK);
  double* RWORK=new double[3*N];
  int INFO=0;
  zheev_( &JOBZ, &UPLO, &N, H, &LDA, W, WORK, &LWORK, RWORK, &INFO);

  cout<<"Eigenvalues : "<<endl;	
  for(int i=0;i<num_Basis;i++)
        cout<<W[i]<<endl;

  cout<<"Eigenvalues : "<<endl;
  for(int i=0;i<num_Basis;i++)
        printf("%lf\n", W[i]);

  for(int i=0;i<num_Basis;i++){
	cout<<endl;
        cout<<"Eigenvector  "<<i<<endl;
        for(int j=0;j<num_Basis;j++){
                cout<<H[i*num_Basis+j].r<<" "<<H[i*num_Basis+j].i<<endl;
        }
        }
//Compute probability densities of each Slater determinant contribution for each state.


  double p_real=0.0;	
  double sum=0.0;
  cout<<"Probability: "<<endl<<endl;
  for(int i=0;i<num_Basis;i++){
        cout<<endl;
        cout<<"Eigenvector  "<<i<<endl;
	sum=0;
        for(int j=0;j<num_Basis;j++){
		p_real=H[i*num_Basis+j].r*H[i*num_Basis+j].r+H[i*num_Basis+j].i*H[i*num_Basis+j].i;
		sum+=p_real;
                cout<<p_real<<endl;
        }
	//cout<<"Sum : "<<sum<<endl;
        }
double singlet_triplet_diff=(W[1]-W[0])*1000;
cout<<"Exchange Energy (2-1 difference)  :   "<<singlet_triplet_diff<<"   (meV)"<<endl;
}


delete [] eList; delete [] hList; delete [] Ee; delete [] Eh;
delete [] myAtoms_begin; delete [] myAtoms_end;

#if (defined MPI3d && !defined FAKE_MPI)
  MPI_Finalize();
#endif
  return 0;

}	
