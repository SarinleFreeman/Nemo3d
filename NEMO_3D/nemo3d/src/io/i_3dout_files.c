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
$Header: /repo/nemo3d/src/io/i_3dout_files.c,v 1.63 2008/09/24 02:09:19 hoonryu Exp $ 
*****************************************************************************/

#ifdef MACOSX
#include <sys/types.h>
#endif
#ifndef DONTBUILDFORMATCONVERTER
#include <netinet/in.h>
#endif
#include "i_3dout_files.h"

#include <sstream>

#define MAX_VARS 10

string replaceSuffix(const char* ch, const char* sfx, int indx) 
{
   stringstream ss(stringstream::in | stringstream::out);
   string fil(ch);
#define GCC_2_96
#ifndef GCC_2_96
   if (indx != -1)
      ss << string(fil, 0, fil.find_last_of(".")) << sfx << "_" << indx;
   else
      ss << string(fil, 0, fil.find_last_of(".")) << sfx;
   ss >> fil;
#else
   // stringstream is broken in old versions of gcc and in PGI
   char indx_str[10];
   
   fil = string(fil, 0, fil.find_last_of(".")) + string(sfx);
   if (indx != -1) {
      sprintf(indx_str,"_%d", indx);
      fil += string(indx_str);
   }
   
#endif   

   return fil;
}


void* malloc_strippedVector(qd_struct d, int sz_datstr, void* dat) 
{
   const int isx = d->cell_s[mpi_n3d_id];
   const int iex = isx + d->cell_ln[mpi_n3d_id];
   const int Natom = d->seg_ln[mpi_n3d_id] / d->NBasisStates;
   const int sz_rel = sz_datstr / sizeof(char);

   char* vIn = (char*) dat;
   char* vOut = (char*) nml_malloc(sz_datstr*Natom);
   char* ptr_v = vOut;
   
   if (!vOut) 
      die("MEMORY ALLOCATION ERROR -- size=%d\n", sz_datstr*Natom);

   for ( int i=isx; i < iex; i++ ) {
   for ( int j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
      if ( d->geo.AtomType[i][j] ) {
	 memcpy(ptr_v, vIn, sz_datstr);
	 ptr_v += sz_rel;
      }
      vIn += sz_rel;
   }}

   return (void*) vOut;
}


/* need to make sure that the allocatec memory is free'ed up with str_free */
char *getfileext(char *fname)
{ 
   char *token,fext[10], *result;
   
   /* Get file extension */
   strtok(fname,".\n");
   while((token=strtok(NULL,".\n")))
      strcpy(fext,token);
   
   result=n3d_strdup_n(fext);
   return result;
}



void writeHeader_siteInfo(qd_struct d, const char* infoType, 
                          const char* filename, const char* label) 
{
   char *usr_comment[MAX_VARS], *dim_lbl[MAX_VARS];
   int Nlbl=0, Ncomment=0, i, dim[MAX_VARS];

   if (mpi_n3d_id != mpi_n3d_masterid) return;

   /* determine number of labels and user comments */
   if (!strcmp(infoType,"basis")) {
      Nlbl=2;
      Ncomment=2;
   }
   else if (!strcmp(infoType,"atom")) {
      Nlbl=1;
      Ncomment=2;
   }
   else if (!strcmp(infoType,"shape")) {
      Nlbl=1;
      Ncomment=2;
   }
   else if (!strcmp(infoType,"cell")) {
      Nlbl=1;
      Ncomment=2;
   }
   else if (!strcmp(infoType,"dxW")) {
      Nlbl=1;
      Ncomment=5;
   }
   else if (!strcmp(infoType,"dxS")) {
      Nlbl=1;
      Ncomment=5;
   }
   else if (!strcmp(infoType,"dxP")) {
      Nlbl=1;
      Ncomment=5;
   }
   else if (!strcmp(infoType,"zcatom")) {
      Nlbl=0;
      Ncomment=7;
   }
   else if (!strcmp(infoType,"singleValue")) {
      Nlbl=1;
      Ncomment=2;
   }
   else if(!strcmp(infoType,"transitionRate")) {
      Nlbl=1;
      Ncomment=2;
   }

   else if(!strcmp(infoType,"overlap")) {
      Nlbl=1;
      Ncomment=2;
   }
 
   /* allocate temporary storage */
   for (i=0; i<Nlbl; i++)
      dim_lbl[i] = (char*) nml_malloc(sizeof(char)*50);
   dim_lbl[Nlbl] = NULL;
   dim[Nlbl] = 0;
   for (i=0; i<Ncomment; i++)
      usr_comment[i] = (char*) nml_malloc(sizeof(char)*50);
   usr_comment[Ncomment]=NULL;
   
   /* map enum type to string for printing */
   if      (d->opt.Dev.band_model==Dev_struct::Bands_20_sp3d5ss_spin)
      strcpy(usr_comment[0], "Bands_20_sp3d5ss_spin # band_model");
   else if (d->opt.Dev.band_model==Dev_struct::Bands_10_sp3ss_spin)
      strcpy(usr_comment[0], "Bands_10_sp3ss_spin  # band_model");
   else if (d->opt.Dev.band_model==Dev_struct::Bands_10_sp3d5ss_nospin)
      strcpy(usr_comment[0], "Bands_10_sp3d5ss_nospin  # band_model");
   else if (d->opt.Dev.band_model==Dev_struct::Bands_1_s_nospin)
      strcpy(usr_comment[0], "Bands_1_s_nospin  # band_model");
   else
      die("ERROR in writeHeader_siteInfo - unknown band model.");
   sprintf(usr_comment[1], "%d", d->opt.Dev.CrystalStruct);
   strcat(usr_comment[1]," # crystal structure");

   /* copy info into temporary storage */
   if (!strcmp(infoType,"basis")) {
      if      (d->opt.Dev.band_model==Dev_struct::Bands_20_sp3d5ss_spin)
	 dim[0] = 20;
      else if (d->opt.Dev.band_model==Dev_struct::Bands_10_sp3ss_spin)
	 dim[0] = 10;
      else if (d->opt.Dev.band_model==Dev_struct::Bands_10_sp3d5ss_nospin)
	 dim[0] = 10;
      else if (d->opt.Dev.band_model==Dev_struct::Bands_1_s_nospin)
	 dim[0] = 1;
      else
	 die("ERROR in writeHeader_siteInfo - unknown band model.");

      dim[1] = d->n_atom_tot;
   
      strcpy(dim_lbl[0],"# order of orbital basis");
      strcpy(dim_lbl[1],"# number of atoms");
      dim_lbl[2]=NULL;
   }
   else if (!strcmp(infoType,"atom")) {
      dim[0] = d->n_atom_tot;
      strcpy(dim_lbl[0],"# number of atoms");
   }
   else if (!strcmp(infoType,"shape")) {
      dim[0] = d->n_atom_tot;
      strcpy(dim_lbl[0],"# number of atoms");
   }
   else if (!strcmp(infoType,"cell")) {
      dim[0] = d->geo.N_Cell;
      strcpy(dim_lbl[0],"# unit cells");
   }
   else if (!strcmp(infoType,"dxP")) {
      char temp[10]; 
      int ev = 0;
      strcpy(temp,"");
      if(strstr(filename,"nd_dx_shapeL") != NULL) sprintf(temp, "shapeL");
      else if(strstr(filename,"nd_dx_shapeS") != NULL) sprintf(temp, "shapeS");            
      sprintf(usr_comment[2],"dim%11.3e%11.3e%11.3e", (real)d->geo.cell_xmax, (real)d->geo.cell_ymax, (real)d->geo.cell_zmax);
      sprintf(usr_comment[3],"%s", temp); 
      sprintf(usr_comment[4],"latticeConsT%11.3e%11.3e%11.3e", (real)d->geo.lattice_x, (real)d->geo.lattice_y, (real)d->geo.lattice_z);
      
      dim[0] = d->geo.N_Cell;
      strcpy(dim_lbl[0],"# unit cells");
   }
   else if (!strcmp(infoType,"dxS")) {
      char temp[5]; 
      int ev = 0;
      strcpy(temp,"");
      if(strstr(filename,"nd_dx_exxS") != NULL) sprintf(temp, "exxS");
      else if(strstr(filename,"nd_dx_exxL") != NULL) sprintf(temp, "exxL");
      else if(strstr(filename,"nd_dx_eyyS") != NULL) sprintf(temp, "eyyS");
      else if(strstr(filename,"nd_dx_eyyL") != NULL) sprintf(temp, "eyyL");
      else if(strstr(filename,"nd_dx_ezzS") != NULL) sprintf(temp, "ezzS");
      else if(strstr(filename,"nd_dx_ezzL") != NULL) sprintf(temp, "ezzL");
      else if(strstr(filename,"nd_dx_exyS") != NULL) sprintf(temp, "exyS");
      else if(strstr(filename,"nd_dx_exyL") != NULL) sprintf(temp, "exyL");
      else if(strstr(filename,"nd_dx_exzS") != NULL) sprintf(temp, "exzS");
      else if(strstr(filename,"nd_dx_exzL") != NULL) sprintf(temp, "exzL");
      else if(strstr(filename,"nd_dx_eyzS") != NULL) sprintf(temp, "eyzS");
      else if(strstr(filename,"nd_dx_eyzL") != NULL) sprintf(temp, "eyzL");                 
      sprintf(usr_comment[2],"dim%11.3e%11.3e%11.3e", (real)d->geo.cell_xmax, (real)d->geo.cell_ymax, (real)d->geo.cell_zmax);
      sprintf(usr_comment[3],"%s", temp); 
      sprintf(usr_comment[4],"latticeConsT%11.3e%11.3e%11.3e", (real)d->geo.lattice_x, (real)d->geo.lattice_y, (real)d->geo.lattice_z);
      
      dim[0] = d->geo.N_Cell;
      strcpy(dim_lbl[0],"# unit cells");
   }
   else if (!strcmp(infoType,"dxW")) {
      char temp[15]; 
      int ev = 0;
      strcpy(temp,"");
      for (int i = 1; i <= d->neigv; i++) {
	      sprintf(temp, "nd_dx_%d_wf", i);
	      if(strstr(filename,temp) != NULL) {
		      ev = i; break;
	      } 
      } 
      sprintf(usr_comment[2],"dim%11.3e%11.3e%11.3e", (real)d->geo.cell_xmax, (real)d->geo.cell_ymax, (real)d->geo.cell_zmax);
      sprintf(usr_comment[3],"ev%11.3e%21.12e", (real)ev, (real)d->WF[ev-1].E*HAMILTONIAN_SCALE_VALUE); 
      sprintf(usr_comment[4],"latticeConsT%11.3e%11.3e%11.3e", (real)d->geo.lattice_x, (real)d->geo.lattice_y, (real)d->geo.lattice_z);
      
      dim[0] = d->geo.N_Cell;
      strcpy(dim_lbl[0],"# unit cells");
   }
   else if (!strcmp(infoType,"zcatom")) {
      
      sprintf(usr_comment[2],"delta %11.3e%11.3e%11.3e", (real)d->geo.lattice_x, (real)d->geo.lattice_y, (real)d->geo.lattice_z);
      sprintf(usr_comment[3],"count %11d%11d%11d", (int)d->geo.cell_xmax, (int)d->geo.cell_ymax, (int)d->geo.cell_zmax);
      sprintf(usr_comment[4],"datatype  double64");
      sprintf(usr_comment[5],"datacount %d", d->geo.N_Cell);
      sprintf(usr_comment[6],"emprymark 1e10");
   
   }   
   else if (!strcmp(infoType,"singleValue")) {

      dim[0] = 4*d->neigv * d->neigv;
      strcpy(dim_lbl[0],"# 4*square of number of eigenvalues");

   }
   else if (!strcmp(infoType,"transitionRate")) {
 
      dim[0] = d->neigv * (d->neigv-1) / 2;
      strcpy(dim_lbl[0],"# eigenvalues*(eigenvalues-1)/2");
  
   }
   else if(!strcmp(infoType,"overlap")) {

      dim[0] = d->neigv * (d->neigv-1) / 2;
      strcpy(dim_lbl[0],"# eigenvalues*(eigenvalues-1)/2");

   }


   writeHeader(filename, label, dim, dim_lbl, usr_comment);

   /* de-allocate temporary storage */
   for (i=0; i<Nlbl; i++)
      nml_free(dim_lbl[i]);
   for (i=0; i<Ncomment; i++)
      nml_free(usr_comment[i]);
}


/* format of label:  {var11,...,var1n} .. {varm1,...,varmn} */
void writeSiteInfo(qd_struct d, const char* infoType, const char* filename, 
                   const char* label, void* dat0, ...) 
{

   FILE* fp;
   int i, j, proc, Ncol_var[MAX_VARS], sz_var[MAX_VARS], Ncol_tot, Nvar;
   int Nreal, Nint, Nfloat;
   char *dat, *datRcv, *ptr_dat[MAX_VARS];
   va_list argptr;
   int sz_datstr;
   int Nrow=0;

   char fieldType[21];
   bzero(fieldType, sizeof(char)*21);

   /* determine Ncol[], Ncol_tot, Nvar and number of reals and ints and floats */
   parseRecordFormat(label, &Nint, &Nreal, &Nfloat, &Ncol_var[0], &sz_var[0], &Nvar, 
                     fieldType);
   Ncol_tot = Nint + Nreal + Nfloat;

   // printf("Numbers of variables detected %d %d %d\n", Nint, Nreal, Nfloat); fflush(stdout);

   sz_datstr = Nreal*sizeof(double) + Nint*sizeof(int) + Nfloat * sizeof(nml_float);

   if (!strcmp(infoType,"basis")) {
      Nrow = d->seg_ln[mpi_n3d_id];
   }
   else if (!strcmp(infoType,"atom")) {
      Nrow = d->seg_ln[mpi_n3d_id] / d->NBasisStates;
   }
   else if (!strcmp(infoType,"shape")) {
      Nrow = d->seg_ln[mpi_n3d_id] / d->NBasisStates;
   }
   else if (!strcmp(infoType,"cell")) {
      Nrow = d->cell_ln[mpi_n3d_id];
   }
   else if (!strcmp(infoType,"dxW")) {
      Nrow = d->cell_ln[mpi_n3d_id];
   }
   else if (!strcmp(infoType,"dxS")) {
      Nrow = d->cell_ln[mpi_n3d_id];
   }
   else if (!strcmp(infoType,"dxP")) {
      Nrow = d->cell_ln[mpi_n3d_id];
   }
   else if (!strcmp(infoType,"zcatom")) {
      Nrow = d->seg_ln[mpi_n3d_id] / d->NBasisStates;
   }
   else if (!strcmp(infoType,"singleValue")) {
      Nrow = 4*d->neigv * d->neigv;
   }
   else if (!strcmp(infoType,"transitionRate")) {
      Nrow = d->neigv * (d->neigv-1) / 2;
   } 
   else if (!strcmp(infoType,"overlap")) {
      Nrow = d->neigv * (d->neigv-1) / 2;
   }

   else {
      die("ERROR -- wrong infoType");
   }
   
   writeHeader_siteInfo(d, infoType, filename, label);

   /* initialize ptr_dat (ptrs to arrays) */
   ptr_dat[0] = (char*) dat0;
   va_start(argptr, dat0);
   for (i=1; i<Nvar; i++) {
      ptr_dat[i] = va_arg(argptr, char*);
   }
   va_end(argptr);
   for (; i<MAX_VARS; i++) {
      ptr_dat[i] = NULL;
   }

   /* copy data to variable 'dat' for convenient message passing and printing*/
   dat = (char*) nml_malloc(sz_datstr*Nrow);
   for (i=0; i<Nrow; i++) {
      char* dat_ij = dat + i*(sz_datstr/sizeof(char));
      for (j=0; j<Nvar; j++) {
	 char* ptr_dat_ij = ptr_dat[j] + (sz_var[j]/sizeof(char))*i;
	 memcpy(dat_ij,	ptr_dat_ij, sz_var[j]);
	 dat_ij += sz_var[j]/sizeof(char);
      }
   }


   /* On a cluster we do not want to append into the same file 
      from different CPUs due to latency issues in the file system.
      We therefore send the data to the master CPU and that CPU 
      takes care of all the writing of the data.  We ASSUME that
      data are ordered in the same way as CPU index so that data
      on CPU i necessarily come before the data on CPU i+1 and
      after data on CPU i-1 */


   if (mpi_n3d_id == mpi_n3d_masterid){
      fp = fopen( filename, "ab" );
      
      /* master CPU writes its data first. */
      fwrite(dat, sz_datstr, Nrow, fp);

      printf("   RECV:  "), fflush(stdout);

      /* receive the data from the other CPUs and write it out */
      for (proc=0; proc<mpi_n3d_numprocs; proc++){
	 if (proc==mpi_n3d_masterid) continue;
	 
	 if (!strcmp(infoType,"basis")) {
	    Nrow = d->seg_ln[proc];
	 }
	 else if (!strcmp(infoType,"atom")) {
	    Nrow = d->seg_ln[proc] / d->NBasisStates;
	 }
	 else if (!strcmp(infoType,"shape")) {
	    Nrow = d->seg_ln[proc] / d->NBasisStates;
	 }
	 else if (!strcmp(infoType,"cell")) {
	    Nrow = d->cell_ln[proc];
	 }
         else if (!strcmp(infoType,"dxW")) {
            Nrow = d->cell_ln[proc]; 
         }
         else if (!strcmp(infoType,"dxS")) {
            Nrow = d->cell_ln[proc]; 
         }
	 else if (!strcmp(infoType,"dxP")) {
            Nrow = d->cell_ln[proc]; 
         }
         else if (!strcmp(infoType,"zcatom")) {
            Nrow = d->seg_ln[proc] / d->NBasisStates;
         }
         else if (!strcmp(infoType,"singleValue")) {
	    Nrow = 4*d->neigv * d->neigv;
	 }	
         else if (!strcmp(infoType,"transitionRate")) {
            Nrow = d->neigv * (d->neigv-1) / 2;
         } 
	 else if (!strcmp(infoType,"overlap")) {
     	    Nrow = d->neigv * (d->neigv-1) / 2;
   	 }


	 datRcv = (char*) nml_malloc(sz_datstr*Nrow);
	 
	 MPI_Recv(&datRcv[0], Nrow*sz_datstr, MPI_BYTE, proc, MPI_ANY_TAG,
		  MPI_COMM_WORLD, &mpi_n3d_status);

	 printf("%d  ", proc), fflush(stdout);

	 fwrite(datRcv, sz_datstr, Nrow, fp);

	 nml_free(datRcv);
      }
      printf("\n\n"), fflush(stdout);

      fclose(fp);
   } else{   /* Send the data to the master cpu. */
      MPI_Send(&dat[0], Nrow*sz_datstr, MPI_BYTE, mpi_n3d_masterid, mpi_n3d_id,
	       MPI_COMM_WORLD);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   nml_free(dat);
}


/* format of label:  {var11,...,var1n} .. {varm1,...,varmn} */
void readSiteInfo(qd_struct d, const char* infoType, const char* filename, void* dat0, ...) 
{
   FILE* fp;
   int i, j, proc, Ncol_var[MAX_VARS], sz_var[MAX_VARS], Ncol_tot, Nvar;
   int Nreal, Nint, Nfloat ;
   char *dat, *datSnd, *ptr_dat[MAX_VARS];
   va_list argptr;
   int sz_datstr;
   int Nrow=0, Nrow_p=0;
   char fieldType[21];
   char band_model[100]; 
   int n_tot, cryst_str, Nbasis;
   int Ndim, dim[MAX_VARS];
   char recFmt[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
   int Nuser;
   
   bzero(fieldType, sizeof(char)*21);

   for (i=0; i<MAX_VARS; i++) {
      dim[i] = 0;
      dim_lbl[i] = (char*) nml_malloc(sizeof(char)*100);
      usr_comment[i] = (char*) nml_malloc(sizeof(char)*100);
   }

   fp = readHeader(filename, recFmt, &Ndim, dim, dim_lbl, &Nuser, usr_comment);

   /* this n_tot is the total number of atoms; this is modified
      later so that n_tot equals (#atoms)*(#basis states) */
   Nbasis = dim[0];
   n_tot = dim[1];

   sscanf(usr_comment[0], "%s", band_model);
   sscanf(usr_comment[1], "%d", &cryst_str);

   /* determine Ncol[], Ncol_tot, Nvar and number of reals and ints and floats*/
   parseRecordFormat(recFmt, &Nint, &Nreal, &Nfloat, &Ncol_var[0], &sz_var[0], &Nvar, 
                     fieldType);
   Ncol_tot = Nint + Nreal + Nfloat;

   sz_datstr = Nreal*sizeof(double) + Nint*sizeof(int) + Nfloat * sizeof(nml_float);

   if (!strcmp(infoType,"basis")) {
      Nrow = d->seg_ln[mpi_n3d_id];
   }
   else if (!strcmp(infoType,"atom")) {
      Nrow = d->seg_ln[mpi_n3d_id] / d->NBasisStates;
   }
   else if (!strcmp(infoType,"shape")) {
      Nrow = d->seg_ln[mpi_n3d_id] / d->NBasisStates;
   }
   else if (!strcmp(infoType,"cell")) {
      Nrow = d->cell_ln[mpi_n3d_id];
   }
   else
      die("ERROR -- wrong infoType");
   
   /* initialize ptr_dat (ptrs to arrays) */
   ptr_dat[0] = (char*) dat0;
   va_start(argptr, dat0);
   for (i=1; i<Nvar; i++)
      ptr_dat[i] = va_arg(argptr, char*);
   va_end(argptr);
   for (; i<MAX_VARS; i++) {
      ptr_dat[i] = NULL;
   }

   // read data into variable 'dat' for convenient message passing and printing
   dat = (char*) nml_malloc(sz_datstr*Nrow);

   /* On a cluster we do not want to append into the same file 
      from different CPUs due to latency issues in the file system.
      We therefore send the data to the master CPU and that CPU 
      takes care of all the writing of the data.  We ASSUME that
      data are ordered in the same way as CPU index so that data
      on CPU i necessarily come before the data on CPU i+1 and
      after data on CPU i-1 */

   if (mpi_n3d_id == mpi_n3d_masterid){

      printf("master reads %d Nrows of size %d\n", Nrow, sz_datstr);

      /* master CPU reads its data first. */
      fread(dat, sz_datstr, Nrow, fp);

      /* read in more data and send to the other CPUs */
      for (proc=0; proc<mpi_n3d_numprocs; proc++){
	 if (proc==mpi_n3d_masterid) continue;
	 
	 if (!strcmp(infoType,"basis")) {
	    Nrow_p = d->seg_ln[proc];
	 }
	 else if (!strcmp(infoType,"atom")) {
	    Nrow_p = d->seg_ln[proc] / d->NBasisStates;
	 }
	 else if (!strcmp(infoType,"shape")) {
	    Nrow_p = d->seg_ln[proc] / d->NBasisStates;
	 }
	 else if (!strcmp(infoType,"cell")) {
	    Nrow_p = d->cell_ln[proc];
	 }
	 datSnd = (char*) nml_malloc(sz_datstr*Nrow_p);
	 fread(datSnd, sz_datstr, Nrow_p, fp);
	 MPI_Send(&datSnd[0], Nrow_p*sz_datstr, MPI_BYTE, proc, mpi_n3d_id,
		  MPI_COMM_WORLD);
	 nml_free(datSnd);
      }
       
   } else{   /* Receive the data from the master cpu. */
      MPI_Recv(&dat[0], Nrow*sz_datstr, MPI_BYTE, mpi_n3d_masterid, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &mpi_n3d_status);
   }

   MPI_Barrier(MPI_COMM_WORLD);

   /* copy back into input data structures */
   for (i=0; i<Nrow; i++) {
      char* dat_ij = dat + i*(sz_datstr/sizeof(char));
      for (j=0; j<Nvar; j++) {
	 char* ptr_dat_ij = ptr_dat[j] + (sz_var[j]/sizeof(char))*i;
	 memcpy(ptr_dat_ij, dat_ij, sz_var[j]);
	 dat_ij += sz_var[j]/sizeof(char);
      }
   }

   nml_free(dat);
   fclose(fp);
}


void writePos(qd_struct d, const char *type)
{
   int i,j;
   int isx = d->cell_s[mpi_n3d_id];
   int iex = isx + d->cell_ln[mpi_n3d_id];
   char label[100];
   char *filename;

   MPI_TIC(mpiTiming.fileio_0);

if ((!strcmp(type,"atom_strain_dx_large") || !strcmp(type,"atom_strain_dx_small"))
     && (d->opt.ExecParam.ElCalc.ElOut.PsiSqr ||
         d->opt.ExecParam.Output.AtomPosAfterStrain ||
         d->opt.ExecParam.Output.AtomPosBeforeStrain ||
         d->opt.ExecParam.ElCalc.ElOut.Eigvect) ) {

      static int counter=0;
      int indx=0;


   if(d->opt.Dev.CrystalStruct == Dev_struct::Zincblende) {

       real* pos_xx = (real*)
           nml_malloc(4*sizeof(real)*d->cell_ln[mpi_n3d_id]);
       real* pos_yy = (real*)
           nml_malloc(4*sizeof(real)*d->cell_ln[mpi_n3d_id]);
       real* pos_zz = (real*)
           nml_malloc(4*sizeof(real)*d->cell_ln[mpi_n3d_id]);
       real* pos_xy = (real*)
           nml_malloc(4*sizeof(real)*d->cell_ln[mpi_n3d_id]);
       real* pos_yz = (real*)
           nml_malloc(4*sizeof(real)*d->cell_ln[mpi_n3d_id]);
       real* pos_xz = (real*)
           nml_malloc(4*sizeof(real)*d->cell_ln[mpi_n3d_id]);

       for ( i=isx; i < iex; i++ ) {

                real exx, eyy, ezz, exy, exz, eyz;
                real unit_exx, unit_eyy, unit_ezz, unit_exy, unit_exz, unit_eyz;

                unit_exx=0;unit_eyy=0;unit_ezz=0;unit_exy=0;unit_exz=0;unit_eyz=0;
                exx=0;eyy=0;ezz=0;exy=0;exz=0;eyz=0; 

                pos_xx[indx] = d->geo.cell__ijk[i][0];
                pos_yy[indx] = d->geo.cell__ijk[i][0];
                pos_zz[indx] = d->geo.cell__ijk[i][0];
                pos_xy[indx] = d->geo.cell__ijk[i][0];
                pos_xz[indx] = d->geo.cell__ijk[i][0];
                pos_yz[indx] = d->geo.cell__ijk[i][0];
                indx++;

               pos_xx[indx] = d->geo.cell__ijk[i][1];
               pos_yy[indx] = d->geo.cell__ijk[i][1];
               pos_zz[indx] = d->geo.cell__ijk[i][1];
               pos_xy[indx] = d->geo.cell__ijk[i][1];
               pos_xz[indx] = d->geo.cell__ijk[i][1];
               pos_yz[indx] = d->geo.cell__ijk[i][1];
               indx++;

               pos_xx[indx] = d->geo.cell__ijk[i][2];
               pos_yy[indx] = d->geo.cell__ijk[i][2];
               pos_zz[indx] = d->geo.cell__ijk[i][2];
               pos_xy[indx] = d->geo.cell__ijk[i][2];
               pos_xz[indx] = d->geo.cell__ijk[i][2];
               pos_yz[indx] = d->geo.cell__ijk[i][2];
               indx++;

         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
            if ( d->geo.AtomType[i][j] ) {


               if (compute_strain_tensor(&exx, &eyy , &ezz, &exy, &exz, &eyz,
               d, i, j) == -1) { exx = 0;
                                 eyy = 0;
                                 ezz = 0;
                                 exy = 0;
                                 exz = 0;
                                 eyz = 0; }
               unit_exx+= exx;
               unit_eyy+= eyy;
               unit_ezz+= ezz;
               unit_exy+= exy;
               unit_exz+= exz;
               unit_eyz+= eyz;
             }
          }
               pos_xx[indx] = unit_exx;
               pos_yy[indx] = unit_eyy;
               pos_zz[indx] = unit_ezz;
               pos_xy[indx] = unit_exy;
               pos_xz[indx] = unit_exz;
               pos_yz[indx] = unit_eyz;
               indx++;
       }

     // string fil;

         if(!strcmp(type,"atom_strain_dx_large")) {
                 sprintf(label, "nd_dx_exxL");
         }
         else {
                 sprintf(label, "nd_dx_exxS");
         }
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
           masterPrint("writing strain profile exx to '%s' for dx", filename);
         writeSiteInfo(d, "dxS", filename, "{real x, real y, real z, real exx }", pos_xx);
         str_free(filename);

         if(!strcmp(type,"atom_strain_dx_large")) {
                 sprintf(label, "nd_dx_eyyL");
         }
         else {
                 sprintf(label, "nd_dx_eyyS");
         }
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile exx to '%s' for dx", filename);
         writeSiteInfo(d, "dxS", filename, "{real x, real y, real z, real eyy }", pos_yy);
         str_free(filename);

         if(!strcmp(type,"atom_strain_dx_large")) {
                 sprintf(label, "nd_dx_ezzL");
         }
         else {
                 sprintf(label, "nd_dx_ezzS");
         }
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile ezz to '%s' for dx", filename);
         writeSiteInfo(d, "dxS", filename, "{real x, real y, real z, real ezz }", pos_zz);
         str_free(filename);

         if(!strcmp(type,"atom_strain_dx_large")) {
                 sprintf(label, "nd_dx_exyL");
         }
         else {
                 sprintf(label, "nd_dx_exyS");
         }
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile exy to '%s' for dx", filename);
         writeSiteInfo(d, "dxS", filename, "{real x, real y, real z, real exy }", pos_xy);
         str_free(filename);

         if(!strcmp(type,"atom_strain_dx_large")) {
                 sprintf(label, "nd_dx_exzL");
         }
         else {
                 sprintf(label, "nd_dx_exzS");
         }
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile exz to '%s' for dx", filename);
         writeSiteInfo(d, "dxS", filename, "{real x, real y, real z, real exz }", pos_xz);
         str_free(filename);

         if(!strcmp(type,"atom_strain_dx_large")) {
                 sprintf(label, "nd_dx_eyzL");
         }
         else {
                 sprintf(label, "nd_dx_eyzS");
         }
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile eyz to '%s' for dx", filename);
         writeSiteInfo(d, "dxS", filename, "{real x, real y, real z, real eyy }", pos_yz);
         str_free(filename);

        nml_free(pos_xx);
        nml_free(pos_yy);
        nml_free(pos_zz);
        nml_free(pos_xy);
        nml_free(pos_yz);
        nml_free(pos_xz);

     }
     // counter++;
  }

  if (!strcmp(type,"shape_dx_large")) {

  	int indx = 0;

     	real *pos_shape = (real*)
        	nml_malloc(4*sizeof(real)*d->cell_ln[mpi_n3d_id]);

     	for(int i=isx; i<iex; i++) {

        	pos_shape[indx++] = d->geo.cell__ijk[i][0];
        	pos_shape[indx++] = d->geo.cell__ijk[i][1];
        	pos_shape[indx++] = d->geo.cell__ijk[i][2];
        	pos_shape[indx++] = 1.0;

     	}

     	sprintf(label, "nd_dx_shapeL");
     	filename = n3d_strdup_n(d->inputfile);
    	n3d_FileTypeSet(&filename, label, TRUE);
     	masterPrint("writing shape profile to '%s' for dx", filename);
     	writeSiteInfo(d, "dxP", filename, "{real x, real y, real z, real shape }", pos_shape);
     	str_free(filename);

     	nml_free(pos_shape);

  }

  if (!strcmp(type,"shape_dx_small")) {
   
	int  indx = 0;

     	real *pos_shape = (real*)
        	nml_malloc(4*sizeof(real)*d->cell_ln[mpi_n3d_id]);

     	for(int i=isx; i<iex; i++) {

        	pos_shape[indx++] = d->geo.cell__ijk[i][0];
        	pos_shape[indx++] = d->geo.cell__ijk[i][1];
        	pos_shape[indx++] = d->geo.cell__ijk[i][2];
        	pos_shape[indx++] = 1.0;

     	}

     	sprintf(label, "nd_dx_shapeS");
     	filename = n3d_strdup_n(d->inputfile);
     	n3d_FileTypeSet(&filename, label, TRUE);
     	masterPrint("writing shape profile to '%s' for dx", filename);
     	writeSiteInfo(d, "dxP", filename, "{real x, real y, real z, real shape }", pos_shape);
     	str_free(filename);

     	nml_free(pos_shape);

  }

  if ( !strcmp(type,"atom_pos_strain") &&
        (d->opt.ExecParam.ElCalc.ElOut.PsiSqr ||
         d->opt.ExecParam.Output.AtomPosAfterStrain ||
         d->opt.ExecParam.Output.AtomPosBeforeStrain ||
         d->opt.ExecParam.ElCalc.ElOut.Eigvect) ) {

      static int counter=0;
      int indx=0;

      real* pos_xx = (real*)
          nml_malloc(4*sizeof(real)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);
      real* pos_yy = (real*)
          nml_malloc(4*sizeof(real)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);
      real* pos_zz = (real*)
          nml_malloc(4*sizeof(real)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);
      real* pos_xy = (real*)
          nml_malloc(4*sizeof(real)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);
      real* pos_yz = (real*)
          nml_malloc(4*sizeof(real)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);
      real* pos_xz = (real*)
          nml_malloc(4*sizeof(real)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);

      for ( i=isx; i < iex; i++ ) {
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
            if ( d->geo.AtomType[i][j] ) {
               d->geo.getPosition_strained(&pos_xx[indx], i, j);
               d->geo.getPosition_strained(&pos_yy[indx], i, j);
               d->geo.getPosition_strained(&pos_zz[indx], i, j);
               d->geo.getPosition_strained(&pos_xy[indx], i, j);
               d->geo.getPosition_strained(&pos_xz[indx], i, j);
               d->geo.getPosition_strained(&pos_yz[indx], i, j);
               indx+=3;

               if (compute_strain_tensor(&pos_xx[indx], &pos_yy[indx] , &pos_zz[indx], &pos_xy[indx], &pos_xz[indx], &pos_yz[indx],
               d, i, j) == -1) { pos_xx[indx] = 0;
                                 pos_yy[indx] = 0;
                                 pos_zz[indx] = 0;
                                 pos_xy[indx] = 0;
                                 pos_xz[indx] = 0;
                                 pos_yz[indx] = 0; }
               indx++;
            }
         }
      }

         sprintf(label, "nd_exx");
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile exx to '%s'", filename);
         writeSiteInfo(d, "atom", filename, "{real x, real y, real z, real exx }", pos_xx);
         str_free(filename);

         sprintf(label, "nd_eyy");
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile eyy to '%s'", filename);
         writeSiteInfo(d, "atom", filename, "{real x, real y, real z, real eyy }", pos_yy);
         str_free(filename);

         sprintf(label, "nd_ezz");
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile ezz to '%s'", filename);
         writeSiteInfo(d, "atom", filename, "{real x, real y, real z, real ezz }", pos_zz);
         str_free(filename);

         sprintf(label, "nd_exy");
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile exy to '%s'", filename);
         writeSiteInfo(d, "atom", filename, "{real x, real y, real z, real exy }", pos_xy);
         str_free(filename);

         sprintf(label, "nd_exz");
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile exz to '%s'", filename);
         writeSiteInfo(d, "atom", filename, "{real x, real y, real z, real exz }", pos_xz);
         str_free(filename);

         sprintf(label, "nd_eyz");
         filename = n3d_strdup_n(d->inputfile);
         n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing strain profile eyz to '%s'", filename);
         writeSiteInfo(d, "atom", filename, "{real x, real y, real z, real eyz }", pos_yz);
         str_free(filename);
                                                                    
            nml_free(pos_xx);
        nml_free(pos_yy);
        nml_free(pos_zz);
        nml_free(pos_xy);
        nml_free(pos_yz);
        nml_free(pos_xz);
   }

   if ( !strcmp(type,"atom_pos") &&
        (d->opt.ExecParam.ElCalc.ElOut.PsiSqr ||
         d->opt.ExecParam.Output.AtomPosAfterStrain ||
         d->opt.ExecParam.Output.AtomPosBeforeStrain ||
         d->opt.ExecParam.ElCalc.ElOut.Eigvect) ) {
      static int counter=0;

      int indx=0;
      real* pos =
         (real*) nml_malloc(3*sizeof(real)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);

      for ( i=isx; i < iex; i++ ) {
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
            if ( d->geo.AtomType[i][j] ) {
               d->geo.getPosition_strained(&pos[indx], i, j);
               indx+=3;
            }
         }
      }

      string fil = replaceSuffix(d->inputfile, ".nd_rAtom", counter);
      // cout_master << "Writing final atomic positions to '" << fil << "'\n";
      writeSiteInfo(d, "atom", fil.c_str(), "{real x, real y, real z}", pos);
      nml_free(pos);
      counter++;

   }

   if (!strcmp(type,"surface_atom_pos")) {
      //int indx=0;
      char* filename = n3d_strdup_n(d->inputfile);
      FILE* f_sa;
      n3d_FileTypeSet(&filename, "nd_rSurfaceAtom", TRUE);
      masterPrint("Writing final surface atomic positions to '%s'", filename);
      f_sa = fopen(filename, "w");
      str_free(filename);
      double pos[3];

      if (!mpi_n3d_id){
         for ( i=isx; i < iex; i++ ) {
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
            if ( d->geo.AtomType[i][j] && d->geo.isOnSurfaceIgnoreBC(i,j) ) {
               d->geo.getPosition_strained(&pos[0], i, j);
               fprintf(f_sa, "%f %f  %f\n", pos[0],pos[1],pos[2]);
            }
         }}
      }
      fclose(f_sa);
   }

   if (d->opt.ExecParam.ElCalc.ElOut.PsiSqr && !strcmp(type,"equil")){
      static int counter=0;
      int indx=0;
      real* pos = (real*)
         nml_malloc(3*sizeof(real)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);

      for ( i=isx; i < iex; i++ ) {
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
             if ( d->geo.AtomType[i][j] ) {
               d->geo.getPosition_equilibrium(&pos[indx], i, j);
               indx+=3;
            }
         }
      }

      string fil = replaceSuffix(d->inputfile, ".nd_rAtomEquil", counter);
      cout_master << "Writing equilibrium atomic positions to '" << fil << "'\n";
      writeSiteInfo(d, "atom", fil.c_str(), "{real x, real y, real z}", pos);
      nml_free(pos);
      counter++;
   }

   if ( d->opt.ExecParam.ElCalc.ElOut.PsiSqr_cell &&
       !strcmp(type,"equil_cell")) {
      static int counter=0;

      int indx=0;
      real* pos = (real*) nml_malloc(3*sizeof(real)*d->cell_ln[mpi_n3d_id]);

      for ( i=isx; i < iex; i++ ) {
         d->geo.getPosition_equilibrium(&pos[indx], i, 0);
         indx+=3;
      }

      string fil = replaceSuffix(d->inputfile, ".nd_rCellEquil", counter);
      cout_master << "Writing equilibrium cellular positions to '" << fil << "'\n";
      writeSiteInfo(d, "cell", fil.c_str(), "{real x, real y, real z}", pos);
      nml_free(pos);
      counter++;
   }

   if (d->opt.ExecParam.Strain.StrainOut.DisplacementData && !strcmp(type,"strain_first")){
      char* filename = n3d_strdup_n(d->inputfile);

      printf("creating malloc_strippedVector"); fflush(stdout);
#ifdef DISPL_SINGLE_PRECISION
      nml_float * disp = (nml_float*)
         malloc_strippedVector(d, 3*sizeof(nml_float),
                               &d->geo.l_m_2_atomdsp[isx][0][0]);
#else /* DISPL_SINGLE_PRECISION */
      real * disp = (real*)
         malloc_strippedVector(d, 3*sizeof(real),
                               &d->geo.l_m_2_atomdsp[isx][0][0]);
#endif /* DISPL_SINGLE_PRECISION */

      n3d_FileTypeSet(&filename, "nd_disp_first", TRUE);
      masterPrint("Writing atomic displacements to '%s' after the first strain calculation.", filename);
#ifdef DISPL_SINGLE_PRECISION
      writeSiteInfo(d, "atom", filename, "{float x, float y, float z}", disp);
#else /* DISPL_SINGLE_PRECISION */
      writeSiteInfo(d, "atom", filename, "{real x, real y, real z}", disp);
#endif /* DISPL_SINGLE_PRECISION */
      str_free(filename);
      nml_free(disp);
   }

   if (d->opt.ExecParam.Strain.StrainOut.DisplacementData && !strcmp(type,"strain_second")){
      char* filename = n3d_strdup_n(d->inputfile);

      printf("creating malloc_strippedVector"); fflush(stdout);

#ifdef DISPL_SINGLE_PRECISION
      nml_float * disp = (nml_float*)
        malloc_strippedVector(d, 3*sizeof(nml_float),
                               &d->geo.l_m_2_atomdsp[isx][0][0]);
#else /* DISPL_SINGLE_PRECISION */
      real * disp = (real*)
         malloc_strippedVector(d, 3*sizeof(real),
                               &d->geo.l_m_2_atomdsp[isx][0][0]);
#endif /* DISPL_SINGLE_PRECISION */

      n3d_FileTypeSet(&filename, "nd_disp_second", TRUE);
      masterPrint("Writing atomic displacements to '%s' after the second strain calculation.", filename);
#ifdef DISPL_SINGLE_PRECISION
      writeSiteInfo(d, "atom", filename, "{float x, float y, float z}", disp);
#else /* DISPL_SINGLE_PRECISION */
      writeSiteInfo(d, "atom", filename, "{real x, real y, real z}", disp);
#endif /* DISPL_SINGLE_PRECISION */
      str_free(filename);
      nml_free(disp);
   }

   if ( d->opt.ExecParam.Strain.StrainOut.DisplacementData && !strcmp(type,"strain-checkpoint")){
      char* filename = n3d_strdup_n(d->inputfile);
#ifdef DISPL_SINGLE_PRECISION
      nml_float * disp = (nml_float*)
         malloc_strippedVector(d, 3*sizeof(nml_float),
                               &d->geo.l_m_2_atomdsp[isx][0][0]);
#else /* DISPL_SINGLE_PRECISION */
      real * disp = (real*)
         malloc_strippedVector(d, 3*sizeof(real),
                               &d->geo.l_m_2_atomdsp[isx][0][0]);
#endif /* DISPL_SINGLE_PRECISION */

      n3d_FileTypeSet(&filename, "nd_disp-checkpoint", TRUE);
      masterPrint("Writing atomic displacements to '%s'", filename);
#ifdef DISPL_SINGLE_PRECISION
      writeSiteInfo(d, "atom", filename, "{float x, float y, float z}", disp);
#else /* DISPL_SINGLE_PRECISION */
      writeSiteInfo(d, "atom", filename, "{real x, real y, real z}", disp);
#endif /* DISPL_SINGLE_PRECISION */
      str_free(filename);
      nml_free(disp);
   }

   if (d->opt.ExecParam.Strain.StrainOut.DisplacementData3D && !strcmp(type,"strain_all")){
      int indx=0;
      char* filename = n3d_strdup_n(d->inputfile);
      real* pos = (real*)
         nml_malloc(7*sizeof(real)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);

      for ( i=isx; i < iex; i++ ) {
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
            if ( d->geo.AtomType[i][j] ) {
               d->geo.getPosition_equilibrium(&pos[indx], i, j);
               indx+=3;
               pos[indx++] = d->geo.l_m_2_atomdsp[i][j][1];
               pos[indx++] = d->geo.l_m_2_atomdsp[i][j][2];
               pos[indx++] = sqrt(pos[indx-1]*pos[indx-1] +
                                  pos[indx-2]*pos[indx-2] +
                                  pos[indx-3]*pos[indx-3]);
            }
         }
      }

      n3d_FileTypeSet(&filename, "nd_disp_3d", TRUE);
      masterPrint("Writing 3D-plottable atomic displacements to '%s'", filename);
      writeSiteInfo(d, "atom", filename,
                    "{real x, real y, real z, real dx, real dy, real dz, real dr2}",
                    pos);
      str_free(filename);
      nml_free(pos);
   }

   if (d->opt.ExecParam.Strain.StrainOut.DisplacementData3D && !strcmp(type,"strain_all_second")){
      int indx=0;
      char* filename = n3d_strdup_n(d->inputfile);
      real* pos = (real*)
         nml_malloc(7*sizeof(real)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);

      for ( i=isx; i < iex; i++ ) {
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
            if ( d->geo.AtomType[i][j] ) {
               d->geo.getPosition_equilibrium(&pos[indx], i, j);
               indx+=3;
               pos[indx++] = d->geo.l_m_2_atomdsp[i][j][0];
               pos[indx++] = d->geo.l_m_2_atomdsp[i][j][1];
               pos[indx++] = d->geo.l_m_2_atomdsp[i][j][2];
               pos[indx++] = sqrt(pos[indx-1]*pos[indx-1] +
                                  pos[indx-2]*pos[indx-2] +
                                  pos[indx-3]*pos[indx-3]);
            }
         }
      }

      n3d_FileTypeSet(&filename, "nd_disp_3d_second", TRUE);
      masterPrint("Writing 3D-plottable atomic displacements to '%s'", filename);
      writeSiteInfo(d, "atom", filename,
                    "{real x, real y, real z, real dx, real dy, real dz, real dr2}",
                    pos);
      str_free(filename);
      nml_free(pos);
   }

   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}

void writeAtomInfo(qd_struct d)
{
   int i,j;
   int isx = d->cell_s[mpi_n3d_id];
   int iex = isx + d->cell_ln[mpi_n3d_id];
   static int counter=0;
   int indx=0;
   int* pos = NULL;

   MPI_TIC(mpiTiming.fileio_0);

   pos = (int*) nml_malloc(sizeof(int)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);
   for ( i=isx; i < iex; i++ ) {
     for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
       if ( d->geo.AtomType[i][j] ) {
           pos[indx++] = (int) d->geo.AtomType[i][j];
       }
     }
   }
   cout_master<< "What is d->inputfile? "<<d->inputfile<<endl; 
   string fil = replaceSuffix(d->inputfile, ".nd_aType", counter);
   cout_master << "Writing atom type to '" << fil << "'\n";
   writeSiteInfo(d, "atom", fil.c_str(), "{int t}", pos);
   nml_free(pos);
   counter++;

   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);
}

void writeShapeInfo(qd_struct d)
{
   int i,j;
   int isx = d->cell_s[mpi_n3d_id];
   int iex = isx + d->cell_ln[mpi_n3d_id];
   static int counter=0;
   int indx=0;
   int* pos = NULL;

   MPI_TIC(mpiTiming.fileio_0);

   pos = (int*) nml_malloc(sizeof(int)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);

   for ( i=isx; i < iex; i++ ) {
     for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
       if ( d->geo.AtomType[i][j] ) {
           pos[indx++] = d->geo.ShapeType[i][j];
       }
     }
   }
   cout_master<< "What is d->inputfile? "<<d->inputfile<<endl; 
   string fil = replaceSuffix(d->inputfile, ".nd_sType", counter);
   cout_master << "Writing shape type to '" << fil << "'\n";
   writeSiteInfo(d, "shape", fil.c_str(), "{int t}", pos);
   nml_free(pos);
   counter++;
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}

void writeNbrIndx(qd_struct d)
{
   int isx = d->cell_s[mpi_n3d_id];
   int iex = isx + d->cell_ln[mpi_n3d_id];
   bool has_periodicity = d->geo.hasPeriodicity();

   int indx=0;

   int neighbors=d->geo.NeighborsMax();

   int* nbrInfo = NULL;

   MPI_TIC(mpiTiming.fileio_0);


   nbrInfo = (int*) nml_malloc(neighbors*sizeof(int)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);
   for (int Zc=isx; Zc < iex; Zc++) {
      int i = d->geo.cell__ijk[Zc][0];
      int j = d->geo.cell__ijk[Zc][1];
      int k = d->geo.cell__ijk[Zc][2];

      for (int Za=0; Za < d->geo.AtomsPerCellMax(); Za++) {
         if ( !d->geo.AtomType[Zc][Za] ) continue;
      
         for (int n = 0; n < d->geo.Neighbors(Za); n++) {
            int i_nbr = i + d->geo.NbrCell(Za,n,0);
            int j_nbr = j + d->geo.NbrCell(Za,n,1);
            int k_nbr = k + d->geo.NbrCell(Za,n,2);

            int Za_nbr = d->geo.NbrCell(Za,n,3);
            int Zc_nbr = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

#ifdef ELIMINATE_SSMAP
            if (Zc_nbr < 0 && has_periodicity && (Zc_nbr + d->geo.N_Cell + 1 ) >= 0)
	      Zc_nbr = Zc_nbr + d->geo.N_Cell + 1 ; 
            if(d->geo.is_tilted)
               Zc_nbr = d->geo.get_cindx_tilted(Zc,Za,n);
#else /* ELIMINATE_SSMAP */
            if (Zc_nbr < 0 && has_periodicity && d->geo.ssmap[Zc][Za][n]!=-1)
               Zc_nbr = d->geo.ssmap[Zc][Za][n];
#endif /* ELIMINATE_SSMAP */
         
            nbrInfo[indx++] = 
               ( Zc_nbr >= 0 )
               ? d->geo.l_m_2_offset_ham_tot[Zc_nbr][Za_nbr] / d->NBasisStates
               : -1;
         }
      }
   }

   char* filename = n3d_strdup_n(d->inputfile);
   n3d_FileTypeSet(&filename, "nd_nbrIndx", TRUE);
   masterPrint("Writing neighbor info to '%s'", filename);
   writeSiteInfo(d, "atom", filename, "{int n1, int n2, int n3, int n4}", nbrInfo);
   str_free(filename);
   nml_free(nbrInfo);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}


void writePhaseInfoForNEMO1D(qd_struct d)
{
   int isx = d->cell_s[mpi_n3d_id];
   int iex = isx + d->cell_ln[mpi_n3d_id];
   bool has_periodicity = d->geo.hasPeriodicity();

   int NeighborsMax = d->geo.NeighborsMax();

   int indx=0;
   int* phaseInfo = (int*) 
      nml_malloc(NeighborsMax*sizeof(int)*d->seg_ln[mpi_n3d_id]/d->NBasisStates);

   for (int Zc=isx; Zc < iex; Zc++) {
      int i = d->geo.cell__ijk[Zc][0];
      int j = d->geo.cell__ijk[Zc][1];
      int k = d->geo.cell__ijk[Zc][2];

      for (int Za=0; Za < d->geo.AtomsPerCellMax(); Za++) {
         if ( !d->geo.AtomType[Zc][Za] ) continue;
      
         for (int n = 0; n < d->geo.Neighbors(Za); n++) {
            int i_nbr = i + d->geo.NbrCell(Za,n,0);
            int j_nbr = j + d->geo.NbrCell(Za,n,1);
            int k_nbr = k + d->geo.NbrCell(Za,n,2);

            int Za_nbr = d->geo.NbrCell(Za,n,3);
            int Zc_nbr = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

#ifdef ELIMINATE_SSMAP
            if (Zc_nbr < 0 && has_periodicity && (Zc_nbr + d->geo.N_Cell + 1 ) >= 0) {
	      Zc_nbr = Zc_nbr + d->geo.N_Cell + 1 ; 
            }
            if(d->geo.is_tilted)
               Zc_nbr = d->geo.get_cindx_tilted(Zc,Za,n);
#else /* ELIMINATE_SSMAP */
            if (Zc_nbr < 0 && has_periodicity && d->geo.ssmap[Zc][Za][n]!=-1)
               Zc_nbr = d->geo.ssmap[Zc][Za][n];
#endif /* ELIMINATE_SSMAP */

            phaseInfo[indx++] = (Zc_nbr >= 0 ? d->geo.phaseInfoForNemo1D(Zc,Za,n) : 0);
         }
      }
   }

   char* filename = n3d_strdup_n(d->inputfile);
   n3d_FileTypeSet(&filename, "nd_phaseInfo", TRUE);
   masterPrint("Writing neighbor info to '%s'", filename);
   writeSiteInfo(d, "atom", filename, "{int ph1, int ph2, int ph3, int ph4}", phaseInfo);
   str_free(filename);
   nml_free(phaseInfo);
}


void writeEigenvector(qd_struct d, complex* wf, int indx, const char* suffix)
{
   if (!d->opt.ExecParam.ElCalc.ElOut.Eigvect)  return;
   
   char label[40];
   sprintf(label,"nd_evec%s_%d", suffix, indx+1);
   char* filename = n3d_strdup_n(d->inputfile);

   MPI_TIC(mpiTiming.fileio_0);

   n3d_FileTypeSet(&filename, label, TRUE);
   masterPrint("writing wave function to '%s'", filename);
   writeSiteInfo(d, "basis", filename, "{real wf.r, real wf.i}", wf);
   str_free(filename);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}

void writeOverLap(qd_struct d, complex *wf, int indx, const char* suffix)
{

   char label[40];
   sprintf(label,"nd_%s", suffix);
   char* filename = n3d_strdup_n(d->inputfile);

   MPI_TIC(mpiTiming.fileio_0);

   n3d_FileTypeSet(&filename, label, TRUE);
   masterPrint("writing Overlap Matrix to '%s'", filename);
   writeSiteInfo(d, "overlap", filename, "{real istate, real fstate, real overlap.r, real overlap.i}", wf);
   str_free(filename);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}

void writeTransitionRate(qd_struct d, real *wf, int indx, const char* suffix)
{

   char label[40];
   sprintf(label,"nd_%s", suffix);
   char* filename = n3d_strdup_n(d->inputfile);
              
   MPI_TIC(mpiTiming.fileio_0);
              
   n3d_FileTypeSet(&filename, label, TRUE);
   masterPrint("writing momentum matrix to '%s'", filename);
   writeSiteInfo(d, "transitionRate", filename, "{real istate, real fstate, real E, real transitionrate}", wf);
   str_free(filename);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}

void writeMomentumMatrix(qd_struct d, complex* wf, int indx, const char* suffix)
{
   char label[40];
   sprintf(label,"nd_momentum_%s", suffix);
   char* filename = n3d_strdup_n(d->inputfile);
              
   MPI_TIC(mpiTiming.fileio_0);
              
   n3d_FileTypeSet(&filename, label, TRUE);
   masterPrint("writing momentum matrix to '%s'", filename);
   writeSiteInfo(d, "singleValue", filename, "{real initState, real finState, real wf.r, real wf.i}", wf);
   str_free(filename);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);
     
}    


void writeOneVector(qd_struct d, complex* x, complex* y) 
{
  char label[40];
  char *filename;

  sprintf(label,"nd_firstx");
  filename = n3d_strdup_n(d->inputfile);
  n3d_FileTypeSet(&filename, label, TRUE);
  writeSiteInfo(d,"basis", filename, "{real x.r, real x.i}", x);

  sprintf(label,"nd_firsty");
  n3d_FileTypeSet(&filename, label, TRUE);
  writeSiteInfo(d,"basis", filename, "{real y.r, real y.i}", y);

  str_free(filename);
  MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}

void writePsiSqr(cmatrix wfc, int n, qd_struct d, int wf_indx)
{
   MPI_TIME_INIT(Tstart);

   char label[100], *filename;
   int i,j,k,l,m, isx, iex, imin, imax;
   real wfsum=0.0, wfmax, wfmax_local;
   rmatrix psisq = NULL;

   MPI_TIC(mpiTiming.fileio_0);

   //psisq = Rmatrix(d->seg_ln[mpi_n3d_id],d->geo.AtomsPerCellMax());
   psisq = Rmatrix(d->cell_ln[mpi_n3d_id],d->geo.AtomsPerCellMax());
   isx = d->cell_s[mpi_n3d_id];
   iex = isx + d->cell_ln[mpi_n3d_id];
   imin = d->seg_s[mpi_n3d_id];
   imax = d->seg_s[mpi_n3d_id] + d->seg_ln[mpi_n3d_id];

   /* Normalize and print the square of the wavefunctions */
   for ( k=0; k < n; k++ ){
      wfmax = -1e100;
      for ( i=isx; i < iex; i++ ){
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ){
            if ( d->geo.AtomType[i][j] ){
               wfsum = 0;
               for ( l=0; l < d->NBasisStates; l++ ){
                  m = d->geo.l_m_2_offset_ham_tot[i][j]+l-imin;
                  wfsum += wfc[k][m].r*wfc[k][m].r +  wfc[k][m].i*wfc[k][m].i;
               }
            }
            psisq[i-isx][j] = wfsum;
            if ( wfsum > wfmax )
               wfmax = wfsum;
         }
      }
      wfmax_local= wfmax;

#if (defined MPI3d && !defined FAKE_MPI)
      MPI_Allreduce(&wfmax_local,&wfmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
      wfmax = wfmax_local;
#endif

      for ( i=isx; i < iex; i++ )
	 for ( j=0; j < d->geo.AtomsPerCellMax(); j++ )
	    if ( d->geo.AtomType[i][j] ) {
	       psisq[i-isx][j] /= wfmax;
	    }

      if (d->opt.ExecParam.ElCalc.ElOut.PsiSqr){
	 real* psi2 = 
            (real*) malloc_strippedVector(d, sizeof(real), &psisq[0][0]);
	 
	 sprintf(label, "nd_wf_%d", wf_indx+k+1);
	 filename = n3d_strdup_n(d->inputfile);
	 n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing probability density to '%s'", filename);
	 writeSiteInfo(d, "atom", filename, "{real psi2}", psi2);
	 str_free(filename);
	 nml_free(psi2);
      }
       
      if (d->opt.ExecParam.ElCalc.ElOut.PsiSqr_cell){

   	 real* psi2 = (real*) nml_malloc(sizeof(real)*d->cell_ln[mpi_n3d_id]); 
    
	 for ( i=isx; i < iex; i++ ) {
	    real dummy=0.0;
	    for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	       if ( d->geo.AtomType[i][j] ) {
		  dummy += psisq[i-isx][j];
	       }
	    }
            psi2[i-isx] = dummy;
	 }
  
	 sprintf(label, "nd_wf_%d_unit", wf_indx+k+1);
	 filename = n3d_strdup_n(d->inputfile);
	 n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing probability density to '%s'", filename);
	 writeSiteInfo(d, "cell", filename, "{real psi2_unit}", psi2);
	 str_free(filename);
	 nml_free(psi2);
      }

     if (d->opt.ExecParam.ElCalc.ElOut.Visualization_3D){
	      
         real* psi2 = (real*) nml_malloc(4*sizeof(real)*d->cell_ln[mpi_n3d_id]);
    
         int index = 0;
    
	 for ( i=isx; i < iex; i++ ) {
	    real dummy=0.0;
	    for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	       if ( d->geo.AtomType[i][j] ) {
		  dummy += psisq[i-isx][j];
	       }
	    }
 
	    psi2[index++] = d->geo.cell__ijk[i][0];
            psi2[index++] = d->geo.cell__ijk[i][1];
            psi2[index++] = d->geo.cell__ijk[i][2];
	    psi2[index++] = dummy;
	 } 

         sprintf(label, "nd_dx_%d_wf", wf_indx+k+1);
	 filename = n3d_strdup_n(d->inputfile);
	 n3d_FileTypeSet(&filename, label, TRUE);
         masterPrint("writing probability density to '%s' for dx", filename);
	 writeSiteInfo(d, "dxW", filename, "{real x, real y, real z, real psi2_unit}", psi2);
	 str_free(filename);
	 nml_free(psi2);

     }

// Newly added routines to generate atomic scale data on which the new visualizer is based.
 #undef ATOM_DATA_DUMP

 //#define ATOM_DATA_DUMP
 #ifdef ATOM_DATA_DUMP

     real* psi2 = (real*) nml_malloc(8*4*sizeof(real)*d->cell_ln[mpi_n3d_id]);
     int index = 0;

     for (i=isx;i<iex;i++) {
          for(j=0;j<d->geo.AtomsPerCellMax();j++) {

              psi2[index++] = d->geo.cell__ijk[i][0];
              psi2[index++] = d->geo.cell__ijk[i][1];
              psi2[index++] = d->geo.cell__ijk[i][2];

              if (d->geo.AtomType[i][j]) psi2[index++] = psisq[i-isx][j];
              else psi2[index++] = 1e10;

          }
     }

     sprintf(label, "nd_zc_%d_wf", wf_indx+k+1);
     filename = n3d_strdup_n(d->inputfile);
     n3d_FileTypeSet(&filename, label, TRUE);
     masterPrint("writing zc data to '%s' for atom-based visualizer", filename);
     writeSiteInfo(d, "zcatom", filename, "{real x, real y, real z, real psi2}", psi2);
     str_free(filename);
     nml_free(psi2);


 #endif //ATOM_DATA_DUMP

 #undef PSI2_INDIVIDUAL_COMPONENT
 
 //#define PSI2_INDIVIDUAL_COMPONENT
 #ifdef PSI2_INDIVIDUAL_COMPONENT
 
 for ( l=0; l < d->NBasisStates; l++ ) {
     
     real* psi2_1 = (real*) nml_malloc(8*4*sizeof(real)*d->cell_ln[mpi_n3d_id]);
     int index2 = 0;
     
     for ( i=isx; i < iex; i++ ){
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ){
             if ( d->geo.AtomType[i][j] ){
                  wfsum = 0;
                  m = d->geo.l_m_2_offset_ham_tot[i][j]+l-imin;
                  wfsum = wfc[k][m].r*wfc[k][m].r +  wfc[k][m].i*wfc[k][m].i;
                  psi2_1[index2++] = d->geo.cell__ijk[i][0];
                  psi2_1[index2++] = d->geo.cell__ijk[i][1];
                  psi2_1[index2++] = d->geo.cell__ijk[i][2];
                  psi2_1[index2++] = wfsum/wfmax;
             }
         }
     }

 
     sprintf(label, "nd_zc_%d_wf_orbital_%d", wf_indx+k+1,l+1);
     filename = n3d_strdup_n(d->inputfile);
     n3d_FileTypeSet(&filename, label, TRUE);
     masterPrint("writing zc data to '%s' for atom-based visualizer", filename);
     writeSiteInfo(d, "atom", filename, "{real x, real y, real z, real psi2_unit}", psi2_1);
     str_free(filename);
     nml_free(psi2_1);
 
 }

 #endif //PSI2_INDIVIDUAL_COMPONENT

 #undef PSI2_SPDSS_SEP
 
 //#define PSI2_SPDSS_SEP
 #ifdef PSI2_SPDSS_SEP
 
 for (int seq = 0; seq < 8; seq++) {
     
     int l_start = 0, l_end = 0;
     switch (seq) {
     
     case 0: l_start = 0; l_end = 1; break;
     case 1: l_start = 1; l_end = 2; break;
     case 2: l_start = 2; l_end = 5; break;
     case 3: l_start = 5; l_end = 10; break;
     case 4: l_start = 10; l_end = 11; break;
     case 5: l_start = 11; l_end = 12; break;
     case 6: l_start = 12; l_end = 15; break;
     case 7: l_start = 15; l_end = 20; break;
     default: printf("\n NOTHING HAPPENS!");
     
     }
     
     real* psi2_2 = (real*) nml_malloc(8*4*sizeof(real)*d->cell_ln[mpi_n3d_id]);
     int index3 = 0;
     
     for ( i=isx; i < iex; i++ ){
        for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ){
             
             wfsum = 0;
             if ( d->geo.AtomType[i][j] ){
                for ( l=l_start; l < l_end; l++ ) {

                  m = d->geo.l_m_2_offset_ham_tot[i][j]+l-imin;
                  wfsum += wfc[k][m].r*wfc[k][m].r +  wfc[k][m].i*wfc[k][m].i;

                }
             }

             psi2_2[index3++] = d->geo.cell__ijk[i][0];
             psi2_2[index3++] = d->geo.cell__ijk[i][1];
             psi2_2[index3++] = d->geo.cell__ijk[i][2];
             psi2_2[index3++] = wfsum/wfmax;

        }
     }

     sprintf(label, "nd_zc_%d_wf_SDPSS_SPINSEP_%d", wf_indx+k+1,seq+1);
     filename = n3d_strdup_n(d->inputfile);
     n3d_FileTypeSet(&filename, label, TRUE);
     masterPrint("writing zc data to '%s' for atom-based visualizer", filename);
     writeSiteInfo(d, "atom", filename, "{real x, real y, real z, real psi2_unit}", psi2_2);
     str_free(filename);
     nml_free(psi2_2);

 }

 #endif //PSI2_SPDSS_SEP 

 #undef PSI2_SPDSS_ALL

 //#define PSI2_SPDSS_ALL
 #ifdef PSI2_SPDSS_ALL

 for (int seq = 0; seq < 4; seq++) {

     int l_start = 0, l_end = 0;
     switch (seq) {

     case 0: l_start = 0; l_end = 1; break;
     case 1: l_start = 1; l_end = 2; break;
     case 2: l_start = 2; l_end = 5; break;
     case 3: l_start = 5; l_end = 10; break;
     default: printf("\n NOTHING HAPPENS!");

     }

     real* psi2_3 = (real*) nml_malloc(8*4*sizeof(real)*d->cell_ln[mpi_n3d_id]);
     int index4 = 0;

     for ( i=isx; i < iex; i++ ){
        for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ){

             wfsum = 0;
             if ( d->geo.AtomType[i][j] ){

               for(int deci = 0; deci < 2; deci++) {
                 for ( l=l_start; l < l_end; l++ )   {

                  m = d->geo.l_m_2_offset_ham_tot[i][j]+ l+deci*10 -imin;
                  wfsum += wfc[k][m].r*wfc[k][m].r +  wfc[k][m].i*wfc[k][m].i;

                 }
               }
             }

             psi2_3[index4++] = d->geo.cell__ijk[i][0];
             psi2_3[index4++] = d->geo.cell__ijk[i][1];
             psi2_3[index4++] = d->geo.cell__ijk[i][2];
             psi2_3[index4++] = wfsum/wfmax;

        }
     }

     sprintf(label, "nd_zc_%d_wf_SDPSS_SPINALL_%d", wf_indx+k+1,seq+1);
     filename = n3d_strdup_n(d->inputfile);
     n3d_FileTypeSet(&filename, label, TRUE);
     masterPrint("writing zc data to '%s' for atom-based visualizer", filename);
     writeSiteInfo(d, "atom", filename, "{real x, real y, real z, real psi2_unit}", psi2_3);
     str_free(filename);
     nml_free(psi2_3);

  }

 #endif //PSI2_SPDSS_ALL
      
 }

   rm_rmatrix(&psisq);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}


void printPsisq(qd_struct d, cvectr p, int n_offset)
{
   cmatrix p_matrix = Cmatrix(1,1);
   cvectr c_dummy = p_matrix[0];
   p_matrix[0] = p;
   writePsiSqr(p_matrix,1,d,n_offset);
   /* print_psisq_gnu_pdb(p_matrix,1,d,n_offset); */
   p_matrix[0] = c_dummy;
   rm_cmatrix(&p_matrix);
}


int print_bond_lengths( qd_struct d )
{
   int i, j, k, l, m, n, p, id, jd, kd, aindx, cindx, nb;
   real x, y, z, bl, blave, disp, dum;
   rvectr nnv;
   FILE *f_bl;
   char *filename=NULL;
   real a_lattice[3];
   
   if (!d->opt.ExecParam.Strain.StrainOut.Bondlength || mpi_n3d_id!=mpi_n3d_masterid)
      return 0;
   
   MPI_TIC(mpiTiming.fileio_0);

   a_lattice[0] = d->geo.lattice_x;
   a_lattice[1] = d->geo.lattice_y;
   a_lattice[2] = d->geo.lattice_z;
   
   filename = n3d_strdup_n(d->inputfile);
   n3d_FileTypeSet(&filename, "nd_bond", TRUE);
   
   f_bl = fopen(filename, "w");
   str_free(filename);
   
   nnv = Rvectr(3);
   
   /* Loop through all unit cells in the structure */
   for ( l=0; l < d->geo.N_Cell; l++ ) {
      
      /* These are real space indexes for the unit cell */
      i = d->geo.cell__ijk[l][0];
      j = d->geo.cell__ijk[l][1];
      k = d->geo.cell__ijk[l][2];
      
      /* Loop through all the atoms in the unit cell */
      for ( m=0; m < d->geo.AtomsPerCellMax(); m++ ) {
	 
	 /*  Check if atom is in simulation domain */ 
	 if ( d->geo.AtomType[l][m] ) {
	    
	    nb = 0;
	    blave = 0;
	    /* Loop through the atom's nearest neighbors */
	    for ( n = 0; n < d->geo.Neighbors(m); n++ ) {
	       
	       /* Get relative position of unit cell where this neighbor lives. */ 
	       id = d->geo.NbrCell(m,n,0);
	       jd = d->geo.NbrCell(m,n,1);
	       kd = d->geo.NbrCell(m,n,2);
	       
	       /* Get the unit cell and atomic index for this neighbor */
	       cindx = d->geo.ijk__cell[i+id][j+jd][k+kd];
	       aindx = d->geo.NbrCell(m,n,3);
	       
	       /* Check if neighbor's unit cell is within simulation domain. */
	       if ( cindx >= 0 ) {
		  
		  /* Check if the neighbor is within the simulation domain */
		  if ( d->geo.AtomType[cindx][aindx] ) {
		     
		     /* Get relative real-space location for this neighbor.  There are two
			terms inside the brackets.  The first gives the relative location
			of the neighbor's unit cell.  The second gives the relative location
			of the neighbor's atom within the unit cell.  The term disp is the
			atom's displacement due to strain */
		     bl = 0;
		     nb++; 
		     for ( p=0; p < 3; p++ ) {
			
			disp = d->geo.l_m_2_atomdsp[cindx][aindx][p] - d->geo.l_m_2_atomdsp[l][m][p];
			dum = a_lattice[p]*( d->geo.NbrCell(m,n,p) + 
                                             (d->geo.PositionInCell(d->geo.NbrCell(m,n,3),p) - 
				      d->geo.PositionInCell(m,p)) ) + disp;
			bl += dum*dum;
			if (d->strain.debugp)
			   printf("\n%1.10f", fabs(dum)*4);
		     }
		     blave += bl;
		     if ( d->strain.debugp )
			printf("\nBL = %1.10f\n", sqrt(bl));
		  }
	       }
	    }	    
	    blave /= nb;
	    x = d->geo.cell__ijk[l][0]*d->geo.lattice_x + d->geo.lattice_x*d->geo.PositionInCell(m,0) ;
	    y = d->geo.cell__ijk[l][1]*d->geo.lattice_y + d->geo.lattice_y*d->geo.PositionInCell(m,1) ;
	    z = d->geo.cell__ijk[l][2]*d->geo.lattice_z + d->geo.lattice_z*d->geo.PositionInCell(m,2) ;
	    fprintf(f_bl,"%f %f %f %f %d %d\n",x,y,z,blave,nb,(int) d->geo.AtomType[l][m]);  
	 }
      }
   }
   
   fclose(f_bl);
   
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

   return 0;
}


/* dump eigenvalues to file.  Possible modes:
   'i':  initialize; create a new file and delete existing contents
   'a':  write eigenvalues and append to current file
*/
void writeEigenvalues(qd_struct d, char mode) 
{
   if (!d->opt.ExecParam.ElCalc.ElOut.EnergyDispersion || mpi_n3d_id)
      return;

   MPI_TIC(mpiTiming.fileio_0);

   char* filename = n3d_strdup_n(d->inputfile);
   n3d_FileTypeSet(&filename, "nd_Ek", TRUE);
   FILE* fp = NULL;

   switch (mode) {
   case 'i':
      if ((fp = fopen(filename, "w"))==NULL)
         die("Could not open the file \"%s\"", filename);
      fclose(fp);
      break;
   case 'a':
      FILE* fpr;
      // first check to see if the file exists; if not create it
      if ((fpr = fopen(filename, "r"))==NULL) {
         FILE* fp_new = fopen(filename, "w");
         if (fp_new==NULL)  die("Could not open the file \"%s\"", filename);
         fclose(fp_new);
      }
      fclose(fpr);

      // file exists -- if no header exists, create one
      if ( (fpr = fopen(filename, "r")) ) {
         int dim[MAX_VARS];
         char* dim_lbl[MAX_VARS];
         char* usr_comment[MAX_VARS];
         char dum[1];
         dim_lbl[0] = new char[100];
         usr_comment[0] = new char[100];
         usr_comment[1] = NULL;

         if ( !(fread(dum, sizeof(char), 1, fpr)) ) {
            // add a header
            dim[0] = d->WF.size();
            dim[1] = 0;

            strcpy(dim_lbl[0], "# number of eigenvalues");
            strcpy(usr_comment[0], "# not used");
            writeHeader(filename, "{real kx, real ky, real kz, real E}", dim,
                        dim_lbl, usr_comment);
            
         }
         else {
            int Ndim, Nuser;
            char recFmt[100];
            fp = readHeader(filename, recFmt, &Ndim, dim, dim_lbl, 
                            &Nuser, usr_comment);
            int sz_dat = 4*dim[0];
            double* copydat = new double[sz_dat];
            fread(copydat, sizeof(double), sz_dat, fp);
            fclose(fp);

            dim[0] += d->WF.size();
            writeHeader(filename, "{real kx, real ky, real kz, real E}", dim,
                        dim_lbl, usr_comment);
            
            fp = fopen(filename, "a");
            fwrite(copydat, sizeof(double), sz_dat, fp);
            fclose(fp);

            delete [] copydat;
         }

         delete [] dim_lbl[0];
         delete [] usr_comment[0];
      }
      
      // append data
      fp = fopen(filename, "a");
      for (unsigned int j=0; j < d->WF.size(); j++) {
         double dum[4];
         dum[0] = d->kxL;
         dum[1] = d->kyL;
         dum[2] = d->kzL;
         dum[3] = d->WF[j].E*HAMILTONIAN_SCALE_VALUE;
         fwrite(dum, sizeof(double), 4, fp);
      }
      fclose(fp);
      break;
   default:
      return;
   }
   
   str_free(filename);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);


}



void multicastDisplacement(qd_struct d)
{
   MPI_TIME_INIT(Tstart);
   int j;
   
   if (mpi_n3d_numprocs>1) {

      if (!mpi_n3d_id) printf("  multicasting displacement data\n"), fflush(stdout);

      /* Broadcast the displacement array to all the other CPUs */
      for (j=0; j<mpi_n3d_numprocs; j++) {
         MPI_Barrier(MPI_COMM_WORLD);

#ifdef DISPL_SINGLE_PRECISION
	 MPI_Bcast(&d->geo.l_m_2_atomdsp[d->cell_s[j]][0][0], 
                   d->cell_ln[j]*d->geo.AtomsPerCellMax()*3, 
                   MPI_FLOAT, j, MPI_COMM_WORLD);
#else /* DISPL_SINGLE_PRECISION */
	 MPI_Bcast(&d->geo.l_m_2_atomdsp[d->cell_s[j]][0][0], 
                   d->cell_ln[j]*d->geo.AtomsPerCellMax()*3, 
                   MPI_DOUBLE, j, MPI_COMM_WORLD);
#endif /* DISPL_SINGLE_PRECISION */

         MPI_Barrier(MPI_COMM_WORLD);
      }
   }
}


void readDisplacement( qd_struct d, rvectr p)
{
   int i, I, j, i_ca, i_a;

   MPI_TIC(mpiTiming.fileio_0);

   /* if no strain computation is to be performed, make sure strain file exists */
   if ( d->opt.ExecParam.Strain.StrainModel == Strain_struct::Read ) {
      if ( d->opt.ExecParam.Strain.Keating.StartFileRead == "NULL") {
         if (mpi_n3d_id==mpi_n3d_masterid)
            printf("ERROR:  There is no strain file specified in the input deck!!\n");
         fflush(stdout);
         die("");
      }
      if ( !fopen(d->opt.ExecParam.Strain.Keating.StartFileRead.c_str(), "r") ) {
         if (mpi_n3d_id==mpi_n3d_masterid)
            printf("ERROR: Unable to read displacement file:  '%s'\n",
                   d->opt.ExecParam.Strain.Keating.StartFileRead.c_str());
         die("");
      }
      else 
         masterPrint("Read displacement data from file: %s\n",
                     d->opt.ExecParam.Strain.Keating.StartFileRead.c_str());
   }
   /* otherwise, be fault tolerant if we can't find the file */
   else {
      if ( d->opt.ExecParam.Strain.Keating.StartFileRead == "NULL" ){
         masterPrint("Strain Calculation without initial guess for displacements");
         return;
      }
      else if ( !fopen(d->opt.ExecParam.Strain.Keating.StartFileRead.c_str(), "r") ) {
         masterPrint("ACHTUNG:  Unable to read displacement file: '%s'\nContinuing anyway",
                     d->opt.ExecParam.Strain.Keating.StartFileRead.c_str());
         return;
      }
      else 
         masterPrint("Read displacement data from file: '%s'",
                     d->opt.ExecParam.Strain.Keating.StartFileRead.c_str());
   }
   
   /* p_tmp is indexed by atom, but we need to fill p and 
      d->geo.l_m_2_atomdsp which are indexed by (cell,atom) */
#ifdef DISPL_SINGLE_PRECISION
   fvectr p_tmp = Fvectr(3 * d->seg_ln[mpi_n3d_id] / d->NBasisStates + 3);
#else /* DISPL_SINGLE_PRECISION */
   rvectr p_tmp = Rvectr(3 * d->seg_ln[mpi_n3d_id] / d->NBasisStates + 3);
#endif /*DISPL_SINGLE_PRECISION */
   readSiteInfo(d, "atom", d->opt.ExecParam.Strain.Keating.StartFileRead.c_str(), 
                &p_tmp[0]);

   I = d->cell_s[mpi_n3d_id];
   i_ca=i_a=0;
   for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, I++ ) {
   for ( j=0; j < d->geo.AtomsPerCellMax(); j++){
      if ( d->geo.AtomType[I][j] ){
	 d->geo.l_m_2_atomdsp[I][j][0] = p[i_ca++] = p_tmp[i_a++];
	 d->geo.l_m_2_atomdsp[I][j][1] = p[i_ca++] = p_tmp[i_a++];
	 d->geo.l_m_2_atomdsp[I][j][2] = p[i_ca++] = p_tmp[i_a++];
      }
      else {
	 i_ca += 3;
         
         printf("proc=%d %d %d -- encountered fake atom", mpi_n3d_id, i, j);
      }
      
   }}
   
//   d->geo.lattice_x = p_tmp[i_a++];
//   d->geo.lattice_y = p_tmp[i_a++];
//   d->geo.lattice_z = p_tmp[i_a++];

#ifdef DISPL_SINGLE_PRECISION
   rm_fvectr(&p_tmp);
#else /* DISPL_SINGLE_PRECISION */
   rm_rvectr(&p_tmp);
#endif /* DISPL_SINGLE_PRECISION */
   /* now that each cpu has its local copy of p and atomdsp, broadcast atomdsp
      so that each proc has a full copy.  Note that this may be overkill */
   multicastDisplacement(d);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}


void printHam ( qd_struct d, char* fmt )
{
   /* run this on one CPU only!!! */
   char *filename;
   char ext[200];
   FILE *fp=NULL;
   int atom, nbr, i, j;

   MPI_TIC(mpiTiming.fileio_0);

   strcat(strcpy(ext,"nd_"), fmt);
   filename = n3d_strdup_n(d->inputfile);
   n3d_FileTypeSet(&filename, ext, TRUE);
   
   if (!strcmp(fmt,"hamB")) {
      if (mpi_n3d_id != mpi_n3d_masterid) {
	 die("printHam() should only be executed when running on a SINGLE cpu!!!");
      }
      fp=fopen(filename,"w");
      for (atom=0;atom<vdim1(d->Ham);atom++) {
      for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++) {
	 fprintf(fp,"%d %d %d \n", d->Ham_offset_y[atom][nbr]/d->NBasisStates, atom, 
		 d->Ham_yc[atom][nbr]);
      }}
   }
   else if (!strcmp(fmt,"surf")) {
      int Zc, Za;
       if (mpi_n3d_id != mpi_n3d_masterid) {
	   die("printHam() should only be executed when running on a SINGLE cpu!!!");
       }
      fp=fopen(filename,"w");
      for (Zc=0; Zc < d->geo.N_Cell; Zc++) {
      for (Za=0; Za < d->geo.AtomsPerCellMax(); Za++) {
	 if (d->geo.isOnSurface(Zc,Za) && d->geo.AtomType[Zc][Za] )
	    fprintf(fp,"%d\n", d->geo.l_m_2_offset_ham_tot[Zc][Za]/d->NBasisStates);
      }}
   }
   else if (!strcmp(fmt,"pos")) {
      int Zc, Za;
      if (mpi_n3d_id != mpi_n3d_masterid) {
	 die("printHam() should only be executed when running on a SINGLE cpu!!!");
      }
      if ( !(fp=fopen(filename,"w")) ) die("unable to open file");

      for (Zc=0; Zc < d->geo.N_Cell; Zc++) {
      for (Za=0; Za < d->geo.AtomsPerCellMax(); Za++) {
	 real X, Y, Z;
	 if (!d->geo.AtomType[Zc][Za] ) continue;

	 X = d->geo.lattice_x*(d->geo.cell__ijk[Zc][0] + d->geo.PositionInCell(Za,0)) + 
	    d->geo.l_m_2_atomdsp[Zc][Za][0];
	 Y = d->geo.lattice_y*(d->geo.cell__ijk[Zc][1] + d->geo.PositionInCell(Za,1)) + 
	    d->geo.l_m_2_atomdsp[Zc][Za][1];
	 Z = d->geo.lattice_z*(d->geo.cell__ijk[Zc][2] + d->geo.PositionInCell(Za,2)) + 
	    d->geo.l_m_2_atomdsp[Zc][Za][2];
	 
	 fprintf(fp, "%18.12f %18.12f %18.12f\n", X, Y, Z);
      }}
   }
   else if (!strcmp(fmt,"ham_diag_nbr")) {
      if (mpi_n3d_id != mpi_n3d_masterid) {
	 die("printHam() should only be executed when running on a SINGLE cpu!!!");
      }
      fp=fopen(filename,"w");
      for (atom=0;atom<vdim1(d->Ham);atom++) {
	 nbr=0;
	 
	 /* for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++) { */
	 fprintf(fp,"%d ",d->Ham_neighbor[atom]);
	 for (i=0; i<d->NBasisStates; i++) {
	    fprintf(fp,"%g ",d->Ham[atom][nbr][i][i].r);
	 }
	 fprintf(fp,"\n");
	 /* } */
      }
   }
   else if (!strcmp(fmt,"ham_full")) {
      if (mpi_n3d_id != mpi_n3d_masterid) {
	 char filename2[255];
	 sprintf(filename2,"ham_full_CPU_%d",mpi_n3d_id);
	 fp=fopen(filename2,"w");
      }else{
	 fp=fopen(filename,"w");
      }
       
      for (atom=0;atom<vdim1(d->Ham);atom++) {
      for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++) {
	 for (i=0; i<d->NBasisStates; i++) {
         for (j=0; j<d->NBasisStates; j++) {
	    if (d->Ham[atom][nbr][i][j].r != 0.0 ||
		d->Ham[atom][nbr][i][j].i != 0.0) {

	       if (d->Ham_transp[atom][nbr]) {
		  fprintf(fp,"%d %d %g %g %g\n", 
			  d->Ham_offset_y[atom][nbr]+i+1, d->NBasisStates*atom+j+1, 
			  d->Ham[atom][nbr][j][i].r,
			  -d->Ham[atom][nbr][j][i].i,
			  sqrt(d->Ham[atom][nbr][j][i].r*d->Ham[atom][nbr][j][i].r + d->Ham[atom][nbr][j][i].i*d->Ham[atom][nbr][j][i].i)
		       );
	       }
	       else {
		  fprintf(fp,"%d %d %g %g %g\n", 
			  d->Ham_offset_y[atom][nbr]+i+1, d->NBasisStates*atom+j+1, 
			  d->Ham[atom][nbr][i][j].r,
			  d->Ham[atom][nbr][i][j].i, 
			  sqrt(d->Ham[atom][nbr][i][j].r*d->Ham[atom][nbr][i][j].r + d->Ham[atom][nbr][i][j].i*d->Ham[atom][nbr][i][j].i)
			  );
	       }
	    }
	 }}
      }}
   }
   else if (!strcmp(fmt,"ham")) {
      if (mpi_n3d_id != mpi_n3d_masterid) {
	 die("printHam() should only be executed when running on a SINGLE cpu!!!");
      }
      fp=fopen(filename,"w");
      for (atom=0;atom<vdim1(d->Ham);atom++) {
      for (nbr=0;nbr<=d->Ham_neighbor[atom];nbr++) {
	 for (i=0; i<d->NBasisStates; i++) {
         for (j=0; j<d->NBasisStates; j++) {
	    if (d->Ham[atom][nbr][i][j].r != 0.0 ||
		d->Ham[atom][nbr][i][j].i != 0.0) {

	       if (d->Ham_transp[atom][nbr]) {
		  fprintf(fp,"%d %d %20.12e %20.12e %d\n", 
			  d->Ham_offset_y[atom][nbr]+i+1, d->NBasisStates*atom+j+1, 
			  d->Ham[atom][nbr][j][i].r,
			  -d->Ham[atom][nbr][j][i].i,
			  d->Ham_neighbor[atom]
		       );
	       }
	       else {
		  fprintf(fp,"%d %d %20.12e %20.12e %d\n", 
			  d->Ham_offset_y[atom][nbr]+i+1, d->NBasisStates*atom+j+1, 
			  d->Ham[atom][nbr][i][j].r,
			  d->Ham[atom][nbr][i][j].i,
			  d->Ham_neighbor[atom]
			  );
	       }
	    }
	 }}
      }}
   }
   else {
      die("Unrecognized format.");
      
   }
   
   str_free(filename);

   fflush(fp);
   fclose(fp);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}


void readWaveFunction(char* filename, qd_struct d, cvectr p) 
{
   char band_model[100]; 
   int cryst_str;

   int Ndim, dim[MAX_VARS];
   char recFmt[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
   int Nuser;
   
   MPI_TIC(mpiTiming.fileio_0);

   for (int i=0; i<MAX_VARS; i++) {
      dim[i] = 0;
      dim_lbl[i] = (char*) nml_malloc(sizeof(char)*100);
      usr_comment[i] = (char*) nml_malloc(sizeof(char)*100);
   }

   FILE* fvec = readHeader(filename, recFmt, &Ndim, dim, dim_lbl, &Nuser, 
                           usr_comment);

   /* this n_tot is the total number of atoms; this is modified
      later so that n_tot equals (#atoms)*(#basis states) */
   //const int Nbasis = dim[0];
   const int n_tot = dim[1] * dim[0];
   sscanf(usr_comment[0], "%s", band_model);
   sscanf(usr_comment[1], "%d", &cryst_str);

   if ( d->opt.Dev.band_model != 
        d->opt.Dev.getband_modelType(string(band_model)) 
        || d->n_ham_tot != n_tot ) {
      masterPrint("size(Ham): %d (current) %d (file)", d->n_ham_tot, n_tot);
      die("Cannot continue -- inconsistency between input file and current parameters.");
   }
   
   // each proc just reads its bit from file
   masterPrint("Reading initial guess from '%s'", filename);
   fseek(fvec, sizeof(complex)/sizeof(char)*d->seg_s[mpi_n3d_id], SEEK_CUR);
   fread(p, sizeof(complex), d->seg_ln[mpi_n3d_id], fvec);

   fclose(fvec);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}


int initStartVect(char *filename, qd_struct d, cvectr p)
{
   char band_model[100]; 
   FILE *fvec;
   int i, ii, n_tot, cryst_str, ip, ip2, Nbasis;
   int print_segments=20;
   int print_lengths;

   int Ndim, dim[MAX_VARS];
   char recFmt[100], *dim_lbl[MAX_VARS], *usr_comment[MAX_VARS];
   int Nuser;
   
   MPI_TIC(mpiTiming.fileio_0);

   for (i=0; i<MAX_VARS; i++) {
      dim[i] = 0;
      dim_lbl[i] = (char*) nml_malloc(sizeof(char)*100);
      usr_comment[i] = (char*) nml_malloc(sizeof(char)*100);
   }

   fvec = readHeader(filename, recFmt, &Ndim, dim, dim_lbl, &Nuser, 
                     usr_comment);

   /* this n_tot is the total number of atoms; this is modified
      later so that n_tot equals (#atoms)*(#basis states) */
   Nbasis = dim[0];
   n_tot = dim[1];
   sscanf(usr_comment[0], "%s", band_model);
   sscanf(usr_comment[1], "%d", &cryst_str);

   if ( d->opt.Dev.band_model == 
        d->opt.Dev.getband_modelType(string(band_model)) ) {
      
      n_tot *= Nbasis;
      
      /* same band model and same number of atoms so ...
	 simply read a starting vector that maps directly into the structure */
      fread(p, 2*sizeof(double), n_tot, fvec);

   } else if ( !strcmp(band_model,"Bands_1_s_nospin") && 
	       d->opt.Dev.band_model == Dev_struct::Bands_10_sp3ss_spin &&
	       cryst_str==2 ) {
      
      /* File has been written using the Bands_1_s_nospin, zincblende models and we are 
	 currently using the Bands_10_sp3ss_spin model, so that the file countains 
	 ten times less orbital information than currently required.  
	 Various mapping algorithms can be tried here to achieve better convergence.
      */
      complex dummy;
      real norm;

      ii=0;
      n_tot *= 10;
	
      print_lengths = n_tot/10/print_segments;
      printf("Length of data: %d lines\n",n_tot),fflush(stdout);
      ip=0, ip2=1;
      for ( i=0; i < n_tot/10; i++ ) {

	 fread(&dummy, 2*sizeof(double), 1, fvec);

	 if (ip++ >print_lengths){
	    printf(" %d/%d",ip2,print_segments),fflush(stdout);
	    ip2++;
	    ip=0;
	 }
	
    
#undef UseForHoles2     /* UseForHoles2 seemed to have worked best in the few cases gekco looked at. */	
#ifdef UseForHoles1	
	 p[ii+1].r = p[ii+6].r = dummy.r/2000000.0;       /* s orbitals have small contribution */
	 p[ii+1].i = p[ii+6].i = dummy.i/2000000.0;      
	 
	 p[ii+0].r = p[ii+5].r = dummy.r / 1000000000.0;  /* s* orbitals have almost no contribution */
	 p[ii+0].i = p[ii+5].i = dummy.i / 1000000000.0;
	    
	 /* px, py, pz orbitals have small contributions */
	 p[ii+2].r = dummy.r / 2.0;
	 p[ii+2].i = dummy.i / 2000.0;
	 p[ii+3].r = p[ii+4].r = p[ii+7].r = p[ii+8].r = p[ii+9].r = dummy.r / 100000.0;
	 p[ii+3].i = p[ii+4].i = p[ii+7].i = p[ii+8].i = p[ii+9].i = dummy.i / 100000.0;
#endif	
#ifdef UseForHoles2	
	 p[ii+1].r = p[ii+6].r = dummy.r/2000000.0;       /* s orbitals have small contribution */
	 p[ii+1].i = p[ii+6].i = dummy.i/2000000.0;      
	    
	 p[ii+0].r = p[ii+5].r = dummy.r / 1000000000.0;  /* s* orbitals have almost no contribution */
	 p[ii+0].i = p[ii+5].i = dummy.i / 1000000000.0;
	    
	 /* px, py, pz orbitals have small contributions */
	 p[ii+3].r = dummy.r / 2.0;
	 p[ii+3].i = dummy.i / 2000.0;
	 p[ii+2].r = p[ii+4].r = p[ii+7].r = p[ii+8].r = p[ii+9].r = dummy.r / 100000.0;
	 p[ii+2].i = p[ii+4].i = p[ii+7].i = p[ii+8].i = p[ii+9].i = dummy.i / 100000.0;
#endif	
#ifdef UseForHoles3	
	 p[ii+1].r = p[ii+6].r = dummy.r/2000000.0;       /* s orbitals have small contribution */
	 p[ii+1].i = p[ii+6].i = dummy.i/2000000.0;      
	 
	 p[ii+0].r = p[ii+5].r = dummy.r / 1000000000.0;  /* s* orbitals have almost no contribution */
	 p[ii+0].i = p[ii+5].i = dummy.i / 1000000000.0;
	 
	 /* px, py, pz orbitals have small contributions */
	 p[ii+4].r = dummy.r / 2.0;
	 p[ii+4].i = dummy.i / 2000.0;
	 p[ii+2].r = p[ii+3].r = p[ii+7].r = p[ii+8].r = p[ii+9].r = dummy.r / 100000.0;
	 p[ii+2].i = p[ii+3].i = p[ii+7].i = p[ii+8].i = p[ii+9].i = dummy.i / 100000.0;
#endif	
#ifdef UseForHoles4	
	 p[ii+1].r = p[ii+6].r = dummy.r/2000000.0;       /* s orbitals have small contribution */
	 p[ii+1].i = p[ii+6].i = dummy.i/2000000.0;      
	 
	 p[ii+0].r = p[ii+5].r = dummy.r / 1000000000.0;  /* s* orbitals have almost no contribution */
	 p[ii+0].i = p[ii+5].i = dummy.i / 1000000000.0;
	 
	 /* px, py, pz orbitals have small contributions */
	 p[ii+2].r = p[ii+3].r = p[ii+4].r = p[ii+7].r = p[ii+8].r = p[ii+9].r = dummy.r / 2.0;
	 p[ii+2].i = p[ii+3].i = p[ii+4].i = p[ii+7].i = p[ii+8].i = p[ii+9].i = dummy.i / 100000.0;
#endif	
#ifdef UseForHoles5	
	 p[ii+1].r = p[ii+6].r = dummy.r/2000000.0;       /* s orbitals have small contribution */
	 p[ii+1].i = p[ii+6].i = dummy.i/2000000.0;      
	 
	 p[ii+0].r = p[ii+5].r = dummy.r / 1000000000.0;  /* s* orbitals have almost no contribution */
	 p[ii+0].i = p[ii+5].i = dummy.i / 1000000000.0;
	 
	 /* px, py, pz orbitals have small contributions */
	 p[ii+2].r = p[ii+3].r = p[ii+4].r = dummy.r / 2.0;
	 p[ii+2].i = p[ii+3].i = p[ii+4].i = dummy.i / 100000.0;
	 p[ii+7].r = p[ii+8].r = p[ii+9].r = dummy.r / 100000.0;
	 p[ii+7].i = p[ii+8].i = p[ii+9].i = dummy.i / 2.0;
#endif	
	 
#define UseForElectrons	
#ifdef UseForElectrons	
	 p[ii+6].r = dummy.r/2.0;       /* s orbitals have the major contribution */
	 p[ii+6].i = dummy.i/2000.0;      
	 p[ii+1].r = -dummy.r/2000.0;       /* s orbitals have the major contribution */
	 p[ii+1].i = -dummy.i/2000000.0;      
	 
	 p[ii+0].r = p[ii+5].r = -dummy.r / 1000000.0;  /* s* orbitals have almost no contribution */
	 p[ii+0].i = p[ii+5].i = -dummy.i / 1000000000.0;
	 
	 /* px, py, pz orbitals have small contributions */
	 p[ii+2].r = p[ii+3].r = p[ii+4].r = p[ii+7].r = p[ii+8].r = p[ii+9].r = -dummy.r / 100000.0;
	 p[ii+2].i = p[ii+3].i = p[ii+4].i = p[ii+7].i = p[ii+8].i = p[ii+9].i = -dummy.i / 100000000.0;
#endif	
#if 1
	 p[ii+1].r = p[ii+6].r = dummy.r/2.0;       /* s orbitals have the major contribution */
	 p[ii+1].i = p[ii+6].i = dummy.i/2.0;      
	 
	 p[ii+0].r = p[ii+5].r = dummy.r / 1000.0;  /* s* orbitals have almost no contribution */
	 p[ii+0].i = p[ii+5].i = dummy.i / 1000.0;
	 
	 /* px, py, pz orbitals have small contributions */
	 p[ii+2].r = p[ii+3].r = p[ii+4].r = p[ii+7].r = p[ii+8].r = p[ii+9].r = dummy.r / 100.0;
	 p[ii+2].i = p[ii+3].i = p[ii+4].i = p[ii+7].i = p[ii+8].i = p[ii+9].i = dummy.i / 100.0;
#endif	

	 ii = ii+10;
      }
      printf("\nNormalize the Data\n"),fflush(stdout);
      /* Normalize the input. */
      norm = 0.0;
      for ( i=0; i < n_tot; i++ ){
	 norm += (p[i].r * p[i].r + p[i].i*p[i].i);
      }
      for ( i=0; i < n_tot; i++ ){
	 p[i].r /= norm;
	 p[i].i /= norm;
      }
	 
      /*		for (iii=1;iii<10;iii++){
			ii++;
			p[ii]=p[ii-1	];
			}	*/
   }  else if ( !strcmp(band_model,"Bands_1_s_nospin") && 
		d->opt.Dev.band_model == Dev_struct::Bands_20_sp3d5ss_spin  &&
		cryst_str==2 ) {
      
      /* File has been written using the Bands_1_s_nospin, zincblende models and we are 
	 currently using the Bands_20_sp3d5ss_spin model, so that the file countains 
	 ten times less orbital information than currently required.  
	 Various mapping algorithms can be tried here to achieve better convergence.
      */
	 
      complex dummy;
      real norm;
      ii=0;
      
      n_tot *= 20;
      
      print_lengths = n_tot/20/print_segments;
      printf("Length of data: %d lines\n",n_tot),fflush(stdout);
      ip=0, ip2=1;
      for ( i=0; i < n_tot/20; i++ ) {
	 fread(&dummy, 2*sizeof(double), 1, fvec);
	 if (ip++ >print_lengths){
	    printf(" %d/%d",ip2,print_segments),fflush(stdout);
	    ip2++;
	    ip=0;
	 }
	 
#undef sp3d5_UseForHoles2     /* UseForHoles2 seemed to have worked best in the few cases gekco looked at. */	
#ifdef sp3d5_UseForHoles2	
	 p[ii+1].r = p[ii+11].r = dummy.r/2000000.0;       /* s orbitals have small contribution */
	 p[ii+1].i = p[ii+11].i = dummy.i/2000000.0;      
	 
	 p[ii+0].r = p[ii+10].r = dummy.r / 1000000000.0;  /* s* orbitals have almost no contribution */
	 p[ii+0].i = p[ii+10].i = dummy.i / 1000000000.0;
	 
	 /* 	px, py, pz orbitals have small contributions */
	 p[ii+3].r = dummy.r / 2.0;
	 p[ii+3].i = dummy.i / 2000.0;
	 p[ii+2].r = p[ii+4].r = p[ii+12].r = p[ii+13].r = p[ii+14].r = dummy.r / 100000.0;
	 p[ii+2].i = p[ii+4].i = p[ii+12].i = p[ii+13].i = p[ii+14].i = dummy.i / 100000.0;
	 
	 /* 	d orbitals assume small contribution as well */
	 p[ii+5].r  = p[ii+6].r  =  p[ii+7].r  =  p[ii+8].r  =  p[ii+9].r  = -dummy.r / 100000.0;
	 p[ii+5].i  = p[ii+6].i  =  p[ii+7].i  =  p[ii+8].i  =  p[ii+9].i  = -dummy.i / 100000000.0;
	 p[ii+15].r = p[ii+16].r =  p[ii+17].r =  p[ii+18].r =  p[ii+19].r  = -dummy.r / 100000.0;
	 p[ii+15].i = p[ii+16].i =  p[ii+17].i =  p[ii+18].i =  p[ii+19].i  = -dummy.i / 100000000.0;
#endif		
	    
	    
#define sp3d5_UseForElectrons	
#ifdef sp3d5_UseForElectrons	
	 p[ii+1].r = dummy.r/2.0;           /* s orbitals have the major contribution */
	 p[ii+1].i = dummy.i/2000.0;      
	 p[ii+11].r = -dummy.r/2000.0;       /* s orbitals have the major contribution */
	 p[ii+11].i = -dummy.i/2000000.0;      
	 
	 p[ii+0].r = p[ii+10].r = -dummy.r / 1000000.0;  /* s* orbitals have almost no contribution */
	 p[ii+0].i = p[ii+10].i = -dummy.i / 1000000000.0;
	 
	 /* 	px, py, pz orbitals have small contributions */
	 p[ii+2].r  = p[ii+3].r  = p[ii+4].r  = -dummy.r / 100000.0;
	 p[ii+2].i  = p[ii+3].i  = p[ii+4].i  = -dummy.i / 100000000.0;
	 p[ii+12].r = p[ii+13].r = p[ii+14].r = -dummy.r / 100000.0;
	 p[ii+12].i = p[ii+13].i = p[ii+14].i = -dummy.i / 100000000.0;
	 
	 /* 	d orbitals assume small contribution as well */
	 p[ii+5].r  = p[ii+6].r  =  p[ii+7].r  =  p[ii+8].r  =  p[ii+9].r  = -dummy.r / 100000.0;
	 p[ii+5].i  = p[ii+6].i  =  p[ii+7].i  =  p[ii+8].i  =  p[ii+9].i  = -dummy.i / 100000000.0;
	 p[ii+15].r = p[ii+16].r =  p[ii+17].r =  p[ii+18].r =  p[ii+19].r  = -dummy.r / 100000.0;
	 p[ii+15].i = p[ii+16].i =  p[ii+17].i =  p[ii+18].i =  p[ii+19].i  = -dummy.i / 100000000.0;
	 
#endif		
	    
	 ii = ii+20;
	 
	    /*			for (iii=1;iii<10;iii++){
				ii++;
				p[ii]=p[	ii-1	];
				}	*/
      }
      /* Normalize the input. */
      printf("\nNormalize Data\n"),fflush(stdout);
      norm = 0.0;
      for ( i=0; i < n_tot; i++ ){
	 norm += (p[i].r * p[i].r + p[i].i*p[i].i);
      }
      for ( i=0; i < n_tot; i++ ){
	 p[i].r /= norm;
	 p[i].i /= norm;
      }
      
   } else {
      die("ERROR -- incompatible input file.");
   }
   
   fclose(fvec);    

   for (i=0; i<MAX_VARS; i++) {
      nml_free(dim_lbl[i]);
      nml_free(usr_comment[i]);
   }

   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);
   return 0;


}


/* Prints the matrix structure of the Hamiltonian */
int print_hmat_struct ( qd_struct d )
{
    int i, j, k, l, idum;
    imatrix h;

    h = Imatrix(d->geo.N_Cell, d->geo.N_Cell);


    for ( l=0; l < d->geo.N_Cell; l++ ){
        i = d->geo.cell__ijk[l][0];
        j = d->geo.cell__ijk[l][1];
        k = d->geo.cell__ijk[l][2];

        /* Hamiltonian diagonal */
        h[l][l] = 1;

        /* Off-diagonal terms */
        idum = d->geo.ijk__cell[i-1][j][k];
        if ( idum >= 0 )
            h[l][idum] = 1;

        idum = d->geo.ijk__cell[i+1][j][k];
	if ( idum >= 0 )
            h[l][idum] = 1;

        idum = d->geo.ijk__cell[i][j-1][k];
        if ( idum >= 0 )
            h[l][idum] = 1;

        idum = d->geo.ijk__cell[i][j+1][k]; 
	if ( idum >= 0 )
            h[l][idum] = 1;

        idum = d->geo.ijk__cell[i][j][k-1];
	if ( idum >= 0 )
            h[l][idum] = 1;

        idum = d->geo.ijk__cell[i][j][k+1];
	if ( idum >= 0 )
            h[l][idum] = 1;

    }

    for ( i=0; i < d->geo.N_Cell; i++ )
    {
        printf("\n");
        for ( j=0; j < d->geo.N_Cell; j++ )
            printf("%i", h[i][j]);
    }
    printf("\n");

    rm_imatrix(&h);
    return(0);
}
   

int print_Psi4_IPR( qd_struct d , rvectr wf4, int n, cvectr eigval)
{
    int i;
    FILE *f_wf;
    char label[50];
    char *filename=NULL;
    if (!d->opt.ExecParam.ElCalc.ElOut.Psi4_IPR){
	return 0;
    }
    
    
    if (mpi_n3d_id == mpi_n3d_masterid){
	/* open the file */
	sprintf(label, "nd_wf4_IPR");
	filename = n3d_strdup_n(d->inputfile);
	n3d_FileTypeSet(&filename, label, TRUE);
	f_wf = fopen( filename, "w" );
	str_free(filename);
	
	/*	 write it's own data. */
	for ( i=0; i < n; i++ ){
	    fprintf(f_wf, "%g %g\n",eigval[i].r,wf4[i]);
	}   
	
	
	fclose(f_wf);
	
    } 
    
    return 0;
}


int print_struct_gmv( qd_struct d )
{
   int i, j, l, na, nb, id, jd, kd, nnbr, ncells; 
   real x, y, z, xmin, xmax, ymin, ymax, zmin, zmax, wf;
   FILE *f_gmv;
   char *filename=NULL;

   if (!d->opt.ExecParam.ElCalc.ElOut.StructureGmv || mpi_n3d_id!=mpi_n3d_masterid)
      return 0;

   MPI_TIC(mpiTiming.fileio_0);

   filename = n3d_strdup_n(d->inputfile);
   n3d_FileTypeSet(&filename, "nd_gmv", TRUE);

   f_gmv = fopen( filename, "w" );
   str_free(filename);
    
   xmin=ymin=zmin=1e100;
   xmax=ymax=zmax=-1e100;

   fprintf(f_gmv,"gmvinput ascii\n\n");
   fprintf(f_gmv,"nodes %i\n",d->n_atom_tot);
   
   /* First print out nodes */
   fprintf(f_gmv,"  ");
   for ( i=0; i < d->geo.N_Cell; i++ ) {
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 
	 /* Check to see if atom is in the structure */ 
	 if ( d->geo.AtomType[i][j] > 0 ) {
	    
	    /* Get location for atom on ith unit cell at
	       jth site on unit cell */
	    x = 10*((real)d->geo.cell__ijk[i][0]*d->geo.lattice_x + 
		    d->geo.lattice_x*d->geo.PositionInCell(j,0)) + d->geo.l_m_2_atomdsp[i][j][0];
	    fprintf(f_gmv," %1.3f",x);
	 }
      }
   }
   fprintf(f_gmv,"\n");
   
   fprintf(f_gmv,"  ");
   for ( i=0; i < d->geo.N_Cell; i++ ) {
      
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 
	 /* Check to see if atom is in the structure */ 
	 if ( d->geo.AtomType[i][j] > 0 ) {
	    
	    /* Get location for atom on ith unit cell at
	       jth site on unit cell */
	    y = 10*((real)d->geo.cell__ijk[i][1]*d->geo.lattice_y + 
		    d->geo.lattice_y*d->geo.PositionInCell(j,1)) + d->geo.l_m_2_atomdsp[i][j][1];
	    fprintf(f_gmv," %1.3f",y);
	 }
      }
   }
   fprintf(f_gmv,"\n");
   
   
   fprintf(f_gmv,"  ");
   for ( i=0; i < d->geo.N_Cell; i++ ) {
      
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 
	 /* Check to see if atom is in the structure */ 
	 if ( d->geo.AtomType[i][j] > 0 ) {
	    
	    /* Get location for atom on ith unit cell at
	       jth site on unit cell */
	    z = 10*((real)d->geo.cell__ijk[i][2]*d->geo.lattice_x + 
		    d->geo.lattice_z*d->geo.PositionInCell(j,2)) + d->geo.l_m_2_atomdsp[i][j][2];
	    fprintf(f_gmv," %1.3f",z);
	 }
      }
   }
   fprintf(f_gmv,"\n\n");
   
#define DO_CELLS	
#ifdef DO_CELLS
   /* Get total number of cells */
   ncells=0;
   for ( i=0; i < d->geo.N_Cell; i++ )
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 
	 if ( d->geo.AtomType[i][j] > 0 ) {
	    
	    na = i*d->geo.AtomsPerCellMax() + j + 1;
	    x = 10*((real)d->geo.cell__ijk[i][0]*d->geo.lattice_x + 
		    d->geo.lattice_x*d->geo.PositionInCell(j,0)) + d->geo.l_m_2_atomdsp[i][j][0];
	    y = 10*((real)d->geo.cell__ijk[i][1]*d->geo.lattice_y + 
		    d->geo.lattice_y*d->geo.PositionInCell(j,1)) + d->geo.l_m_2_atomdsp[i][j][1];
	    z = 10*((real)d->geo.cell__ijk[i][2]*d->geo.lattice_z + 
		    d->geo.lattice_z*d->geo.PositionInCell(j,2)) + d->geo.l_m_2_atomdsp[i][j][2];
	    if ( x > xmax )
	       xmax = x;
	    if ( x < xmin )
	       xmin = x;
	    if ( y > ymax )
	       ymax = y;
	    if ( y < ymin )
	       ymin = y;
	    if ( z > zmax )
	       zmax = z;
	    if ( z < zmin )
	       zmin = z;
	    nnbr=0;
	    /* First get number of bonded atoms */
	    for ( l=0; l < d->geo.Neighbors(j); l++ ) {
	       
	       id = d->geo.cell__ijk[i][0] + d->geo.NbrCell(j,l,0);
	       jd = d->geo.cell__ijk[i][1] + d->geo.NbrCell(j,l,1);
	       kd = d->geo.cell__ijk[i][2] + d->geo.NbrCell(j,l,2);
	       if ( d->geo.ijk__cell[id][jd][kd] >= 0 ) {
		  
		  if ( d->geo.AtomType[d->geo.ijk__cell[id][jd][kd]][d->geo.NbrCell(j,l,3)] > 0 )
		     nnbr++;
	       }
	    }
	    if ( nnbr >= 2 && nnbr <= 4 )
	       ncells++;
	 }
      }
   
   fprintf(f_gmv,"\ncells %i\n",ncells);
   
   /* Now print cell types and vertices */    
   for ( i=0; i < d->geo.N_Cell; i++ )
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 
	 if ( d->geo.AtomType[i][j] > 0 ) {
	    
	    na = i*d->geo.AtomsPerCellMax() + j + 1;
	    x = 10*((real)d->geo.cell__ijk[i][0]*d->geo.lattice_x + 
		    d->geo.lattice_x*d->geo.PositionInCell(j,0)) + d->geo.l_m_2_atomdsp[i][j][0];
	    y = 10*((real)d->geo.cell__ijk[i][1]*d->geo.lattice_y + 
		    d->geo.lattice_y*d->geo.PositionInCell(j,1)) + d->geo.l_m_2_atomdsp[i][j][1];
	    z = 10*((real)d->geo.cell__ijk[i][2]*d->geo.lattice_z + 
		    d->geo.lattice_z*d->geo.PositionInCell(j,2)) + d->geo.l_m_2_atomdsp[i][j][2];
	    nnbr=0;
	    /* First get number of bonded atoms */
	    for ( l=0; l < d->geo.Neighbors(j); l++ ) {
	       
	       id = d->geo.cell__ijk[i][0] + d->geo.NbrCell(j,l,0);
	       jd = d->geo.cell__ijk[i][1] + d->geo.NbrCell(j,l,1);
	       kd = d->geo.cell__ijk[i][2] + d->geo.NbrCell(j,l,2);
	       if ( d->geo.ijk__cell[id][jd][kd] >= 0 ) {
		  
		  if ( d->geo.AtomType[d->geo.ijk__cell[id][jd][kd]][d->geo.NbrCell(j,l,3)] > 0 )
		     nnbr++;
	       }
	    }
	    /* Print out cell indices */
	    switch(nnbr) {
	    case(2): fprintf(f_gmv,"line 2\n");
	       break;
	    case(3): fprintf(f_gmv,"tri 3\n");
	       break;
	    case(4): fprintf(f_gmv,"tet 4\n");
	    }
	    if ( nnbr >= 2 && nnbr <= 4 ) {
	       
	       for ( l=0; l < d->geo.Neighbors(j); l++ ) {
                  
		  id = d->geo.cell__ijk[i][0] + d->geo.NbrCell(j,l,0);
		  jd = d->geo.cell__ijk[i][1] + d->geo.NbrCell(j,l,1);
		  kd = d->geo.cell__ijk[i][2] + d->geo.NbrCell(j,l,2);
		  if ( d->geo.ijk__cell[id][jd][kd] >= 0 ) {
		     
		     if ( d->geo.AtomType[d->geo.ijk__cell[id][jd][kd]][d->geo.NbrCell(j,l,3)] > 0 ) {
			
			nb = d->geo.ijk__cell[id][jd][kd]*d->geo.AtomsPerCellMax() + 
			     d->geo.NbrCell(j,l,3) + 1;
			fprintf(f_gmv, " %i",nb);
		     }
		  }
	       }
	       fprintf(f_gmv,"\n");
	    }
	 }
      }
#else
   fprintf(f_gmv,"\ncells 0\n");
#endif
   
   fprintf(f_gmv,"\nvariable\n");
   fprintf(f_gmv,"\nwavefunc 1\n");
   fprintf(f_gmv,"  ");
   for ( i=0; i < d->geo.N_Cell; i++ ) {
      
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 
	 /* Check to see if atom is in the structure */ 
	 if ( d->geo.AtomType[i][j] > 0 ) {
	    
	    /* Get location for atom on ith unit cell at
	       jth site on unit cell */
	    x = 10*((real)d->geo.cell__ijk[i][0]*d->geo.lattice_x + 
		    d->geo.lattice_x*d->geo.PositionInCell(j,0)) + d->geo.l_m_2_atomdsp[i][j][0];
	    y = 10*((real)d->geo.cell__ijk[i][1]*d->geo.lattice_y + 
		    d->geo.lattice_y*d->geo.PositionInCell(j,1)) + d->geo.l_m_2_atomdsp[i][j][1];
	    z = 10*((real)d->geo.cell__ijk[i][2]*d->geo.lattice_z + 
		    d->geo.lattice_z*d->geo.PositionInCell(j,2)) + d->geo.l_m_2_atomdsp[i][j][2];
	    wf = sin(3.14159*(x-xmin)/(xmax-xmin))*
	       sin(3.14159*(y-ymin)/(ymax-ymin))*
	       sin(3.14159*(z-zmin)/(zmax-zmin));
	    fprintf(f_gmv," %1.3f",wf);
	 }
      }
   }
   fprintf(f_gmv,"\n");     
   fprintf(f_gmv,"\nendvars\n");
   fprintf(f_gmv,"\nendgmv\n");
   fclose(f_gmv);      
   
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);


   return 0;
}

int print_struct_pdb( qd_struct d )
{
   int i, j, l, na=-1, nb, id, jd, kd; 
   real x, y, z;
   /* enum {None,B,C,N,O,F,Al,Si,P,Zn,Ga,Ge,As,Cd,In,Sn,Sb,H,S}; */
   FILE *f_pdb;
   char *filename=NULL;
   
   if (!d->opt.ExecParam.ElCalc.ElOut.StructurePdb || mpi_n3d_id!=mpi_n3d_masterid)
      return 0;
   
   MPI_TIC(mpiTiming.fileio_0);

   filename = n3d_strdup_n(d->inputfile);
   n3d_FileTypeSet(&filename, "nd_pdb", TRUE);
   f_pdb = fopen( filename, "w" );
   str_free(filename);
   
   fprintf(f_pdb,"ORIGX1      %1.6f  %1.6f  %1.6f        %1.5f\n",1.0,0.0,0.0,0.0);
   fprintf(f_pdb,"ORIGX2      %1.6f  %1.6f  %1.6f        %1.5f\n",0.0,1.0,0.0,0.0);
   fprintf(f_pdb,"ORIGX3      %1.6f  %1.6f  %1.6f        %1.5f\n",0.0,0.0,1.0,0.0);
   fprintf(f_pdb,"SCALE1      %1.6f  %1.6f  %1.6f        %1.5f\n",.1,0.0,0.0,0.0);
   fprintf(f_pdb,"SCALE2      %1.6f  %1.6f  %1.6f        %1.5f\n",0.0,.1,0.0,0.0);
   fprintf(f_pdb,"SCALE2      %1.6f  %1.6f  %1.6f        %1.5f\n",0.0,0.0,.1,0.0);
   
   /* First print out atom locations */
   for ( i=0; i < d->geo.N_Cell; i++ ) {
      
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 
	 /* Check to see if atom is in the structure */ 
	 if ( d->geo.AtomType[i][j] > 0 ) {
	    
	    /* Get location for atom on ith unit cell at
	       jth site on unit cell */
	    na = i*d->geo.AtomsPerCellMax() + j +1;
	    x = 10*((real)d->geo.cell__ijk[i][0]*d->geo.lattice_x + 
		    d->geo.lattice_x*d->geo.PositionInCell(j,0)) + d->geo.l_m_2_atomdsp[i][j][0];
	    y = 10*((real)d->geo.cell__ijk[i][1]*d->geo.lattice_y + 
		    d->geo.lattice_y*d->geo.PositionInCell(j,1)) + d->geo.l_m_2_atomdsp[i][j][1];
	    z = 10*((real)d->geo.cell__ijk[i][2]*d->geo.lattice_z + 
		    d->geo.lattice_z*d->geo.PositionInCell(j,2)) + d->geo.l_m_2_atomdsp[i][j][2];	  
	    /* the output of multimillion atom simulations requires the ability to look at 
	     * more than 5 digits of integer numbers.
	     * The original formatting required these limitations to 5 digits.
	     *
	     * Now we have a macro that can be used to change the formatting length. */

	    switch(d->geo.AtomType[i][j]) {
	    case(EL_None): 
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM"              1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_B):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" B            1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_C):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" C            1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;		  
	    case(EL_N):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" N            1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_O):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" O            1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_F):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" F            1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_Al):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" Al           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_Si):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" Si           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_P):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" P            1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;		    
	    case(EL_Zn):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" Zn           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_Ga):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" Ga           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_Ge):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" Ge           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_As):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" As           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_Cd):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" Cd           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;		    
	    case(EL_In):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" In           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_Sn):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" Sn           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_Sb):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" Sb           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	       break;
	    case(EL_H):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" H            1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);	
	       break;
	    case(EL_S):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" S            1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	    case(EL_Se):
	       fprintf(f_pdb, "ATOM  "PDBFORMATATOM" Se           1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	    case(EL_EXPLICIT):
	    case(EL_EXPLICIT+1):
	    case(EL_EXPLICIT+2):
	    case(EL_EXPLICIT+3):
	    case(EL_EXPLICIT+4):
	       fprintf(f_pdb, "ATOM  "FORMATOUT" S            1    %8.3f%8.3f%8.3f"COORDINATEWRITE"\n",na,x,y,z);
	    }
	    /* if ( d->geo.AtomType[i][j] == 12 )
	       fprintf(f_pdb, "ATOM  "FORMATOUT" As           1    %8.3f%8.3f%8.3f  "COORDINATEWRITE"\n",na,x,y,z);
	       else if ( d->geo.AtomType[i][j] == 10 )
	       fprintf(f_pdb, "ATOM  "FORMATOUT" Ga           1    %8.3f%8.3f%8.3f  "COORDINATEWRITE"\n",na,x,y,z);
	       else if ( d->geo.AtomType[i][j] == 14 )
	       fprintf(f_pdb, "ATOM  "FORMATOUT" In           1    %8.3f%8.3f%8.3f  "COORDINATEWRITE"\n",na,x,y,z); */
	 }
      }
   }
/*   fprintf(f_pdb, "TER   "FORMATOUT"\n",na+1); */ /* Commented to make outputfile readable in rasmol */ 
   
   /* Next print out connections */
   for ( i=0; i < d->geo.N_Cell; i++ )
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 
	 if ( d->geo.AtomType[i][j] > 0 ) {
	    
	    na = i*d->geo.AtomsPerCellMax() + j + 1;
	    x = 10*((real)d->geo.cell__ijk[i][0]*d->geo.lattice_x + 
		    d->geo.lattice_x*d->geo.PositionInCell(j,0)) + d->geo.l_m_2_atomdsp[i][j][0];
	    y = 10*((real)d->geo.cell__ijk[i][1]*d->geo.lattice_y + 
		    d->geo.lattice_y*d->geo.PositionInCell(j,1)) + d->geo.l_m_2_atomdsp[i][j][1];
	    z = 10*((real)d->geo.cell__ijk[i][2]*d->geo.lattice_z + 
		    d->geo.lattice_z*d->geo.PositionInCell(j,2)) + d->geo.l_m_2_atomdsp[i][j][2];
	    for ( l=0; l < d->geo.Neighbors(j); l++ ) {
	       
	       id = d->geo.cell__ijk[i][0] + d->geo.NbrCell(j,l,0);
	       jd = d->geo.cell__ijk[i][1] + d->geo.NbrCell(j,l,1);
	       kd = d->geo.cell__ijk[i][2] + d->geo.NbrCell(j,l,2);
	       if ( d->geo.ijk__cell[id][jd][kd] >= 0 ) {
		  
		  if ( d->geo.AtomType[d->geo.ijk__cell[id][jd][kd]][d->geo.NbrCell(j,l,3)] > 0 ) {
		     
		     nb = d->geo.ijk__cell[id][jd][kd]*d->geo.AtomsPerCellMax() + 
			d->geo.NbrCell(j,l,3) + 1;
		     fprintf(f_pdb, "CONECT "PDBFORMATCONCT" "PDBFORMATCONCT"\n",na,nb);
		  }
	       }
	    }
	 }
      }
   fprintf(f_pdb,"END\n");
   fclose(f_pdb);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

   return(0);
}


int print_isosurf_pdb( qd_struct d, rmatrix data, real isoval, int eigvallabel )
{
   int i, j, l, na=-1, nb, id, jd, kd; 
   real x, y, z;
   /* enum {None,B,C,N,O,F,Al,Si,P,Zn,Ga,Ge,As,Cd,In,Sn,Sb,H,S}; */
   FILE *f_pdb;
   char label[30];
   char *filename=NULL;
   
   if (!d->opt.ExecParam.ElCalc.ElOut.PsiSqrPdb || mpi_n3d_id!=mpi_n3d_masterid)
      return 0;
   
   MPI_TIC(mpiTiming.fileio_0);

   filename = n3d_strdup_n(d->inputfile);
   sprintf(label,"nd_wf_%d_pdb",eigvallabel+1);
   n3d_FileTypeSet(&filename, label, TRUE);
   f_pdb = fopen(  filename, "w" );
   str_free(filename);
   
   fprintf(f_pdb,"ORIGX1      %1.6f  %1.6f  %1.6f        %1.5f\n",1.0,0.0,0.0,0.0);
   fprintf(f_pdb,"ORIGX2      %1.6f  %1.6f  %1.6f        %1.5f\n",0.0,1.0,0.0,0.0);
   fprintf(f_pdb,"ORIGX3      %1.6f  %1.6f  %1.6f        %1.5f\n",0.0,0.0,1.0,0.0);
   fprintf(f_pdb,"SCALE1      %1.6f  %1.6f  %1.6f        %1.5f\n",.1,0.0,0.0,0.0);
   fprintf(f_pdb,"SCALE2      %1.6f  %1.6f  %1.6f        %1.5f\n",0.0,.1,0.0,0.0);
   fprintf(f_pdb,"SCALE2      %1.6f  %1.6f  %1.6f        %1.5f\n",0.0,0.0,.1,0.0);
   
   /* First print out atom locations */
   for ( i=0; i < d->geo.N_Cell; i++ ) {
      
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 
	 /* Check to see if atom is in the structure */ 
	 if ( d->geo.AtomType[i][j] > 0 ) {
	    
	    /* Get location for atom on ith unit cell at
	       jth site on unit cell */
	    if ( data[i][j] >= isoval ) {
	       
	       na = i*d->geo.AtomsPerCellMax() + j +1;
	       x = 10*((real)d->geo.cell__ijk[i][0]*d->geo.lattice_x + 
		       d->geo.lattice_x*d->geo.PositionInCell(j,0)) + d->geo.l_m_2_atomdsp[i][j][0];
	       y = 10*((real)d->geo.cell__ijk[i][1]*d->geo.lattice_y + 
		       d->geo.lattice_y*d->geo.PositionInCell(j,1)) + d->geo.l_m_2_atomdsp[i][j][1];
	       z = 10*((real)d->geo.cell__ijk[i][2]*d->geo.lattice_z + 
		       d->geo.lattice_z*d->geo.PositionInCell(j,2)) + d->geo.l_m_2_atomdsp[i][j][2];	  
	       switch(d->geo.AtomType[i][j]) {
	       case(EL_None): 
		  fprintf(f_pdb, "ATOM  "FORMATOUT"              1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_B):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" B            1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_C):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" C            1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;		  
	       case(EL_N):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" N            1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_O):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" O            1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_F):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" F            1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_Al):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" Al           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_Si):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" Si           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_P):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" P            1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;		    
	       case(EL_Zn):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" Zn           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_Ga):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" Ga           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_Ge):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" Ge           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_As):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" As           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_Cd):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" Cd           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;		    
	       case(EL_In):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" In           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_Sn):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" Sn           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_Sb): 
		  fprintf(f_pdb, "ATOM  "FORMATOUT" Sb           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
		  break;
	       case(EL_H):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" H            1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);	
		  break;
	       case(EL_S):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" S            1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
	       case(EL_Se):
		  fprintf(f_pdb, "ATOM  "FORMATOUT" Se           1    %8.3f%8.3f%8.3f  1.00  0.00\n",na,x,y,z);
	       }
	    }
	 }
      }
   }
   fprintf(f_pdb, "TER   "FORMATOUT"\n",na+1);
   
   /* Next print out connections */
   for ( i=0; i < d->geo.N_Cell; i++ )
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 
	 if ( d->geo.AtomType[i][j] > 0 ) {
	    
	    /* if ( data[i][j] >= isoval )
	       { */
	    na = i*d->geo.AtomsPerCellMax() + j + 1;
	    x = 10*((real)d->geo.cell__ijk[i][0]*d->geo.lattice_x + 
		    d->geo.lattice_x*d->geo.PositionInCell(j,0)) + d->geo.l_m_2_atomdsp[i][j][0];
	    y = 10*((real)d->geo.cell__ijk[i][1]*d->geo.lattice_y + 
		    d->geo.lattice_y*d->geo.PositionInCell(j,1)) + d->geo.l_m_2_atomdsp[i][j][1];
	    z = 10*((real)d->geo.cell__ijk[i][2]*d->geo.lattice_z + 
		    d->geo.lattice_z*d->geo.PositionInCell(j,2)) + d->geo.l_m_2_atomdsp[i][j][2];
	    for ( l=0; l < d->geo.Neighbors(j); l++ ) {
	       
	       id = d->geo.cell__ijk[i][0] + d->geo.NbrCell(j,l,0);
	       jd = d->geo.cell__ijk[i][1] + d->geo.NbrCell(j,l,1);
	       kd = d->geo.cell__ijk[i][2] + d->geo.NbrCell(j,l,2);
	       if ( d->geo.ijk__cell[id][jd][kd] >= 0 ) {
		  
		  if ( d->geo.AtomType[d->geo.ijk__cell[id][jd][kd]][d->geo.NbrCell(j,l,3)] > 0 ) {
		     
		     /* if ( data[d->geo.ijk__cell[id][jd][kd]][d->geo.NbrCell(j,l,3)] >= isoval )
			{ */
		     nb = d->geo.ijk__cell[id][jd][kd]*d->geo.AtomsPerCellMax() + 
			d->geo.NbrCell(j,l,3) + 1;
		     fprintf(f_pdb, "CONECT "FORMATOUT" "FORMATOUT"\n",na,nb);
		     /* } */
		  }
	       }
	    }
	    /*  } */
	 }
      }
   fprintf(f_pdb,"END\n");
   fclose(f_pdb);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

   return(0);
}


#define DELIMITER_HC " \n"


int scatter_print_all(qd_struct d, rmatrix psisq, FILE *f_out,int isx, int iex)
{
   int i,j;
   real xa=0.0,ya=0.0,za=0.0;
   
   for ( i=isx; i < iex; i++ ) {
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 if ( d->geo.AtomType[i][j] ) {
	    xa = d->geo.lattice_x*d->geo.cell__ijk[i][0] + d->geo.lattice_x*d->geo.PositionInCell(j,0);
	    ya = d->geo.lattice_y*d->geo.cell__ijk[i][1] + d->geo.lattice_y*d->geo.PositionInCell(j,1);
	    za = d->geo.lattice_z*d->geo.cell__ijk[i][2] + d->geo.lattice_z*d->geo.PositionInCell(j,2);
	    fprintf(f_out, "%21.12e %21.12e %21.12e", xa, ya, za );
	    fprintf(f_out, "%21.12e", psisq[i-isx][j]); 
	    fprintf(f_out,"\n");
	 }
      }
   }
   
   return 0;
}


int scatter_print_unitcell(qd_struct d, rmatrix psisq, FILE *f_out,int isx, int iex)
{
   int i,j;
   real xa=0.0,ya=0.0,za=0.0;
   real dummy;
   for ( i=isx; i < iex; i++ ) {
      j=0;
      xa = d->geo.lattice_x*d->geo.cell__ijk[i][0] + d->geo.lattice_x*d->geo.PositionInCell(j,0);
      ya = d->geo.lattice_y*d->geo.cell__ijk[i][1] + d->geo.lattice_y*d->geo.PositionInCell(j,1);
      za = d->geo.lattice_z*d->geo.cell__ijk[i][2] + d->geo.lattice_z*d->geo.PositionInCell(j,2);
      dummy=0.0;
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
	 if ( d->geo.AtomType[i][j] ) {
	    dummy += psisq[i-isx][j];
	 }
      }
      fprintf(f_out, "%21.12e %21.12e %21.12e %21.12e", xa, ya, za , dummy);
      /* fprintf(f_out, "%21.12e\n", dummy); */
      fprintf(f_out,"\n");
   }
   
   return 0;
}


int print_psisqr_gnuplot( qd_struct d, rmatrix psisq, int k)
{
   MPI_TIME_INIT(Tstart);

   FILE *f_all=NULL, *f_unit=NULL;
   int isx = d->cell_s[mpi_n3d_id];
   int iex = isx + d->cell_ln[mpi_n3d_id];
   int proc;
   rmatrix psisqr_others=NULL;
   char *filename=NULL;
   char label[40];
   
   if (!d->opt.ExecParam.ElCalc.ElOut.PsiSqr && !d->opt.ExecParam.ElCalc.ElOut.PsiSqr_cell)
      return 0;
   
   MPI_TIC(mpiTiming.fileio_0);

   /* On a cluster we do not want to append into the same file 
      from different CPUs due to latency issues in the file system.
      We therefore send the data to the master CPU and that CPU 
      takes care of all the writing of the data.  */
   if (mpi_n3d_id==mpi_n3d_masterid){
       if (d->opt.ExecParam.ElCalc.ElOut.PsiSqr){
	   sprintf(label, "nd_wf_%d_gnu", k);
	   filename = n3d_strdup_n(d->inputfile);
	   n3d_FileTypeSet(&filename, label, TRUE);
	   f_all = fopen( filename, "w" );
	   str_free(filename);
	   scatter_print_all(d, psisq, f_all,isx, iex);
       }
       
	
       if (d->opt.ExecParam.ElCalc.ElOut.PsiSqr_cell){
	   sprintf(label, "nd_wf_%d_gnu_unit", k);
	   filename = n3d_strdup_n(d->inputfile);
	   n3d_FileTypeSet(&filename, label, TRUE);
	   f_unit = fopen( filename, "w" );
	   str_free(filename);
	   scatter_print_unitcell(d, psisq, f_unit,isx, iex);
       }
       
	
      /* receive the data from the other CPUs and write it out */
      for (proc=0;proc<mpi_n3d_numprocs;proc++){
	 if (proc!=mpi_n3d_masterid){
	    isx = d->cell_s[proc];
	    iex = isx + d->cell_ln[proc];

	    /* check that we have enough memory allocated in the matrix */
	    psisqr_others = check_rmatrix(psisqr_others,iex-isx,d->geo.AtomsPerCellMax());

	    MPI_Recv(&psisqr_others[0][0],(iex-isx)*d->geo.AtomsPerCellMax(),
                     MPI_DOUBLE, proc,MPI_ANY_TAG,MPI_COMM_WORLD, &mpi_n3d_status);

            if (d->opt.ExecParam.ElCalc.ElOut.PsiSqr) scatter_print_all(d, psisqr_others, f_all,isx, iex);
	    if (d->opt.ExecParam.ElCalc.ElOut.PsiSqr_cell) scatter_print_unitcell(d, psisqr_others, f_unit,isx, iex);
	 }
      }
      if (d->opt.ExecParam.ElCalc.ElOut.PsiSqr) fclose(f_all);
      if (d->opt.ExecParam.ElCalc.ElOut.PsiSqr_cell) fclose(f_unit);
      rm_rmatrix(&psisqr_others);
   } else {

      /* Send the data to the master cpu. */
      MPI_Send(&psisq[0][0],(iex-isx)*d->geo.AtomsPerCellMax(),MPI_DOUBLE,
	       mpi_n3d_masterid,mpi_n3d_id,MPI_COMM_WORLD);
   }

   MPI_Barrier(MPI_COMM_WORLD);

   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

   return 0;
}

void print_psisq_gnu_pdb_k(cmatrix wfc,int k,qd_struct d,int n_offset)
{
   MPI_TIME_INIT(Tstart);

   int i,j,l,m, isx, iex, imin, imax;
   real wfsum=0.0, wfmax, wfmax_local;
   //rmatrix psisq = Rmatrix(d->seg_ln[mpi_n3d_id],d->geo.AtomsPerCellMax());
   rmatrix psisq = Rmatrix(d->cell_ln[mpi_n3d_id],d->geo.AtomsPerCellMax());

   isx = d->cell_s[mpi_n3d_id];
   iex = isx + d->cell_ln[mpi_n3d_id];
   imin = d->seg_s[mpi_n3d_id];
   imax = d->seg_s[mpi_n3d_id] + d->seg_ln[mpi_n3d_id];

   printf("Generate psisqr\n"),fflush(stdout);
   /* Normalize and print the square of the wavefunctions */
      wfmax = -1e100;
      for ( i=isx; i < iex; i++ ){
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ){
            if ( d->geo.AtomType[i][j] ){
               wfsum = 0;
               for ( l=0; l < d->NBasisStates; l++ ){
                  m = d->geo.l_m_2_offset_ham_tot[i][j]+l-imin;
                  wfsum += wfc[k][m].r*wfc[k][m].r +  wfc[k][m].i*wfc[k][m].i;
               }
            }
            psisq[i-isx][j] = wfsum;
            if ( wfsum > wfmax )
               wfmax = wfsum;
         }
      }
      wfmax_local= wfmax;

#if (defined MPI3d && !defined FAKE_MPI)
    MPI_Allreduce(&wfmax_local,&wfmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
    wfmax = wfmax_local;
#endif

    for ( i=isx; i < iex; i++ )
       for ( j=0; j < d->geo.AtomsPerCellMax(); j++ )
          if ( d->geo.AtomType[i][j] ) {
             psisq[i-isx][j] /= wfmax;
          }

    /* print_psisqr_gnuplot( d, psisq, k+n_offset+1); */
    /* print_isosurf_pdb( d, psisq, d->opt.ExecParam.ElCalc.ElOut.PsiSqrCutoffPdb , k+n_offset+1);  */

   rm_rmatrix(&psisq);

   return;
}


void print_psisq_gnu_pdb(cmatrix wfc,int n,qd_struct d,int n_offset)
{
   MPI_TIME_INIT(Tstart);

   int i,j,k,l,m, isx, iex, imin, imax;
   real wfsum=0.0, wfmax, wfmax_local;
   //rmatrix psisq = Rmatrix(d->seg_ln[mpi_n3d_id],d->geo.AtomsPerCellMax());
   rmatrix psisq = Rmatrix(d->cell_ln[mpi_n3d_id],d->geo.AtomsPerCellMax());

   isx = d->cell_s[mpi_n3d_id];
   iex = isx + d->cell_ln[mpi_n3d_id];
   imin = d->seg_s[mpi_n3d_id];
   imax = d->seg_s[mpi_n3d_id] + d->seg_ln[mpi_n3d_id];

   if (mpi_n3d_id==mpi_n3d_masterid) printf("Generate psisqr\n"),fflush(stdout);
   /* Normalize and print the square of the wavefunctions */
   for ( k=0; k < n; k++ ){
      wfmax = -1e100;
      for ( i=isx; i < iex; i++ ){
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ){
            if ( d->geo.AtomType[i][j] ){
               wfsum = 0;
               for ( l=0; l < d->NBasisStates; l++ ){
                  m = d->geo.l_m_2_offset_ham_tot[i][j]+l-imin;
                  wfsum += wfc[k][m].r*wfc[k][m].r +  wfc[k][m].i*wfc[k][m].i;
               }
            }
            psisq[i-isx][j] = wfsum;
            if ( wfsum > wfmax )
               wfmax = wfsum;
         }
      }
      wfmax_local= wfmax;

#if (defined MPI3d && !defined FAKE_MPI)
      MPI_Allreduce(&wfmax_local,&wfmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
      wfmax = wfmax_local;
#endif

      for ( i=isx; i < iex; i++ )
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ )
            if ( d->geo.AtomType[i][j] ) {
               psisq[i-isx][j] /= wfmax;
            }

     /*  print_psisqr_gnuplot( d, psisq, k+n_offset+1); */
    /* print_isosurf_pdb( d, psisq, d->opt.ExecParam.ElCalc.ElOut.PsiSqrCutoffPdb , k+n_offset+1);  */
   }

   rm_rmatrix(&psisq);

   return;
}


#undef CHECK_COMMUNICATION
void write_displacement( qd_struct d)
{
   int i, j, k;
   FILE *fp=NULL;
   FILE *fp2=NULL;
   char *filename=NULL;
   real xa, ya, za;
    
#ifdef CHECK_COMMUNICATION
   char cpu[4];
   sprintf(cpu,"_%d",mpi_n3d_id);
   if ((!d->opt.ExecParam.Strain.StrainOut.DisplacementData && !d->opt.ExecParam.Strain.StrainOut.DisplacementData3D) || 
       mpi_n3d_id!=mpi_n3d_masterid)
      return;
#else
   if ((!d->opt.ExecParam.Strain.StrainOut.DisplacementData && !d->opt.ExecParam.Strain.StrainOut.DisplacementData3D) || 
       mpi_n3d_id!=mpi_n3d_masterid)
      return;
#endif

   MPI_TIC(mpiTiming.fileio_0);

    
   if (d->opt.ExecParam.Strain.StrainOut.DisplacementData){
      filename = n3d_strdup_n(d->inputfile);
#ifdef CHECK_COMMUNICATION
      n3d_FileTypeSet(&filename, "nd_disp"cpu, TRUE);
#else
      n3d_FileTypeSet(&filename, "nd_disp", TRUE);
#endif
      printf("Dump out the bond distortion to 	filename: %s\n",filename),fflush(stdout);
       
      fp = fopen( filename, "w" );
      for ( i=0; i < d->geo.N_Cell; i++ ){
	 for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ){
	    if ( d->geo.AtomType[i][j] ){
	       for ( k=0; k < 3; k++ ){
		  fprintf(fp, "%21.12e  ",d->geo.l_m_2_atomdsp[i][j][k]);
	       }
	       fprintf(fp,"\n");
	    }
	 }
      }
      fclose(fp);

      str_free(filename);
   }

   if (d->opt.ExecParam.Strain.StrainOut.DisplacementData3D){
      filename = n3d_strdup_n(d->inputfile);
#ifdef CHECK_COMMUNICATION
      n3d_FileTypeSet(&filename, "nd_disp_3d"cpu, TRUE);
#else
      n3d_FileTypeSet(&filename, "nd_disp_3d", TRUE);
#endif
      printf("Dump out 3d-plottable the bond distortion to filename: %s\n",filename),fflush(stdout);
      fp2 = fopen( filename, "w" );
      str_free(filename);
      for ( i=0; i < d->geo.N_Cell; i++ ){
	 for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ){
	    if ( d->geo.AtomType[i][j] ){
	       xa = d->geo.lattice_x*d->geo.cell__ijk[i][0] + 
                    d->geo.lattice_x*d->geo.PositionInCell(j,0);
	       ya = d->geo.lattice_y*d->geo.cell__ijk[i][1] + 
                    d->geo.lattice_y*d->geo.PositionInCell(j,1);
	       za = d->geo.lattice_z*d->geo.cell__ijk[i][2] + 
                    d->geo.lattice_z*d->geo.PositionInCell(j,2);

	       fprintf(fp2, "%e %e %e", xa, ya, za );
	       fprintf(fp2, " %e",d->geo.l_m_2_atomdsp[i][j][0]);
	       fprintf(fp2, " %e",d->geo.l_m_2_atomdsp[i][j][1]);
	       fprintf(fp2, " %e",d->geo.l_m_2_atomdsp[i][j][2]);
	       fprintf(fp2, " %e",sqrt(d->geo.l_m_2_atomdsp[i][j][0]*d->geo.l_m_2_atomdsp[i][j][0]+
				       d->geo.l_m_2_atomdsp[i][j][1]*d->geo.l_m_2_atomdsp[i][j][1]+
				       d->geo.l_m_2_atomdsp[i][j][2]*d->geo.l_m_2_atomdsp[i][j][2]));
	       fprintf(fp2,"\n");
	    }
	 }
      }
      fclose(fp2);
   }
   
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

   return;
}


void read_displacement( qd_struct d, rvectr p)
{
   int i, I, j, id, proc;
   FILE *fp=NULL;
   char line[100], item[30], *token, *tail;
   
   MPI_TIC(mpiTiming.fileio_0);

   if (mpi_n3d_id == mpi_n3d_masterid) {
      printf("Read strain distribution\n"),fflush(stdout);
      if ( !(fp = fopen( d->opt.ExecParam.Strain.Keating.StartFileRead.c_str(), "r" )) ) {
	 die("File \"%s\" not found in working directory",
	       d->opt.ExecParam.Strain.Keating.StartFileRead.c_str());
      }

      for  (proc=0; proc<mpi_n3d_numprocs; proc++) {
	 rvectr p_snd = Rvectr(3*d->cell_ln[proc]*d->geo.AtomsPerCellMax());
	 I = d->cell_s[proc];
	 for ( i=0; i < d->cell_ln[proc]; i++, I++ ) {
	    for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ){
	       if ( d->geo.AtomType[I][j] ){
		  fgets(line,100,fp);
		  token = strtok(line,DELIMITERS);
		  strcpy(item,token);
		  id = i*d->geo.AtomsPerCellMax()*3 + j*3;
		  p_snd[id++] = d->geo.l_m_2_atomdsp[I][j][0] = strtod(item,&tail);
		  token = strtok(NULL,DELIMITERS);
		  strcpy(item,token);
		  p_snd[id++] = d->geo.l_m_2_atomdsp[I][j][1] = strtod(item,&tail);
		  token = strtok(NULL,DELIMITERS);
		  strcpy(item,token);
		  p_snd[id++] = d->geo.l_m_2_atomdsp[I][j][2] = strtod(item,&tail);
	       }
	    }
	 }

	 if (proc) {
	    /* might use non-blocking send? */
	    MPI_Send(p_snd, 3*d->geo.AtomsPerCellMax()*d->cell_ln[proc], 
		     MPI_DOUBLE, proc, mpi_n3d_masterid, MPI_COMM_WORLD);
	 }
	 else {
	    for (i=0; i<3*d->cell_ln[proc]*d->geo.AtomsPerCellMax(); i++)
	       p[i] = p_snd[i];
	 }
	 rm_rvectr(&p_snd);
      }
      fclose(fp);
   } 
   else {
      MPI_Recv(p, 3*d->geo.AtomsPerCellMax()*d->cell_ln[mpi_n3d_id], MPI_DOUBLE, 
	       mpi_n3d_masterid, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_n3d_status);

      I = d->cell_s[mpi_n3d_id];
      for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, I++ ) {
	 for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ){
	    if ( d->geo.AtomType[I][j] ){
	       id = i*d->geo.AtomsPerCellMax()*3 + j*3;
	       d->geo.l_m_2_atomdsp[I][j][0] = p[id++];
	       d->geo.l_m_2_atomdsp[I][j][1] = p[id++];
	       d->geo.l_m_2_atomdsp[I][j][2] = p[id++];
	    }
	 }
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}


/* some code in this function assumes zinc-blende construction */
void write_nbr( qd_struct d, int first_or_second ) 
//The variable first_or_second tells whether this function is called the first or second time.
{
   real locbs_xmin, locbs_xmax, locbs_ymin, locbs_ymax, locbs_zmin, locbs_zmax;

   bool has_periodicity = d->geo.hasPeriodicity();
   int b, l, m, n, p;
   double dl, x, y, z, alat;
   ivectr cell_nbr=NULL;
   FILE* fp;
   double dr[4][3], dr0[4][3];
   char label[100];
   
   int isx = d->cell_s[mpi_n3d_id];
   int iex = isx + d->cell_ln[mpi_n3d_id];

   int map_mm_to_n[8][8];
   int map_mb_to_m[8][4];  /* map (cation m, block label of nbr) to m of nbr */ 
   int map_n_to_l[4];

   real a_lattice[3];
   a_lattice[0]=d->geo.lattice_x;
   a_lattice[1]=d->geo.lattice_y;
   a_lattice[2]=d->geo.lattice_z;

   if (!d->opt.ExecParam.Strain.StrainOut.Neighborhood ) return;

   MPI_TIC(mpiTiming.fileio_0);

   locbs_xmin = d->opt.ExecParam.Subdomain.xMin;
   locbs_xmax = d->opt.ExecParam.Subdomain.xMax;
   locbs_ymin = d->opt.ExecParam.Subdomain.yMin;
   locbs_ymax = d->opt.ExecParam.Subdomain.yMax;
   locbs_zmin = d->opt.ExecParam.Subdomain.zMin;
   locbs_zmax = d->opt.ExecParam.Subdomain.zMax;

   if (!mpi_n3d_id) {
      printf("\nEach processor is currently writing out neighborhood info ... ");
      fflush(stdout);
   }
   
   /* need map from nbr indx (n) to atom indx (m) */
   for ( m=0; m < d->geo.AtomsPerCellMax(); m++ ) {
      for ( n = 0; n < d->geo.Neighbors(m); n++ ) {
         int m_nbr = d->geo.NbrCell(m,n,3);
         map_mm_to_n[m][m_nbr] = n;
      }
   }
   
   for ( m=0; m < 4; m++ ) {
      map_mb_to_m[m][0] = 4 + m;
      map_mb_to_m[m][1] = 4 + (7-m)%4;
      map_mb_to_m[m][2] = 4 + (2+m)%4;
      map_mb_to_m[m][3] = 4 + (5-m)%4;
   }
   for ( m=4; m < 8; m++ ) {
      map_mb_to_m[m][0] = (m-2)%4;
      map_mb_to_m[m][1] = (9-m)%4;
      map_mb_to_m[m][2] =     m%4;
      map_mb_to_m[m][3] = (7-m)%4;
   }
   
   
   {
      char *filename=NULL;
      if(first_or_second == 1) sprintf(label,"nd_bondLen_first_%d", mpi_n3d_id);
      else  sprintf(label,"nd_bondLen_second_%d", mpi_n3d_id);
      filename = n3d_strdup_n(d->inputfile);
      n3d_FileTypeSet(&filename, label, TRUE);
      fp = fopen(filename,"w");
      str_free(filename);
   }
   
   for ( l=isx; l < iex; l++ ) {
      
      int flag;
      int i = d->geo.cell__ijk[l][0];
      int j = d->geo.cell__ijk[l][1];
      int k = d->geo.cell__ijk[l][2];
      
      for ( m=0; m < d->geo.AtomsPerCellMax(); m++ ) {
         if ( ! d->geo.AtomType[l][m] ) continue;
         
         x = d->geo.lattice_x*(d->geo.cell__ijk[l][0] + d->geo.PositionInCell(m,0));
         y = d->geo.lattice_y*(d->geo.cell__ijk[l][1] + d->geo.PositionInCell(m,1));
         z = d->geo.lattice_z*(d->geo.cell__ijk[l][2] + d->geo.PositionInCell(m,2));
         
         
         /* only look at those cations within this box */
         if (x < locbs_xmin || x > locbs_xmax) continue;
         if (y < locbs_ymin || y > locbs_ymax) continue;
         if (z < locbs_zmin || z > locbs_zmax) continue;
         
         /* check that (l,m) has all its neighbors */
         flag = 0;
         for ( n = 0; n < d->geo.Neighbors(m); n++ ) {
            int i_nbr = i + d->geo.NbrCell(m,n,0);
            int j_nbr = j + d->geo.NbrCell(m,n,1);
            int k_nbr = k + d->geo.NbrCell(m,n,2);
            
            int aindx = d->geo.NbrCell(m,n,3);
            int cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

#ifdef ELIMINATE_SSMAP
            if (cindx < 0 && has_periodicity && (cindx + d->geo.N_Cell + 1 ) >= 0) {
	      cindx = cindx + d->geo.N_Cell + 1 ; 
            }
            if(d->geo.is_tilted)
               cindx = d->geo.get_cindx_tilted(l,m,n);
#else /* ELIMINATE_SSMAP */
            if (cindx < 0 && has_periodicity && d->geo.ssmap[l][m][n]!=-1) {
               cindx = d->geo.ssmap[l][m][n];
            }
#endif /* ELIMINATE_SSMAP */
            
            map_n_to_l[n] = cindx;
            
            if (cindx<0 || !d->geo.AtomType[cindx][aindx]) { flag=1; break; } 
         }
         
         if (flag) continue;  /* incomplete set of neighbors found; throw out */

#ifdef DISPL_SINGLE_PRECISION         
         fvectr disp_this = d->geo.l_m_2_atomdsp[l][m];
#else /* DISPL_SINGLE_PRECISION */
         rvectr disp_this = d->geo.l_m_2_atomdsp[l][m];
#endif /* DISPL_SINGLE_PRECISION */         
         for ( b = 0; b < d->geo.Neighbors(m); b++ ) {
            int cindx;
            int aindx = map_mb_to_m[m][b];
            n = map_mm_to_n[m][aindx];
            cindx = map_n_to_l[n];
            
            alat = d->strain.amat[ (int) d->geo.AtomType[l][m] ][ (int) d->geo.AtomType[cindx][aindx] ];
#ifdef DISPL_SINGLE_PRECISION
            fvectr disp_nbr = d->geo.l_m_2_atomdsp[cindx][aindx];
#else /* DISPL_SINGLE_PRECISION */
            rvectr disp_nbr = d->geo.l_m_2_atomdsp[cindx][aindx];
#endif /* DISPL_SINGLE_PRECISION */
            rvectr rAtom_nbr = d->geo.PositionInCell(aindx);
            rvectr rAtom_this = d->geo.PositionInCell(m);
            cell_nbr = d->geo.NbrCell(m,n);
            
            p=0;
            for (p=0; p<3; p++) {
               dl = cell_nbr[p] + rAtom_nbr[p] - rAtom_this[p];
               dr[b][p] = a_lattice[p]*dl + (disp_nbr[p] - disp_this[p]);
               dr0[b][p] = alat*dl;
            }
         }
         
         fprintf(fp, "%3d   %11.8f %11.8f %11.8f   %11.8f %11.8f %11.8f   %11.8f %11.8f %11.8f   %11.8f %11.8f %11.8f   %11.8f %11.8f %11.8f   %11.8f %11.8f %11.8f   %11.8f %11.8f %11.8f   %11.8f %11.8f %11.8f   %11.8f %11.8f %11.8f\n",
                 (int) d->geo.AtomType[l][m],
                 x+disp_this[0], y+disp_this[1], z+disp_this[2],
                 dr[0][0], dr[0][1], dr[0][2],
                 dr[1][0], dr[1][1], dr[1][2],
                 dr[2][0], dr[2][1], dr[2][2],
                 dr[3][0], dr[3][1], dr[3][2],
                 dr0[0][0], dr0[0][1], dr0[0][2],
                 dr0[1][0], dr0[1][1], dr0[1][2],
                 dr0[2][0], dr0[2][1], dr0[2][2],
                 dr0[3][0], dr0[3][1], dr0[3][2]
                 );
         
         
      }
   }
   fclose(fp);
   
   printf(" %d ", mpi_n3d_id);  fflush(stdout);
  
   MPI_Barrier(MPI_COMM_WORLD);

   if (!mpi_n3d_id) {
      printf("done!\n");
      fflush(stdout);
   }
   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);

}


/* This is primarily intended for the simple singel band case where the 
   overall structure may be a free-standing pyramid and some atoms might have been 
   out of the original cubic domain.  We need to fill in these missing atoms
   for visualization on a regular mesh */
int write_cubic_struct_wavefunctions( qd_struct d )
{
  
   int i, j, cell, ev;
   int xi, yi, zi;

   i3tensor geo_xi_yi_zi_2_cell=NULL;
   b3tensor geo_xi_yi_zi_inside=NULL;
   cvectr wfc=NULL;
   FILE *f_dx;
   char label[100], *filename=NULL;

   if (!d->opt.ExecParam.ElCalc.ElOut.PsiSqr || mpi_n3d_id!=mpi_n3d_masterid)
      return 0;

   MPI_TIC(mpiTiming.fileio_0);

    
   geo_xi_yi_zi_2_cell = I3tensor(0,d->geo.Nx_Cell,0,d->geo.Ny_Cell,0,d->geo.Nz_Cell);
   geo_xi_yi_zi_inside = B3tensor(0,d->geo.Nx_Cell,0,d->geo.Ny_Cell,0,d->geo.Nz_Cell);
   

   /* for the app-qdot application on nanoHUB it would be nice to outline the shapes of the 
      embedded quantum dots.  In the simplest of all QD calculations there is only one material 
      and the confinement is "air", therefore the outer edge of all included atoms / sites is 
      really forming the confinement area.
      However if it is an embedded quantum dot, then we are looking to put the shape into the transition 
      of the confining material to the confinement material.  Here at the output level we have really 
      lost all shape information already.  Really this processing should happen at the input scanning routine!
   */
   /* Loop over the whole structure and figure out the minimum an the maximum 
      indices in the three coordinate axes */
   for ( i=0; i < d->geo.N_Cell; i++ ){
      /* Only loop over the cells here! 
         Do not loop over the atoms in the cell. 
         Assume that the overall cell has a cubic topology.  */
      j = 0;  
      /* set the atom in cell index to the first one.
         This will work exactly for the cubic case. */
      if ( d->geo.AtomType[i][j] > 0 ) {
	    
         xi = d->geo.cell__ijk[i][0];
         yi = d->geo.cell__ijk[i][1];
         zi = d->geo.cell__ijk[i][2];
         
         geo_xi_yi_zi_2_cell[xi][yi][zi] = i;
         geo_xi_yi_zi_inside[xi][yi][zi] = 1;
      }
   }
   
   
   for (ev=0;ev<d->neigv;ev++){
      
      sprintf(label, "nd_dx_%d", ev+1);
      filename = n3d_strdup_n(d->inputfile);
      n3d_FileTypeSet(&filename, label, TRUE);


      f_dx = fopen( filename, "w" );
      str_free(filename);
      fprintf(f_dx,"#### Generated directly from NEMO3D\n");
      fprintf(f_dx,"#### eigenvalue # %d\n", ev+1);
      fprintf(f_dx,"#### Eigenvalue %g\n", d->WF[ev].E*HAMILTONIAN_SCALE_VALUE);
      fprintf(f_dx,"object 1 class gridpositions counts %d %d %d\n",d->geo.Nx_Cell+1,d->geo.Ny_Cell+1,d->geo.Nz_Cell+1);
      fprintf(f_dx,"origin 0 0 0\n");
      fprintf(f_dx,"delta  %g 0 0\n", d->geo.lattice_x);
      fprintf(f_dx,"delta  0 %g 0\n", d->geo.lattice_y);
      fprintf(f_dx,"delta  0 0 %g\n", d->geo.lattice_z);
      fprintf(f_dx,"object 2 class gridconnections counts %d %d %d\n",d->geo.Nx_Cell+1,d->geo.Ny_Cell+1,d->geo.Nz_Cell+1);
      fprintf(f_dx,"attribute \"element type\" string \"cubes\"\n");
      fprintf(f_dx,"attribute \"ref\" string \"positions\"\n");
      fprintf(f_dx,"object 3 class array type float rank 0 items %d data follows\n",(d->geo.Nx_Cell+1)*(d->geo.Ny_Cell+1)*(d->geo.Nz_Cell+1));
      if (1){
         wfc = d->wfc[ev];
         for (xi=0;xi<=d->geo.Nx_Cell;xi++){
            for (yi=0;yi<=d->geo.Ny_Cell;yi++){
               for (zi=0;zi<=d->geo.Nz_Cell;zi++){
                  if (geo_xi_yi_zi_inside[xi][yi][zi]){
                     cell = geo_xi_yi_zi_2_cell[xi][yi][zi]; 
                     fprintf(f_dx,"%g\n",wfc[cell].r*wfc[cell].r+wfc[cell].i*wfc[cell].i);
                  }else{
                     fprintf(f_dx,"0.00\n");
                  }
               }
            }
         }
      }else{
         for (xi=0;xi<=d->geo.Nx_Cell;xi++){
            for (yi=0;yi<=d->geo.Ny_Cell;yi++){
               for (zi=0;zi<=d->geo.Nz_Cell;zi++){
                  if (geo_xi_yi_zi_inside[xi][yi][zi]){
                     cell = geo_xi_yi_zi_2_cell[xi][yi][zi]; 
                     fprintf(f_dx,"%g\n",wfc[cell].r*wfc[cell].r+wfc[cell].i*wfc[cell].i);
                  }else{
                     fprintf(f_dx,"0.00\n");
                  }
               }
            }
         }
      }
      fprintf(f_dx,"attribute \"dep\" string \"positions\"\n");
      fprintf(f_dx,"object \"regular positions regular connections\" class field\n");
      fprintf(f_dx,"component \"positions\" value 1\n");
      fprintf(f_dx,"component \"connections\" value 2\n");
      fprintf(f_dx,"component \"data\" value 3\n");
      fprintf(f_dx,"end\n");
      fclose(f_dx);      
   }
   
      sprintf(label, "nd_dx_shape");
      filename = n3d_strdup_n(d->inputfile);
      n3d_FileTypeSet(&filename, label, TRUE);


      f_dx = fopen( filename, "w" );
      str_free(filename);
      fprintf(f_dx,"#### Generated directly from NEMO3D - confinement domain\n");
      fprintf(f_dx,"object 1 class gridpositions counts %d %d %d\n",d->geo.Nx_Cell+1,d->geo.Ny_Cell+1,d->geo.Nz_Cell+1);
      fprintf(f_dx,"origin 0 0 0\n");
      fprintf(f_dx,"delta  %g 0 0\n", d->geo.lattice_x);
      fprintf(f_dx,"delta  0 %g 0\n", d->geo.lattice_y);
      fprintf(f_dx,"delta  0 0 %g\n", d->geo.lattice_z);
      fprintf(f_dx,"object 2 class gridconnections counts %d %d %d\n",d->geo.Nx_Cell+1,d->geo.Ny_Cell+1,d->geo.Nz_Cell+1);
      fprintf(f_dx,"attribute \"element type\" string \"cubes\"\n");
      fprintf(f_dx,"attribute \"ref\" string \"positions\"\n");
      fprintf(f_dx,"object 3 class array type float rank 0 items %d data follows\n",(d->geo.Nx_Cell+1)*(d->geo.Ny_Cell+1)*(d->geo.Nz_Cell+1));
      wfc = d->wfc[ev];
      for (xi=0;xi<=d->geo.Nx_Cell;xi++){
         for (yi=0;yi<=d->geo.Ny_Cell;yi++){
            for (zi=0;zi<=d->geo.Nz_Cell;zi++){
               fprintf(f_dx,"%d\n",(int) geo_xi_yi_zi_inside[xi][yi][zi]);
            }
         }
      }
      fprintf(f_dx,"attribute \"dep\" string \"positions\"\n");
      fprintf(f_dx,"object \"regular positions regular connections\" class field\n");
      fprintf(f_dx,"component \"positions\" value 1\n");
      fprintf(f_dx,"component \"connections\" value 2\n");
      fprintf(f_dx,"component \"data\" value 3\n");
      fprintf(f_dx,"end\n");
      fclose(f_dx);      


   rm_i3tensor(&geo_xi_yi_zi_2_cell,0);
   rm_b3tensor(&geo_xi_yi_zi_inside,0);
   

   MPI_TOC(mpiTiming.fileio, mpiTiming.fileio_0);
   

   return 0;
}
