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
$Header: /repo/nemo3d/src/base/nemo3d_entry.c,v 1.25 2007/08/10 15:03:33 baeh Exp $ 
*****************************************************************************/

#include "nemo3d_entry.h"

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>

//#include <fstream.h>
#include <fstream>

// for getpid
#include <sys/types.h>
#include <unistd.h>

#ifdef MPI_TIMING 
void TimeReport(char* d_filename){
   {
      // gather timing data from all involved CPUs and store the data into 
      // some common files.
      if (mpi_n3d_id==mpi_n3d_masterid){
         int cpu=0;
         int idata=0;
         rmatrix out_matrix = Rmatrix(mpi_n3d_numprocs,mpiTiming.out_data_length);
         rvectr  rec_data = Rvectr(mpiTiming.out_data_length);
         rvectr  max_data = Rvectr(mpiTiming.out_data_length);
         rvectr  min_data = Rvectr(mpiTiming.out_data_length);
         for (idata=0;idata<mpiTiming.out_data_length;idata++){
             min_data[idata] = 1e20;
         }
         

         for (cpu=0;cpu< mpi_n3d_numprocs; cpu++){
             if (cpu==mpi_n3d_masterid){
                 rm_rvectr(&rec_data);
                 rec_data = mpiTiming.GetData();
             }else{
                 MPI_Recv(&rec_data[0], mpiTiming.out_data_length, MPI_DOUBLE, cpu, cpu, MPI_COMM_WORLD, &mpi_n3d_status);
             }

             for (idata=0;idata<mpiTiming.out_data_length;idata++){
                 out_matrix[cpu][idata] = rec_data[idata];
                 max_data[idata] = max(max_data[idata],rec_data[idata]);
                 min_data[idata] = min(max_data[idata],rec_data[idata]);
             }
         }
         {
             char* filename = n3d_strdup_n(d_filename);
             FILE *ftime=NULL;

             if (mpi_n3d_numprocs>1){
                n3d_FileTypeSet(&filename, "nd_time_summary", TRUE);
                ftime = fopen(filename,"a");

                fprintf(ftime,"# CPU(min/max), total (min/max), setup (min/max), strain (min/max), electr (min/max), fileio (min/max), sendrecev (min/max), allreduce (min/max), barrier (min/max), bcast (min/max), send(min/max), receive (min/max), matmul (min/max), strain_sendrecv (min/max), strain_allreduce (min/max), strain_barrier (min/max), strain_bcast (min/max), strain_send (min/max), strain_recv (min/max), electr_sendrecv (min/max), electr_allreduce (min/max), electr_barrier (min/max), electr_bcast (min/max), electr_send (min/max), electr_recv (min/max)\n");
                fprintf(ftime,"%3d ",mpi_n3d_numprocs);
                for (idata=0;idata<mpiTiming.out_data_length;idata++){
                   fprintf(ftime,"%g %g ",min_data[idata],max_data[idata]);
                }
                fprintf(ftime,"\n");
                fclose(ftime);
             }
   
             n3d_FileTypeSet(&filename, "nd_time_detail", TRUE);
             ftime = fopen(filename,"a");
             fprintf(ftime,"# cpu total, setup, strain, electr, fileio, sendrecev, allreduce, barrier, bcast, send, receive, matmul, strain_sendrecv, strain_allreduce, strain_barrier, strain_bcast, strain_send, strain_recv, electr_sendrecv, electr_allreduce, electr_barrier, electr_bcast, electr_send, electr_recv\n");

             if (mpi_n3d_numprocs>1){
                fprintf(ftime,"%3d ",mpi_n3d_numprocs);
                fprintf(ftime,"%g ",-2.0);
                for (idata=1;idata<mpiTiming.out_data_length;idata++){
                   fprintf(ftime,"%g ",min_data[idata]);
                }
                fprintf(ftime,"\n");

                fprintf(ftime,"%3d ",mpi_n3d_numprocs);
                fprintf(ftime,"%g ",-1.0);
                for (idata=1;idata<mpiTiming.out_data_length;idata++){
                   fprintf(ftime,"%g ",max_data[idata]);
                }
                fprintf(ftime,"\n");
 
             }
             for (cpu=0;cpu< mpi_n3d_numprocs; cpu++){
                fprintf(ftime,"%3d ",mpi_n3d_numprocs);
                for (idata=0;idata<mpiTiming.out_data_length;idata++){
                   fprintf(ftime,"%g ",out_matrix[cpu][idata]);
                }
                fprintf(ftime,"\n");
             }
             fclose(ftime);

             str_free(filename);
         }

         rm_rvectr(&rec_data);
         rm_rvectr(&max_data);
         rm_rvectr(&min_data);
         rm_rmatrix(&out_matrix);
 
      }else if (mpi_n3d_numprocs>1){
         rvectr  out_data = mpiTiming.GetData();
         MPI_Send(&out_data[0],mpiTiming.out_data_length,MPI_DOUBLE,mpi_n3d_masterid, mpi_n3d_id, MPI_COMM_WORLD);
         rm_rvectr(&out_data);
      }
   }
}
#endif /* MPI_TIMING */


int nemo3d(const char* inputFile, const char *inFileName)
{
#ifdef MPIP
		MPI_Pcontrol(0); // disable whole-program profiling
#endif
   cout_master << "===========================================\n"
               << "   NEMO-3D\n"
               << "===========================================\n";

   // Create memory for the central simulation structure.
   qd_struct d = new QD_struct;
   
   // Assign the input filename from the command line options.
   d->inputfile = n3d_copy_str(inFileName);

   // Read in data from input file into the top_struct data structure
   d->setTopStruct(inputFile, inFileName);
   
   // check for (and modify if necessary) any incompatiblities in input
   d->checkInput();

   // this only sets up sMatList
   d->set_sMatList();
   
   set_random_seed(d->opt.ExecParam.RandomSeed);

   MPI_TOC(mpiTiming.setup, mpiTiming.setup_0);

   setup_hamiltonian(d);

   if(d->opt.ExecParam.Output.ShapeInfo) writeShapeInfo(d); 
   if(d->opt.ExecParam.Output.AtomInfo) writeAtomInfo(d); 

   if(d->opt.ExecParam.Output.NbrIndx)  {
      writeNbrIndx(d);
      writePhaseInfoForNEMO1D(d);
   }
   
   MPI_TIC(mpiTiming.work_0);

   cout_master << "Finished with set-up of data structures.\n\n\n";
   //Compute electron eigenstates
   if (d->opt.ExecParam.CalcElDisp)
     {
       
       if (d->opt.ExecParam.Phys.HbarOverTau > 0.0) {
         switch (d->opt.ExecParam.ElCalc.ElAlg.ResFind) {
         case ElAlg_struct::Lancz:
           ns_lanc_driver( d ); 
           break;
         case ElAlg_struct::Rayleigh:
           die("Rayleigh quotient method not implemented\n"
	       "for a non-hermitian system.\n");
           break;
         case ElAlg_struct::Direct:
           h_ceigval( d, 0 ); 
           break;
         default:
           die("unknown resonance finding specification.\n");
         }
       } else {
         d->Ham_use_hermiticity=1;
         compute_electronic_structure(d);
       }
     }//If d->opt.ExecParam.CalcElDisp
      //Compute phonon eigenstates
   if (d->opt.ExecParam.CalcPhonDisp)
     {
       // include Atomic mass in qd_stuct
       d->Amass = Rvectr(EL_TOTAL);
       for (vector<MaterialHandle>::iterator v=d->sMatList.begin(); 
	    v != d->sMatList.end(); v++) {

         const Material_struct * const smh = v->hndl;
         if (!smh) continue;
         
         for (int i=0; i < (int) smh->M_atom.size(); i++) {
           int atom = (i < (int) v->cation_id.size())
             ? v->cation_id[i]
             : v->anion_id[i - v->cation_id.size()];
           
           d->Amass[atom] = smh->M_atom[i];
         }
       }
       masterPrint("Start phonon calculation\n\n");
       phon_disp(d);
     }//If d->opt.ExecParam.CalcPhonDisp
   //Postprocessing:
   if (d->opt.ExecParam.ElCalc.ElOut.OptMatrixEl) {
      calc_opt_mat_cmplx(d,
                   d->opt.ExecParam.Abs_Ramp.Ephi,
                   d->opt.ExecParam.Abs_Ramp.Etheta);
   }
   
   if (d->opt.ExecParam.ElCalc.ElOut.Absorption) {
      real min=1.0e100, max=1.0e100, *target;
      int num=-1;

      /*
      if (d->opt.ExecParam.Absorp.AbsCoeff) {
	 calc_opt_rates(d,d->opt.ExecParam.Abs_Ramp.Estrength,
                        d->opt.ExecParam.Abs_Ramp.Ephi,
                        d->opt.ExecParam.Abs_Ramp.Etheta);
      }
      */

      if (d->opt.ExecParam.Absorp.IntegrAbs || d->opt.ExecParam.Absorp.AbsCoeff) {
	 calc_opt_rates_int(d,d->opt.ExecParam.Abs_Ramp.Estrength,
                            d->opt.ExecParam.Abs_Ramp.Ephi,
                            d->opt.ExecParam.Abs_Ramp.Etheta);
      }
      
      target = NULL;
      switch (d->opt.ExecParam.Absorp.ramp){
      case Absorp_struct::fermi:    
	 target = &(d->opt.ExecParam.Abs_Ramp.Efermi);
	 min = d->opt.ExecParam.Abs_Ramp.Efermi_min;
	 max = d->opt.ExecParam.Abs_Ramp.Efermi_max;
	 num = d->opt.ExecParam.Abs_Ramp.Efermi_n;
	 break;
      case Absorp_struct::temp:
	 target = &(d->opt.ExecParam.Abs_Ramp.Etemp);
	 min = d->opt.ExecParam.Abs_Ramp.Etemp_min;
	 max = d->opt.ExecParam.Abs_Ramp.Etemp_max;
	 num = d->opt.ExecParam.Abs_Ramp.Etemp_n;
	 break;
      case Absorp_struct::strength:
	 target = &(d->opt.ExecParam.Abs_Ramp.Estrength);
	 min = d->opt.ExecParam.Abs_Ramp.Estrength_min;
	 max = d->opt.ExecParam.Abs_Ramp.Estrength_max;
	 num = d->opt.ExecParam.Abs_Ramp.Estrength_n;
	 break;
      case Absorp_struct::phi:
	 target = &(d->opt.ExecParam.Abs_Ramp.Ephi);
	 min = d->opt.ExecParam.Abs_Ramp.Ephi_min;
	 max = d->opt.ExecParam.Abs_Ramp.Ephi_max;
	 num = d->opt.ExecParam.Abs_Ramp.Ephi_n;
	 break;
      case Absorp_struct::theta:
	 target = &(d->opt.ExecParam.Abs_Ramp.Etheta);
	 min = d->opt.ExecParam.Abs_Ramp.Etheta_min;
	 max = d->opt.ExecParam.Abs_Ramp.Etheta_max;  
	 num = d->opt.ExecParam.Abs_Ramp.Etheta_n;
         break;
      default:
	 die("ERROR in Absorption Ramp initialization\n");
	 break;
      }
      
      if (target != NULL)
	 calc_opt_rates_int_ramp(d, target, min, max, num); 
   }

   MPI_TOC(mpiTiming.total, mpiTiming.total_0);
   MPI_TOC(mpiTiming.work, mpiTiming.work_0);

   cout_master << "\n\n";

#ifdef MPI_TIMING
   mpiTiming.print();
   TimeReport(d->inputfile);
#endif /* MPI_TIMING */

    nml_memory_report("before clean-up");
   // d->cleanup();
   if (d->wfc)  rm_cmatrix(&(d->wfc));
    rm_i3tensor(&d->mat,0);
    nml_memory_report("after clean-up");

   cout_master << "\nExiting NEMO-3D\n"
               << "===========================================\n";
   return 0;
}

int nemo3d_wrapper(const char* str, char *filename) 
{
   string str_str(str);
   int str_sz = str_str.size();

   char procID[20];
   pid_t pid = getpid();
   sprintf(procID,"_%d_%d", pid, mpi_n3d_id);

#ifdef FIFOMETHOD
   string fifo = string("/tmp/fifo_nemo3d_") + string(procID);
      // string(n3d_FileNameGet(filename)) + string(procID);

   unlink(fifo.c_str());                               // remove old version

   if (mkfifo(fifo.c_str(), 0666) < 0) {               // create fifo
      perror("mkfifo");
      exit(1);
   }
  
   pid_t pid = fork ();
   if (pid < (pid_t) 0) {                          // fork is borked
      fprintf (stderr, "Fork failed.\n");
      return EXIT_FAILURE;
   }
   else if (pid == (pid_t) 0) {                    // child writes to fifo
      /* for now allow only master to write to fifo.  However, if processes
         have access to a local filesystem (e..g /tmp), it may be preferable
         for each process to write to a different fifo */

      //if (!mpi_n3d_id) {
         int fd;
         if ((fd = open(fifo.c_str(), O_WRONLY)) < 0) {
            perror("open");
            exit(1);
         }
         write(fd, str, str_sz);
         close(fd);
         //}
      
      exit(0);
   } else {                                        // parent reads from fifo
      nemo3d(fifo.c_str(), filename);

      unlink(fifo.c_str());
      return EXIT_SUCCESS;
   }
#else

   
#if 0
   //string fifo = string("/tmp/fifo_nemo3d_") + 
   string fifo = string("./fifo_nemo3d_") + 
      string(n3d_FileNameGet(filename)) + string(procID);
      int fd;
         if ((fd = open(fifo.c_str(), O_WRONLY)) < 0) {
            perror("open");
            exit(1);
         }
         write(fd, str, str_sz);
         close(fd);
#else
   //string fifo = string("/tmp/fifo_nemo3d_") + 
//   string fifo = string("/tmp/fifo_nemo3d_") + 
//      string(n3d_FileNameGet(filename)) + string(procID);

   string fifo = string("/tmp/fifo_nemo3d") + string(procID);
      FILE *fd;
         fd = fopen(fifo.c_str(), "w");
         fprintf(fd,"%s\n",str);
         fclose(fd);
#endif
      nemo3d(fifo.c_str(), filename);

      return EXIT_SUCCESS;

#endif

}

// The definition below establishes an unmangled C function compiled with C++
// such that this function can be called from a C-linked executable
#ifdef __cplusplus
extern "C" {
#endif
int nemo3d_wrapper_c(const char* str,  char *filename)
{
    return nemo3d_wrapper(str, filename);
}
#ifdef __cplusplus
}
#endif

