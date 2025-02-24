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
$Header: /repo/nemo3d/src/base/strain3d.c,v 1.47 2008/07/23 16:49:37 nkharche Exp $ 
*****************************************************************************/

#include "strain3d.h"
#include <fstream>

#ifdef INTEL_MODIFIED
    #include <emmintrin.h>
    #include <pmmintrin.h>

    #define _PTR2OFFS(p) (((char*)(p)) - ((char*)(0)))
    #define _OFFS2PTR(n) (((char*)(0)) + (n))
    #define _MAKE_ALIGNED_PTR(type, p, a) \
        (type*)_OFFS2PTR((_PTR2OFFS(a-1)+_PTR2OFFS(p))&~_PTR2OFFS(a-1))
#endif

short _new_minimize = 1;

#define swapInt(g, h)  { int y=(g); (g)=(h); (h)=y; }

static void quick_sort(int* x, int beg, int end) 
{
    int p, i;
   
    if (beg>=end) return;
   
    swapInt(x[beg], x[(beg+end)/2]);
   
    p = beg;

    for (i = beg+1; i<=end; i++)
        if(x[i] < x[beg]){
	        ++p;
	        swapInt(x[i], x[p]);
        }
   
    swapInt(x[beg], x[p]);
   
    quick_sort(x, beg, p-1);
    quick_sort(x, p+1, end);
}



/*! ========================================================================= */
/*!
   vectors p, xi are padded with 3 additional dimensions only on master proc
   and if d->strain.doPeriodMinimization()
   and if _new_minimize is set

   here is the basic structure:

      min_estrain:  create msg passing data struct
                    call macoptII_par which minimizes calc_estrain_macopt

      calc_estrain_macopt: distribute atomic positions
                               call calc_estrain
                               set xi[id] <- d->strain.dp

      calc_estrain:
      
*/
/*! ========================================================================= */

void errorCheck(qd_struct d, char* msg) {
    int isx = d->cell_s[mpi_n3d_id];
    int iex = isx + d->cell_ln[mpi_n3d_id];
    bool has_periodicity = d->geo.hasPeriodicity();

    for (int Zc=isx; Zc < iex; Zc++ ) {
        int i = d->geo.cell__ijk[Zc][0];
        int j = d->geo.cell__ijk[Zc][1];
        int k = d->geo.cell__ijk[Zc][2];

        for (int  Za=0; Za < d->geo.AtomsPerCellMax(); Za++ ) {
            int n, Zc_nbr, Za_nbr;
            int AtomType_this = (int) d->geo.AtomType[Zc][Za];

            if ( !AtomType_this ) continue;
      
            for ( n = 0; n < d->geo.Neighbors(Za); n++ ) {
                int i_nbr = i + d->geo.NbrCell(Za,n,0);
                int j_nbr = j + d->geo.NbrCell(Za,n,1);
                int k_nbr = k + d->geo.NbrCell(Za,n,2);

                Za_nbr = d->geo.NbrCell(Za,n,3);
                Zc_nbr = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

                #ifdef ELIMINATE_SSMAP
                if (Zc_nbr < 0 && has_periodicity && (Zc_nbr + d->geo.N_Cell + 1 ) >= 0)
	                Zc_nbr = Zc_nbr + d->geo.N_Cell + 1 ; 
                if(d->geo.is_tilted)
                    Zc_nbr = d->geo.get_cindx_tilted(Zc,Za,n);

                #else /*! ELIMINATE_SSMAP */
                if (/*Zc_nbr < 0 &&*/ has_periodicity && d->geo.ssmap[Zc][Za][n]!=-1) /*! When crystal is tilted ssmap should be used for connecting atoms in cell with index = 0 as well */
                    Zc_nbr = d->geo.ssmap[Zc][Za][n];
                #endif /*! ELIMINATE_SSMAP */
         
                if ( Zc_nbr >= 0 ) {
                    int AtomType_nbr = (int) d->geo.AtomType[Zc_nbr][Za_nbr];
               
                    if (AtomType_nbr != 12 && AtomType_nbr != 10 && AtomType_nbr != 14) {
                        printf("ACHTUNG(%d)!!!: nbr atom type = %d  atom type = %d  (indx=%d %d)\n", 
                         mpi_n3d_id, AtomType_nbr, AtomType_this, Zc_nbr, Za_nbr);
    }}}}}
}
   

void min_estrain(qd_struct d) 
/*! Comment_N: This is like a main() function in strain calculation. It is called from 'setup_hamiltonian' function in ham3d_final.c */
//!first_or_second must be either 1 or 2.
//!If first_or_second is 1, min_estrain is called the first time.
//!If first_or_second is 2, min_estrain is called the second time.
{
    if (d->opt.ExecParam.Strain.StrainModel == Strain_struct::Read) {
        /*! Comment_N: Read displacements from a file. File name is provided from the input deck. Displacements read from this file are stored in l_m_2_atomdsp[cell_indx][atom_indx][coord_indx]. coord_indx: 0,1,2 for x,y,z respectively. If strain calculation is turned off (from the input deck) then these will be the final displacements of relaxed structure from original structure. If Strain calculation is on then these displacements will be used as initial guess. */
        const bool do_min = d->strain.doPeriodMinimization();
        d->strain.debugp = 0;
        int Npad = (!mpi_n3d_id && do_min ? 3 : 0);
        int n_loc = 3*d->cell_ln[mpi_n3d_id]*d->geo.AtomsPerCellMax() + Npad;
        rvectr p = Rvectr(n_loc);

        readDisplacement(d,p);

        rm_rvectr(&p);
        write_nbr(d, 1);

        return;
    }
    else if ((d->opt.ExecParam.Strain.StrainModel==Strain_struct::VFF_keating) ||
             (d->opt.ExecParam.Strain.StrainModel==Strain_struct::VFF_keating_strained)){
        /*! Comment_N: VFF_Keating minimization from scratch i.e. no initial guess is provided. */  
        MPI_TIME_INIT(Tstart);

        nml_memory_report("Strain_minimization- beginning");

        double E, E_min, E_max, latt_min=0.0, latt_max, latt;
      
        //! Are we minimizing the strain energy 
        //! with respect to the lattice constant in at least one direction?
        const bool do_min =d->strain.doPeriodMinimization();
        /*! Comment_N: doPeriodMinimization() is coded in file StrainParam.h. Relevent portion of the input deck is as follows:
                            <name>Minimize wrt lattice const</name>
                            <cTag>MinWrtLatt</cTag>
                            <desc>allow the periodicity to relax</desc>
                            <value>z</value>
                            <option>none</option>
                            <option>x</option>
                            <option>y</option>
                            <option>z</option>
                            <option>xy</option>
                            <option>xz</option>
                            <option>yz</option>
                            <option>xyz</option>
In the above calculation a boolean called d->strain.MinWrtZ will be true and lattice constant in z-direction will be allowed to vary. a_x and a_y will remain fixed. Initial values of a_x, a_y and a_z are specified in the begining of the input deck (Just below the band model). */
        int Npad = (!mpi_n3d_id && do_min ? 3 : 0);
        int n_loc = 3*d->cell_ln[mpi_n3d_id]*d->geo.AtomsPerCellMax() + Npad;
        int n_tot = 3*d->n_atom_tot + Npad;
        rvectr p = Rvectr(n_loc);
/*! Comment_N: p - is the vector of independent variables which will be interatively modified in the Conjugate Gradient (CG) routine. Total length of this vecot is 
  n_tot = n_loc[0] + n_loc[1] + ... + n_loc[#CPUs-1]. 
  For all CPUs except the mater CPU (CPU #0):
  n_loc = 3*d->cell_ln[mpi_n3d_id]*d->geo.AtomsPerCellMax(). 3 - (x,y,z). d->cell_ln[mpi_n3d_id] - # of unit cells on proc #mpi_n3d_id and AtomsPerCellMax() = # atoms in each unit cell 
  For master CPU:
  n_loc = 3*d->cell_ln[mpi_n3d_id]*d->geo.AtomsPerCellMax() + Npad. 
  Npad = 3 if lattice constant minimization is turned on. Lattice constants are stored on the master CPU and they are broadcasted to all other CPUs for calculation of gradient for next iteration. */
        macopt_args amac;
 
        macopt_defaults ( &amac, d->opt.ExecParam.Strain.Keating.MaxIter,
                          d->opt.ExecParam.Strain.Keating.verbose,
                          d->opt.ExecParam.Strain.Keating.tol,
                          int(d->opt.ExecParam.Strain.Keating.ConvCrit));
/*! Comment_N: macoptII_par is a CG routine. And above are the arguments passed to this routine. Their values are set from the input deck. */
        amac.filename_its=NULL;
        amac.filename_last=NULL;

        if (d->opt.ExecParam.Strain.StrainOut.EnergyCGiter){
            char* filename = n3d_strdup_n(d->inputfile);

            if(d->strain.first_or_second==1) n3d_FileTypeSet(&filename, "nd_estrain_iter_first", TRUE);
            else n3d_FileTypeSet(&filename, "nd_estrain_iter_second", TRUE);

            amac.filename_its = filename;
        }

        if (d->opt.ExecParam.Strain.StrainOut.EnergyFinal){
            char* filename = n3d_strdup_n(d->inputfile);

            if(d->strain.first_or_second==1) n3d_FileTypeSet(&filename, "nd_estrain_first", TRUE);
            else n3d_FileTypeSet(&filename, "nd_estrain_second", TRUE);

            amac.filename_last = filename;
        }
      
        d->strain.dp = Rvectr(n_loc);
/*! Comment_N: strain.dp is the gradient of vector p. It is passed to macoptII_par for calculation of next approximation to p. Note that the length of strain.dp on every processor are same as that of p which is equal to n_loc for that processor. */
        d->strain.debugp = 0;
/*! Comment_N: If debugp == 1. Then code runs in debug mode. In this mode it prints out some important values which could ease the dibugging process. */ 
     
        if(d->strain.first_or_second == 1) readDisplacement(d,p);
        //! In the first strain calculation, initialize both p and d->geo.l_m_2_atomdsp from a file,
        //! if the input file exists.
        else {
            /*! In the second strain calculation, */
            /*! read d->geo.l_m_2_atomdsp and initialize p vector */
            int Icell = d->cell_s[mpi_n3d_id];
            int i_ca=0;

            for ( int i=0; i < d->cell_ln[mpi_n3d_id]; i++, Icell++ ) {
                for ( int j=0; j < d->geo.AtomsPerCellMax(); j++){
                    p[i_ca++] = d->geo.l_m_2_atomdsp[Icell][j][0];
                    p[i_ca++] = d->geo.l_m_2_atomdsp[Icell][j][1];
                    p[i_ca++] = d->geo.l_m_2_atomdsp[Icell][j][2];
            }}
        }

        #ifdef DISPL_SINGLE_PRECISION
        /*
        {
	        d->geo.l_m_2_atomdsp_onedomain = F3tensor(0,d->cell_ln[mpi_n3d_id]-1,
						  0,d->geo.AtomsPerCellMax() - 1, 0, 2) ;
            int Icell = d->cell_s[mpi_n3d_id];

            for ( int i=0; i < d->cell_ln[mpi_n3d_id]; i++, Icell++ ) {
	            for ( int j=0; j < d->geo.AtomsPerCellMax(); j++){
	                d->geo.l_m_2_atomdsp_onedomain[i][j][0] = d->geo.l_m_2_atomdsp[Icell][j][0] ;
	                d->geo.l_m_2_atomdsp_onedomain[i][j][1] = d->geo.l_m_2_atomdsp[Icell][j][1] ;
	                d->geo.l_m_2_atomdsp_onedomain[i][j][2] = d->geo.l_m_2_atomdsp[Icell][j][2] ;
	        }}

	        rm_f3tensor(&d->geo.l_m_2_atomdsp,0) ;
	        nml_memory_report("REPLACED THE STRAIN ARRAY WITH A TEMPORARY COPY - FOR MACOPT") ;
        }

        NOTE: THIS CODE CAUSES AN INCORRECT OPERATION OF THE PROGRAM
        IF THE USER CHOOSES TO CHECKPOINT THE DISPLACEMENTS

        */
        #endif /*! DISPL_SINGLE_PRECISION */

        nml_memory_report("Strain_minimization- basic vectors initialized");
        initStrainComm(d);
/*! Comment_N: This function initializes the communication arrays to communicate calculated strain results among participating CPUs. */
        nml_memory_report("Strain_minimization- strain communicator intitalized ");
      
        if (!do_min) {
            /*! don't try to minimize energy with respect to lattice constant */
/*! Comment_N: To remove minimization wrt lattice constant this part should be executed always. */
            macoptII_par(p, n_loc, n_tot, (calc_estrain_macopt), (void*) d, &amac);
/*! Comment_N: macoptII_par is the actual CG routine. It is coded in the file NEMO_3D/eigen/src/macopt_par.c */

            E = amac.track_one_value;
        }
        else {   /*! do try to minimize energy with respect to lattice constant */
/*! Comment_N: We want to remove the following part which minimizes energy wrt lattice constant. */
            /*         
            latt_min = (d->geo.lattice_x + d->geo.lattice_y + d->geo.lattice_z) / 3.0;
            d->geo.lattice_x = latt_min;
            d->geo.lattice_y = latt_min;
            d->geo.lattice_z = latt_min;
            */         

            if (Npad) {
/*! Comment_N: Npad is non-zero only for the master proc. Last 3 elements of p on the master proc are used to store (variable) lattice constants. SCALE_FCTR is provided from the input deck. Relevant portion of the input deck is as follows: 
                        <param type="real">
                            <name>What is the scaling factor for dE/da?</name>
                            <cTag>SCALE_FCTR</cTag>
                            <desc>What is the scaling factor (S) for dE/da?
                                  Use number between 1e-2 to 1e-5 depending
                                  on the number of atoms (N) in strain calculation.
                                  An ideal S is about 1.0/sqrt(N).
                            </desc>
                            <value>1e-1</value>
                        </param>
I don't know the use of this scaling factor. */
                double SCALE_FCTR  = d->strain.SCALE_FCTR;
                p[n_loc-3] = d->geo.lattice_x / SCALE_FCTR;
                p[n_loc-2] = d->geo.lattice_y / SCALE_FCTR;
                p[n_loc-1] = d->geo.lattice_z / SCALE_FCTR;
/*! Comment_N: Above three lines assign initial values of lattice constants to appropriate elements of p. */
            }
        
            macoptII_par(p, n_loc, n_tot, (calc_estrain_macopt), (void*) d, &amac); 
/*! Comment_N: macoptII_par is the actual CG routine. It is coded in the file NEMO_3D/eigen/src/macopt_par.c */
            E_min = amac.track_one_value;
         
            if (!_new_minimize) {
/*! Comment_N: _new_minimize appears to be always one. Its value is set on line #52 in the same (this) file. So te following code is not executed. */
                latt_max = 0.62;
                d->geo.lattice_x = latt_max;
                d->geo.lattice_y = latt_max;
                d->geo.lattice_z = latt_max;
            
                macoptII_par(p, n_loc, n_tot, (calc_estrain_macopt), (void*) d, &amac); 

                E_max = amac.track_one_value;
            
                while (fabs(latt_max-latt_min) > 1e-8) {
                    latt = .5*(latt_min + latt_max);
                    d->geo.lattice_x = latt;
                    d->geo.lattice_y = latt;
                    d->geo.lattice_z = latt;
               
                    macoptII_par(p, n_loc, n_tot, (calc_estrain_macopt), (void*) d, &amac ); 
               
                    E = amac.track_one_value;
               
                    if (fabs(E-E_min) < fabs(E-E_max)) {
                        E_max = E;
                        latt_max = latt;
                    }
                    else {
                        E_min = E;
                        latt_min = latt;
                    }
               
                    /* latt_max=latt_min; */
                }

                printf("LATT in (%e %e)\n", latt_min, latt_max);
            }
         
/*! Comment_N: Minimization is done. Print out new lattice constants on the master CPU. */
            cout_master << "a=(" 
                     << d->geo.lattice_x << ","
                     << d->geo.lattice_y << ","
                     << d->geo.lattice_z << ")\n";
        }
      
        d->strain.debugp = 0;

        #ifdef DISPL_SINGLE_PRECISION
        /*
        {
	        d->geo.l_m_2_atomdsp = F3tensor(0,d->geo.N_Cell - 1,
					0,d->geo.AtomsPerCellMax() - 1 ,0, 2) ;
            int Icell = d->cell_s[mpi_n3d_id];
            int i_ca=0;

            for ( int i=0; i < d->cell_ln[mpi_n3d_id]; i++, Icell++ ) {
	            for ( int j=0; j < d->geo.AtomsPerCellMax(); j++){
	                d->geo.l_m_2_atomdsp[Icell][j][0] = d->geo.l_m_2_atomdsp_onedomain[i][j][0] ;
	                d->geo.l_m_2_atomdsp[Icell][j][1] = d->geo.l_m_2_atomdsp_onedomain[i][j][1] ;
	                d->geo.l_m_2_atomdsp[Icell][j][2] = d->geo.l_m_2_atomdsp_onedomain[i][j][2] ;
	        }}

	        rm_f3tensor(&d->geo.l_m_2_atomdsp_onedomain,0) ;
	        nml_memory_report("REPLACED THE TEMPORARY COPY WITH THE STRAIN ARRAY - AFTER MACOPT") ;
        }

        NOTE: THIS CODE CAUSES AN INCORRECT OPERATION OF THE PROGRAM
        IF THE USER CHOOSES TO CHECKPOINT THE DISPLACEMENTS
        */
        #endif /* DISPL_SINGLE_PRECISION */
      
        nml_memory_report("Strain_minimization- before multicast");
        multicastDisplacement(d);
        if(d->strain.first_or_second==1){
            write_nbr(d,1);  /*! write if d->opt.ExecParam.Strain.StrainOut.Neighborhood is on */
            writePos(d, "strain_first"); /*! write if d->opt.ExecParam.Strain.StrainOut.DisplacementData is on */
            //! cout_master <<" Printing Surface Atomic Position. "<<d->strain.first_or_second<<endl; 
            //! writePos(d, "surface_atom_pos");
        }
        else {
            write_nbr(d,2);  /*! write if d->opt.ExecParam.Strain.StrainOut.Neighborhood is on */
            writePos(d, "strain_second"); /*! write if d->opt.ExecParam.Strain.StrainOut.DisplacementData is on */
        }

        if (amac.filename_its)
            str_free(amac.filename_its), amac.filename_its=NULL;
        if (amac.filename_last)
            str_free(amac.filename_last), amac.filename_last=NULL;

        rm_rvectr(&p);
        rm_rvectr(&d->strain.dp);
      
        /*! release memory used for strain comm */
        nml_free(d->strain.sndL.data);
        nml_free(d->strain.sndR.data);
        nml_free(d->strain.rcvL.data);
        nml_free(d->strain.rcvR.data);
        nml_memory_report("Strain_minimization- ending");
    } /*! Comment_N: if VFF_Keating loop ends. */
    else { 
        cout_master << "No strain calculation requested "<<endl;
        d->strain.debugp = 0;
        return; 
    }
}


/*! set up d->lm_{snd,rcv}{L,R} data structure for doing sends and receives */
void initStrainComm(qd_struct d) 
{
/*! Comment_N: Hopefully the same communication arrays will also work after we eliminate lattice constant minimization. We will still be using same data structures, so this function should work just fine. I didn't go through this function carefully. */
    int countL=0, countR=0;
    int sz_bndryL=0, sz_bndryR=0;
    int sz_msg;
    int l, m, n, ii;
    int isx = d->cell_s[mpi_n3d_id];
    bool has_periodicity = d->geo.hasPeriodicity();
   
    //! if we are on edge, processor to past edge is NULL (-1) unless we
    //! have periodicity in x.
	//! N CPU ring for Periodic structure
	//! CPU(0) CPU(1) ... CPU(N-1)
	//! CPU(-1) = CPU(N-1) / CPU(N) = CPU(0)
    const int maxNbr = d->geo.NeighborsMax();
    int mpi_left = ( mpi_n3d_id 
                    ? mpi_n3d_id-1 
                    : (d->geo.isPeriodicX() ? mpi_n3d_numprocs-1 : -1) );
    int mpi_right = ( mpi_n3d_id < mpi_n3d_numprocs-1 
                     ? mpi_n3d_id+1 
                     : (d->geo.isPeriodicX() ? 0 : -1) );
    int left_min = (mpi_left>=0 ? d->seg_s[mpi_left] : 9999);
    int left_max = (mpi_left>=0 ? d->seg_s[mpi_left] + d->seg_ln[mpi_left] : -9999);
    int right_min = (mpi_right>=0 ? d->seg_s[mpi_right] : 9999);
    int right_max = (mpi_right>=0 ? d->seg_s[mpi_right] + d->seg_ln[mpi_right] : -9999);
   
    int szDatStr = sizeof(int) + sizeof(real);
    int szCommL = (mpi_left>=0 ? 6*maxNbr*d->commsize[mpi_n3d_id][mpi_left] : 1);
    int szCommR = (mpi_right>=0 ? 6*maxNbr*d->commsize[mpi_n3d_id][mpi_right] : 1);
   
    ivectr tmp0_sndL, tmp0_sndR, tmp_sndL, tmp_sndR, tmp_rcvL, tmp_rcvR;
   
    /*! FIX:  replace this with a better estimate!!! */
    d->strain.sndL.data = (void*) nml_malloc(szDatStr*szCommL);
    d->strain.sndL.N = 0;
    d->strain.sndR.data = (void*) nml_malloc(szDatStr*szCommR);
    d->strain.sndR.N = 0;
    d->strain.rcvL.data = (void*) nml_malloc(szDatStr*szCommL);
    d->strain.rcvL.N = 0;
    d->strain.rcvR.data = (void*) nml_malloc(szDatStr*szCommR);
    d->strain.rcvR.N = 0;

    #define MEMCHECK
    #ifdef MEMCHECK
        double strain_memory=0.0, tl_memory=0.0;
        strain_memory = 2* (szDatStr*szCommL + szDatStr*szCommR);
        tl_memory += strain_memory;
        printf("\n strain memory allocation\nmemory in MB strain_memory=%g tl_memory=%g\n",strain_memory,  tl_memory),fflush(stdout);
    #endif
   
    if (mpi_n3d_numprocs==1)  return; //! only set up if there will be communication
   
    /*! it's simplest to send equal sized-messages; set sz_msg to be max among procs */
    sz_msg=1;

    for (ii=0; ii<mpi_n3d_numprocs; ii++) {
        int sz_ii = d->cell_ln[ii] * d->geo.AtomsPerCellMax();
        if (sz_msg < sz_ii)  sz_msg = sz_ii;
    }

    sz_msg++; /*! allow for size of data in msg to be included in msg */
   
    tmp0_sndL = Ivectr(sz_msg);
    tmp0_sndR = Ivectr(sz_msg);
   
    tmp_sndL = Ivectr(sz_msg);
    tmp_sndR = Ivectr(sz_msg);
    tmp_rcvL = Ivectr(sz_msg);
    tmp_rcvR = Ivectr(sz_msg);
         
    for ( l=isx; l < isx + d->cell_ln[mpi_n3d_id]; l++ ) {
        /*! These are real space indexes for the unit cell */
        int i = d->geo.cell__ijk[l][0];
        int j = d->geo.cell__ijk[l][1];
        int k = d->geo.cell__ijk[l][2];

        /*! Loop through all the atoms in the unit cell within simulation domain */
        for ( m=0; m < d->geo.AtomsPerCellMax(); m++ ) {
            if ( !d->geo.AtomType[l][m] ) continue;
	 
            /*! Loop through the atom's nearest neighbors */
            for ( n = 0; n < d->geo.Neighbors(m); n++ ) {
                /*! Get relative position of unit cell where this neighbor lives. */ 
                int i_nbr = i + d->geo.NbrCell(m,n,0);
                int j_nbr = j + d->geo.NbrCell(m,n,1);
                int k_nbr = k + d->geo.NbrCell(m,n,2);
	            int aindx = d->geo.NbrCell(m,n,3);

                /*! Get the unit cell and atomic index for this neighbor */
                int cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

                #ifdef ELIMINATE_SSMAP
                if (cindx < 0 && has_periodicity && (cindx + d->geo.N_Cell + 1 ) >= 0 ) 
	                cindx = cindx + d->geo.N_Cell + 1 ; 
                if(d->geo.is_tilted)
                   cindx = d->geo.get_cindx_tilted(l,m,n);
                #else /*! ELIMINATE_SSMAP */
                if (/*cindx < 0 &&*/ has_periodicity && d->geo.ssmap[l][m][n]!=-1) /*! When crystal is tilted ssmap should be used for connecting atoms in cell with index = 0 as well */
                    cindx = d->geo.ssmap[l][m][n];
                #endif /*! ELIMINATE_SSMAP */

                /*! Check if neighbor's unit cell is within simulation domain. */
                if ( cindx >= 0 && d->geo.AtomType[cindx][aindx] ) {
                    int isn = d->geo.l_m_2_offset_ham_tot[cindx][aindx];

                    /*! Set nstat to {+1,-1} if neighbor is {upper,lower} triangular */ 
                    if (isn >= left_min && isn < left_max) 
                       tmp0_sndL[countL++] = m + l*d->geo.AtomsPerCellMax();

                    if (isn >= right_min && isn < right_max) 
                        tmp0_sndR[countR++] = m + l*d->geo.AtomsPerCellMax();
                }
            }
        }
    }

    if (countL==0) 
      mpi_left = -1;  /*! don't even bother sending to the left */
    else {
        /*! weed out duplicates -- left sends */
        quick_sort(&tmp0_sndL[0], 0, countL-1);
        sz_bndryL = 2;
        tmp_sndL[1] = tmp0_sndL[0];

        for (ii=1; ii<countL; ii++)
            if (tmp0_sndL[ii] != tmp0_sndL[ii-1])  
                tmp_sndL[sz_bndryL++] = tmp0_sndL[ii];
        tmp_sndL[0] = sz_bndryL-1;
    }
         
    if (countR==0) 
      mpi_right = -1;  /*! don't even bother sending to the right */
    else {
        /*! weed out duplicates -- right sends */
        quick_sort(&tmp0_sndR[0], 0, countR-1);
        sz_bndryR = 2;
        tmp_sndR[1] = tmp0_sndR[0];

        for (ii=1; ii<countR; ii++)
            if (tmp0_sndR[ii] != tmp0_sndR[ii-1])  
                tmp_sndR[sz_bndryR++] = tmp0_sndR[ii];
        tmp_sndR[0] = sz_bndryR-1;
    }
   
    /*! at this point, we have a temporary index mapping files {l,m}_{tmp_snd,tmp_rcv}{L,R}.
       we need to send this mapping to neighboring processors so that they know
       what to do with the actual displacement data they will receive. */
    if (mpi_n3d_id % 2) {
        if (mpi_left >= 0) {
            MPI_Sendrecv(&tmp_sndL[0], sz_msg, MPI_INT, mpi_left, 1,
                         &tmp_rcvL[0], sz_msg, MPI_INT, mpi_left, 0, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
        }

        if (mpi_right >= 0) {
            MPI_Sendrecv(&tmp_sndR[0], sz_msg, MPI_INT, mpi_right, 0,
                         &tmp_rcvR[0], sz_msg, MPI_INT, mpi_right, 1, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
        }
    }
    else {
        if (mpi_right >= 0) {
            MPI_Sendrecv(&tmp_sndR[0], sz_msg, MPI_INT, mpi_right, 0,
                         &tmp_rcvR[0], sz_msg, MPI_INT, mpi_right, 1, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
        }

        if (mpi_left >= 0) {
            MPI_Sendrecv(&tmp_sndL[0], sz_msg, MPI_INT, mpi_left, 1,
                         &tmp_rcvL[0], sz_msg, MPI_INT, mpi_left, 0, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
        }
    }
   
    /*! commsize should be same size as tmp_{rcv,snd}{R,L} -- we really didn't
       need to recompute this, but for now, we redo the calculation and then
       verify below.  Note that commsize need not be symmetric! */
    if ( (mpi_left>=0 ? d->commsize[mpi_n3d_id][mpi_left] : 0) != tmp_rcvL[0] ||
         (mpi_right >= 0 ? d->commsize[mpi_n3d_id][mpi_right] : 0) != tmp_rcvR[0] ||
         (mpi_left>=0 ? d->commsize[mpi_left][mpi_n3d_id] : 0) != tmp_sndL[0] ||
         (mpi_right >= 0 ? d->commsize[mpi_right][mpi_n3d_id] : 0) != tmp_sndR[0] ) {
        printf("%d => %d <= %d      %d <= %d => %d   (%d,%d)  (%d,%d)\n", 
             (mpi_left>=0 ? d->commsize[mpi_n3d_id][mpi_left] : 0), 
             mpi_n3d_id, 
             (mpi_right >= 0 ? d->commsize[mpi_n3d_id][mpi_right] : 0), 
             (mpi_left>=0 ? d->commsize[mpi_left][mpi_n3d_id] : 0), 
             mpi_n3d_id, 
             (mpi_right >= 0 ? d->commsize[mpi_right][mpi_n3d_id] : 0), 
             tmp_sndL[0], tmp_sndR[0],
             tmp_rcvL[0], tmp_rcvR[0]
        );
    }

    /*! now that we've done the sends/receives, we can fill up the permanent
       data structures in qd_struct d */
    d->strain.lm_sndL = Ivectr(tmp_sndL[0]);
    memcpy(&d->strain.lm_sndL[0], tmp_sndL + 1, sizeof(int)*tmp_sndL[0]);
   
    d->strain.lm_sndR = Ivectr(tmp_sndR[0]);
    memcpy(&d->strain.lm_sndR[0], tmp_sndR + 1, sizeof(int)*tmp_sndR[0]);
   
    d->strain.lm_rcvL = Ivectr(tmp_rcvL[0]);
    memcpy(&d->strain.lm_rcvL[0], tmp_rcvL + 1, sizeof(int)*tmp_rcvL[0]);
   
    d->strain.lm_rcvR = Ivectr(tmp_rcvR[0]);
    memcpy(&d->strain.lm_rcvR[0], tmp_rcvR + 1, sizeof(int)*tmp_rcvR[0]);
   
    #ifdef DISPL_SINGLE_PRECISION
    /*! analyze the structure of lm_rcvL and rm_rcvR :
     * - find the max and min cell indices
     * - find the overall max and min atom indices
     * - initialize the imatrices rcvL_lm and rcvR_lm
     * - fill those matrices with indices referring to
     *    the position of the corresponding entry in lm_rcvL and lm_rcvR, respectively
     */
   
    /*! first rcvL */
    if ( countL != 0 )
    {
        d->strain.cellLMin = d->strain.cellLMax = d->strain.lm_rcvL[0] / d->geo.AtomsPerCellMax() ;
        d->strain.atomLMin = d->strain.atomLMax = d->strain.lm_rcvL[0] % d->geo.AtomsPerCellMax() ;
       
        for ( ii = 0 ; ii < nml_iv_extent(d->strain.lm_rcvL) ; ii++ )
	    {
	        l = d->strain.lm_rcvL[ii] / d->geo.AtomsPerCellMax() ;
	        m = d->strain.lm_rcvL[ii] % d->geo.AtomsPerCellMax() ;

	        if ( l  < d->strain.cellLMin ) d->strain.cellLMin = l ;
	        if ( l  > d->strain.cellLMax ) d->strain.cellLMax = l ;
	        if ( m  < d->strain.atomLMin ) d->strain.atomLMin = m ;
	        if ( m  > d->strain.atomLMax ) d->strain.atomLMax = m ;
	    }
       
        if ( d->strain.cellLMax != d->strain.cellLMin )
	    {
	        d->strain.rcvL_lm = Imatrix( (d->strain.cellLMax - d->strain.cellLMin + 1) , 
					(d->strain.atomLMax - d->strain.atomLMin + 1) ) ;
	    }
        else
	        d->strain.rcvL_lm = NULL ;
       
        for ( l = 0 ; l < (d->strain.cellLMax - d->strain.cellLMin + 1) ; l++ ) {
	        for ( m = 0 ; m < (d->strain.atomLMax - d->strain.atomLMin + 1) ; m++ ) 
	            d->strain.rcvL_lm[l][m] = -1 ;
        }
       
        for ( ii = 0 ; ii < nml_iv_extent(d->strain.lm_rcvL) ; ii++ )
	    {
	        l = d->strain.lm_rcvL[ii] / d->geo.AtomsPerCellMax() ;
	        m = d->strain.lm_rcvL[ii] % d->geo.AtomsPerCellMax() ;
	   
	        l = l - d->strain.cellLMin ;
	        m = m - d->strain.atomLMin ;
	   
	        d->strain.rcvL_lm[l][m] = ii ;
	    }
    } else /*! i.e., if countL == 0 */
        d->strain.rcvL_lm = NULL ;
     
    /*! now rcvR */
    if ( countR != 0 )
    {
        d->strain.cellRMin = d->strain.cellRMax = d->strain.lm_rcvR[0] / d->geo.AtomsPerCellMax() ;
        d->strain.atomRMin = d->strain.atomRMax = d->strain.lm_rcvR[0] % d->geo.AtomsPerCellMax() ;
       
        for ( ii = 0 ; ii < nml_iv_extent(d->strain.lm_rcvR) ; ii++ )
	    {
	        l = d->strain.lm_rcvR[ii] / d->geo.AtomsPerCellMax() ;
	        m = d->strain.lm_rcvR[ii] % d->geo.AtomsPerCellMax() ;

	        if ( l  < d->strain.cellRMin ) d->strain.cellRMin = l ;
	        if ( l  > d->strain.cellRMax ) d->strain.cellRMax = l ;
	        if ( m  < d->strain.atomRMin ) d->strain.atomRMin = m ;
	        if ( m  > d->strain.atomRMax ) d->strain.atomRMax = m ;
	    }
       
        if ( d->strain.cellRMax != d->strain.cellRMin )
	    {
	        d->strain.rcvR_lm = Imatrix( (d->strain.cellRMax - d->strain.cellRMin + 1) , 
					(d->strain.atomRMax - d->strain.atomRMin + 1) ) ;
	    }
        else
	        d->strain.rcvR_lm = NULL ;
       
        for ( l = 0 ; l < (d->strain.cellRMax - d->strain.cellRMin + 1) ; l++ ) {
	        for ( m = 0 ; m < (d->strain.atomRMax - d->strain.atomRMin + 1) ; m++ ) 
	            d->strain.rcvR_lm[l][m] = -1 ;
        }
       
        for ( ii = 0 ; ii < nml_iv_extent(d->strain.lm_rcvR) ; ii++ )
	    {
	        l = d->strain.lm_rcvR[ii] / d->geo.AtomsPerCellMax() ;
	        m = d->strain.lm_rcvR[ii] % d->geo.AtomsPerCellMax() ;
	   
	        l = l - d->strain.cellRMin ;
	        m = m - d->strain.atomRMin ;
	   
	        d->strain.rcvR_lm[l][m] = ii ;
	    }
    } else /*! i.e., if countR == 0 */
        d->strain.rcvR_lm = NULL ;
    #endif /*! DISPL_SINGLE_PRECISION */
   
    rm_ivectr(&tmp0_sndL);
    rm_ivectr(&tmp0_sndR);
   
    rm_ivectr(&tmp_sndL);
    rm_ivectr(&tmp_sndR);
    rm_ivectr(&tmp_rcvL);
    rm_ivectr(&tmp_rcvR);
}

void calc_estrain_macopt( rvectr p, int sz, rvectr xi, void *data, int itercount, real *onevalue ) 
{
/*! Comment_N: This function is called from macoptII_par(). macoptII_par is called as follows:             macoptII_par(p, n_loc, n_tot, (calc_estrain_macopt), (void*) d, &amac);
Corresponding argument in macoptII_par() is 
      void (*dfunc)(nml_dvector*, int, nml_dvector*, void*, int, double*).
calc_estrain_macopt() calculates gradient. Actually this function calls function calc_estrain() which does final gradient calculation and stores it in d->strain.dp.
calc_estrain_macopt() takes care of boundary conditions.
*/
    int i, ioff, j, k, id, first, thislength ;
    real es;
    qd_struct d = (qd_struct)(data);
    double SCALE_FCTR = d->strain.SCALE_FCTR;

    //! are we minimizing the lattice constant in at least one direction?
    const bool do_min =d->strain.doPeriodMinimization();

    rvectr sndL ;
    rvectr sndR ;
    rvectr rcvL ;
    rvectr rcvR ;

    ioff = d->cell_s[mpi_n3d_id]; 

/*! Comment_N: Apply boundary conditions to vector p. Vector p has values calculated in previous iteration. Only those components which are allowed by the boundary conditions are used. All remaining elements are not used for updating displacements and/or lattice constants. */
    switch (d->opt.ExecParam.Strain.Keating.BoundCond) {
        case Keating_struct::FixOne:
            first = (mpi_n3d_id==mpi_n3d_masterid ? 1 : 0);
            for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
                for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
	                if ( d->geo.AtomType[ioff][j] ){
	                    if ( first == 0 ){
	                        for ( k=0; k < 3; k++ ) {
		                        id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
		                        d->geo.l_m_2_atomdsp[ioff][j][k] = p[id];
	                        }
	                    } else {
	                        first = 0;

                            #ifdef DISPL_SINGLE_PRECISION
	                        for ( k=0; k < 3; k++ ) {
		                        id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;

		                        /* p[id] = d->geo.l_m_2_atomdsp_onedomain[i][j][k] ; */
		                        p[id] = d->geo.l_m_2_atomdsp[ioff][j][k] ;
                            }
                            #endif /*! DISPL_SINGLE_PRECISION */
                        }
	                }

                    #ifdef DISPL_SINGLE_PRECISION
	                else {
	                    for ( k=0; k < 3; k++ ) {
	                        id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;

	                        /* p[id] = d->geo.l_m_2_atomdsp_onedomain[i][j][k] ; */
	                        p[id] = d->geo.l_m_2_atomdsp[ioff][j][k] ;
	                    }
	                }
                    #endif /*! DISPL_SINGLE_PRECISION */
                }
			}
        break;

        case Keating_struct::AllFree:
/*! Comment_N: All atoms are allowed to displace from their original positions. i.e. all d->geo.l_m_2_atomdsp[ioff][j][k] are allowed to change as per p[id] calculated from CG routine macoptII_par. */
            for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
                for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
	                if ( d->geo.AtomType[ioff][j] ){
	                    for ( k=0; k < 3; k++ ){
	                        id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;

	                        /* d->geo.l_m_2_atomdsp_onedomain[i][j][k] = p[id]; */
	                        d->geo.l_m_2_atomdsp[ioff][j][k] = p[id];
	                    }
	                }
                }
			}
        break;

        case Keating_struct::FixSurfaceAll:
/*! Comment_N: All atoms except those on the surface are allowed to displace. */
            for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
                for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
	                if ( d->geo.AtomType[ioff][j] && !d->geo.isOnSurface(ioff,j) ) {
	                    for ( k=0; k < 3; k++ ){
	                        id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;

	                        /* d->geo.l_m_2_atomdsp_onedomain[i][j][k] = p[id];*/
	                        d->geo.l_m_2_atomdsp[ioff][j][k] = p[id];
	                    }
	                }

                    #ifdef DISPL_SINGLE_PRECISION
	                else {
	                    for ( k=0; k < 3; k++ ){
	                        id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;

	                        /* p[id] = d->geo.l_m_2_atomdsp_onedomain[i][j][k] ;*/
	                        p[id] = d->geo.l_m_2_atomdsp[ioff][j][k] ;
	                    }
	                }
                    #endif /*! DISPL_SINGLE_PRECISION */
                }
			}
        break;

   case Keating_struct::FixZmin:
/*! Comment_N: All atoms except those in the bottom-most (in zidirection) unit cell are allowed to displace. */
      for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
	 if ( d->geo.AtomType[ioff][j] && !d->geo.isInZ_min(ioff,j) ) {
	    for ( k=0; k < 3; k++ ){
	       id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
#ifdef DISPL_SINGLE_PRECISION
	       /*	       d->geo.l_m_2_atomdsp_onedomain[i][j][k] = p[id]; */
	       d->geo.l_m_2_atomdsp[ioff][j][k] = p[id];
#else /*! DISPL_SINGLE_PRECISION */
	       d->geo.l_m_2_atomdsp[ioff][j][k] = p[id];
#endif /*! DISPL_SINGLE_PRECISION */
	    }
	 }
#ifdef DISPL_SINGLE_PRECISION
	 else {
	    for ( k=0; k < 3; k++ ){
	       id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
	       /*	       p[id] = d->geo.l_m_2_atomdsp_onedomain[i][j][k] ;*/
	       p[id] = d->geo.l_m_2_atomdsp[ioff][j][k] ;
	    }
	 }
#endif /*! DISPL_SINGLE_PRECISION */
      }}
      break;
   case Keating_struct::FixInterior:
/*! Comment_N: Only surface atoms are allowed to displace. */
      for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
	 if ( d->geo.AtomType[ioff][j] && d->geo.isOnSurface(ioff,j) ) {
	    for ( k=0; k < 3; k++ ){
	       id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
#ifdef DISPL_SINGLE_PRECISION
	       /*	       d->geo.l_m_2_atomdsp_onedomain[i][j][k] = p[id]; */
	       d->geo.l_m_2_atomdsp[ioff][j][k] = p[id];
#else /*! DISPL_SINGLE_PRECISION */
	       d->geo.l_m_2_atomdsp[ioff][j][k] = p[id];
#endif /*! DISPL_SINGLE_PRECISION */
	    }
	 }
#ifdef DISPL_SINGLE_PRECISION
	 else {
	    for ( k=0; k < 3; k++ ){
	       id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
	       /* p[id] = d->geo.l_m_2_atomdsp_onedomain[i][j][k] ;*/
	       p[id] = d->geo.l_m_2_atomdsp[ioff][j][k] ;
	    }
	 }
#endif /*! DISPL_SINGLE_PRECISION */
      }}
      break;
   default:
      die("Unidentified Boundary Condition in calc_estrain_macopt\n");
      break;
   }


#ifdef Test_Displacement_Data
  if(itercount==0){
      ioff = d->cell_s[mpi_n3d_id];
      ofstream fout("displacement.dat");
      for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) {
         if ( d->geo.AtomType[ioff][j] ) {
            for ( k=0; k < 3; k++ ){
               id = ioff*d->geo.AtomsPerCellMax()*3 + j*3 + k;
#ifdef DISPL_SINGLE_PRECISION
               fout<< id<< " "<<d->geo.l_m_2_atomdsp[ioff][j][k]<<endl;
	       /*               fout<< id<< " "<<d->geo.l_m_2_atomdsp_onedomain[i][j][k]<<endl; */
#else /*! DISPL_SINGLE_PRECISION */
               fout<< id<< " "<<d->geo.l_m_2_atomdsp[ioff][j][k]<<endl;
#endif /*! DISPL_SINGLE_PRECISION */

            }
         }
      }}
  } 
#endif

   //! master cpu holds new value of lattice constant; broadcast to other cpus
   if (_new_minimize && do_min) {
/*! Comment_N: If do_min is true (i.e. lattice const minimization is tured on) then lattice constant latt[0-2] is allowed to vary. It may appear that all three lattice constants are allowed to vary. latt[0-2] has values of new lattice constants but only those lattice constants which are permitted to vary according to boundary conditions are updated. This is done below, please scroll down. */
      real latt[3];
	 
      if (!mpi_n3d_id) {
	 real* p_ptr = &p[3*d->cell_ln[mpi_n3d_id]*d->geo.AtomsPerCellMax()];
	 latt[0] = *(p_ptr++) * SCALE_FCTR;
	 latt[1] = *(p_ptr++) * SCALE_FCTR;
	 latt[2] = *(p_ptr++) * SCALE_FCTR;
      }
      
      MPI_Bcast(&latt[0], 3, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);
      
      /*! if we are running the special case of Fixed interior pts, we 
         want to make sure that interior atoms are not moved by a change 
         of lattice constant */
      if (d->opt.ExecParam.Strain.Keating.BoundCond==Keating_struct::FixInterior) {
         ioff=d->cell_s[mpi_n3d_id];
         for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
         for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
            if ( !d->geo.AtomType[ioff][j] && !d->geo.isOnSurface(ioff,j) )  continue;
#ifdef DISPL_SINGLE_PRECISION
	    /*            d->geo.l_m_2_atomdsp_onedomain[i][j][0] += 
               (d->geo.lattice_x - latt[0]) * 
               ( d->geo.cell__ijk[ioff][0] + d->geo.PositionInCell(j,0));
            d->geo.l_m_2_atomdsp_onedomain[i][j][1] += 
               (d->geo.lattice_y - latt[1]) * 
               ( d->geo.cell__ijk[ioff][1] + d->geo.PositionInCell(j,1));
            d->geo.l_m_2_atomdsp_onedomain[i][j][2] += 
               (d->geo.lattice_z - latt[2]) * 
               ( d->geo.cell__ijk[ioff][2] + d->geo.PositionInCell(j,2)); */

            d->geo.l_m_2_atomdsp[ioff][j][0] += 
               (d->geo.lattice_x - latt[0]) * 
               ( d->geo.cell__ijk[ioff][0] + d->geo.PositionInCell(j,0));
            d->geo.l_m_2_atomdsp[ioff][j][1] += 
               (d->geo.lattice_y - latt[1]) * 
               ( d->geo.cell__ijk[ioff][1] + d->geo.PositionInCell(j,1));
            d->geo.l_m_2_atomdsp[ioff][j][2] += 
               (d->geo.lattice_z - latt[2]) * 
               ( d->geo.cell__ijk[ioff][2] + d->geo.PositionInCell(j,2));

#else /*! DISPL_SINGLE_PRECISION */
            d->geo.l_m_2_atomdsp[ioff][j][0] += 
               (d->geo.lattice_x - latt[0]) * 
               ( d->geo.cell__ijk[ioff][0] + d->geo.PositionInCell(j,0));
            d->geo.l_m_2_atomdsp[ioff][j][1] += 
               (d->geo.lattice_y - latt[1]) * 
               ( d->geo.cell__ijk[ioff][1] + d->geo.PositionInCell(j,1));
            d->geo.l_m_2_atomdsp[ioff][j][2] += 
               (d->geo.lattice_z - latt[2]) * 
               ( d->geo.cell__ijk[ioff][2] + d->geo.PositionInCell(j,2));
#endif /*! DISPL_SINGLE_PRECISION */
         }}
      }

/*! Comment_N: Impose boundary conditions on lattice constant. Only those lattice constants which are permitted to vary accoring to the boundary condtions are updated. Variables MinErtX, Y or Z are set from the input deck. */
      if (d->strain.MinWrtX) d->geo.lattice_x = latt[0];
      if (d->strain.MinWrtY) d->geo.lattice_y = latt[1];
      if (d->strain.MinWrtZ) d->geo.lattice_z = latt[2];
      
   }

   /*! send displacement data to neighbors (NB: this algorithm suffers 
      if there are an odd number of procs */
   if (mpi_n3d_numprocs > 1) {
/*! Comment_N: Assign values to the communication arrays. */
      //! if we are on edge, processor to past edge is NULL (-1) unless we
      //! have periodicity in x.
      int mpi_left = ( mpi_n3d_id 
                       ? mpi_n3d_id-1 
                       : (d->geo.isPeriodicX() ? mpi_n3d_numprocs-1 : -1) );
      int mpi_right = ( mpi_n3d_id < mpi_n3d_numprocs-1 
                        ? mpi_n3d_id+1 
                        : (d->geo.isPeriodicX() ? 0 : -1) );

      sndL = Rvectr(3*nml_iv_extent(d->strain.lm_sndL));
      sndR = Rvectr(3*nml_iv_extent(d->strain.lm_sndR));
      rcvL = Rvectr(3*nml_iv_extent(d->strain.lm_rcvL));
      rcvR = Rvectr(3*nml_iv_extent(d->strain.lm_rcvR));
      real* ptr;


      /*! copy data to be sent into a flat array */
      ptr = &sndL[0];
      thislength = nml_iv_extent(d->strain.lm_sndL) ; 
      for (i=0; i<thislength; i++) {
         int l = d->strain.lm_sndL[i] / d->geo.AtomsPerCellMax();
         int m = d->strain.lm_sndL[i] % d->geo.AtomsPerCellMax();

#ifdef DISPL_SINGLE_PRECISION
	 int cellInThisDomain = l - d->cell_s[mpi_n3d_id] ;
	 id = cellInThisDomain * d->geo.AtomsPerCellMax() * 3 + m*3 ;
         real* dsp = &(p[id]);
#else /*! DISPL_SINGLE_PRECISION */
         real* dsp = &(d->geo.l_m_2_atomdsp[l][m][0]);
#endif /*! DISPL_SINGLE_PRECISION */

         *ptr++ = *dsp++;
         *ptr++ = *dsp++;
         *ptr++ = *dsp++;
      }

      ptr = &sndR[0];
      thislength = nml_iv_extent(d->strain.lm_sndR) ; 
      for (i=0; i<thislength; i++) {
         int l = d->strain.lm_sndR[i] / d->geo.AtomsPerCellMax();
         int m = d->strain.lm_sndR[i] % d->geo.AtomsPerCellMax();
#ifdef DISPL_SINGLE_PRECISION
	 int cellInThisDomain = l - d->cell_s[mpi_n3d_id] ;
	 id = cellInThisDomain * d->geo.AtomsPerCellMax() * 3 + m*3 ;
         real* dsp = &(p[id]);
#else /*! DISPL_SINGLE_PRECISION */
         real* dsp = &(d->geo.l_m_2_atomdsp[l][m][0]);
#endif /*! DISPL_SINGLE_PRECISION */
         *ptr++ = *dsp++;
         *ptr++ = *dsp++;
         *ptr++ = *dsp++;
      }

      if (mpi_n3d_id % 2) {
         if (mpi_left >= 0) {
            MPI_Sendrecv(&sndL[0], nml_dv_extent(sndL), MPI_DOUBLE, mpi_left, 1,
                         &rcvL[0], nml_dv_extent(rcvL), MPI_DOUBLE, mpi_left, 0, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
         }
         if (mpi_right >= 0) {
            MPI_Sendrecv(&sndR[0], nml_dv_extent(sndR), MPI_DOUBLE, mpi_right, 0,
                         &rcvR[0], nml_dv_extent(rcvR), MPI_DOUBLE, mpi_right, 1, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
         }
      }
      else {
         if (mpi_right >= 0) {
            MPI_Sendrecv(&sndR[0], nml_dv_extent(sndR), MPI_DOUBLE, mpi_right, 0,
                         &rcvR[0], nml_dv_extent(rcvR), MPI_DOUBLE, mpi_right, 1, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
         }
         if (mpi_left >= 0) {
            MPI_Sendrecv(&sndL[0], nml_dv_extent(sndL), MPI_DOUBLE, mpi_left, 1,
                         &rcvL[0], nml_dv_extent(rcvL), MPI_DOUBLE, mpi_left, 0, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
         }
      }
   
#ifndef DISPL_SINGLE_PRECISION
      ptr = &rcvL[0];
      thislength = nml_iv_extent(d->strain.lm_rcvL) ; 
      for (i=0; i<thislength; i++) {
         int l = d->strain.lm_rcvL[i] / d->geo.AtomsPerCellMax();
         int m = d->strain.lm_rcvL[i] % d->geo.AtomsPerCellMax();
         real* dsp = &(d->geo.l_m_2_atomdsp[l][m][0]);
         *dsp++ = *ptr++;
         *dsp++ = *ptr++;
         *dsp++ = *ptr++;
      }

      ptr = &rcvR[0];
      thislength = nml_iv_extent(d->strain.lm_rcvR) ; 
      for (i=0; i<thislength; i++) {
         int l = d->strain.lm_rcvR[i] / d->geo.AtomsPerCellMax();
         int m = d->strain.lm_rcvR[i] % d->geo.AtomsPerCellMax();
         real* dsp = &(d->geo.l_m_2_atomdsp[l][m][0]);
         *dsp++ = *ptr++;
         *dsp++ = *ptr++;
         *dsp++ = *ptr++;
      }
#endif /*! DISPL_SINGLE_PRECISION */
      /*! Note that if defined DISPL_SINGLE_PRECISION, displacements are only stored
       * in an array encompassing only the atoms on current CPU, so it is not possible
       * to store values from rcvL and rcvR */

   }
   else /*!  there is only 1 CPU */
     {
       rcvL = NULL ;
       rcvR = NULL ;
     }
   
   *onevalue = es = calc_estrain(d, p, rcvL, rcvR);
/*! Comment_N: Function calc_estrain() does the actual gradient and strain energy calculation. */

   if (mpi_n3d_numprocs > 1) {
      rm_rvectr(&sndL);
      rm_rvectr(&sndR);
      rm_rvectr(&rcvL);
      rm_rvectr(&rcvR);
   }

   if (d->opt.ExecParam.Strain.StrainOut.EnergyCGall && mpi_n3d_id==mpi_n3d_masterid){
      FILE *fp;
      char *filename = n3d_strdup_n(d->inputfile);
      n3d_FileTypeSet(&filename, "nd_estrain_all", TRUE);
      fp = fopen( filename, "a" );
      str_free(filename);
      fprintf(fp,"%d %g\n",itercount, es);
      fclose(fp);
   }

   if(mpi_n3d_id==mpi_n3d_masterid){
       FILE *fp = NULL;
       char* filename = n3d_strdup_n(d->inputfile);
       n3d_FileTypeSet(&filename, "nd_lattice_iter", TRUE);
       fp  = fopen(filename, "a" );
       str_free(filename);
       fprintf(fp,"%d %e %e %e \n",
                   itercount, d->geo.lattice_x, d->geo.lattice_y, d->geo.lattice_z);
       fclose(fp);
   }


/*! Comment_N: Now gradient and strain energy calculations are done. We update only those components of vector xi (which is the gradient vector supplied to CG routine) with calculated components of strain.dp vector which are allowed to be used in minimization (as per the boundary condtions). */
   ioff = d->cell_s[mpi_n3d_id]; 
   switch (d->opt.ExecParam.Strain.Keating.BoundCond) {
   case Keating_struct::FixOne:
      first = (mpi_n3d_id==mpi_n3d_masterid ? 1 : 0);
      for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
	 if ( d->geo.AtomType[ioff][j] ){
	    if ( first == 0 ){
	       for ( k=0; k < 3; k++ ){
		  id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
		  xi[id] = d->strain.dp[id];
	       }
	    } else{
	       for ( k=0; k < 3; k++ ){
		  id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
		  xi[id] = 0.0;
	       }
	       first = 0;
	    }
	 }
      }}
      break;
   case Keating_struct::AllFree:
      for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
	 if ( d->geo.AtomType[ioff][j] ){
	    for ( k=0; k < 3; k++ ){
	       id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
	       xi[id] = d->strain.dp[id];
	    }
	 }
      }}
      break;
   case Keating_struct::FixSurfaceAll:
      for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
	 if ( d->geo.AtomType[ioff][j]){
	    if (!d->geo.isOnSurface(ioff,j)) {
	       for ( k=0; k < 3; k++ ) {
		  id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
		  xi[id] = d->strain.dp[id];
	       }
	    } else {
	       for ( k=0; k < 3; k++ ){
		  id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
		  xi[id] = 0.0;
	       }
	    }
	 }
      }}
      break;
   case Keating_struct::FixZmin:
      for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
	 if ( d->geo.AtomType[ioff][j]){
	    if (!d->geo.isInZ_min(ioff,j)) {
	       for ( k=0; k < 3; k++ ) {
		  id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
		  xi[id] = d->strain.dp[id];
	       }
	    } else {
	       for ( k=0; k < 3; k++ ){
		  id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
		  xi[id] = 0.0;
	       }
	    }
	 }
      }}
      break;
   case Keating_struct::FixInterior:
      for ( i=0; i < d->cell_ln[mpi_n3d_id]; i++, ioff++ ) {
      for ( j=0; j < d->geo.AtomsPerCellMax(); j++ ) { 
	 if ( d->geo.AtomType[ioff][j]){
	    if ( d->geo.isOnSurface(ioff,j)) {
	       for ( k=0; k < 3; k++ ) {
		  id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
		  xi[id] = d->strain.dp[id];
	       }
	    } else {
	       for ( k=0; k < 3; k++ ){
		  id = i*d->geo.AtomsPerCellMax()*3 + j*3 + k;
		  xi[id] = 0.0;
	       }
	    }
	 }
      }}
      break;
   default:
      die("Unidentified Boundary Condition in calc_estrain_macopt\n");
      break;
   }

   if (_new_minimize && !mpi_n3d_id && do_min) {
      int offset = 3*d->cell_ln[mpi_n3d_id]*d->geo.AtomsPerCellMax();
      real* xi_0 = &xi[offset];
      real* dp_0 = &(d->strain.dp[offset]);
      
      if (d->strain.MinWrtX) xi_0[0] = dp_0[0];
      if (d->strain.MinWrtY) xi_0[1] = dp_0[1];
      if (d->strain.MinWrtZ) xi_0[2] = dp_0[2];
   }

   fflush(stdout);
    
   if (itercount && d->opt.ExecParam.Strain.Keating.strainCheckPointInterval && 
       (itercount % d->opt.ExecParam.Strain.Keating.strainCheckPointInterval) == 0) {
      cout_master << "Checkpointing strain calculation:\n";
      multicastDisplacement(d);
      writePos(d, "strain-checkpoint");
   }

}


real calc_estrain( qd_struct d , rvectr displacements, rvectr rcvL, rvectr rcvR )
{
/*! Comment_N: This function does the actual gradient and strain energy calculation. */
#if (defined MPI3d && !defined FAKE_MPI)
   real Einfo_tot[4];
#endif

   double SCALE_FCTR = d->strain.SCALE_FCTR;
   int AtomType_cindx_aindx;
   int i, j, k, l, m, n, p, q, aindx=-1, cindx=-1, is0, isn;
   int cn_L, cn_R, cn_this, cp_L, cp_R, cp_this , wrapped ;
   real alat, d0dot, dsdot, den, dssq, d0sq, bldl, bldc, bbdl, bbdc;
   real Eloc, Etot;
   real a_lattice[3];

   real bld=0.0, bbd=0.0;

   real* nnv_p=NULL;
   real* nnv_n=NULL;
   real* nnv0_p=NULL;
   real* nnv0_n=NULL;
   
   int sz_loc = 3*d->cell_ln[mpi_n3d_id]*d->geo.AtomsPerCellMax() +
      ( _new_minimize && !mpi_n3d_id && d->strain.doPeriodMinimization()
        ? 3 : 0);
#ifdef INTEL_MODIFIED
__declspec(align(16)) real dE_dlatt[3], dE_dlatt_tot[3];
#else   
   real dE_dlatt[3], dE_dlatt_tot[3];
#endif
   real Einfo_loc[4];

   int Nsleft=0, Nsright=0, Nrleft=0, Nrright=0;

   real bbdc_bbdl_2;
   real r1,r2;
   int i1, i3;
   
   bool has_periodicity = d->geo.hasPeriodicity();

   //! if we are on edge, processor to past edge is NULL (-1) unless we
   //! have periodicity in x.
   const int maxNbr = d->geo.NeighborsMax();
   int mpi_left = ( mpi_n3d_id 
                    ? mpi_n3d_id-1 
                    : (d->geo.isPeriodicX() ? mpi_n3d_numprocs-1 : -1) );
   int mpi_right = ( mpi_n3d_id < mpi_n3d_numprocs-1 
                     ? mpi_n3d_id+1 
                     : (d->geo.isPeriodicX() ? 0 : -1) );

   int cell0 = d->cell_s[mpi_n3d_id];
   int cellF = cell0 + d->cell_ln[mpi_n3d_id];
   int cell0_L = (mpi_left>=0 ? d->cell_s[mpi_left] : -1);
   int cell0_R = (mpi_right>=0 -1 ? d->cell_s[mpi_right] : -1);

   int imin = d->seg_s[mpi_n3d_id];
   int imax = d->seg_s[mpi_n3d_id] + d->seg_ln[mpi_n3d_id];

   int szDatStr = sizeof(int) + sizeof(real);
   //! pad with one extra element; probably not necessary
   int szCommL = 1 + (mpi_left>=0 ? 6*maxNbr*d->commsize[mpi_n3d_id][mpi_left] : 0);
   int szCommR = 1 + (mpi_right>=0 ? 6*maxNbr*d->commsize[mpi_n3d_id][mpi_right] : 0);
 
   rvectr sleft       = Rvectr(szCommL);
   rvectr sright      = Rvectr(szCommR);
   ivectr indx_sleft  = Ivectr(szCommL);
   ivectr indx_sright = Ivectr(szCommR);

   rvectr  alpha = Rvectr(maxNbr);
   rvectr  beta = Rvectr(maxNbr);
   real Beta;//! effective beta for deformed doublet of bonds
   ivectr  nstat = Ivectr(maxNbr);
#ifdef INTEL_MODIFIED
#pragma message ("INTEL_MODIFIED: Strain nnv[] allign ON")
   rmatrix nnv;
   double* tmp[maxNbr];
   for (int i=0; i<maxNbr; i++){
   tmp[i] = (double*)_mm_malloc(3*sizeof(real), 16);
   }
   nnv = (double**) &tmp;
#else
   rmatrix nnv = Rmatrix(maxNbr,3);
#endif
   rmatrix nnv0 = Rmatrix(maxNbr,3);
   imatrix indxmat = Imatrix(maxNbr,2);
   rmatrix eta = Rmatrix(maxNbr,3);

   dE_dlatt[0] = dE_dlatt[1] = dE_dlatt[2] = 0.0;
   dE_dlatt_tot[0] = dE_dlatt_tot[1] = dE_dlatt_tot[2] = 0.0;

/*! Comment_N: Current value of the lattice constant. Note that we updated d->geo.lattice_x etc in calc_estrain_macopt() */
   a_lattice[0]=d->geo.lattice_x;
   a_lattice[1]=d->geo.lattice_y;
   a_lattice[2]=d->geo.lattice_z;
   
   for ( i=0; i < sz_loc; i++ )
      d->strain.dp[i] = 0.0;

/*! Comment_N: Following loop performs gradient and strain energy calculation. */   
   /*! Loop through all unit cells associated with this processor */
   for ( l=cell0; l < cellF; l++ ) {
      /*! These are real space indexes for the unit cell */
      i = d->geo.cell__ijk[l][0];
      j = d->geo.cell__ijk[l][1];
      k = d->geo.cell__ijk[l][2];

#ifdef _DBUG_SEGFLT
      if (!mpi_n3d_id)
         printf("%d...", l), fflush(stdout);
      errorCheck(d, "A");
#endif

      
      /*! Loop through all the atoms in the unit cell within simulation domain */
      for ( m=0; m < d->geo.AtomsPerCellMax(); m++ ) {
	 int AtomType_l_m = (int) d->geo.AtomType[l][m];

	 if ( !AtomType_l_m ) continue;
	 
	 is0 = d->geo.l_m_2_offset_ham_tot[l][m];
	 /*! Loop through the atom's nearest neighbors */
	 for ( n = 0; n < d->geo.Neighbors(m); n++ ) {
/*! Comment_N: This loop calculated nearest neighbor vectors (nnv's). nnv's calculated in this loop are used in later loops over nearest neighbors. */	    
	    /*! Get relative position of unit cell where this neighbor lives. */ 
	    int i_nbr = i + d->geo.NbrCell(m,n,0);
	    int j_nbr = j + d->geo.NbrCell(m,n,1);
	    int k_nbr = k + d->geo.NbrCell(m,n,2);
	    
	    /*! Initialize neighbor status to zero.  This vector is used to insure
	       that double counting does not occur.  It is set to 0 if the neigbor
	       doesn't exist, one of the neighbor's hamiltonian matrix element is
	       upper triangular and  -1 if the matrix element is lower triangular.  By
	       including only upper or lower triangular contributions we insure that
	       bond length distortion contributions are only counted once in the strain
	       energy sum. */
	    nstat[n] = 0;
	    
	    /*! Get the unit cell and atomic index for this neighbor */
	    cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];

#ifdef ELIMINATE_SSMAP
            if (cindx < 0 && has_periodicity && (cindx + d->geo.N_Cell + 1 ) >= 0 ) {
	      cindx = cindx + d->geo.N_Cell + 1 ; 
            }
            if(d->geo.is_tilted)
               cindx = d->geo.get_cindx_tilted(l,m,n);
#else /*! ELIMINATE_SSMAP */
            if (/*cindx < 0 &&*/ has_periodicity && d->geo.ssmap[l][m][n]!=-1) {/*! When crystal is tilted ssmap should be used for connecting atoms in cell with index = 0 as well */
               cindx = d->geo.ssmap[l][m][n];
            }
#endif /*! ELIMINATE_SSMAP */

	    aindx = d->geo.NbrCell(m,n,3);

	    /*! Check if neighbor's unit cell is within simulation domain. */
	    if ( cindx >= 0 ) {

	       AtomType_cindx_aindx = (int) d->geo.AtomType[cindx][aindx];
	       
	       /*! Check if the neighbor is within the simulation domain */
	       if ( AtomType_cindx_aindx ) {
#ifdef DISPL_SINGLE_PRECISION
/*! Comment_N: DISPL_SINGLE_PRECISION is defined in geo_struct.c. Therefore following code will be executed. */
		 /*! Here displacements are transmitted via displacements, rcvL , and rcvR.
		  * We need to assign the atoms appropriately
		  * We know that the cell number l is within the domain on this processor
		  * but we do not know about the cell number cindx
		  * and we need to check it
		  */

		 /*! First deal with the atom [l][m] */
		 int cellIndexInThisDomain = l - cell0 ;
		 int entryNumber = cellIndexInThisDomain * d->geo.AtomsPerCellMax()*3 + m*3 ;
		 real * rv2 = &(displacements[entryNumber]) ;
/*! Comment_N: *rv2 = l_m_2_atomdsp for atom [l][m] */
		 /*! There. Now find the atom [cindx][aindx].
		  * It can be within this domain (on this CPU) or it can be
		  * inside rcvL or rcvR */
		 real * rv1 ;
		 int theEntry ;
		 if ( (cindx >= cell0) && (cindx < cellF)  )
		   {  /*! this neighbor is on this CPU and we find it in displacements */

		     int cellIndexInThisDomain = cindx - cell0 ;
		     int entryNumber = cellIndexInThisDomain * d->geo.AtomsPerCellMax()*3 + aindx*3 ;
		     rv1 = &(displacements[entryNumber]) ;
/*! Comment_N: *rv1 = l_m_2_atomdsp for nbr [cindx][aindx] */

		   }
		 else if ( ( cindx >= d->strain.cellLMin ) && ( cindx <= d->strain.cellLMax ) ) 
		   { /*! this neighbor is on the left neighboring CPU and we find it in rcvL */

		     theEntry =  d->strain.rcvL_lm[cindx - d->strain.cellLMin][aindx - d->strain.atomLMin] ;
		     if ( ( theEntry >= 0 ) && ( d->strain.rcvL_lm != NULL ) ) 
		       {
			 rv1 = &(rcvL[theEntry * 3]) ;		       
/*! Comment_N: *rv1 = l_m_2_atomdsp for nbr [cindx][aindx] */
		       }
		     else
		       { 
			 die("Problem with finding neighbors when calculating estrain, rcvL") ;
		       }
		   }
		 else if ( ( cindx >= d->strain.cellRMin ) && ( cindx <= d->strain.cellRMax ) ) 
		   { /*! this neighbor is on the right neighboring CPU and we find it in rcvR */

		     theEntry =  d->strain.rcvR_lm[cindx - d->strain.cellRMin][aindx - d->strain.atomRMin] ;
		     if ( ( theEntry >= 0 ) && (d->strain.rcvR_lm != NULL )  ) 
		       {
			 rv1 = &(rcvR[theEntry * 3]) ;
/*! Comment_N: *rv1 = l_m_2_atomdsp for nbr [cindx][aindx] */
		       }
		     else
		       { 
			 die("Problem with finding neighbors when calculating estrain, rcvR") ;
		       }
		   }

#else /*! DISPL_SINGLE_PRECISION */
/*! Comment_N: DISPL_SINGLE_PRECISION is defined in geo_struct.c. Therefore these two statements will not be executed. */
		  real* rv1 = d->geo.l_m_2_atomdsp[cindx][aindx];
		  real* rv2 = d->geo.l_m_2_atomdsp[l][m];
#endif /*! DISPL_SINGLE_PRECISION */
		  real* rv3 = d->geo.PositionInCell(aindx);
/*! Comment_N: *rv3 = Position of nbr in its unit cell (i.e. cell with index cindx.) PositionInCell have values like (0,0,0), (0,-0.5,-0.5), (-0.25,-0.25,-0.25) etc */
		  real* rv4 = d->geo.PositionInCell(m);
/*! Comment_N: Same as rv3 but for atom m in cell l. */

		  int* iv1 = d->geo.NbrCell(m,n);
/*! Comment_N: *iv1 = relative position of unit cell of nbr (cindx,aindx) wrt unit cell of atom (l,m). NbrCell is a vecor of length 3 with elements as either 0 or 1. e.g. 
1. If nbr lies in the same unit cell as that of an atom then NbrCell = (0,0,0) 
2. If nbr lies in the unit cell which is just above (in z-direction) the unit cell of that atom then NbrCell = (0,0,1). */

		  isn = d->geo.l_m_2_offset_ham_tot[cindx][aindx];
		  /*! Set nstat to {+1,-1} if neighbor is {upper,lower} triangular */ 
		  if ( isn >= imax || isn < imin ) {
		     nstat[n] = 99;
		  }
		  else if ( isn > is0 )
		     nstat[n] = 1;
		  else
		     nstat[n] = -1;
		  
		  indxmat[n][0] = cindx;
		  indxmat[n][1] = aindx;

		  alat = d->strain.amat[AtomType_l_m][AtomType_cindx_aindx];
/*! Comment_N: Unstrained lattice constant of a atom-nbr combination e.g. GaAs, SiGe, SiSi etc. These values are provided through the input deck. If material is pure Si then alat = Si lattice constant. */

		  /*! Get relative real-space location for this neighbor.  There are two
		     terms inside the brackets.  The first gives the relative location
		     of the neighbor's unit cell.  The second gives the relative location
		     of the neighbor's atom within the unit cell.  The term disp is the
		     atom's displacement due to strain */
		  nnv_n = nnv[n];
		  nnv0_n = nnv0[n];

		  for ( p=0; p < 3; p++, rv1++, rv2++, rv3++, rv4++, iv1++ ) {
		     real disp = *rv1 - *rv2;
		     r1 = (double)(*iv1) + (*rv3 - *rv4);
		     eta[n][p] = r1;
		     nnv_n[p] = a_lattice[p] * r1 + disp;
/* Addition for ZB_110 */
                     real alat_temp = alat;
                     if(d->geo.isCellZB_110() || d->geo.isCellZB_111() || d->geo.isCellUserDef()) {
                        alat_temp *= d->geo.get_a_latt_rot_norm(p);
                     }
/* End Addition for ZB_110 */

/*! Comment_N: R_ij's in Keating energy expression. Note that a_lattice has current values of lattice_x,y,z which are allowed to vary depending on boundary conditions. */
                     nnv0_n[p] = alat_temp*r1; /* Added for ZB_110 and ZB_111 */
//		     nnv0_n[p] = alat*r1;
/*! Comment_N: d_ij's in Keating ekenrgy expression. */
		  }

                  if (nnv0[n][0]==0.0 && nnv0[n][1]==0.0 && nnv0[n][2]==0.0) {
                     printf("WARNING:  nnv0 vanishes!!! on proc %d -- Try to analyse:\n", mpi_n3d_id); 
                     printf("W   (l,m,n) = (%d,%d,%d)\n", l, m, n);
                     printf("W   (l_nbr,m_nbr) = (%d,%d)\n", cindx, aindx);
                     printf("W   alat = (%g %g)\n", alat, d->strain.amat[AtomType_l_m][AtomType_cindx_aindx]);
                     printf("W   (AtomType_l_m,AtomType_cindx_aindx) = (%d,%d)\n",
                            AtomType_l_m, AtomType_cindx_aindx);
                     printf("W   d->nnmap[m][n]=(%d,%d,%d)\n",
                            d->geo.NbrCell(m,n,0),
                            d->geo.NbrCell(m,n,1),
                            d->geo.NbrCell(m,n,2));
                     printf("W CellAtom_2_ijk[aindx]=(%e %e %e)\n", 
                            d->geo.PositionInCell(aindx,0), 
                            d->geo.PositionInCell(aindx,1), 
                            d->geo.PositionInCell(aindx,2));
                     printf("W CellAtom_2_ijk[m]=(%e %e %e)\n", 
                            d->geo.PositionInCell(m,0),
                            d->geo.PositionInCell(m,1), 
                            d->geo.PositionInCell(m,2));
                  }
                  /////////////////////////////////////////////////////////////////////////////
                  if((d->opt.ExecParam.Strain.StrainModel==Strain_struct::VFF_keating_strained)
                     && 1) {
                    dssq =  nnv_n[0]* nnv_n[0] +  nnv_n[1]* nnv_n[1] +  nnv_n[2]* nnv_n[2];
                    d0sq = nnv0_n[0]*nnv0_n[0] + nnv0_n[1]*nnv0_n[1] + nnv0_n[2]*nnv0_n[2];
                    
                    alpha[n]=0;
                    beta[n]=0;
                    if (d0sq==0.0) printf("ERROR:  d0sq vanishes!!!\n"), fflush(stdout);
                    alpha[n] = d->strain.alpha_eff(AtomType_l_m,AtomType_cindx_aindx,dssq,d0sq,1,false);
                    /*printf("alpha_eff[%d]=%g   alpha[%d]=%g\n",
                           n,alpha[n],
                           n,d->strain.alpha[AtomType_l_m][AtomType_cindx_aindx]);*/
                  }
                  else if((d->opt.ExecParam.Strain.StrainModel==Strain_struct::VFF_keating)
                          && 1)
                    alpha[n] = d->strain.alpha[AtomType_l_m][AtomType_cindx_aindx];

		  beta[n]  =  d->strain.beta[AtomType_l_m][AtomType_cindx_aindx];

	       }
	    }
	 }

         /*! try to catch Nan errors -- seems to happen when nnv0[nbr][xyz] = 0;
            this leads to vanishing d0sq and den */
	 for ( n=0; n < d->geo.Neighbors(m); n++ ) {
	    if ( nstat[n] == 0 ) continue;
            
            if (nnv0[n][0]==0.0 && nnv0[n][1]==0.0 && nnv0[n][2]==0.0) {
               printf("Error!  nnv0 vanishes!!! on proc %d -- Try to analyse:\n", mpi_n3d_id); 
               printf("-> nstat[%d] = %d\n", n, nstat[n]);
               printf("   (l,m,n) = (%d,%d,%d)\n", l, m, n);
               printf("   (l_nbr,m_nbr) = (%d,%d)\n", cindx, aindx);
               printf("   d->nnmap[m][n]=(%d,%d,%d)\n",
                      d->geo.NbrCell(m,n,0),
                      d->geo.NbrCell(m,n,1),
                      d->geo.NbrCell(m,n,2));
               printf("d->CellAtom_2_ijk[aindx]=(%e %e %e)\n", 
                      d->geo.PositionInCell(aindx,0), 
                      d->geo.PositionInCell(aindx,1), 
                      d->geo.PositionInCell(aindx,2));
               printf("d->CellAtom_2_ijk[m]=(%e %e %e)\n", 
                      d->geo.PositionInCell(m,0),
                      d->geo.PositionInCell(m,1), 
                      d->geo.PositionInCell(m,2));
               fflush(stdout);
            }
         }

	 /*! Calculate strain energy */
#ifdef INTEL_MODIFIED
#pragma message ("INTEL_MODIFIED: Strain, move invariant (part1) is ON")
        i3 = (l-cell0)*d->geo.AtomsPerCellMax()*3 + m*3;
        i1 = i3;
#else
//! move invariant

#endif

	 for ( n=0; n < d->geo.Neighbors(m); n++ ) {
	    if ( nstat[n] == 0 ) continue;     //!  Check if neighbor exists

	    if (d->strain.debugp) printf("\nStrain loop");
            
            /*! Bond length distortion */
            nnv_n = nnv[n];
            nnv0_n = nnv0[n];
            
            dssq =  nnv_n[0]* nnv_n[0] +  nnv_n[1]* nnv_n[1] +  nnv_n[2]* nnv_n[2];
            d0sq = nnv0_n[0]*nnv0_n[0] + nnv0_n[1]*nnv0_n[1] + nnv0_n[2]*nnv0_n[2];
            
            if (d0sq==0.0) printf("ERROR:  d0sq vanishes!!!\n"), fflush(stdout);
            
            /*! Strain Energy due to bond length distortion */
#define _CORRECTION_BOYKIN
#ifdef _CORRECTION_BOYKIN
            bldc = 0.5*alpha[n]/d0sq;
#else
            bldc = alpha[n]/d0sq;
#endif
            bldl = (dssq-d0sq);
            bld += bldc*bldl*bldl;
/*! Comment_N: bld is bond length distortion contribution to the strain energy. */
            /*! Derivative for minimization */
#ifdef INTEL_MODIFIED
#pragma message ("INTEL_MODIFIED: Strain, move invariant, [n] prefetch (part2) is ON")
//! move invariant
#else
            i1 = (l-cell0)*d->geo.AtomsPerCellMax()*3 + m*3;
#endif
            nnv_n = nnv[n];
            r1 = 4*bldc*bldl;
#ifdef INTEL_MODIFIED
/*! Comment_N: INTEL_MODIFIED is defined. So following code will be executed. */
#pragma message ("INTEL_MODIFIED: Strain loop unroll is ON")
                        real tmp1 = r1*nnv_n[0];
                        real tmp2 = r1*nnv_n[1];
                        d->strain.dp[i3] -= tmp1;
                        d->strain.dp[i3+1] -= tmp2;
                        dE_dlatt[0] += eta[n][0]*tmp1;
                        dE_dlatt[1] += eta[n][1]*tmp2;
                        real tmp3 = r1*nnv_n[2];
                        d->strain.dp[i3+2] -= tmp3;
            dE_dlatt[2] += eta[n][2]*tmp3;
/*! Comment_N: strain.dp = dE/dR_ij, dE_dlatt = dE/da */
#else
/*! Comment_N: INTEL_MODIFIED is defined. So following code will not be executed. */
            for ( p=0; p < 3; p++ )
               d->strain.dp[ i1 + p] -= r1*nnv_n[p]; 
            dE_dlatt[0] += r1*eta[n][0]*nnv_n[0];
            dE_dlatt[1] += r1*eta[n][1]*nnv_n[1];
            dE_dlatt[2] += r1*eta[n][2]*nnv_n[2];
#endif            
            /*! Bond bending distortion */
            for ( p=n+1; p < d->geo.Neighbors(m); p++ ) {
               if ( nstat[p] == 0 ) continue;
               
               /*! d0dot == dot product between vectors pointing to undistorted atoms n and p 
                  dsdot == dot product between vectors pointing to strained atoms n and p
                  den == product of undistorted distances to atoms n and p */
#ifdef INTEL_MODIFIED
#pragma message ("INTEL_MODIFIED: Strain move invariant (part4) is ON")
#else
               nnv_n = nnv[n];
               nnv0_n = nnv0[n];
#endif
               nnv_p = nnv[p];
               nnv0_p = nnv0[p];
               d0dot = nnv0_n[0]*nnv0_p[0] + nnv0_n[1]*nnv0_p[1] + nnv0_n[2]*nnv0_p[2];
               dsdot =  nnv_n[0]* nnv_p[0] +  nnv_n[1]* nnv_p[1] +  nnv_n[2]* nnv_p[2];
               den   = nnv0_p[0]*nnv0_p[0] + nnv0_p[1]*nnv0_p[1] + nnv0_p[2]*nnv0_p[2];
               
               if (den==0.0) printf("ERROR:  den vanishes!!!\n"), fflush(stdout);

               Beta=0.0;
               if((d->opt.ExecParam.Strain.StrainModel==Strain_struct::VFF_keating_strained)
                  && 1){
                 
                 Beta=d->strain.beta_eff((int) d->geo.AtomType[l][m],//! type of the m-th atom
                                         (int) d->geo.AtomType[indxmat[n][0]][indxmat[n][1]],//! type of the n1-th atom
                                         (int) d->geo.AtomType[indxmat[p][0]][indxmat[p][1]],//! type of the n2-th atom
                                         nnv0[n],nnv0[p],//! vectors of the relaxed bonds of the V doublet
                                         nnv[n],nnv[p],//! vectors of the strained bonds of the V doublet
                                         1,//! use anharmonicity corrections for strain and acoustical phonons
                                         false);
                 /*printf("%d->%d : beta_eff=%g   beta=%g\n",
                        n,p,
                        Beta,
                        sqrt(beta[n]*beta[p]));*/
               }
               else if((d->opt.ExecParam.Strain.StrainModel==Strain_struct::VFF_keating)
                       && 1)
#ifdef INTEL_MODIFIED
#pragma message ("INTEL_MODIFIED: Strain arithmetic opt's is ON")
               Beta=sqrt(beta[n]*beta[p]);
               bbdc = Beta/*sqrt(beta[n]*beta[p])*//sqrt(den*d0sq);

                bbdl = (dsdot-d0dot);
                real bb_mul = bbdc*bbdl;
                bbd += bb_mul*bbdl;
/*! Comment_N: bbd is bond bending contribution to the strain energy. */
                bbdc_bbdl_2 = 2*bb_mul;

//Already
//               nnv_p = nnv[p];
//               nnv_n = nnv[n];
               /*! Factor of 2 added in the following 3 lines. Without factor of two these 
                  expressions are not correct. refer eqn 14 in C. Pryor et al. JAP vol 83, 
                  p 2548, Mar 1 1998. Factor of 3/8 is ignored in the code because it will 
                  not affect minimization procedure. */
               dE_dlatt[0] += 2*bb_mul * (eta[n][0]*nnv_p[0] + eta[p][0]*nnv_n[0]);
               dE_dlatt[1] += 2*bb_mul * (eta[n][1]*nnv_p[1] + eta[p][1]*nnv_n[1]);
               dE_dlatt[2] += 2*bb_mul * (eta[n][2]*nnv_p[2] + eta[p][2]*nnv_n[2]);
#else
		Beta=sqrt(beta[n]*beta[p]);
               	bbdc = Beta/*sqrt(beta[n]*beta[p])*//sqrt(den*d0sq);
               	if (d->strain.debugp)
                 printf("\ndot_prod %e %e", dsdot, d0dot);
               	bbdl = (dsdot-d0dot);
				/* bbdl = (dsdot+0.60583*.60583/16); */
               bbd += bbdc*bbdl*bbdl;
				/*! Derivatives for minimization */
               nnv_p = nnv[p];
               nnv_n = nnv[n];
               bbdc_bbdl_2 = 2*bbdc*bbdl;
               //i3 = (l-cell0)*d->geo.AtomsPerCellMax()*3 + m*3;
               
//#define KLUDGE1  0.5 /*! needed for proper minimization */
#define KLUDGE1  1.0 /*! Above line was a bug.
                        For expressions of strain energy, refer eqn 14 in C. Pryor et al. JAP vol 83, p 2548, Mar 1 1998.
                        Factor of 3/8 is ignored in the code because it will not affect minimization procedure. */

               dE_dlatt[0] += KLUDGE1*bbdc_bbdl_2 * (eta[n][0]*nnv_p[0] +
                                                     eta[p][0]*nnv_n[0]);
               dE_dlatt[1] += KLUDGE1*bbdc_bbdl_2 * (eta[n][1]*nnv_p[1] +
                                                     eta[p][1]*nnv_n[1]);
               dE_dlatt[2] += KLUDGE1*bbdc_bbdl_2 * (eta[n][2]*nnv_p[2] +
                                                     eta[p][2]*nnv_n[2]);
#endif

#ifdef INTEL_MODIFIED
#pragma message ("INTEL_MODIFIED: Strain, move invariant (part3) is ON")
#else
               i3 = (l-cell0)*d->geo.AtomsPerCellMax()*3 + m*3;
#endif

               if (mpi_left>=0) {
                  cn_L = indxmat[n][0] - cell0_L;
                  cp_L = indxmat[p][0] - cell0_L;
               }
               else {
                  cn_L = cp_L = -1; /*! ie undefined */
               }
               if (mpi_right>=0) {
                  cn_R = indxmat[n][0] - cell0_R;
                  cp_R = indxmat[p][0] - cell0_R;
               }
               else {
                  cn_R = cp_R = -1; /*! ie undefined */
               }
               cn_this = indxmat[n][0] - cell0;
               cp_this = indxmat[p][0] - cell0;
               
               if (cn_this>=0 && cn_this<d->cell_ln[mpi_n3d_id]) {
                  int offset = 3*(cn_this*d->geo.AtomsPerCellMax() + indxmat[n][1]);
#ifdef INTEL_MODIFIED
/*! Comment_N: INTEL_MODIFIED is defined. So following code will be executed. Note that following code and code in #else is same. The following code is optimized to run faster. Same is true for all #ifdef INTEL_MODIFIED loops below. */
#pragma message ("INTEL_MODIFIED: Strain, small loop (1) vectorization is ON")
                                  __m128d __r, __b, __d, __nnv_p;
// Unaligned nnv_n, d->strain.dp
                                  // r2[0] = nnv_n[0] *bbdc_bbdl_2;
                                  // r2[1] = nnv_n[1] *bbdc_bbdl_2;
                                  __b = _mm_load_sd((double*)&bbdc_bbdl_2);
                                  __b = _mm_unpacklo_pd(__b, __b);
                                  __nnv_p = _mm_load_pd((double*)nnv_p);
                                  __r = _mm_mul_pd(__nnv_p, __b);

                                  __d = _mm_loadu_pd(&d->strain.dp[offset]);
                                  __d = _mm_add_pd(__d, __r);
                                  _mm_storeu_pd(&d->strain.dp[offset], __d);

                                  __d = _mm_loadu_pd((double*)&d->strain.dp[i3]);
                                  __d = _mm_sub_pd(__d, __r);
                                  _mm_storeu_pd(&d->strain.dp[i3], __d);

                                  r1 = nnv_p[2]*bbdc_bbdl_2;
                                  d->strain.dp[offset+2] += r1;
                                  d->strain.dp[i3+2] -= r1;
                  
#else
/*! Comment_N: INTEL_MODIFIED is defined. So following code will not be executed. */
		     for ( q=0; q < 3; q++ ) {
                     r1 = nnv_p[q]*bbdc_bbdl_2;
                     d->strain.dp[ q + offset ] += r1;
                     d->strain.dp[ q + i3 ] -= r1;
                  }
#endif
               }
               else if (cn_L>=0 && cn_L<d->cell_ln[mpi_left]) {
                  int offset = 3*(cn_L*d->geo.AtomsPerCellMax() + indxmat[n][1]);
                  for ( q=0; q < 3; q++ ) {
                     r1 = nnv_p[q]*bbdc_bbdl_2;
                     if (r1 != 0.0) {
                        sleft[Nsleft] = r1;
                        indx_sleft[Nsleft++] = q + offset;
                     }
                     d->strain.dp[ i3 + q] -= r1;
                  }
               }
               else if (cn_R>=0 && cn_R<d->cell_ln[mpi_right]) {
                  int offset = 3*(cn_R*d->geo.AtomsPerCellMax() + indxmat[n][1]);
                  for ( q=0; q < 3; q++ ) {
                     r1 = nnv_p[q]*bbdc_bbdl_2;
                     if (r1 != 0.0) {
                        sright[Nsright] = r1;
                        indx_sright[Nsright++] = q + offset;
                     }
                     d->strain.dp[ i3 + q] -= r1;
                  }
               }
               else {
                  die("invalid range: %d (%d %d) %d\n", mpi_n3d_id, 
                        mpi_left, mpi_right, indxmat[n][0]);
               }
               
               if (cp_this>=0 && cp_this<d->cell_ln[mpi_n3d_id]) {
                  int offset = 3*(cp_this*d->geo.AtomsPerCellMax() + indxmat[p][1]);
#ifdef INTEL_MODIFIED
#pragma message ("INTEL_MODIFIED: Strain, small loop (2) vectorization is ON")
                                  __m128d __r, __b, __d, __nnv_n;
// Unaligned nnv_n, d->strain.dp
                                  // r2[0] = nnv_n[0] *bbdc_bbdl_2;
                                  // r2[1] = nnv_n[1] *bbdc_bbdl_2;
                                  __b = _mm_load_sd((double*)&bbdc_bbdl_2);
                                  __b = _mm_unpacklo_pd(__b, __b);
                                  __nnv_n = _mm_load_pd((double*)nnv_n);
                                  __r = _mm_mul_pd(__nnv_n, __b);

                                  __d = _mm_loadu_pd(&d->strain.dp[offset]);
                                  __d = _mm_add_pd(__d, __r);
                                  _mm_storeu_pd(&d->strain.dp[offset], __d);

                                  __d = _mm_loadu_pd((double*)&d->strain.dp[i3]);
                                  __d = _mm_sub_pd(__d, __r);
                                  _mm_storeu_pd(&d->strain.dp[i3], __d);

                                  real r2_2 = nnv_n[2]*bbdc_bbdl_2;
                                  d->strain.dp[offset+2] += r2_2;
                                  d->strain.dp[i3+2] -= r2_2;
#else
       		     for ( q=0; q < 3; q++ ) {
                     r2 = nnv_n[q]*bbdc_bbdl_2;
                     d->strain.dp[ q + offset] += r2;
                     d->strain.dp[ q + i3] -= r2;
                  }
#endif
               }
               else if (cp_L>=0 && cp_L<d->cell_ln[mpi_left]) {
                  int offset = 3*(cp_L*d->geo.AtomsPerCellMax() + indxmat[p][1]);
                  for ( q=0; q < 3; q++ ) {
                     r2 = nnv_n[q]*bbdc_bbdl_2;
                     if (r2 != 0.0) {
                        sleft[Nsleft] = r2;
                        indx_sleft[Nsleft++] = q + offset;
                     }
                     d->strain.dp[ i3 + q] -= r2;
                  }
               }
               else if (cp_R>=0 && cp_R<d->cell_ln[mpi_right]) {
                  int offset = 3*(cp_R*d->geo.AtomsPerCellMax() + indxmat[p][1]);
                  for ( q=0; q < 3; q++ ) {
                     r2 = nnv_n[q]*bbdc_bbdl_2;
                     if (r2 != 0.0) {
                        sright[Nsright] = r2;
                        indx_sright[Nsright++] = q + offset;
                     }
                     d->strain.dp[ i3 + q] -= r2;
                  }
               }
               else {
                  die("invalid range: %d (%d %d) %d\n", mpi_n3d_id, 
                        mpi_left, mpi_right, indxmat[p][0]);
               }
            }
	 }
      }
   }
   
   Eloc = bld+bbd;
/*! Comment_N: Strain energy of atoms on this processor. */
   
   /*! test E_loc (i.e. make sure it's not nan) */
   if ( !(Eloc <=0) && !(Eloc>0)) {
      printf("ERROR:  Eloc is Nan!!!! (ID=%d) . . .   ", mpi_n3d_id), fflush(stdout);
      if ( !(bld <=0) && !(bld>0)) {
         printf("bld is Nan!!!! "), fflush(stdout);
      }
      if ( !(bbd <=0) && !(bbd>0)) {
         printf("bbd is Nan!!!! "), fflush(stdout);
      }
      printf("\n"); fflush(stdout);
      die("");
   }
   
   Einfo_loc[0] = dE_dlatt[0];
   Einfo_loc[1] = dE_dlatt[1];
   Einfo_loc[2] = dE_dlatt[2];
   Einfo_loc[3] = Eloc;

#if (defined MPI3d && !defined FAKE_MPI)
   MPI_Allreduce(&Einfo_loc[0], &Einfo_tot[0], 4, MPI_DOUBLE, MPI_SUM, 
                 MPI_COMM_WORLD);
   dE_dlatt_tot[0] = Einfo_tot[0];
   dE_dlatt_tot[1] = Einfo_tot[1];
   dE_dlatt_tot[2] = Einfo_tot[2];
   Etot            = Einfo_tot[3];
#else
   dE_dlatt_tot[0] = Einfo_loc[0];
   dE_dlatt_tot[1] = Einfo_loc[1];
   dE_dlatt_tot[2] = Einfo_loc[2];
   Etot            = Einfo_loc[3];
#endif

   if (!mpi_n3d_id) {
      int ii;
      double nrm2=0;
      
      for (ii=0; ii<3*d->cell_ln[mpi_n3d_id]*d->geo.AtomsPerCellMax(); ii++) 
	 nrm2 += d->strain.dp[ii] * d->strain.dp[ii];

      if ( _new_minimize && d->strain.doPeriodMinimization() ) {
/*! Comment_N: last 3 entries of strain.dp on the master processor are derivative dE/da. */
	 d->strain.dp[sz_loc-3] = dE_dlatt_tot[0] * SCALE_FCTR;
	 d->strain.dp[sz_loc-2] = dE_dlatt_tot[1] * SCALE_FCTR;
	 d->strain.dp[sz_loc-1] = dE_dlatt_tot[2] * SCALE_FCTR;

#ifdef _SAME_LATT_CONST
         /*! this enforces same lattice constant in x,y,z */
	 d->strain.dp[sz_loc-3] = (dE_dlatt_tot[0]+dE_dlatt_tot[1]+dE_dlatt_tot[2])/3.0* SCALE_FCTR;
	 d->strain.dp[sz_loc-2] = (dE_dlatt_tot[0]+dE_dlatt_tot[1]+dE_dlatt_tot[2])/3.0* SCALE_FCTR;
	 d->strain.dp[sz_loc-1] = (dE_dlatt_tot[0]+dE_dlatt_tot[1]+dE_dlatt_tot[2])/3.0* SCALE_FCTR;
#endif
      }
      printf("\n  STRAIN:   E=%6.3e  a=(%f,%f,%f)", 
             Etot, a_lattice[0], a_lattice[1], a_lattice[2]);
      printf("\n     dE/da=(%6.3e,%6.3e,%6.3e)     |dE/dr|^2=%6.3e", 
             dE_dlatt_tot[0], dE_dlatt_tot[1], dE_dlatt_tot[2], nrm2);
   }
      
   rm_rmatrix(&eta);
#ifdef INTEL_MODIFIED
   for (int i=0; i<maxNbr; i++)
   _mm_free(nnv[i]);
#else
   rm_rmatrix(&nnv);
#endif
   rm_rmatrix(&nnv0);
   rm_rvectr(&alpha);
   rm_rvectr(&beta);
   rm_ivectr(&nstat);
   rm_imatrix(&indxmat);

   /*! sanity check */
   if (Nsleft >= szCommL) {
      printf("EXCEEDED ARRAY BOUNDS!!!(%d) -- Nsleft=%d\n", mpi_n3d_id, Nsleft), fflush(stdout);
   }
   
   if (Nsright >= szCommR) {
      
      printf("EXCEEDED ARRAY BOUNDS!!!(%d) -- Nsright=%d\n", mpi_n3d_id, Nsright), fflush(stdout);
   }
   
      
   if (mpi_n3d_numprocs > 1) {
      /*! send size info to neighbors */
      if (mpi_n3d_id % 2) {
	 if (mpi_left >= 0) {
	    MPI_Sendrecv(&Nsleft, 1, MPI_INT, mpi_left, 1,
			 &Nrleft, 1, MPI_INT, mpi_left, 0, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
	 }
	 if (mpi_right >= 0) {
	     MPI_Sendrecv(&Nsright, 1, MPI_INT, mpi_right, 0,
			 &Nrright, 1, MPI_INT, mpi_right, 1, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
	 }
      }
      else {
	 if (mpi_right >= 0) {
	    MPI_Sendrecv(&Nsright, 1, MPI_INT, mpi_right, 0,
			 &Nrright, 1, MPI_INT, mpi_right, 1, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
	 }
	 if (mpi_left >= 0) {
	    MPI_Sendrecv(&Nsleft, 1, MPI_INT, mpi_left, 1,
			 &Nrleft, 1, MPI_INT, mpi_left, 0, 
                         MPI_COMM_WORLD, &mpi_n3d_status);
	 }
      }

#ifdef _PRINTSTUFF1
      printf("%d <= [%d] => %d   %d => [%d] <= %d   (%d)\n", 
	     Nsleft, mpi_n3d_id, Nsright, Nrleft, mpi_n3d_id, Nrright, 
	     !mpi_n3d_id ? d->commsize[mpi_n3d_id][mpi_right] : 
             d->commsize[mpi_n3d_id][mpi_left] );
      fflush(stdout);
#endif	

      /*! copy send data */
      d->strain.sndL.indx = (int*) d->strain.sndL.data;
      d->strain.sndL.rdata = (real*) ((int*) d->strain.sndL.data + Nsleft);
      d->strain.sndR.indx = (int*) d->strain.sndR.data;
      d->strain.sndR.rdata = (real*) ((int*) d->strain.sndR.data + Nsright);
      for (j=0; j<Nsleft; j++) {
	 d->strain.sndL.indx[j] = indx_sleft[j];
	 d->strain.sndL.rdata[j] = sleft[j];
      }
      for (j=0; j<Nsright; j++) {
	 d->strain.sndR.indx[j] = indx_sright[j];
	 d->strain.sndR.rdata[j] = sright[j];
      }
      
      /*! set receive data pointers */
      d->strain.rcvL.indx = (int*) d->strain.rcvL.data;
      d->strain.rcvL.rdata = (real*) ((int*) d->strain.rcvL.data + Nrleft);
      d->strain.rcvR.indx = (int*) d->strain.rcvR.data;
      d->strain.rcvR.rdata = (real*) ((int*) d->strain.rcvR.data + Nrright);

      /*! send index data to neighbors */
      if (mpi_n3d_id % 2) {
	 if (mpi_left >= 0) {
	    MPI_Sendrecv(d->strain.sndL.data, Nsleft*szDatStr, MPI_BYTE, 
                         mpi_left, 1, d->strain.rcvL.data, Nrleft*szDatStr, 
                         MPI_BYTE, mpi_left, 0, MPI_COMM_WORLD, &mpi_n3d_status);
	 }
	 if (mpi_right >= 0) {
	    MPI_Sendrecv(d->strain.sndR.data, Nsright*szDatStr, MPI_BYTE, 
                         mpi_right, 0, d->strain.rcvR.data, Nrright*szDatStr,
                         MPI_BYTE, mpi_right, 1, MPI_COMM_WORLD, &mpi_n3d_status);
	 }
      }
      else {
	 if (mpi_right >= 0) {
	    MPI_Sendrecv(d->strain.sndR.data, Nsright*szDatStr, MPI_BYTE, 
                         mpi_right, 0, d->strain.rcvR.data, Nrright*szDatStr,
                         MPI_BYTE, mpi_right, 1, MPI_COMM_WORLD, &mpi_n3d_status);
	 }
	 if (mpi_left >= 0) {
	    MPI_Sendrecv(d->strain.sndL.data, Nsleft*szDatStr, MPI_BYTE, 
                         mpi_left, 1, d->strain.rcvL.data, Nrleft*szDatStr, 
                         MPI_BYTE, mpi_left, 0, MPI_COMM_WORLD, &mpi_n3d_status);
	 }
      }

      for (j=0; j<Nrleft; j++)
	 d->strain.dp[ d->strain.rcvL.indx[j] ] += d->strain.rcvL.rdata[j];
      for (j=0; j<Nrright; j++)
	 d->strain.dp[ d->strain.rcvR.indx[j] ] += d->strain.rcvR.rdata[j];
      
   }
   
#ifdef _PRINT_STUFF2
      for (j=0; j<mpi_n3d_numprocs; j++) {
	 int ii, jj;
	 MPI_Barrier(MPI_COMM_WORLD);
	 if (mpi_n3d_id != j) continue;
	 for ( ii=0; ii < d->cell_ln[mpi_n3d_id]; ii++) {
	    int II = ii + cell0;
	    for ( jj=0; jj < d->geo.AtomsPerCellMax(); jj++ ) { 
	       if ( d->geo.AtomType[II][jj] ) {
		  id = 3*(ii*d->geo.AtomsPerCellMax() + jj);
		  if (fabs(d->strain.dp[id]) > 1e-10) {
		     printf("%d %d    ",II,jj);
		     printf(" %e %e %e\n", d->strain.dp[id], d->strain.dp[id+1], d->strain.dp[id+2]);
		     fflush(stdout);
		  }
	       }
	    }
	 }
	 fflush(stdout);
      }
#endif   
 
      rm_rvectr(&sleft);
      rm_rvectr(&sright);
      rm_ivectr(&indx_sleft);
      rm_ivectr(&indx_sright);
      
      return Etot;
}

/**************************************************************************/

int entryInRcvWithThisAtom(int cell , int atom , ivectr rcvEntries , int noAtomsPerCell )
{
  /*!   Find the entry corresponding to the cellNumber in  rcvEntries */

  int lengthOfRcv = nml_iv_extent(rcvEntries) ; 
  int i , l , a , found = 0  ;
      for (i=0; i<lengthOfRcv; i++) {
	l = rcvEntries[i] / noAtomsPerCell ;
	a = rcvEntries[i] % noAtomsPerCell ;
	if ( ( l == cell ) && ( a == atom )  ) 
	  {
	    found = 1 ;
	    break ;
	  }
      }

      /*! here we have or have not found the entry... analyze! */
      
      if ( 0 == found  )
	{
/*	  printf("\nIn entryInRcvWithThisAtom: The cell %d atom %d NOT FOUND! \n", cell, atom) ;
	  fflush(stdout) ; */
	  i = -1 ;
	}

      return( i ) ;
}

