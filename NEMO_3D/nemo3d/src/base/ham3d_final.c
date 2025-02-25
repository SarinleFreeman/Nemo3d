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

*****************************************************************************/

#include <fstream>

#include "eigsys3d_par.h"
#include "ham3d_final.h"

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

// #undef READ_ATYPE
#define READ_ATYPE

void adjust_topology(qd_struct d) {
  get_total_hamiltonian_dimension_and_offsets(d);
  /* if the mpi structures were initialized before,
     we have to redo it here, since the number of active
     atoms in the Hamiltonian changed. */
  mpi_segment_init(d, mpi_n3d_numprocs);

  d->geo.SupercellReInit();
  d->geo.fill_is_DiffNeighbor();
  d->geo.FindSurfaceAtoms();
  d->geo.isInvalidPeriodicity();
}

void adjust_crystal_for_electronic_calculation(qd_struct d) {
  //{ ofstream outi("geo_poststrain"); outi << d->geo; outi.close(); }

  d->geo.setPeriodicityToEl();

  // Figure out first if there are atoms to be trimmed away.
  // If there are atoms which are treated for strain only and need to
  // be trimmed away, that is is computational and memory intensive
  // task that may not have to be conducted.
  vector<Shape_struct>::iterator sh = d->opt.Dev.ShapeList.begin();
  int all_strain_and_electronic = 1;
  for (; sh != d->opt.Dev.ShapeList.end(); sh++) {
    all_strain_and_electronic =
        all_strain_and_electronic && sh->use_for_strain_and_electronic;
  }

  if (all_strain_and_electronic) {
    cout_master << "\nno removal of strain-only atoms needed";
  } else {

    cout_master << "\nRemoving strain-only atoms.";
    bool ReducedBuffer = d->geo.InvalidateStrainOnly();

    /* To make sure that all the CPU's see the same atomtype
       in the case of alloy disorder or other random
       placements in the strucutre (maybe the random variables were distributed
       differently on different CPUs), we broadcast the
       AtomType matrix to all the CPU's */
#ifdef MEMORY_SAVE_NEWSTRUCT
    MPI_Bcast(&d->geo.AtomType[0][0],
              d->geo.N_Cell * d->geo.cellgeom.AtomsPerCellMax(), MPI_CHAR,
              mpi_n3d_masterid, MPI_COMM_WORLD);
#else  /* MEMORY_SAVE_NEWSTRUCT  */
    MPI_Bcast(&d->geo.AtomType[0][0],
              d->geo.N_Cell * d->geo.cellgeom.AtomsPerCellMax(), MPI_INT,
              mpi_n3d_masterid, MPI_COMM_WORLD);
#endif /* MEMORY_SAVE_NEWSTRUCT  */

    adjust_topology(d);
  }

  //{ ofstream outi("geo_postsurfacerelax"); outi << d->geo; outi.close(); }

  print_bond_lengths(d);

  /* If we are not using periodic boundary conditions we can
     eliminate single bonded atoms on the edge of the outer cube,
     which tend to cause edge states, rather than delocalized states
     due to their localized state behavior. */
  int Nremoved;
  if (!d->geo.hasPeriodicity()) { /* The block below should be executed only if
        the structure is not periodic. But the periodicity was not checked!
        This bug was fixed while adding tilt */
    /* remove loose, unbonded atoms from the list */
    if ((Nremoved = (d->geo.remove_few_bonded_atoms(0))) > 0) {
      cout_master << "Removed " << Nremoved
                  << " unbonded atoms from the structure\n";
      adjust_topology(d);
    }
  } /* if hasPeriodicity. */
  //   int removed_atoms = d->geo.remove_user_defined_atoms();
  //   cout<<"the number of atoms removed by user define option
  //   "<<removed_atoms<<endl;
  if (d->opt.ExecParam.ElCalc.ElGeo.RemoveSingleBonds &&
      !d->geo.hasPeriodicity() &&
      (Nremoved = (d->geo.remove_few_bonded_atoms(1))) > 0) {
    cout_master << "Removed " << Nremoved
                << " single bonded atoms from the structure\n";

    adjust_topology(d);
  }

  /* Allow for atoms of one species only to occur on the
     open surface of the simulation domain.  */
  if (d->opt.ExecParam.ElCalc.ElGeo.HomogenizeSurface &&
      !d->geo.hasPeriodicity() &&
      (Nremoved = (d->geo.remove_other_surface_species())) > 0) {

    cout_master << "Removed " << Nremoved
                << " single bonded atoms from the structure\n";

    adjust_topology(d);
  }

  print_struct_gmv(d);
  print_struct_pdb(d);

  /* take care of mapping of atom to (l,m) */
  int l, m, indx = 0;
  int isx = d->cell_s[mpi_n3d_id];
  int iex = isx + d->cell_ln[mpi_n3d_id];

  d->Atom_to_lm = Imatrix(d->cell_ln[mpi_n3d_id] * d->geo.AtomsPerCellMax(), 2);
  for (l = isx; l < iex; l++) {
    for (m = 0; m < d->geo.AtomsPerCellMax(); m++) {
      if (!d->geo.AtomType[l][m])
        continue;

      d->Atom_to_lm[indx][0] = l;
      d->Atom_to_lm[indx][1] = m;
      indx++;
    }
  }
}

/* This initializes the communication map for matrix-vector multiplication */
int vecmap_init(qd_struct d, int proc, ivectr vecmap_proc) {
  /* ssmap should be calculated before vecmap when crystal is tilted */
  if (d->geo.is_tilted)
    d->geo.SupercellReInit();
  /* End ssmap should be calculated */

  int i, j, k, l, m, n, p, isn, target_i;
  int isx, iex;
  int imin, imax, isyc, yccnt, check;
  rvectr nnv;

  nnv = Rvectr(3);

  /* Initialize segment boundary indices for processor proc.
     isx and iex contain first and last atom assigned to processor proc.
     imin and imax contain the first and last row assigned to processor
     proc.  There are d->NBasisStates rows for each atom. */
  isx = d->cell_s[proc];
  iex = isx + d->cell_ln[proc];
  isyc = 0;
  yccnt = 0;
  imin = d->seg_s[proc];
  imax = d->seg_s[proc] + d->seg_ln[proc];

  /* Loop through all unit cells associated with this processor. */
  for (l = isx; l < iex; l++) {
    /* These are real space indexes for the unit cell */
    i = d->geo.cell__ijk[l][0];
    j = d->geo.cell__ijk[l][1];
    k = d->geo.cell__ijk[l][2];

    // Loop through all the atoms in the unit cell within simulation domain
    for (m = 0; m < d->geo.AtomsPerCellMax(); m++) {
      if (!d->geo.AtomType[l][m])
        continue;

      /* Loop through atom's neighbors that are within simulation domain */
      for (n = 0; n < d->geo.Neighbors(m); n++) {
        /* Get virtual neighbor's unit cell position */
        int i_nbr = i + d->geo.NbrCell(m, n, 0);
        int j_nbr = j + d->geo.NbrCell(m, n, 1);
        int k_nbr = k + d->geo.NbrCell(m, n, 2);

        /* Get location indices of the unit cell that we will couple to.
           This is the opposite face, edge or corner of the rectangular
           solid.  The if statements check to see if we are at a
           terminating boundary in that dimension.  If so, we couple to
           the opposite boundary in that dimension. Otherwise we couple
           to our current location plus the neighbor displacement */
        int cindx, aindx;
        if (!d->geo.is_tilted) { /* Connect atoms as follows when lattice is NOT
                                    tilted */
          int in = (d->geo.isPeriodicX() &&
                    (i_nbr > d->geo.cell_xmax || i_nbr < d->geo.cell_xmin))
                       ? i - (i_nbr - i) * (d->geo.cell_xmax - d->geo.cell_xmin)
                       : i_nbr;

          int jn = (d->geo.isPeriodicY() &&
                    (j_nbr > d->geo.cell_ymax || j_nbr < d->geo.cell_ymin))
                       ? j - (j_nbr - j) * (d->geo.cell_ymax - d->geo.cell_ymin)
                       : j_nbr;

          int kn = (d->geo.isPeriodicZ() &&
                    (k_nbr > d->geo.cell_zmax || k_nbr < d->geo.cell_zmin))
                       ? k - (k_nbr - k) * (d->geo.cell_zmax - d->geo.cell_zmin)
                       : k_nbr;

          /* Get the unit cell and atomic index for this neighbor */
          cindx = d->geo.ijk__cell[in][jn][kn];
          aindx = d->geo.NbrCell(m, n, 3);
        } else { /* Use already calculated ssmap when lattice is tilted */
          aindx = d->geo.NbrCell(m, n, 3);
          cindx = d->geo.ijk__cell[i_nbr][j_nbr][k_nbr];
#ifdef ELIMINATE_SSMAP
          if (cindx < 0 && (cindx + d->geo.N_Cell + 1) >= 0) {
            cindx = cindx + d->geo.N_Cell + 1;
          }
          if (d->geo.is_tilted)
            cindx = d->geo.get_cindx_tilted(l, m, n);

#else  /* ELIMINATE_SSMAP */
          if (/*cindx < 0 &&*/ d->geo.ssmap[l][m][n] !=
              -1) { /* When crystal is tilted ssmap should be used for
                       connecting atoms in cell with index = 0 as well */
            cindx = d->geo.ssmap[l][m][n];
          }
#endif /* ELIMINATE_SSMAP */
        } /* ! is_tilted  */

        if (cindx < 0 || !d->geo.AtomType[cindx][aindx])
          continue;

        /* Get the location on the global vector for the proper atom
           in the unit cell */
        isn = d->geo.l_m_2_offset_ham_tot[cindx][aindx];

        /* Check if neighbor is outside the x vector range */
        if (isn >= imax || isn < imin) {

          /* Check to see if this section has been previously visited.
             If so, set the offset isyc to the appropriate place */
          check = 1;
          for (p = 0; p < yccnt; p++)
            if (isn == vecmap_proc[p]) {
              check = 0;
              isyc = p * d->NBasisStates;
            }

          /* If this section has not been previously visited, add it
             to the mapping vector */
          if (check) {
            vecmap_proc[yccnt] = isn;
            isyc = yccnt * d->NBasisStates;
            yccnt++;
          }
        }
      }
    }
  }

  for (target_i = 0; target_i < d->nproc; target_i++) {
    if (proc == target_i)
      continue;

    for (i = 0; i < yccnt; i++) {
      if (vecmap_proc[i] >= d->seg_s[target_i] &&
          vecmap_proc[i] < d->seg_s[target_i] + d->seg_ln[target_i]) {
        d->commsize[proc][target_i] += 1;
      }
    }
  }

  d->nvmap[proc] = yccnt;

  rm_rvectr(&nnv);

  return yccnt;
}

/* This initializes segmentation of structure for MPI */
int mpi_segment_init(qd_struct d, int nproc) {
  int i, j, cnt = 0, source_i, target_i;
#ifdef NEW_VECMAP
  ivectr vecmap_temp = NULL;
#endif

  static short initialized_before = 0;

  cout_master << "Perform MPI segment Inititalization\n";

  d->nproc = nproc;
  if (d->nproc > d->geo.N_Cell)
    d->nproc = d->geo.N_Cell;

  if (initialized_before) {
    rm_ivectr(&d->seg_s);
    rm_ivectr(&d->seg_ln);
    rm_ivectr(&d->cell_s);
    rm_ivectr(&d->cell_ln);
    rm_ivectr(&d->nvmap);

#ifdef NEW_VECMAP
    rm_ivectrvectr(&d->vecmap);
#else
    rm_imatrix(&d->vecmap);
#endif

    rm_imatrix(&d->commsize);
    rm_imatrix(&d->commorder);
  }

  initialized_before = 1;

  d->seg_s = Ivectr(d->nproc);
  d->seg_ln = Ivectr(d->nproc);
  d->cell_s = Ivectr(d->nproc);
  d->cell_ln = Ivectr(d->nproc);
  d->nvmap = Ivectr(d->nproc);

#ifdef NEW_VECMAP
  d->vecmap = Ivectrvectr(d->nproc);
#else
  d->vecmap = Imatrix(d->nproc, d->n_atom_tot);
#endif

  d->commsize = Imatrix(d->nproc, d->nproc);
  d->commorder = Imatrix(d->nproc, 2);

  d->maxseg = 0;

  /* there are d->geo.N_cell unit cells in the structure.
   * The cells are evenly distributed accross all the CPUs available.
   * d->cell_ln[i] contains the length (number) of cells inside a single
   * processor number "i"
   * */
  /* set d->cell_ln */
  do {
    for (i = 0; i < d->nproc; i++) {
      d->cell_ln[i]++;
      cnt++;
      if (cnt == d->geo.N_Cell)
        break;
    }
  } while (cnt < d->geo.N_Cell);

  /*
   * * d->cell_s[i]    starting cell index for CPU #i
   * 			CPU 0 starts with index 0
   * 			CPU 1 starts after cell_ln[0]
   * 			CPU 2 starts after cell_ln[0]+cell_ln[1]
   * 			etc....
   * 			*/
  /* set d->cell_s */
  d->cell_s[0] = 0;
  for (i = 1; i < d->nproc; i++)
    d->cell_s[i] += d->cell_s[i - 1] + d->cell_ln[i - 1];

  /* ******************************************/
  /* loop through all processors */
  /* determine and set d->seg_s and d->seg_ln */
  for (i = 0; i < d->nproc; i++) {
    /* set the first atom index first */
    d->seg_s[i] = -1;
    /* loop through all the cells that are sitting on the particular CPU #i */
    for (int c = d->cell_s[i]; c < d->cell_s[i] + d->cell_ln[i]; c++) {
      /* loop through the atoms in each cell, break out if there is a initiated
       * atom found */
      for (j = 0; j < d->geo.AtomsPerCellMax(); j++) {
        if (d->geo.AtomType[c][j]) {
          d->seg_s[i] = d->geo.l_m_2_offset_ham_tot[c][j];
          break;
        }
      }
      /* break out of the loop if there has been a valid atom found */
      if (d->seg_s[i] != -1)
        break;
    }

    /* determine the length of the vector segment for CPU #i */
    d->seg_ln[i] = -1;
    /* loop through all the cells that live on CPU #i
     * begin from the back end of the included cells and cound back down to the
     * beginning*/
    for (int c = d->cell_s[i] + d->cell_ln[i] - 1; c >= d->cell_s[i]; c--) {
      /* loop through the atoms in each cell starting from the top of the list,
       * check if each individual atom is initialized/included
       * determine the length of the overall vector on CPU #i based the
       * Hamiltonian resulution */
      for (j = d->geo.AtomsPerCellMax() - 1; j >= 0; j--) {
        if (d->geo.AtomType[c][j]) {
          d->seg_ln[i] =
              d->geo.l_m_2_offset_ham_tot[c][j] + d->NBasisStates - d->seg_s[i];
          break;
        }
      }
      if (d->seg_ln[i] != -1)
        break;
    }
    d->maxseg = max(d->maxseg, d->seg_ln[i]);
  }
  /* have determined and set d->seg_s and d->seg_ln */
  /* end of loop through all processors */
  /* ******************************************/

  /*
   * d->geo.N_cell    	total number of unit cells in the overall structure
   *
   * d->cell_ln[i]	number of unit cells inside a single CPU #i
   * 			total number of cells are evenly ditributed across all
   * CPUs d->cell_s[i]      starting cell index for CPU #i CPU 0 starts with
   * index 0 CPU 1 starts after cell_ln[0] CPU 2 starts after
   * cell_ln[0]+cell_ln[1] d->seg_s[i]       starting atom index for CPU #i
   *                   starting index of the first orbital element of the
   * Hamltonian on CPU #i starting index of the first vector element in the
   * Hamiltonian d->seg_ln[i]      length of the vector in the Hamiltonian basis
   *                   the length determination explicitly depende on the
   * existence of uniform unit cells and a constant number of orbitals on each
   * atom d->maxseg         maximum segment length of all seg_ln[i]
   *
   * */
#ifdef NEW_VECMAP
  for (i = 0; i < d->nproc; i++) {
    if (!vecmap_temp)
      vecmap_temp = Ivectr(d->n_atom_tot);
    /* initialize the vecmap in a temporary array */
    vecmap_init(d, i, vecmap_temp);

    /* If there are entries in the array,
       attach the temporary array into the holding structure in the
       appropriate processor slot and swet the pointer to the temporary
       array to null, such that it can be used again. */
    if (d->commsize[mpi_n3d_id][i] > 0 || d->commsize[i][mpi_n3d_id] > 0 ||
        mpi_n3d_id == i) {
      add_elem_ivectrvectr(vecmap_temp, i, d->vecmap);
      vecmap_temp = NULL;
    }
  }
  if (vecmap_temp)
    rm_ivectr(&vecmap_temp);
#else
  for (i = 0; i < d->nproc; i++)
    printf("\n%-5.5i", d->nvmap[i]);
  vecmap_init(d, i, d->vecmap[i]);
#endif

#ifdef PRINT_COMMUNICATOR
  for (i = 0; i < d->nproc; i++) {
    printf("\n");
    for (j = 0; j < d->nproc; j++)
      printf("%-5.5i ", d->commsize[i][j]);
  }
#endif /* PRINT_COMMUNICATOR */

  for (i = 0; i < d->nproc; i++) {
    for (j = 0; j < d->nproc; j++) {
      if (d->commsize[i][j] != 0) {
        if ((abs(i - j) != 1) && (abs(i - j) != d->nproc - 1)) {

          cout_master << "\nToo few atoms per processor..."
                      << "Run on fewer processors.\n\n";
          die("");
        }
      }
    }
  }

  for (source_i = 0; source_i < d->nproc; source_i++)
    for (cnt = 0; cnt < 2; cnt++)
      d->commorder[source_i][cnt] = -1;

  for (cnt = 0; cnt < 2; cnt++) {
    for (source_i = 0; source_i < d->nproc; source_i++) {
      if (cnt == 0) {
        if ((source_i % 2) == 0) {
          target_i = source_i + 1;
          if (target_i >= d->nproc) {
            target_i -= d->nproc;
          }
        } else {
          target_i = source_i - 1;
          if (target_i <= -1) {
            target_i += d->nproc;
          }
        }
      } else {
        if ((source_i % 2) == 0) {
          target_i = source_i - 1;
          if (target_i <= -1) {
            target_i += d->nproc;
          }
        } else {
          target_i = source_i + 1;
          if (target_i >= d->nproc) {
            target_i -= d->nproc;
          }
        }
      }
      if (d->commsize[source_i][target_i]) {
        d->commorder[source_i][cnt] = target_i;
      }
    }
    if (d->nproc == 2)
      break;
  }

#ifdef PRINT_COMMUNICATOR
  for (source_i = 0; source_i < d->nproc; source_i++)
    printf("\nsource=%d targets=%2i,%2i", source_i, d->commorder[source_i][0],
           d->commorder[source_i][1]);
#endif /* PRINT_COMMUNICATOR */
  return (0);
}

/* Plow through all the atoms in the cells, count them and establish
   an overall Hamiltonian map of the offsets.  */
int get_total_hamiltonian_dimension_and_offsets(qd_struct d) {
  d->n_atom_tot = 0;
  for (int l = 0; l < d->geo.N_Cell; l++) {
    int n = 0;
    for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
      if (d->geo.AtomType[l][m]) {
        d->geo.l_m_2_offset_ham_tot[l][m] =
            (d->n_atom_tot + n) * d->NBasisStates;
        n++;
      }
    }
    d->n_atom_tot += n;
  }

  d->n_ham_tot = d->n_atom_tot * d->NBasisStates;

  cout_master << "\nHamiltonian offsets and structure are now set:";
  cout_master << "   * total number of atoms = " << d->n_atom_tot;
  cout_master << "   * Hamiltonian dimension = " << d->n_ham_tot << "\n";

  return 0;
}

/*  This function sets members of d that are related to the choice of basis.
    where U=[0...Nb/2-1],  D=[Nb/2...Nb-1]
    Basis_Spin[U] = 1         Basis_Spin[D] = 0
    Basis_Anion[U] = i        Basis_Anion[D] = i-10
    Basis_Cation[U] = i+10    Basis_Anion[D] = i+10
*/
void basis_init(qd_struct d) {
  if (d->opt.Dev.band_model == Dev_struct::Bands_1_s_nospin) {
    d->NBasisStates = 1;
    d->BandModel = BM_1_s_nospin;
  } else if (d->opt.Dev.band_model == Dev_struct::Bands_10_sp3ss_spin) {
    d->NBasisStates = 10;
    d->BandModel = BM_10_sp3ss_spin;
  } else if (d->opt.Dev.band_model == Dev_struct::Bands_20_sp3d5ss_spin) {
    d->NBasisStates = 20;
    d->BandModel = BM_20_sp3d5ss_spin;
  } else if (d->opt.Dev.band_model == Dev_struct::Bands_10_sp3d5ss_nospin) {
    d->NBasisStates = 10;
    d->BandModel = BM_10_sp3d5ss_nospin;
  } else
    die("Unsupported bandstructure model specified (Bands_1_s_nospin, "
        "Bands_10_sp3ss_spin, Bands_20_sp3d5ss_spin, "
        "Bands_10_sp3d5ss_nospin)\n");

  d->Basis_Spin = Ivectr(d->NBasisStates);
  d->BasisAnion = Ivectr(d->NBasisStates);
  d->BasisCation = Ivectr(d->NBasisStates);
  if (d->BandModel == BM_1_s_nospin) {
    d->Basis_Spin[0] = 1;
    d->BasisAnion[0] = Sa;  /* Sa=1 */
    d->BasisCation[0] = Sc; /* Sc=11  */
  } else if (d->BandModel == BM_20_sp3d5ss_spin ||
             d->BandModel == BM_10_sp3ss_spin) {
    for (int i = 0; i < d->NBasisStates / 2; i++) {
      d->Basis_Spin[i] = 1;
      d->Basis_Spin[i + d->NBasisStates / 2] = 0;
      d->BasisAnion[i] = i;
      d->BasisAnion[i + d->NBasisStates / 2] = i;
      d->BasisCation[i] = i + 10;
      d->BasisCation[i + d->NBasisStates / 2] = i + 10;
    }
  } else if (d->BandModel == BM_10_sp3d5ss_nospin) {
    for (int i = 0; i < d->NBasisStates; i++) {
      d->Basis_Spin[i] = -1;
      d->BasisAnion[i] = i;
      d->BasisCation[i] = i + 10;
    }
  } else
    die("Unsupported bandstructure model specified "
        "(Bands_1,Bands_10_so_s,Bands_20_s,Bands_10_sp3d5s)\n");
}

void setup_crystal_structure(qd_struct d) {
  geo_struct di;
  di.InitializeFromInputFile(d->opt);

  /* Additions for tilt */
  vector<Shape_struct>::iterator sh = d->opt.Dev.ShapeList.begin();
  int Nxx = (int)(sh->x / di.lattice_x);
  int Nyy = (int)(sh->y / di.lattice_y);
  FILE *fp;
  /*   fp = fopen("GeoInfo.in","r");
     if(fp != NULL) {
        fclose(fp);
  //      di.readTiltParams(Nxx,Nyy);  /* this function reads attitional /*
  tilt_fix */
  //                           information needed for defining tilted geometry
  //                           */ /* tilt_fix */
  di.readStepRoughness(Nxx, Nyy); /* this function reads attitional
                       information needed for defining tilted geometry */
                                  /*     }
                                     else {
                                        printf("Could not open file 'GeoInfo.in' for reading.");
                                        printf("\nLattice not tilted");
                                        di.is_tilted = 0;
                                       }
                                  /* End Additions for tilt */

  if (d->opt.Dev.StructFile != "NULL") {
    ham_file_init_spds(d, &di);
  } else {
    cout_master << "\nSet up final tight-binding data structures:";
    basis_init(d);
    param_init(d);

    d->strain.SetMatParm(d->sMatList);
    d->strain.SetOptions(d->opt);

    cout_master << "Processing physical structures:\n";
    vector<Shape_struct>::iterator sh = d->opt.Dev.ShapeList.begin();
    // special case:  do some extra initializations for the first shape
    di.FillFirstShape(*sh, d->opt, d->sMatList);

    d->ho = Coperator(d->NBasisStates);
    d->hd = Coperator(d->NBasisStates);
    d->kxL = d->kyL = d->kzL = 0;

    const int gridbuf = 2;
    int Nx = (int)(sh->x / di.lattice_x) + gridbuf;
    int Ny = (int)(sh->y / di.lattice_y) + gridbuf;
    int Nz = (int)(sh->z / di.lattice_z) + gridbuf;

    d->mat = I3tensor(0, Nx - 1, 0, Ny - 1, 0, Nz - 1);
    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nz; k++) {
          d->mat[i][j][k] = sh->matid3d;
        }
      }
    }

    sh++; // point to the second shape
    for (; sh != d->opt.Dev.ShapeList.end(); sh++) {
      if (sh->has_cell_granularity)
        di.FillShape(*sh, d->sMatList, true);
      else
        di.FillShape(*sh, d->sMatList, false);
    }
  }

  d->geo.Reduce(&di);

  // remove unbonded atoms from the list
  if (!d->geo.is_tilted) // Do not remove these atoms when lattice is tilted
    d->geo.remove_few_bonded_atoms(0);

  /* To make sure that all the CPU's see the same atomtype
     in the case of alloy disorder or other random
     placements in the strucutre (maybe the random variables were distributed
     differently on different CPUs), we broadcast the
     AtomType matrix to all the CPU's */
#ifdef MEMORY_SAVE_NEWSTRUCT
  MPI_Bcast(&d->geo.AtomType[0][0],
            d->geo.N_Cell * d->geo.cellgeom.AtomsPerCellMax(), MPI_CHAR,
            mpi_n3d_masterid, MPI_COMM_WORLD);
#else  /* MEMORY_SAVE_NEWSTRUCT  */
  MPI_Bcast(&d->geo.AtomType[0][0],
            d->geo.N_Cell * d->geo.cellgeom.AtomsPerCellMax(), MPI_INT,
            mpi_n3d_masterid, MPI_COMM_WORLD);
#endif /* MEMORY_SAVE_NEWSTRUCT  */

  d->geo.setPeriodicityToStrain();
  adjust_topology(d);

  di.Deallocate();

  if (d->opt.ExecParam.Output.ShapeInfo)
    writeShapeInfo(d);

#ifdef DISPL_SINGLE_PRECISION
  {
    d->geo.l_m_2_atomdsp =
        F3tensor(0, d->geo.N_Cell - 1, 0, d->geo.AtomsPerCellMax() - 1, 0, 2);
    if (0 == mpi_n3d_id) {
      // nml_memory_report("RESTORED THE STRAIN ARRAY AFTER
      // SETUP_OF_CRYSTAL_STRUCTURE") ;
    }
  }
#endif /* DISPL_SINGLE_PRECISION */
}

void setup_hamiltonian(qd_struct d) {
  MPI_TIC(mpiTiming.setup_0);

  // set up d->mstar[material ID] with some default values
  trivial_mat_init(d);

  // nml_memory_report("setup_hamiltonian: After setup of trivial_mat_init");

  setup_crystal_structure(d);

  // adjust atom types
#ifdef READ_ATYPE

  cout_master << "Reading Atom type file." << endl;
  FILE *fp5 = fopen("Steps.txt", "r");
  int at = 0;
  int count = 0;
  int skip = 0;

  if (fp5) {
    while (skip < 6) {
      fscanf(fp5, "%*[^\n]\n");
      skip++;
    }
    for (int l = 0; l < d->geo.N_Cell; l++) {
      for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
        if (!d->geo.AtomType[l][m])
          continue;
        fscanf(fp5, "%d", &at);

        count++;
        if (count < 5)
          cout << at << endl;
        d->geo.AtomType[l][m] = at;
      }
    }
    cout_master << "ok\n";
    fclose(fp5);
  } else {
    cout_master << "NO Atom Type FILE to read.\n";
  }

#endif

  nml_memory_report(
      "setup_hamiltonian: After setup of setup_crystal_structure");
  // print atomic positions <--> d->opt.ExecParam.ElCalc.ElOut.PsiSqr is on

  // cout << "FAO:  uncommented writePos() in setup_hamiltonian()" << endl;
  // cout << "FAO:  the following file contains the PRE-STRAIN atomic positions"
  // << endl;
  if (d->opt.ExecParam.Output.AtomPosBeforeStrain)
    writePos(d, "atom_pos");

  // cout_master << d->geo << endl;
  // exit(1);

  if (d->opt.ExecParam.Strain.Keating.UpdateLatt) {
    double tmp_x = d->opt.ExecParam.Strain.Keating.UpdateLattConst[0];
    double tmp_y = d->opt.ExecParam.Strain.Keating.UpdateLattConst[1];
    double tmp_z = d->opt.ExecParam.Strain.Keating.UpdateLattConst[2];
    d->geo.updateOnlyLattConst(tmp_x, tmp_y, tmp_z);
    cout_master << "Updating lattice constants before strain "
                   "calculaton/reading displacements\n";
    cout_master << " New lattice Constants: [" << d->geo.lattice_x << ", "
                << d->geo.lattice_y << ", " << d->geo.lattice_z << "]\n";
  }

  MPI_TOC(mpiTiming.setup, mpiTiming.setup_0);

  if (d->opt.ExecParam.CalcStr) {
#undef FILEIO_BARRIER_TEST
#ifdef FILEIO_BARRIER_TEST
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    MPI_TIME_INIT(Tstart);
    MPI_TIC(Tstart);
#ifdef MPIP
    MPI_Pcontrol(2); // discard previous profile data
    MPI_Pcontrol(1); // turn on profiling
#endif
#ifdef MPI_TIMING
    mpiTiming.section = mpiTiming.strain_mpi;
#endif
    cout_master
        << "Begin the First Strain Calculation or reading dispacements\n";

    d->strain.first_or_second = 1;
    min_estrain(d); // the first strain calculation
    cout_master << "End the First Strain Calculation or reading dispacements\n";

#ifdef MPI_TIMING
    mpiTiming.section = mpiTiming.other_mpi;
#endif
#ifdef MPIP
    MPI_Pcontrol(3); // generate verbose report
    MPI_Pcontrol(0); // disable profiling
#endif
    MPI_TOC(mpiTiming.strain, Tstart);
    // print atomic positions <--> d->opt.ExecParam.ElCalc.ElOut.PsiSqr is on
    // print species type <--> d->opt.ExecParam.Output.AtomInfo is on
    //     writePos(d, "atom_pos");

    if (d->opt.ExecParam.Strain.StrainOut.TensionOutput) {
      if (d->opt.ExecParam.ElCalc.ElOut.Visualization_3D) {

        writePos(d, "atom_strain_dx_large");
        writePos(d, "shape_dx_large");
      }
    }

#ifdef FILEIO_BARRIER_TEST
    MPI_Barrier(MPI_COMM_WORLD);
#endif // FILEIO_BARRIER_TEST
  }

  MPI_TIC(mpiTiming.setup_0);
  nml_memory_report("setup_hamiltonian: After strain calculation ");

  adjust_crystal_for_electronic_calculation(d);
  /********************************************************/
  // Dump strain for atoms in electronic domain ONLY to .nd_disp_small file:
  if (d->opt.ExecParam.Strain.StrainOut.StrainElOnly) {
    int isx = d->cell_s[mpi_n3d_id];

    char *filename = n3d_strdup_n(d->inputfile);
    printf("creating malloc_strippedVector");
    fflush(stdout);

#ifdef DISPL_SINGLE_PRECISION
    nml_float *disp = (nml_float *)malloc_strippedVector(
        d, 3 * sizeof(nml_float), &d->geo.l_m_2_atomdsp[isx][0][0]);
#else  /* DISPL_SINGLE_PRECISION */
    real *disp = (real *)malloc_strippedVector(
        d, 3 * sizeof(real), &d->geo.l_m_2_atomdsp[isx][0][0]);
#endif /* DISPL_SINGLE_PRECISION */

    n3d_FileTypeSet(&filename, "nd_disp_small", TRUE);
    masterPrint(
        "Writing atomic displacements to '%s' for the electronic domain.",
        filename);
#ifdef DISPL_SINGLE_PRECISION
    writeSiteInfo(d, "atom", filename, "{float x, float y, float z}", disp);
#else  /* DISPL_SINGLE_PRECISION */
    writeSiteInfo(d, "atom", filename, "{real x, real y, real z}", disp);
#endif /* DISPL_SINGLE_PRECISION */
    str_free(filename);
    nml_free(disp);
  }
  /********************************************************/

  if (d->opt.ExecParam.Strain.StrainOut.TensionOutput) {

    // writePos(d,"atom_pos_strain");
    if (d->opt.ExecParam.ElCalc.ElOut.Visualization_3D) {
      writePos(d, "atom_strain_dx_small");
    }
  }

  nml_memory_report("setup_hamiltonian: adjusting memory for electronic "
                    "structure calculation ");
  // print atomic positions <--> d->opt.ExecParam.ElCalc.ElOut.PsiSqr is on
  MPI_TOC(mpiTiming.setup, mpiTiming.setup_0);

  // cout << "FAO:  the following file contains the POST-STRAIN atomic
  // positions" << endl;
  if (d->opt.ExecParam.Output.AtomPosAfterStrain)
    writePos(d, "atom_pos");

  MPI_TIC(mpiTiming.setup_0);
  // initialize electrostatic potentil if necessary
  if (d->opt.ExecParam.Phys.PotentialFile != "NULL" ||
      d->opt.ExecParam.Phys.ImpNumber > 0 || d->opt.ExecParam.Phys.EfieldON) {
    d->phi = R2tensor(0, d->geo.N_Cell - 1, 0, d->geo.AtomsPerCellMax() - 1);
    for (int l = 0; l < d->geo.N_Cell; l++) {
      for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
        d->phi[l][m] = 0.0;
      }
    }
  }

  if (d->opt.ExecParam.Phys.PotentialFile != "NULL") {
    // KLUDGE:  add in field here
    cout_master << "Reading " << d->opt.ExecParam.Phys.PotentialFile << " ... ";
    FILE *fp = fopen(d->opt.ExecParam.Phys.PotentialFile.c_str(), "r");
    if (fp) {
      for (int l = 0; l < d->geo.N_Cell; l++) {
        for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
          if (!d->geo.AtomType[l][m])
            continue;
          fscanf(fp, "%lf", &(d->phi[l][m]));
        }
      }
      cout_master << "ok\n";
      fclose(fp);
    } else {
      cout_master
          << "NO SUCH FILE -- setting electrostatic potential to zero.\n";
    }
  }

  /* Apply electric field when no impurity is present */
  if (d->opt.ExecParam.Phys.EfieldON && d->opt.ExecParam.Phys.ImpNumber <= 0) {
    double pos[3];
    for (int l = 0; l < d->geo.N_Cell; l++) {
      for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
        if (!d->geo.AtomType[l][m])
          continue;
        d->geo.getPosition_strained(pos, l, m);
        real V = -(d->opt.ExecParam.Phys.Ex *
                   (pos[0] - d->opt.Dev.ShapeList[0].xorigin)) -
                 (d->opt.ExecParam.Phys.Ey *
                  (pos[1] - d->opt.Dev.ShapeList[0].yorigin)) -
                 (d->opt.ExecParam.Phys.Ez *
                  (pos[2] - d->opt.Dev.ShapeList[0].zorigin));
        d->phi[l][m] += V;
      } // m
    } // l
  } /* if */

  if (d->opt.ExecParam.Phys.ImpNumber > 0) {
    const real U0 = d->opt.ExecParam.Phys.U0; // eV
    const real inv_eps = 1.0 / 11.9;
    const real esqr = 1.4399766; // eV-nm
    int numImpurities = d->opt.ExecParam.Phys.ImpNumber;
    int *Zcell = (int *)malloc(sizeof(int) * numImpurities);
    int *Zatom = (int *)malloc(sizeof(int) * numImpurities);
    real dist_min = 999999;
    float *xyz_desired = (float *)malloc(sizeof(float) * numImpurities * 3);
    bool flag = false;
    int imp = 0;
    int k = 0;
    int i, j;

    // read impurity coords.
    for (i = 0; i < numImpurities; i++) {
      for (j = 0; j < 3; j++) {
        Zcell[i] = -1;
        Zatom[i] = -1;
        xyz_desired[i * 3 + j] = d->opt.ExecParam.Phys.ImpsXYZ[k];
        k++;
      }
    }

    for (int i = 0; i < numImpurities; i++) {
      if ((xyz_desired[3 * i] < d->opt.Dev.ShapeList[0].xorigin) ||
          (xyz_desired[3 * i] >
           (d->opt.Dev.ShapeList[0].xorigin + d->opt.Dev.ShapeList[0].x)) ||
          (xyz_desired[3 * i + 1] < d->opt.Dev.ShapeList[0].yorigin) ||
          (xyz_desired[3 * i + 1] >
           (d->opt.Dev.ShapeList[0].yorigin + d->opt.Dev.ShapeList[0].y)) ||
          (xyz_desired[3 * i + 2] < d->opt.Dev.ShapeList[0].zorigin) ||
          (xyz_desired[3 * i + 2] >
           (d->opt.Dev.ShapeList[0].zorigin + d->opt.Dev.ShapeList[0].z))) {
        flag = true;
        imp = i + 1;
        break;
      }
    }

    if (flag)
      cout_master << "Error in specifying coordinates of impurity " << imp
                  << endl;

    real pos[3];
    real *pos0 = (real *)malloc(sizeof(real) * numImpurities * 3);

    // find closest atom to desired location
    for (int i = 0; i < numImpurities; i++) {
      dist_min = 999999;
      for (int l = 0; l < d->geo.N_Cell; l++) {
        for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
          if (!d->geo.AtomType[l][m])
            continue;

          d->geo.getPosition_strained(pos, l, m);
          real R2 =
              (pos[0] - xyz_desired[3 * i]) * (pos[0] - xyz_desired[3 * i]) +
              (pos[1] - xyz_desired[3 * i + 1]) *
                  (pos[1] - xyz_desired[3 * i + 1]) +
              (pos[2] - xyz_desired[3 * i + 2]) *
                  (pos[2] - xyz_desired[3 * i + 2]);

          if (R2 < dist_min) {
            dist_min = R2;
            Zcell[i] = l;
            Zatom[i] = m;
          }
        } // m
      } // l

      d->geo.getPosition_strained(&pos0[3 * i], Zcell[i], Zatom[i]);

      // a more sophisticated Phosphorus impurity model
      if (d->opt.ExecParam.Phys.imp_model == 2)
        d->geo.AtomType[Zcell[i]][Zatom[i]] = EL_P;

      cout_master << "Placing a donor at r = (" << pos0[3 * i] << ","
                  << pos0[3 * i + 1] << "," << pos0[3 * i + 2]
                  << ") with U0=" << U0 << " eV ... " << endl;

    } // end of impurity loop

    real *R = (real *)malloc(sizeof(real) * numImpurities);
    // real *Rxy = (real *) malloc(sizeof(real)*numImpurities);// YuWang
    // real *Rz = (real *) malloc(sizeof(real)*numImpurities);// YuWang
    real V_nuc = 0.0;
    if (d->opt.ExecParam.Phys.include_nuclear_term && (numImpurities > 1)) {
      real dd = (pos0[0] - pos0[3]) * (pos0[0] - pos0[3]) +
                (pos0[1] - pos0[4]) * (pos0[1] - pos0[4]) +
                (pos0[2] - pos0[5]) * (pos0[2] - pos0[5]);
      V_nuc = inv_eps * esqr / sqrt(dd);
      cout << "Repulsive potential of the Nucleus: " << V_nuc << endl;
    }

    FILE *fp1 = fopen("donor_potential", "w");
    FILE *fp2 = fopen("Applied_potential", "w");

    // double Ld = 0.317; // YuWang
    // double tempu;
    for (int l = 0; l < d->geo.N_Cell; l++) {
      for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
        if (!d->geo.AtomType[l][m])
          continue;

        for (int i = 0; i < numImpurities; i++) {

          d->geo.getPosition_strained(pos, l, m);

          R[i] = sqrt((pos[0] - pos0[3 * i]) * (pos[0] - pos0[3 * i]) +
                      (pos[1] - pos0[3 * i + 1]) * (pos[1] - pos0[3 * i + 1]) +
                      (pos[2] - pos0[3 * i + 2]) * (pos[2] - pos0[3 * i + 2]));
          // Rxy[i]=sqrt((pos[1]-pos0[3*i+1])*(pos[1]-pos0[3*i+1]) +
          // (pos[2]-pos0[3*i+2])*(pos[2]-pos0[3*i+2])); // YuWang Rz[i] =
          // sqrt((pos[0]-pos0[3*i])*(pos[0]-pos0[3*i]));
          // tempu=inv_eps*esqr*(exp(-Rxy[i]/Ld)/Rxy[i] + 1.0/Rz[i]);
          // tempu = inv_eps*esqr/fabs(pos[2]-pos0[3*i+2]);
          // tempu = inv_eps * esqr / R[i];
          // d->phi[l][m] -= ( l==Zcell[i] && m==Zatom[i] ? U0 : tempu);
          d->phi[l][m] +=
              (l == Zcell[i] && m == Zatom[i] ? U0 : inv_eps * esqr / R[i]);
          // d->phi[l][m] = 0;
          //  sunnylee - for impurity donor potential used in Si:P 1D
          //  bandstructure it will be commented out for this is only applied to
          //  specific application. Again, This should NOT be uncommented...
#if 0 
					// pos0 is original donor position fix a vertual P donor by displacing multiples of a0*2
					// pos is silicon atom position
					{
					    //cout<<"single impurity is repeated virtually in 2 unit cell span"<<endl;
						for(int k=1;k<=100;k++)
						{
                            R[i]=sqrt((pos[0]-(pos0[3*i]+2*numImpurities*k*0.543095))*(pos[0]-(pos0[3*i]+2*numImpurities*k*0.543095))+
                                      (pos[1]-pos0[3*i+1])*(pos[1]-pos0[3*i+1]) +
                                      (pos[2]-pos0[3*i+2])*(pos[2]-pos0[3*i+2]));
							Rxy[i]=sqrt((pos[0]-pos0[3*i])*(pos[0]-pos0[3*i]) +
                                   (pos[1]-pos0[3*i+1])*(pos[1]-pos0[3*i+1])); // YuWang
                            d->phi[l][m] += inv_eps * esqr / R[i];
							//d->phi[l][m] -= ( l==Zcell[i] && m==Zatom[i] ? U0 : inv_eps * esqr*(exp(-Rxy[i]/Ld)/Rxy[i] + 1.0/(pos[2]-pos0[3*i+2])));
                            R[i]=sqrt((pos[0]-(pos0[3*i]-2*numImpurities*k*0.543095))*(pos[0]-(pos0[3*i]-2*numImpurities*k*0.543095))+
                                      (pos[1]-pos0[3*i+1])*(pos[1]-pos0[3*i+1]) +
                                      (pos[2]-pos0[3*i+2])*(pos[2]-pos0[3*i+2]));
							   
                            Rxy[i]=sqrt((pos[0]-pos0[3*i])*(pos[0]-pos0[3*i]) +
                                       (pos[1]-pos0[3*i+1])*(pos[1]-pos0[3*i+1])); // YuWang
                           d->phi[l][m] += inv_eps * esqr / R[i];
					       //d->phi[l][m] -= ( l==Zcell[i] && m==Zatom[i] ? U0 : inv_eps * esqr*(exp(-Rxy[i]/Ld)/Rxy[i] + 1.0/(pos[2]-pos0[3*i+2])));
						}
					}
#endif
          // end of sunnylee

        } // end of impurity loop

        if (d->opt.ExecParam.Phys.include_nuclear_term)
          d->phi[l][m] -= V_nuc;

        if (d->opt.ExecParam.Phys.EfieldON) {
          real V = -(d->opt.ExecParam.Phys.Ex *
                     (pos[0] - d->opt.Dev.ShapeList[0].xorigin)) -
                   (d->opt.ExecParam.Phys.Ey *
                    (pos[1] - d->opt.Dev.ShapeList[0].yorigin)) -
                   (d->opt.ExecParam.Phys.Ez *
                    (pos[2] - d->opt.Dev.ShapeList[0].zorigin));
          d->phi[l][m] += V;
        }

        if (!mpi_n3d_id) {
          if (d->opt.ExecParam.Phys.donor_file)
            fprintf(fp1, "%f %f %f %e\n", pos[0], pos[1], pos[2], d->phi[l][m]);
          // if(d->opt.ExecParam.Phys.electric_file)
          // fprintf(fp2,"%f %f %f %e\n", pos[0], pos[1], pos[2], V);
        }
      }
    }

    free(R);
    // free(Rxy);
    // free(Rz);
    free(Zatom);
    free(Zcell);
    free(pos0);
    fclose(fp1);
    fclose(fp2);

    if (d->opt.ExecParam.Phys.EfieldON && d->opt.ExecParam.Phys.ImpNumber > 0)
      cout_master << "done - donor potential and applied E-field" << endl;
  }

  MPI_TOC(mpiTiming.setup, mpiTiming.setup_0);

  // commented out by Shaikh for piezoelectric calculations
  /*
  #ifdef COMPUTE_STRAIN_TENSOR
     MPI_TIC(mpiTiming.strain_total_0);
     for (int l=0; l < d->geo.N_Cell; l++ ) {
     for (int m=0; m < d->geo.AtomsPerCellMax(); m++ ) {
        if ( !d->geo.AtomType[l][m] )  continue;

        real exx, eyy, ezz, exy, exz, eyz;
        if (compute_strain_tensor(&exx, &eyy, &ezz, &exy, &exz, &eyz, d, l, m)
  == -1) continue;
     }}
     MPI_TOC(mpiTiming.strain_total,mpiTiming.strain_total_0);
  #endif
  */
  //

  // The following portion is added by Shaikh

  cout << "breakpoint" << endl;
  // #ifdef COMPUTE_STRAIN_TENSOR
  cout << "second check" << endl;

  int isx = d->cell_s[mpi_n3d_id];
  int iex = isx + d->cell_ln[mpi_n3d_id];
  FILE *fp, *fpc;

  printf("mpi_n3d_id=%d\n", mpi_n3d_id);
  printf("mpi_n3d_masterid=%d\n", mpi_n3d_masterid);

  if (mpi_n3d_id == mpi_n3d_masterid) {

    fpc = fopen("strain3d_cation", "w");
    fp = fopen("strain3d", "w");
    real pos[3];

    if (d->geo.isCellZB_110() || d->geo.isCellZB_110_small() ||
        d->geo.isCellZB_111() || d->geo.isCellUserDef()) {
      fprintf(fp, "'compute_strain_tensor' not implemented for unit cells "
                  "different than cubic Zincblende unit cell");
      fprintf(fpc, "'compute_strain_tensor' not implemented for unit cells "
                   "different than cubic Zincblende unit cell");
    } else {
      // MPI_TIC(mpiTiming.strain_total_0);
      // for ( l=isx; l < iex; l++ ) {
      for (int l = 0; l < d->geo.N_Cell; l++) {
        for (int m = 0; m < d->geo.AtomsPerCellMax(); m++) {
          d->geo.getPosition_strained(pos, l, m);
          fprintf(fp, "%e %e %e", pos[0], pos[1], pos[2]);
          if (!d->geo.AtomType[l][m])
            continue;

          real exx, eyy, ezz, exy, exz, eyz;

          if (compute_strain_tensor(&exx, &eyy, &ezz, &exy, &exz, &eyz, d, l,
                                    m) == -1) {
            fprintf(fp, " 0 0 0 0 0 0 \n");
            continue;
          }
          fprintf(fp, " %e %e %e %e %e %e \n", exx, eyy, ezz, exy, exz, eyz);

          // For Cations Only
          if (d->geo.isCation(l, m)) {
            d->geo.getPosition_strained(pos, l, m);
            fprintf(fpc, "%e %e %e", pos[0], pos[1], pos[2]);
            if (!d->geo.AtomType[l][m])
              continue;
            real exx, eyy, ezz, exy, exz, eyz;
            if (compute_strain_tensor(&exx, &eyy, &ezz, &exy, &exz, &eyz, d, l,
                                      m) == -1) {
              fprintf(fpc, " 0 0 0 0 0 0 \n");
              continue;
            }
            fprintf(fpc, " %e %e %e %e %e %e \n", exx, eyy, ezz, exy, exz, eyz);
          }
        }
      }
      fclose(fpc);
      fclose(fp);
    }
  }
  // MPI_TOC(mpiTiming.strain_total,mpiTiming.strain_total_0);
  // #endif

  // Shaikh's addition ends here

  nml_memory_report("setup_hamiltonian: end of function ");
}

void setup_hamiltonian_for_postprocessing(qd_struct d) {
  trivial_mat_init(d);
  setup_crystal_structure(d);
  if (d->opt.ExecParam.Strain.Keating.UpdateLatt) {
    double tmp_x = d->opt.ExecParam.Strain.Keating.UpdateLattConst[0];
    double tmp_y = d->opt.ExecParam.Strain.Keating.UpdateLattConst[1];
    double tmp_z = d->opt.ExecParam.Strain.Keating.UpdateLattConst[2];
    d->geo.updateOnlyLattConst(tmp_x, tmp_y, tmp_z);
    cout_master
        << "Lattice constants are updated before reading displacements\n";
    cout_master << " New lattice Constants: [" << d->geo.lattice_x << ", "
                << d->geo.lattice_y << ", " << d->geo.lattice_z << "]\n";
  }

  bool do_min = (d->strain.MinWrtX || d->strain.MinWrtY || d->strain.MinWrtZ);
  int Npad = (!mpi_n3d_id && do_min ? 3 : 0);
  int n_loc = 3 * d->cell_ln[mpi_n3d_id] * d->geo.AtomsPerCellMax() + Npad;
  rvectr p = Rvectr(n_loc);

  readDisplacement(d, p);
  adjust_crystal_for_electronic_calculation(d);

  d->Ham_use_hermiticity = 1;
  d->_primitiveCell_granularity = 0;
  d->_bandstruct_on_full_domain = 1;
  d->_bandstruct_on_local_domain = 0;
  rmatrix kr = NULL;
  int Nk = 0;
  generate_wave_vectors(d, &kr, &Nk);
}

void setup_crystal_structure_for_postprocessing(qd_struct d) {
  trivial_mat_init(d);
  setup_crystal_structure(d);
  adjust_crystal_for_electronic_calculation(d);
}

void ham_file_init_spds(qd_struct d, geo_struct *di) {
  char *token, ftype[10];

  if (d->opt.Dev.StructFile == "NULL")
    die("No structure filename specified\n");
  // fix this kludge
  char strFile[1000];
  strcpy(strFile, d->opt.Dev.StructFile.c_str());
  strtok(strFile, ".\n");
  while ((token = strtok(NULL, ".\n")))
    strcpy(ftype, token);
  if (strcasecmp(ftype, "HIN") == 0) {
    di->InitializeForHyperChem(d->inputfile);
    di->lattice_x = di->lattice_y = di->lattice_z = 1;
  } else
    die("File format .%s is not currently supported, convert to .hin\n", ftype);

  basis_init(d);
  d->ho = Coperator(d->NBasisStates);
  d->hd = Coperator(d->NBasisStates);
}
