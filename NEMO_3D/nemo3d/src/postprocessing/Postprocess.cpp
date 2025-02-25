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
$Header: /repo/nemo3d/src/postprocessing/Postprocess.cpp,v 1.1 2004/08/26
21:28:07 swlee Exp $
*****************************************************************************/

#include "Postprocess.h"
#include "AngularMomentum.h"
#include "Coulomb.h"
#include "Decoherence.h"
#include "Dipole.h"
#include "MPI_Timing.h"
#include "Projection.h"
#include "ReadData.h"
#include "SimpleParser.h"
#include "TimeReversal.h"
#include "eigsys3d.h"
#include "eigsys3d_par.h"
#include "ham3d.h"
#include "ham3d_final.h"
#include "i_mat.h"
#include "io_utils.h"
#include "mat3d.h"
#include "math3d.h"
#include "opt3d.h"
#include "qd_struct.h"
#include "run3d_mpi.h"
#include "top_struct.H"

using namespace std;

int nemo3d_postprocessing(int argc, char *argv[]) {

  int _is_mpi_initialized = 0;
#if (defined MPI3d && !defined FAKE_MPI)
  // Initialize the MPI message passing interface for parallel execution.
  MPI_Initialized(&_is_mpi_initialized);
  if (_is_mpi_initialized)
    printf("One process already initialized!\n"), fflush(stdout);
  printf("x"), fflush(stdout);
  MPI_Init(&argc, &argv);
#else  /* not (defined MPI3d && !defined FAKE_MPI) */
  printf("Initializing Fake MPI.\n");
#endif /* (defined MPI3d && !defined FAKE_MPI) */
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_n3d_numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_n3d_id);
  printf("\nProcess [%i]/%i:  activated\n", mpi_n3d_id + 1, mpi_n3d_numprocs);
  fflush(stdout);

#ifdef MPI_TIMING
  mpiTiming.Initialize();
#endif /* MPI_TIMING */

  setup_global_streams();

  // Processing command options =======================================
  char data[123];
  char inputXML[123];
  char evec_file[123];
  int NumElectrons = 0;
  int NumHoles = 0;
  int NumAtoms = 0;
  int NumOrbitals = 0;
  int *eList;
  int *hList;
  bool ComputeCoulomb = false;
  bool ComputeDipole = false;
  bool ComputeDipoleApprox = false;
  bool ComputeProjection = false;
  bool ComputeAngularMomentum = false;
  bool ComputeDecoherence = false;
  bool ComputeTimeReversal = false;
  bool HasEvecFile = false;
  vector<string> options;
  for (int i = 0; i < argc; i++)
    options.push_back(argv[i]);

  size_t iopt = 0;
  while (iopt < options.size()) {
    if (options[iopt] == "--data") {
      strcpy(data, options[++iopt].c_str());
      sprintf(inputXML, "%s.xml", data);
    } else if (options[iopt] == "--electron") {
      NumElectrons = atoi(options[++iopt].c_str());
      if (NumElectrons > 0)
        eList = new int[NumElectrons];
      for (int ne = 0; ne < NumElectrons; ne++) {
        eList[ne] = atoi(options[++iopt].c_str());
      }
    } else if (options[iopt] == "--hole") {
      NumHoles = atoi(options[++iopt].c_str());
      if (NumHoles > 0)
        hList = new int[NumHoles];
      for (int nh = 0; nh < NumHoles; nh++)
        hList[nh] = atoi(options[++iopt].c_str());
    } else if (options[iopt] == "--coulomb") {
      ComputeCoulomb = true;
    } else if (options[iopt] == "--dipole") {
      ComputeDipole = true;
    } else if (options[iopt] == "--dipole_approximation") {
      ComputeDipoleApprox = true;
    } else if (options[iopt] == "--projection") {
      ComputeProjection = true;
    } else if (options[iopt] == "--angularmomentum") {
      ComputeAngularMomentum = true;
    } else if (options[iopt] == "--decoherence") {
      ComputeDecoherence = true;
    } else if (options[iopt] == "--time_reversal") {
      ComputeTimeReversal = true;
    } else if (options[iopt] == "--evec_file") {
      HasEvecFile = true;
      strcpy(evec_file, options[++iopt].c_str());
    }
    iopt++;
  }

  if (ComputeTimeReversal && !HasEvecFile)
    die("Option '--time_reversal' should be used with option '--evec_file "
        "filename'.");

  // Prepare Hamiltonian =============================================
  qd_struct d = new QD_struct;
  d->inputfile = n3d_copy_str(inputXML);
  d->setTopStruct(inputXML, inputXML);
  d->set_sMatList();
  masterPrint("Input file has been processed. \n");
  set_random_seed(d->opt.ExecParam.RandomSeed);
  d->opt.ExecParam.Phys.MagneticFieldOn = false;
  if (ComputeAngularMomentum)
    setup_hamiltonian_for_postprocessing(d);
  else
    setup_crystal_structure_for_postprocessing(d);
  masterPrint("Finished with set-up of data structures.\n\n\n");

  // Read NEMO3D Output Files ==========================================
  complex *wf_e;
  wf_e = 0;
  complex *wf_h;
  wf_h = 0;
  complex *wf;
  wf = 0;
  int *atom_id;
  atom_id = 0;
  int *shape_rank;
  shape_rank = 0;
  int *neighbor;
  neighbor = 0;
  double *lattice;
  lattice = 0;

  read_natoms_norbitals(data, &NumAtoms, &NumOrbitals);
  cout << "NumAtoms = " << NumAtoms << endl;
  cout << "NumOrbitals = " << NumOrbitals << endl << endl;

  if (ComputeTimeReversal)
    read_wavefunction_only(evec_file, NumAtoms, NumOrbitals, &wf);
  else
    read_data_for_postprocessing(data, NumElectrons, eList, NumHoles, hList,
                                 &wf_e, &wf_h, &atom_id, &shape_rank, &neighbor,
                                 &lattice);

  cout << "Finished reading NEMO3D outputfiles" << endl;

  // Compute Time Reversal Operator ========================================
  if (ComputeTimeReversal) {
    TimeReversal timereversal = TimeReversal(d, mpi_n3d_id);
    timereversal.get_data(wf);
    timereversal.apply_time_reversal_operator(evec_file);
  }

  // Compute Angular Momentum ==============================================
  if (ComputeAngularMomentum) {
    AngularMomentum angularmomentum =
        AngularMomentum(d, NumElectrons, NumHoles, mpi_n3d_id);
    angularmomentum.get_data(wf_e, wf_h, atom_id, shape_rank, neighbor,
                             lattice);
    angularmomentum.prepare_MatrixL();
    angularmomentum.set_origin(9.33, 9.33, 3.39);
    angularmomentum.compute_angular_momentum_spin('x');
  }

  // Compute Dipole Momentum ===============================================
  if (ComputeDipole) {
    ifstream file_input("dipole.table");
    if (!file_input.good())
      die("File 'dipole.table' does not exist.");
    Dipole dipole = Dipole(d, NumElectrons, NumHoles, mpi_n3d_id);
    dipole.get_data(wf_e, wf_h, atom_id, shape_rank, neighbor, lattice);
    dipole.assign_table(file_input);
    dipole.compute_dipole();
    dipole.compute_dipole_with_spin_suppressed();
  }

  // Compute Dipole Momentum Approximation =================================
  if (ComputeDipoleApprox) {
    Dipole dipole = Dipole(d, NumElectrons, NumHoles, mpi_n3d_id);
    dipole.get_data(wf_e, wf_h, lattice);
    dipole.compute_dipole_approximation();
  }

  // Compute Coulomb and Exchange
  // ===============================================
  if (ComputeCoulomb) {
    ifstream file_input("coulomb.table");
    if (!file_input.good())
      die("File 'coulomb.table' does not exist.");
    Coulomb coulomb = Coulomb(d, NumElectrons, NumHoles, mpi_n3d_id);
    coulomb.get_data(wf_e, wf_h, atom_id, neighbor, lattice);
    coulomb.assign_table(file_input);
    //    culomb.compute_coulomb_exchange();
    //    culomb.print_phase();
    cout_master << "Exchange energy along y direction" << endl;
    cout_master << "Overlaped Unit Cells, Exchange Energy (eV)" << endl;
    for (int overlap = 1; overlap <= 72; overlap++) {
      coulomb.compute_overlap_vs_dot_separation(overlap, 2);
    }
    for (int overlap = 1; overlap <= 30; overlap++) {
      coulomb.compute_exchange_vs_dot_separation(overlap, 2);
    }
  }

  // Project wavefunctions onto bulk-band Bloch functions ==================
  if (ComputeProjection) {
    Projection projection = Projection(d, NumElectrons, NumHoles, mpi_n3d_id);
    projection.get_data(wf_e, wf_h, atom_id, shape_rank, neighbor, lattice);
    projection.compute_projection();
    projection.compute_extent(0);
    /*
        projection.compute_overlap(0,0);
        projection.compute_overlap(0,1);
        projection.compute_extent(0);
        projection.compute_weight_on_shape(0);

        if(NumElectrons>0 && NumHoles>0)
          projection.compute_electron_hole_alignment(0,0);
        if(NumElectrons >1) {
           projection.compute_overlap(0,1);
           projection.compute_overlap_without_spin(0,1);
        }
    */
  }

  // Compute Decoherence Time ===========================================
  if (ComputeDecoherence) {
    Decoherence decoherence =
        Decoherence(d, NumElectrons, NumHoles, mpi_n3d_id);
    decoherence.get_data(wf_e, wf_h, atom_id, shape_rank, neighbor, lattice);
    decoherence.get_geometry(3.0, 3.0, 3.0, 15.0, 15.0, 6.0, 1.0);
    decoherence.compute_decoherence_time();
    //    decoherence.MagneticField=0.0;
    //    decoherence.RandomOrientation=true;
    //    decoherence.compute_spin_correlator('e', 0);
    //    decoherence.RandomOrientation=false;
    //    decoherence.compute_spin_correlator('e', 0);
  }

  // Deleting allocated memory =========================================
  if (wf_e)
    rm_cvectr(&wf_e);
  if (wf_h)
    rm_cvectr(&wf_h);
  if (wf)
    rm_cvectr(&wf);
  if (lattice)
    delete[] lattice;
  if (atom_id)
    delete[] atom_id;
  if (shape_rank)
    delete[] shape_rank;
  if (neighbor)
    delete[] neighbor;

  // Finalizing MPI ====================================================

#ifdef MPI3d
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif /* MPI3d */

  return 0;
}
