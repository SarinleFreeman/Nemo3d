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
$Header: /repo/nemo3d/src/base/geo_struct.c,v 1.38 2008/07/23 16:47:14 nkharche
Exp $
*****************************************************************************/

#include "geo_struct.h"
#include "io_utils.h"
#include "mat_def.h"
#include "util_stringops.h"

// private helper method used by copy-init constructor and assignment operator
void geo_struct::set_geo_struct(
    int periodicity, int N_Cell, int cell_xmin, int cell_ymin, int cell_zmin,
    int cell_xmax, int cell_ymax, int cell_zmax,
    /* Additional variables to tilt the lattice */
    // ivectr tiltdir,
    int dx, int dy, int dz, int is_tilted,
    //                       int NaddALY,
    int NaddALZ,
    //                       double adjThicknessY,
    double adjThicknessZ, i2tensor delta_step,
    //                       double adjustedLenX,
    /* End: Additional variables to tilt the lattice */
    real lattice_x, real lattice_y, real lattice_z, int unitcelldivide,
    const Cell &cellgeom, i3tensor ijk__cell, i2tensor cell__ijk,
#ifdef MEMORY_SAVE_NEWSTRUCT
    b2tensor AtomType, s2tensor ShapeType, b2tensor isDiffNbr,
    b2tensor surfaceInfo, b2tensor atom_strain_only,
#else  /* MEMORY_SAVE_NEWSTRUCT  */
    i2tensor AtomType, i2tensor ShapeType, i2tensor isDiffNbr,
    i2tensor surfaceInfo, i2tensor atom_strain_only,
#endif /* MEMORY_SAVE_NEWSTRUCT  */
    i2tensor l_m_2_offset_ham_tot,
#ifdef DISPL_SINGLE_PRECISION
    f3tensor l_m_2_atomdsp
#else  /* DISPL_SINGLE_PRECISION */
    r3tensor l_m_2_atomdsp
#endif /* DISPL_SINGLE_PRECISION */
#ifndef ELIMINATE_SSMAP
    ,
    i3matrix ssmap /* The comma is not a bug! */
#endif             /* ELIMINATE_SSMAP */
) {
  this->periodicity = periodicity;
  this->N_Cell = N_Cell;
  this->Nx_Cell = Nx_Cell;
  this->Ny_Cell = Ny_Cell;
  this->Nz_Cell = Nz_Cell;
  this->cell_xmin = cell_xmin;
  this->cell_ymin = cell_ymin;
  this->cell_zmin = cell_zmin;
  this->cell_xmax = cell_xmax;
  this->cell_ymax = cell_ymax;
  this->cell_zmax = cell_zmax;
  /* Additional variables to tilt the lattice */
  //   this->tiltdir = Ivectr(3);
  //   this->tiltdir[0] = tiltdir[0]; this->tiltdir[1] = tiltdir[1];
  //   this->tiltdir[2] = tiltdir[2];
  this->dx = dx;
  this->dy = dy;
  this->dz = dz;
  this->is_tilted = is_tilted;
  //   this->NaddALY = NaddALY;
  this->NaddALZ = NaddALZ;
  //   this->adjThicknessY = adjThicknessY;
  this->adjThicknessZ = adjThicknessZ;
  if (this->is_tilted) {
    int Nxx = d_1h(ijk__cell) - d_1l(ijk__cell) + 1 - 2;
    int Nyy = d_2h(ijk__cell) - d_2l(ijk__cell) + 1 - 2;
    int No_steps2 = double(Nxx) / (double(this->dx) / 2) + 1;
    this->delta_step = I2tensor(0, int(4 * Nyy), 0, No_steps2 - 1);
    for (int i = 0; i < int(4 * Nyy); i++)
      for (int j = 0; j < No_steps2; j++)
        this->delta_step[i][j] = delta_step[i][j];
  } /* if is_tilted */
  //   this->adjustedLenX = adjustedLenX;
  /* End: Additional variables to tilt the lattice */
  this->lattice_x = lattice_x;
  this->lattice_y = lattice_y;
  this->lattice_z = lattice_z;
  this->unitcelldivide = unitcelldivide;
  this->cellgeom = cellgeom;

  const int N_MaxAtomCell = cellgeom.AtomsPerCellMax();

  const int Nx = d_1h(ijk__cell) - d_1l(ijk__cell) + 1;
  const int Ny = d_2h(ijk__cell) - d_2l(ijk__cell) + 1;
  const int Nz = d_3h(ijk__cell) - d_3l(ijk__cell) + 1;

  // nml_memory_report("before geo_struct::set_geo_struct");

  this->ijk__cell = I3tensor(0, Nx - 1, 0, Ny - 1, 0, Nz - 1);
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        this->ijk__cell[i][j][k] = ijk__cell[i][j][k];
      }
    }
  }

  this->cell__ijk = I2tensor(0, N_Cell - 1, 0, 2);
  for (int c = 0; c < N_Cell; c++) {
    this->cell__ijk[c][0] = cell__ijk[c][0];
    this->cell__ijk[c][1] = cell__ijk[c][1];
    this->cell__ijk[c][2] = cell__ijk[c][2];
  }

#ifndef ELIMINATE_SSMAP

#define MEMCHECK
#ifdef MEMCHECK
  double int2_memory = 0.0, ijk_cell_memory = 0.0, dspl_memory = 0.0,
         ssmap_memory = 0.0, tl_memory = 0.0;
  ;
  ijk_cell_memory = Nx * Ny * Nz * 4 / 1.0e6;
  tl_memory += ijk_cell_memory;
  int2_memory =
      N_Cell * N_MaxAtomCell * sizeof(int) * 5 / 1.0e6; // 5 tensors of I2tensor
  tl_memory += int2_memory;
  ssmap_memory += N_Cell * N_MaxAtomCell * cellgeom.N_MaxNeighbor *
                  sizeof(int) / 1.0e6; // 1 tensor of I3tensor
  tl_memory += ssmap_memory;
  dspl_memory = N_Cell * N_MaxAtomCell * 3 * 8 / 1.0e6; // 1 r3tensor
  tl_memory += dspl_memory;
  printf("\n set_geo_struct memory allocation in MB ijk_cell_memory=%g "
         "int2_memory=%g dspl_memory=%g ssmap_mem=%g tl_memory=%g\n",
         ijk_cell_memory, int2_memory, dspl_memory, ssmap_memory, tl_memory),
      fflush(stdout);
#endif

#endif /* ELIMINATE_SSMAP */

#ifdef MEMORY_SAVE_NEWSTRUCT
  this->ShapeType = S2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  // NoNeed_isDiffNbr this->isDiffNbr             =
  // B2tensor(0,N_Cell-1,0,N_MaxAtomCell-1);

  this->AtomType = B2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  this->surfaceInfo = B2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  this->atom_strain_only = B2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);

#else  /* MEMORY_SAVE_NEWSTRUCT  */
  this->ShapeType = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  // NoNeed_isDiffNbr this->isDiffNbr             =
  // I2tensor(0,N_Cell-1,0,N_MaxAtomCell-1);

  this->AtomType = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  this->surfaceInfo = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  this->atom_strain_only = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
#endif /* MEMORY_SAVE_NEWSTRUCT  */

  this->l_m_2_offset_ham_tot = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);

#ifdef DISPL_SINGLE_PRECISION
  this->l_m_2_atomdsp = F3tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1, 0, 2);
#else  /* DISPL_SINGLE_PRECISION */
  this->l_m_2_atomdsp = R3tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1, 0, 2);
#endif /* DISPL_SINGLE_PRECISION */

  for (int c = 0; c < N_Cell; c++) {
    for (int a = 0; a < N_MaxAtomCell; a++) {
      this->ShapeType[c][a] = ShapeType[c][a];
      // NoNeed_isDiffNbr this->isDiffNbr[c][a] = isDiffNbr[c][a];
      this->AtomType[c][a] = AtomType[c][a];
      this->l_m_2_offset_ham_tot[c][a] = l_m_2_offset_ham_tot[c][a];
      this->surfaceInfo[c][a] = surfaceInfo[c][a];
      this->atom_strain_only[c][a] = atom_strain_only[c][a];
      this->l_m_2_atomdsp[c][a][0] = l_m_2_atomdsp[c][a][0];
      this->l_m_2_atomdsp[c][a][1] = l_m_2_atomdsp[c][a][1];
      this->l_m_2_atomdsp[c][a][2] = l_m_2_atomdsp[c][a][2];
    }
  }

#ifndef ELIMINATE_SSMAP
  this->ssmap = I3matrix(N_Cell, N_MaxAtomCell, cellgeom.N_MaxNeighbor);
  for (int c = 0; c < N_Cell; c++) {
    for (int a = 0; a < N_MaxAtomCell; a++) {
      for (int n = 0; n < cellgeom.N_MaxNeighbor; n++) {
        this->ssmap[c][a][n] = ssmap[c][a][n];
      }
    }
  }
#endif /* ELIMINATE_SSMAP */

  // nml_memory_report("after geo_struct::set_geo_struct");
}

/* This function reads step roughness profile */
// void geo_struct::readTiltParams(int Nxx, int Nyy) /* tilt_fix */
void geo_struct::readStepRoughness(int Nxx, int Nyy) /* tilt_fix */
{
  FILE *fp;
  double az;

  printf("\n tilt_dir = [%d %d %d] NaddALZ = %d\n", this->dx, this->dy,
         this->dz, this->NaddALZ);
  printf("\n\nNxx = %d Nyy = %d\n\n", Nxx, Nyy);

  int No_steps2 = int(double(Nxx) / (double(this->dx) / 2)) + 1;
  //   printf("\n\nNo_steps2 = %d\n\n",No_steps2);
  this->delta_step =
      I2tensor(0, 4 * Nyy, 0,
               No_steps2 - 1); /* This tensor has step roughness information */

  fp = fopen("step_roughness.in", "r");
  if (fp == NULL) {
    printf("No 'step_roughness' file in the current directory.");
    printf("\nNO step roughness introduced\n");
    for (int i = 0; i < 4 * Nyy; i++)
      for (int j = 0; j < No_steps2; j++)
        this->delta_step[i][j] = 0;
  } /* if fp */
  else {
    printf("\nStep roughness introduced\n");
    int rmax, cmax;
    fscanf(fp, "%d%d", &rmax, &cmax);
    int **temp_delta_step;
    temp_delta_step = I2tensor(0, rmax - 1, 0, cmax - 1);
    for (int i = 0; i < rmax; i++) {
      for (int j = 0; j < cmax; j++)
        fscanf(fp, "%d", &temp_delta_step[i][j]);
    }
    fclose(fp);

    for (int i = 0; i < 4 * Nyy; i++) {
      for (int j = 0; j < No_steps2; j++)
        delta_step[i][j] = temp_delta_step[i][j];
    }
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < No_steps2; j++)
        delta_step[i][j] = 0;
    }
    for (int i = 4 * Nyy - 2; i < 4 * Nyy; i++) {
      for (int j = 0; j < No_steps2; j++)
        delta_step[i][j] = 0;
    }
  } /* else fp */

  /* Test Print */
  /*   for(int i=0;i<4*Nyy;i++) {
      for(int j=0;j<No_steps2;j++)
          printf("%d ",delta_step[i][j]);
      printf("\n");
     } */

} /* readTiltParams*/

/* This function determines adjThicknessZ and is_tilted */
void geo_struct::calc_tilt_params() {
  int NoTiltedUnitCellsInX = int((this->cell_xmax - this->cell_xmin + 1) / dx);
  //   this->adjThicknessY = this->lattice_y*(this->dy*NoTiltedUnitCellsInX/* +
  //   double(4-this->NaddALY)/4*/);
  this->adjThicknessZ = this->lattice_z * (this->dz * NoTiltedUnitCellsInX +
                                           double(4 - this->NaddALZ) / 4);
  /* subtract these numbers from Ly and Lz entered in the input deck */
  // printf("\n\nadjThicknessZ = %12.10le\n\n",this->adjThicknessZ);

  /* tilt_fix
     if (this->dx==1 && this->dy==0 && this->dz==0 && this->NaddALZ==0)
  this->is_tilted = 1; else this->is_tilted = 1; end tilt_fix */
  this->is_tilted = 1; /* TiltedBox => is_tilted=1 */

} /* calc_tilt_params */

void geo_struct::Allocate(int N_Cell, int Nx, int Ny, int Nz, int veryfirst,
                          const Cell &cellgeom) {
  const int LARGE_INT = 9999999;
  const int N_MaxAtomCell = cellgeom.AtomsPerCellMax();
  const int N_MaxNbr = cellgeom.N_MaxNeighbor;

  this->N_Cell = N_Cell;

  this->cell_xmin = this->cell_ymin = this->cell_zmin = LARGE_INT;
  this->cell_xmax = this->cell_ymax = this->cell_zmax = -LARGE_INT;

  this->cellgeom = cellgeom;

  this->Nx_Cell = Nx;
  this->Ny_Cell = Ny;
  this->Nz_Cell = Nz;

#ifndef ELIMINATE_SSMAP

#define MEMCHECK
#ifdef MEMCHECK
  double int2_memory = 0.0, ijk_cell_memory = 0.0, dspl_memory = 0.0,
         ssmap_memory = 0.0, tl_memory = 0.0;
  ;
  ijk_cell_memory = Nx * Ny * Nz * sizeof(int) / 1.0e6;
  tl_memory += ijk_cell_memory;
  int2_memory =
      N_Cell * N_MaxAtomCell * sizeof(int) * 5 / 1.0e6; // 6 tensors of I2tensor
  int2_memory += N_Cell * 2 * sizeof(int) / 1.0e6;      // 1 tensor of I2tensor
  tl_memory += int2_memory;
  dspl_memory =
      N_Cell * N_MaxAtomCell * 3 * sizeof(double) / 1.0e6; // 1 r3tensor
  tl_memory += dspl_memory;
  ssmap_memory =
      N_Cell * N_MaxAtomCell * N_MaxNbr * sizeof(int) / 1.0e6; // 1 i3matrix
  tl_memory += ssmap_memory;
  printf("\n geometry memory allocation in MB ijk_cell_memory=%g "
         "int2_memory=%g dspl_memory=%g ssmap_memory=%g tl_memory=%g\n",
         ijk_cell_memory, int2_memory, dspl_memory, ssmap_memory, tl_memory),
      fflush(stdout);
  printf("\n memory datasize  sizeof(long)=%d sizeof(double)=%d, "
         "sizeof(float)=%d sizeof(int)=%d, sizeof(short)=%d, sizeof(char)=%d\n",
         (int)sizeof(long), (int)sizeof(double), (int)sizeof(float),
         (int)sizeof(int), (int)sizeof(short), (int)sizeof(char));
  fflush(stdout);

  int inttestvalue = 3;
  char charstore = '\0';
  int inttestresult = 0;

  charstore = (char)inttestvalue;
  inttestresult = charstore;
  printf("memory inttestvalue=%d, inttestresult=%d \n", inttestvalue,
         inttestresult);
  inttestresult = (int)charstore;
  printf("memory inttestvalue=%d, inttestresult=%d \n", inttestvalue,
         inttestresult);

  int test_in128 = 128;
  int test_in64 = 64;

  char c_store128, c_store64, c_store128_64;

  c_store128 = (char)test_in128;
  c_store64 = (char)test_in64;

  c_store128_64 = c_store128 + c_store64;
  c_store128_64 = (char)((int)c_store128 + (int)c_store64);

  printf("c_store128=%c c_store64=%c  testequal1=%d testnonequal1=%d "
         "testequal2=%d\n",
         c_store128, c_store64, (((int)c_store128_64) & 128),
         (((int)c_store128_64) & 64), (((int)c_store64) & 64));
  fflush(stdout);
#endif

#endif /* ELIMINATE_SSMAP */

  // nml_memory_report("before geo_struct Allocate");

  this->cell__ijk = I2tensor(0, N_Cell - 1, 0, 2);

#ifdef MEMORY_SAVE_NEWSTRUCT
  this->AtomType = B2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  this->atom_strain_only = B2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  this->ShapeType = S2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
#else  /* MEMORY_SAVE_NEWSTRUCT  */
  this->AtomType = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  this->atom_strain_only = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  this->ShapeType = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
#endif /* MEMORY_SAVE_NEWSTRUCT  */

  if (!veryfirst) {
    this->ijk__cell = I3tensor(0, Nx - 1, 0, Ny - 1, 0, Nz - 1);

    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nz; k++) {
#ifdef ELIMINATE_SSMAP
          this->ijk__cell[i][j][k] = -1 - (this->N_Cell) - 1;
#else  /* ELIMINATE_SSMAP */
          this->ijk__cell[i][j][k] = -1;
#endif /* ELIMINATE_SSMAP */
        }
      }
    }

    this->l_m_2_offset_ham_tot = I2tensor(
        0, N_Cell - 1, 0,
        N_MaxAtomCell -
            1); // This array is only used in the second geometry assignment

#ifdef MEMORY_SAVE_NEWSTRUCT
    // NoNeed_isDiffNbr this->isDiffNbr             = B2tensor(0,N_Cell-1,
    // 0,N_MaxAtomCell-1);  // This array is only used in the second geometry
    // assignment
    this->surfaceInfo = B2tensor(
        0, N_Cell - 1, 0,
        N_MaxAtomCell -
            1); // This array is only used in the second geometry assignment ---
                // will be assigned values of 4,8,16,32,64,128 to indicate the
                // position of the surface position

#else  /* MEMORY_SAVE_NEWSTRUCT  */

    // NoNeed_isDiffNbr this->isDiffNbr             = I2tensor(0,N_Cell-1,
    // 0,N_MaxAtomCell-1);  // This array is only used in the second geometry
    // assignment
    this->surfaceInfo = I2tensor(
        0, N_Cell - 1, 0,
        N_MaxAtomCell -
            1); // This array is only used in the second geometry assignment ---
                // will be assigned values of 4,8,16,32,64,128 to indicate the
                // position of the surface position
#endif /* MEMORY_SAVE_NEWSTRUCT  */

#ifdef DISPL_SINGLE_PRECISION
    this->l_m_2_atomdsp = F3tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1, 0, 2);
#else  /* DISPL_SINGLE_PRECISION */
    this->l_m_2_atomdsp = R3tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1, 0, 2);
#endif /* DISPL_SINGLE_PRECISION */

    // MEMCHECK_noneed this->ssmap = I3matrix(N_Cell, N_MaxAtomCell, N_MaxNbr);
    // MEMCHECK_noneed for (int c=0; c < N_Cell; c++) {
    // MEMCHECK_noneed for (int a=0; a < cellgeom.AtomsPerCellMax(); a++) {
    // MEMCHECK_noneed for (int n=0; n < cellgeom.N_MaxNeighbor; n++) {
    // MEMCHECK_noneed this->ssmap[c][a][n] = -1;
    // MEMCHECK_noneed }}}
  }
  nml_memory_report("after geo_struct Allocate");
}

void geo_struct::Allocate_NOSTRAIN(int N_Cell, int Nx, int Ny, int Nz,
                                   int veryfirst, const Cell &cellgeom) {
  const int LARGE_INT = 9999999;
  const int N_MaxAtomCell = cellgeom.AtomsPerCellMax();
  const int N_MaxNbr = cellgeom.N_MaxNeighbor;

  this->N_Cell = N_Cell;

  this->cell_xmin = this->cell_ymin = this->cell_zmin = LARGE_INT;
  this->cell_xmax = this->cell_ymax = this->cell_zmax = -LARGE_INT;

  this->cellgeom = cellgeom;

  this->Nx_Cell = Nx;
  this->Ny_Cell = Ny;
  this->Nz_Cell = Nz;

#ifndef ELIMINATE_SSMAP

#undef MEMCHECK
#ifdef MEMCHECK
  double int2_memory = 0.0, ijk_cell_memory = 0.0, dspl_memory = 0.0,
         ssmap_memory = 0.0, tl_memory = 0.0;
  ;
  ijk_cell_memory = Nx * Ny * Nz * sizeof(int) / 1.0e6;
  tl_memory += ijk_cell_memory;
  int2_memory =
      N_Cell * N_MaxAtomCell * sizeof(int) * 5 / 1.0e6; // 6 tensors of I2tensor
  int2_memory += N_Cell * 2 * sizeof(int) / 1.0e6;      // 1 tensor of I2tensor
  tl_memory += int2_memory;
  dspl_memory =
      N_Cell * N_MaxAtomCell * 3 * sizeof(double) / 1.0e6; // 1 r3tensor
  tl_memory += dspl_memory;
  ssmap_memory =
      N_Cell * N_MaxAtomCell * N_MaxNbr * sizeof(int) / 1.0e6; // 1 i3matrix
  tl_memory += ssmap_memory;
  printf("\n geometry memory allocation in MB ijk_cell_memory=%g "
         "int2_memory=%g dspl_memory=%g ssmap_memory=%g tl_memory=%g\n",
         ijk_cell_memory, int2_memory, dspl_memory, ssmap_memory, tl_memory),
      fflush(stdout);
  printf("\n memory datasize  sizeof(long)=%d sizeof(double)=%d, "
         "sizeof(float)=%d sizeof(int)=%d, sizeof(short)=%d, sizeof(char)=%d\n",
         (int)sizeof(long), (int)sizeof(double), (int)sizeof(float),
         (int)sizeof(int), (int)sizeof(short), (int)sizeof(char));
  fflush(stdout);

  int inttestvalue = 3;
  char charstore = '\0';
  int inttestresult = 0;

  charstore = (char)inttestvalue;
  inttestresult = charstore;
  printf("memory inttestvalue=%d, inttestresult=%d \n", inttestvalue,
         inttestresult);
  inttestresult = (int)charstore;
  printf("memory inttestvalue=%d, inttestresult=%d \n", inttestvalue,
         inttestresult);

  int test_in128 = 128;
  int test_in64 = 64;

  char c_store128, c_store64, c_store128_64;

  c_store128 = (char)test_in128;
  c_store64 = (char)test_in64;

  c_store128_64 = c_store128 + c_store64;
  c_store128_64 = (char)((int)c_store128 + (int)c_store64);

  printf("c_store128=%c c_store64=%c  testequal1=%d testnonequal1=%d "
         "testequal2=%d\n",
         c_store128, c_store64, (((int)c_store128_64) & 128),
         (((int)c_store128_64) & 64), (((int)c_store64) & 64));
  fflush(stdout);
#endif

#endif /* ELIMINATE_SSMAP */

  // nml_memory_report("before geo_struct Allocate_NOSTRAIN");
  this->cell__ijk = I2tensor(0, N_Cell - 1, 0, 2);

#ifdef MEMORY_SAVE_NEWSTRUCT

  this->AtomType = B2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);

  this->atom_strain_only = B2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);

#else  /* MEMORY_SAVE_NEWSTRUCT  */
  this->AtomType = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
  this->atom_strain_only = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
#endif /* MEMORY_SAVE_NEWSTRUCT  */

  if (!veryfirst) {
    this->ijk__cell = I3tensor(0, Nx - 1, 0, Ny - 1, 0, Nz - 1);

    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < Nz; k++) {
#ifdef ELIMINATE_SSMAP
          this->ijk__cell[i][j][k] = -1 - (this->N_Cell) - 1;
#else  /* ELIMINATE_SSMAP */
          this->ijk__cell[i][j][k] = -1;
#endif /* ELIMINATE_SSMAP */
        }
      }
    }

    this->l_m_2_offset_ham_tot = I2tensor(
        0, N_Cell - 1, 0,
        N_MaxAtomCell -
            1); // This array is only used in the second geometry assignment

#ifdef MEMORY_SAVE_NEWSTRUCT
    this->ShapeType = S2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
    // NoNeed_isDiffNbr this->isDiffNbr             = B2tensor(0,N_Cell-1,
    // 0,N_MaxAtomCell-1);  // This array is only used in the second geometry
    // assignment
    this->surfaceInfo = B2tensor(
        0, N_Cell - 1, 0,
        N_MaxAtomCell -
            1); // This array is only used in the second geometry assignment ---
                // will be assigned values of 4,8,16,32,64,128 to indicate the
                // position of the surface position

#else  /* MEMORY_SAVE_NEWSTRUCT  */
    this->ShapeType = I2tensor(0, N_Cell - 1, 0, N_MaxAtomCell - 1);
    // NoNeed_isDiffNbr this->isDiffNbr             = I2tensor(0,N_Cell-1,
    // 0,N_MaxAtomCell-1);  // This array is only used in the second geometry
    // assignment
    this->surfaceInfo = I2tensor(
        0, N_Cell - 1, 0,
        N_MaxAtomCell -
            1); // This array is only used in the second geometry assignment ---
                // will be assigned values of 4,8,16,32,64,128 to indicate the
                // position of the surface position
#endif /* MEMORY_SAVE_NEWSTRUCT  */

    /* STRAIN IS NOT INITIALIZED!!!!! */

#ifdef DISPL_SINGLE_PRECISION
    //      this->l_m_2_atomdsp         =
    //      F3tensor(0,N_Cell-1,0,N_MaxAtomCell-1,0,2);
#else  /* DISPL_SINGLE_PRECISION */
    //   this->l_m_2_atomdsp         =
    //   R3tensor(0,N_Cell-1,0,N_MaxAtomCell-1,0,2);
#endif /* DISPL_SINGLE_PRECISION */

    // MEMCHECK_noneed this->ssmap = I3matrix(N_Cell, N_MaxAtomCell, N_MaxNbr);
    // MEMCHECK_noneed for (int c=0; c < N_Cell; c++) {
    // MEMCHECK_noneed for (int a=0; a < cellgeom.AtomsPerCellMax(); a++) {
    // MEMCHECK_noneed for (int n=0; n < cellgeom.N_MaxNeighbor; n++) {
    // MEMCHECK_noneed this->ssmap[c][a][n] = -1;
    // MEMCHECK_noneed }}}
  }
  nml_memory_report("after geo_struct Allocate_NOSTRAIN");
}

void geo_struct::InitializeFromInputFile(const top_struct &opt) {
  // set periodicity
  this->periodicity = 0;
  switch (opt.ExecParam.ElCalc.ElGeo.Periodicity) {
  case ElGeo_struct::Closed:
    break;
  case ElGeo_struct::Periodic_x: // set 7st bit
    this->periodicity += 64;
    break;
  case ElGeo_struct::Periodic_y: // set 8nd bit
    this->periodicity += 128;
    break;
  case ElGeo_struct::Periodic_z: // set 9rd bit
    this->periodicity += 256;
    break;
  case ElGeo_struct::Periodic_xy: // set 7st, 8nd bits
    this->periodicity += 192;
    break;
  case ElGeo_struct::Periodic_yz: // set 8nd, 9rd bits
    this->periodicity += 384;
    break;
  case ElGeo_struct::Periodic_xz: // set 7st, 94d bits
    this->periodicity += 320;
    break;
  case ElGeo_struct::Periodic_xyz: // set bits 7,8,9
    this->periodicity += 448;
    break;
  default:
    break;
  }

  switch (opt.ExecParam.Strain.Keating.Periodicity_strain) {
  case Keating_struct::SameAsGlobal:
    this->periodicity += (this->periodicity) >> 3;
    break;
  case Keating_struct::Closed:
    break;
  case Keating_struct::Periodic_x: // set bit 4
    this->periodicity += 8;
    break;
  case Keating_struct::Periodic_y: // set bit 5
    this->periodicity += 16;
    break;
  case Keating_struct::Periodic_z: // set bit 6
    this->periodicity += 32;
    break;
  case Keating_struct::Periodic_xy: // set bits 4,5
    this->periodicity += 24;
    break;
  case Keating_struct::Periodic_yz: // set bits 5,6
    this->periodicity += 48;
    break;
  case Keating_struct::Periodic_xz: // set bits 4,6
    this->periodicity += 40;
    break;
  case Keating_struct::Periodic_xyz: // set bits 4,5,6
    this->periodicity += 56;
    break;
  default:
    break;
  }

  this->lattice_x = opt.Dev.a_lattice_x;
  this->lattice_y = opt.Dev.a_lattice_y;
  this->lattice_z = opt.Dev.a_lattice_z;
  /* tilt_fix */
  this->dx = opt.Dev.tilt_h;
  this->dy = opt.Dev.tilt_k;
  this->dz = opt.Dev.tilt_l;
  //   this->NaddALY = opt.Dev.NaddALY; /* Not active at this stage */
  this->NaddALZ = opt.Dev.NaddALZ;
  printf("\nTilt_dir = [%d %d %d]\n", this->dx, this->dy, this->dz);
  printf("NaddALZ = %d\n", this->NaddALZ);
  /* end tilt_fix */
  if (opt.Dev.CrystalStruct == Dev_struct::Diamond_FCC ||
      opt.Dev.CrystalStruct == Dev_struct::Zincblende_FCC) {
    this->unitcelldivide = 2;
    this->lattice_x /= 2.0;
    this->lattice_y /= 2.0;
    this->lattice_z /= 2.0;
  } else {
    this->unitcelldivide = 1;
  }
}

void geo_struct::Deallocate() {
  // nml_memory_report("geo_struct::Dealloc   before the action");
  if (this->ijk__cell) {
    rm_i3tensor(&this->ijk__cell, 0);
    this->ijk__cell = NULL;
  }
  if (this->cell__ijk) {
    rm_i2tensor(&this->cell__ijk, 0);
    this->cell__ijk = NULL;
  }

#ifdef MEMORY_SAVE_NEWSTRUCT
  if (this->AtomType) {
    rm_b2tensor(&this->AtomType, 0);
    this->AtomType = NULL;
  }
  if (this->ShapeType) {
    rm_s2tensor(&this->ShapeType, 0);
    this->ShapeType = NULL;
  }
  // NoNeed_isDiffNbr if (this->isDiffNbr){
  // NoNeed_isDiffNbr    rm_b2tensor(&this->isDiffNbr,0);
  // NoNeed_isDiffNbr    this->isDiffNbr=NULL;
  // NoNeed_isDiffNbr }
  if (this->surfaceInfo) {
    rm_b2tensor(&this->surfaceInfo, 0);
    this->surfaceInfo = NULL;
  }
  if (this->atom_strain_only) {
    rm_b2tensor(&this->atom_strain_only, 0);
    this->atom_strain_only = NULL;
  }
#else  /* MEMORY_SAVE_NEWSTRUCT  */
  if (this->AtomType) {
    rm_i2tensor(&this->AtomType, 0);
    this->AtomType = NULL;
  }
  if (this->ShapeType) {
    rm_i2tensor(&this->ShapeType, 0);
    this->ShapeType = NULL;
  }
  // NoNeed_isDiffNbr if (this->isDiffNbr){
  // NoNeed_isDiffNbr    rm_i2tensor(&this->isDiffNbr,0);
  // NoNeed_isDiffNbr    this->isDiffNbr=NULL;
  // NoNeed_isDiffNbr }
  if (this->surfaceInfo) {
    rm_i2tensor(&this->surfaceInfo, 0);
    this->surfaceInfo = NULL;
  }
  if (this->atom_strain_only) {
    rm_i2tensor(&this->atom_strain_only, 0);
    this->atom_strain_only = NULL;
  }
#endif /* MEMORY_SAVE_NEWSTRUCT  */
  if (this->l_m_2_offset_ham_tot) {
    rm_i2tensor(&this->l_m_2_offset_ham_tot, 0);
    this->l_m_2_offset_ham_tot = NULL;
  }
  if (this->l_m_2_atomdsp) {
#ifdef DISPL_SINGLE_PRECISION
    rm_f3tensor(&this->l_m_2_atomdsp, 0);
#else  /* DISPL_SINGLE_PRECISION */
    rm_r3tensor(&this->l_m_2_atomdsp, 0);
#endif /* DISPL_SINGLE_PRECISION */

    this->l_m_2_atomdsp = NULL;
  }

#ifndef ELIMINATE_SSMAP
  if (this->ssmap) {
    rm_i3matrix(&this->ssmap);
    this->ssmap = NULL;
  }
#endif /* ELIMINATE_SSMAP */
       // nml_memory_report("geo_struct::Dealloc  after the action");
}

void geo_struct::zero() {
  this->N_Cell = 0;

  this->Deallocate();
  this->cellgeom.Deallocate();
}

bool geo_struct::InvalidateStrainOnly() {

  geo_struct g = *this; // copy this to a temporary object
  this->zero();

  ivectr indx_map = Ivectr(g.N_Cell);

  // initial values -- set

  int xMax = d_1l(g.ijk__cell);
  int yMax = d_2l(g.ijk__cell);
  int zMax = d_3l(g.ijk__cell);
  int xMin = d_1h(g.ijk__cell);
  int yMin = d_2h(g.ijk__cell);
  int zMin = d_3h(g.ijk__cell);

  /* set indx_map, a list of cell indices containing at least one
     valid atom */
  int Zcell = 0;

  for (int l = 0; l < g.N_Cell; l++) {
    for (int m = 0; m < g.cellgeom.AtomsPerCellMax(); m++)
      if ((g.AtomType[l][m] != 0) && (!g.atom_strain_only[l][m])) {
        // add this cell to list
        indx_map[Zcell] = l;
        Zcell++;
        if (g.cell__ijk[l][0] > xMax)
          xMax = g.cell__ijk[l][0];
        if (g.cell__ijk[l][1] > yMax)
          yMax = g.cell__ijk[l][1];
        if (g.cell__ijk[l][2] > zMax)
          zMax = g.cell__ijk[l][2];
        if (g.cell__ijk[l][0] < xMin)
          xMin = g.cell__ijk[l][0];
        if (g.cell__ijk[l][1] < yMin)
          yMin = g.cell__ijk[l][1];
        if (g.cell__ijk[l][2] < zMin)
          zMin = g.cell__ijk[l][2];
        break;
      }
  }
  bool ReducedBuffer = false;
  if (Zcell != g.N_Cell) {
    ReducedBuffer = true;
    cout_master << "The number of removed cells \n"
                << "from strain to electronic structure calculations is "
                << g.N_Cell - Zcell << " \n";
  }

  // actual N is Max-Min + 1, but we pad with zeros on both sides.
  int Nx = xMax - xMin + 3;
  int Ny = yMax - yMin + 3;
  int Nz = zMax - zMin + 3;

  /* Also, initialize ijk__cell to -1 */

  this->Allocate(Zcell, Nx, Ny, Nz, 0, g.cellgeom);

  // this->cell_xmin = xMin;
  // this->cell_ymin = yMin;
  // this->cell_zmin = zMin;
  // this->cell_xmax = xMax;
  // this->cell_ymax = yMax;
  // this->cell_zmax = zMax;

  // NOTE:  structure is shifted back to origin here.
  this->cell_xmin = 1;
  this->cell_xmax = d_1h(this->ijk__cell) - 1;
  this->cell_ymin = 1;
  this->cell_ymax = d_2h(this->ijk__cell) - 1;
  this->cell_zmin = 1;
  this->cell_zmax = d_3h(this->ijk__cell) - 1;

  // Copy scratch-pad into final hamiltonian data-structures

  for (int l = 0; l < this->N_Cell; l++) {
    int *Zr = this->cell__ijk[l];

    Zr[0] = g.cell__ijk[indx_map[l]][0] - xMin + 1;
    Zr[1] = g.cell__ijk[indx_map[l]][1] - yMin + 1;
    Zr[2] = g.cell__ijk[indx_map[l]][2] - zMin + 1;
    this->ijk__cell[Zr[0]][Zr[1]][Zr[2]] = l;

    for (int m = 0; m < g.cellgeom.AtomsPerCellMax(); m++) {
      this->ShapeType[l][m] = g.ShapeType[indx_map[l]][m];
      this->atom_strain_only[l][m] = g.atom_strain_only[indx_map[l]][m];
      // Exclude strain_only_atoms from electronic structure calculation
      if (this->atom_strain_only[l][m])
        this->AtomType[l][m] = 0;
      else
        this->AtomType[l][m] = g.AtomType[indx_map[l]][m];
      this->l_m_2_atomdsp[l][m][0] = g.l_m_2_atomdsp[indx_map[l]][m][0];
      this->l_m_2_atomdsp[l][m][1] = g.l_m_2_atomdsp[indx_map[l]][m][1];
      this->l_m_2_atomdsp[l][m][2] = g.l_m_2_atomdsp[indx_map[l]][m][2];
    }

    /* These are set later:
       i2tensor l_m_2_offset_ham_tot;
       i2tensor isDiffNbr;
       i2tensor surfaceInfo;
    */
  }

  rm_ivectr(&indx_map);
  return ReducedBuffer;
}

/* Allocates space for geo_struct and copies ijk__cell, AtomType,
   atom_strain_only and cell_{xyz}{min,max} from geo_struct* di */
void geo_struct::Reduce(const geo_struct *const di) {
  ivectr indx_map = Ivectr(di->N_Cell);

  /* set indx_map, a list of cell indices containing at least one
     valid atom */
  int Zcell = 0;
  for (int l = 0; l < di->N_Cell; l++)
    for (int m = 0; m < di->AtomsPerCellMax(); m++)
      if (di->AtomType[l][m] != 0) {
        indx_map[Zcell] = l;
        Zcell++;
        break;
      }

#define MEMCHECKmods
#ifdef MEMCHECKmods
  int Nx = di->Nx_Cell;
  int Ny = di->Ny_Cell;
  int Nz = di->Nz_Cell;
#else
  int Nx = d_1h(di->ijk__cell) - d_1l(di->ijk__cell) + 1;
  int Ny = d_2h(di->ijk__cell) - d_2l(di->ijk__cell) + 1;
  int Nz = d_3h(di->ijk__cell) - d_3l(di->ijk__cell) + 1;
#endif

  /* Also, initialize ijk__cell to -1 */
#ifdef DISPL_SINGLE_PRECISION
  this->Allocate_NOSTRAIN(Zcell, Nx, Ny, Nz, 0, di->cellgeom);
#else  /* DISPL_SINGLE_PRECISION */
  this->Allocate(Zcell, Nx, Ny, Nz, 0, di->cellgeom);
#endif /* DISPL_SINGLE_PRECISION */

  this->periodicity = di->periodicity;
  this->cell_xmin = di->cell_xmin;
  this->cell_ymin = di->cell_ymin;
  this->cell_zmin = di->cell_zmin;
  this->cell_xmax = di->cell_xmax;
  this->cell_ymax = di->cell_ymax;
  this->cell_zmax = di->cell_zmax;

  /* Additional variables for tilting the lattice */
  this->dx = di->dx;
  this->dy = di->dy;
  this->dz = di->dz;
  this->is_tilted = di->is_tilted;
  //   this->NaddALY = di->NaddALY;
  this->NaddALZ = di->NaddALZ;
  //   this->adjThicknessY = di->adjThicknessY;
  this->adjThicknessZ = di->adjThicknessZ;

  if (this->is_tilted) {
    int No_steps2 = (double(di->Nx_Cell) - 2) / (double(this->dx) / 2) + 1;
    this->delta_step = I2tensor(0, 4 * (di->Ny_Cell - 2), 0, No_steps2 - 1);
    for (int i = 0; i < int(4 * (di->Ny_Cell - 2)); i++)
      for (int j = 0; j < No_steps2; j++)
        this->delta_step[i][j] = di->delta_step[i][j];
  } /* if is_tilted */
  //   this->adjustedLenX = di->adjustedLenX;
  /* End: Additional variables for tilting the lattice */

  this->lattice_x = di->lattice_x;
  this->lattice_y = di->lattice_y;
  this->lattice_z = di->lattice_z;
  this->unitcelldivide = di->unitcelldivide;

  /* Copy scratch-pad into final hamiltonian data-structures */
  for (int l = 0; l < this->N_Cell; l++) {
    int *Zr = this->cell__ijk[l];

    Zr[0] = di->cell__ijk[indx_map[l]][0];
    Zr[1] = di->cell__ijk[indx_map[l]][1];
    Zr[2] = di->cell__ijk[indx_map[l]][2];
    this->ijk__cell[Zr[0]][Zr[1]][Zr[2]] = l;

    for (int m = 0; m < di->AtomsPerCellMax(); m++) {
      this->AtomType[l][m] = di->AtomType[indx_map[l]][m];

      this->ShapeType[l][m] = di->ShapeType[indx_map[l]][m];

      /* MEMCHECKclip  atom_strain_only has not been assigned in earlier
       * initializations! */
      this->atom_strain_only[l][m] = di->atom_strain_only[indx_map[l]][m];
    }
  }

  rm_ivectr(&indx_map);
}

/* Addition for tilt */
void geo_struct::isNbrOut(int l, int m, int n, bool *_isNotInXnegtv,
                          bool *_isNotInXpostv, bool *_isNotInYnegtv,
                          bool *_isNotInYpostv, bool *_isNotInZnegtv,
                          bool *_isNotInZpostv) {
  /* These are real space indexes for the unit cell */
  int i = this->cell__ijk[l][0];
  int j = this->cell__ijk[l][1];
  int k = this->cell__ijk[l][2];

  // Get relative position of unit cell of virtual neighbor
  int i_nbr = i + this->NbrCell(m, n, 0);
  int j_nbr = j + this->NbrCell(m, n, 1);
  int k_nbr = k + this->NbrCell(m, n, 2);

  // Get the unit cell and atomic index for this neighbor
  int cindx = this->ijk__cell[i_nbr][j_nbr][k_nbr];
  int aindx = this->NbrCell(m, n, 3);

  /* get position of this neighbor */
  double posnbr[3];
  double *Ratom = this->PositionInCell(
      aindx); /* position of this neighbor in the unit cell */
  posnbr[0] = this->lattice_x * (i_nbr + Ratom[0]);
  posnbr[1] = this->lattice_y * (j_nbr + Ratom[1]);
  posnbr[2] = this->lattice_z * (k_nbr + Ratom[2]);

  double delta = 1e-6;

  /* calculate minimum and maximum x, y coordinates of atoms */
  /* Comments: Crystal starts from (0 0 0) before strain calculation.
               After removing strain only atoms crystal is shifted
               back to the origin. So, following boundary definitions
               will be always valid */

  double xminn, xmaxx;
  /* Neerav's addition for tilting ZB_111 */
  if (this->isCellZB_111()) {
    xminn =
        this->lattice_x * (cell_xmin + this->PositionInCell(/*aindx=*/3, 0));
    xmaxx =
        this->lattice_x * (cell_xmax + this->PositionInCell(/*aindx=*/0, 0));
  } else {
    xminn =
        this->lattice_x * (cell_xmin + this->PositionInCell(/*aindx=*/6, 0));
    xmaxx =
        this->lattice_x * (cell_xmax + this->PositionInCell(/*aindx=*/0, 0));
  }

  /* End Neerav's addition for tilting ZB_111 */

  double yminn = ((posnbr[0] / this->lattice_x) * this->lattice_y) *
                     (double(this->dy) / double(this->dx)) +
                 delta;
  /* exclude all atoms on lower boundary */
  double ymaxx;
  ymaxx = yminn + this->lattice_y * (this->cell_ymax - this->cell_ymin + 1 +
                                     this->dy); // - this->adjThicknessY;
  /* include all atoms on upper boundary */

  double zminn = ((posnbr[0] / this->lattice_x) * this->lattice_z) *
                     (double(this->dz) / double(this->dx)) +
                 delta;
  /* exclude all atoms on lower boundary */

  //   if(shape_rank==2) { // If this line is uncommented roughness is
  //   introduced in Si QW only
  // adjust step edges for roughness
  /*      double stepLen =
  (this->lattice_x/4)/(double(this->dz)/double(this->dx)); int
  delta_indx=int(posnbr[1]/(this->lattice_y/4)); // row index in delta_step
        double Delta;

        int stepNo = int((posnbr[0]-delta)/stepLen);
        int stepNo2 = int((posnbr[0]-delta)/(2*stepLen));

  //if (delta_indx>4*(this->Ny_Cell-2)) delta_indx--;
  if (delta_indx==0) delta_indx++;

  //printf("\n\nstepNo = %d stepNo2 = %d delta_indx = %d this->Ny_Cell =
  %d\n\n",stepNo,stepNo2,delta_indx,this->Ny_Cell);

        Delta = delta_step[delta_indx-1][stepNo2]*this->lattice_x/4;

  //printf("\ndelta_indx = %d stepNo = %d stepNo2 = %d Delta = %12.10le stepLen
  = %12.10le",delta_indx,stepNo,stepNo2,Delta,stepLen);

        if(posnbr[0]>(double(stepNo2))*2*stepLen &&
  posnbr[0]<(double(stepNo2)+1)*2*stepLen) { if(Delta>0 &&
  stepNo==(2*stepNo2+1)) { if(posnbr[0]>(double(stepNo))*stepLen &&
  posnbr[0]<(double(stepNo))*stepLen+Delta) zminn = (stepNo)*this->lattice_z/4 -
  delta;
             }
           if(Delta<0 && stepNo==2*stepNo2) {
              if(posnbr[0]>(stepNo+1)*stepLen+Delta &&
  posnbr[0]<(stepNo+1)*stepLen) zminn = (stepNo+2)*this->lattice_z/4 - delta;
             }
          }
  //     } /* if shape_rank */

  double zmaxx;
  zmaxx = zminn +
          this->lattice_z * (this->cell_zmax - this->cell_zmin + 1 + this->dz) -
          this->adjThicknessZ;
  /* include all atoms on upper boundary */

  /* define booleans */
  *_isNotInXnegtv = (posnbr[0] < xminn); /* neighbor is out in -ve x-dir */
  *_isNotInXpostv = (posnbr[0] > xmaxx); /* neighbor is out in +ve x-dir */
  *_isNotInYnegtv = (posnbr[1] < yminn); /* neighbor is out in -ve y-dir */
  *_isNotInYpostv = (posnbr[1] > ymaxx); /* neighbor is out in +ve y-dir */
  *_isNotInZnegtv = (posnbr[2] < zminn); /* neighbor is out in -ve z-dir */
  *_isNotInZpostv = (posnbr[2] > zmaxx); /* neighbor is out in +ve z-dir */
}

/* Initializes the supercell coupling map, which contains three indices
   ssmap[i][j][k]. i and j index an atom's unit cell and atom
   (within the unit cell) respectively. and k is a nearest neigbor index
   for the atom indexed by i and j.  The mapping is to a unit cell  */
void geo_struct::SupercellReInit() {

  /* Readjust adjThicknessZ */
  if (this->is_tilted) {
    static int count_SupercellReInit = 0;
    //     if(count_SupercellReInit>0 && this->dy==0)
    //     {  this->adjThicknessY = 0.0;  }
    if (count_SupercellReInit > 0 && this->dz == 0) {
      this->adjThicknessZ = 0.0;
    }
    count_SupercellReInit++;
  } /*if is_tilted */
  /* End Readjust adjThicknessZ */

  // nml_memory_report("SupercellReInit - beginning");
  // deallocate memory

#ifndef ELIMINATE_SSMAP
  if (this->ssmap) {
    rm_i3matrix(&this->ssmap);
  }

#undef MEMCHECK
#ifdef MEMCHECK
  double ssmap_memory = 0.0, tl_memory = 0.0;
  ;
  ssmap_memory =
      this->N_Cell * this->AtomsPerCellMax() * this->NeighborsMax() * 4 / 1.0e6;
  tl_memory += ssmap_memory;
  printf("\n SupercellReinit  allocation in MB ssmap_memory=%g tl_memory=%g\n",
         ssmap_memory, tl_memory),
      fflush(stdout);
#endif

  // reallocate memory and reinitialize to -1
  this->ssmap =
      I3matrix(this->N_Cell, this->AtomsPerCellMax(), this->NeighborsMax());

  for (int l = 0; l < this->N_Cell; l++) {
    for (int m = 0; m < this->AtomsPerCellMax(); m++) {
      for (int n = 0; n < this->Neighbors(m); n++) {
        this->ssmap[l][m][n] = -1; // initialize
      }
    }
  }
#endif /* ELIMINATE_SSMAP */

  if (!this->hasPeriodicity())
    return;

  /* Loop through all unit cells in the structure */
  for (int l = 0; l < this->N_Cell; l++) {
    /* These are real space indexes for the unit cell */
    int i = this->cell__ijk[l][0];
    int j = this->cell__ijk[l][1];
    int k = this->cell__ijk[l][2];

    // Loop through all the atoms in the unit cell within simulation domain
    for (int m = 0; m < this->AtomsPerCellMax(); m++) {
      if (!this->AtomType[l][m])
        continue;

      // Loop through nearest neighbors.
      for (int n = 0; n < this->Neighbors(m); n++) {

#ifndef ELIMINATE_SSMAP
        this->ssmap[l][m][n] = -1; // ******* add this line!!!!
#endif                             /* ELIMINATE_SSMAP */

        // Get relative position of unit cell of virtual neighbor
        int i_nbr = i + this->NbrCell(m, n, 0);
        int j_nbr = j + this->NbrCell(m, n, 1);
        int k_nbr = k + this->NbrCell(m, n, 2);

        // Get the unit cell and atomic index for this neighbor
        int cindx = this->ijk__cell[i_nbr][j_nbr][k_nbr];
        // int aindx = this->NbrCell(m,n,3);

        int in, jn, kn;
        // if virtual cell is in fact real, no need to look further
        if (!this->is_tilted) { /* Connect atoms as follows when lattice is NOT
                                   tilted */
          if (cindx >= 0)
            continue;

          bool _isNotInX = (i_nbr > this->cell_xmax || i_nbr < this->cell_xmin);
          bool _isNotInY = (j_nbr > this->cell_ymax || j_nbr < this->cell_ymin);
          bool _isNotInZ = (k_nbr > this->cell_zmax || k_nbr < this->cell_zmin);

          /* if neighboring cell doesn't exist and we are not periodic,
             throw out */
          if ((_isNotInX && !this->isPeriodicX()) ||
              (_isNotInY && !this->isPeriodicY()) ||
              (_isNotInZ && !this->isPeriodicZ())) {
            continue;
          }

          /* Cell (i_nbr, j_nbr, k_nbr) is outside of simulation domain.
             Therefore, we record the periodic boundary coupling information
             into the super cell map ssmap.  */
          in = _isNotInX
                   ? i - (i_nbr - i) * (this->cell_xmax - this->cell_xmin) /
                             this->unitcelldivide
                   : i_nbr;
          jn = _isNotInY
                   ? j - (j_nbr - j) * (this->cell_ymax - this->cell_ymin) /
                             this->unitcelldivide
                   : j_nbr;
          kn = _isNotInZ
                   ? k - (k_nbr - k) * (this->cell_zmax - this->cell_zmin) /
                             this->unitcelldivide
                   : k_nbr;
        } /*if !is_tilted */

        if (this->is_tilted) { /* Connect atoms as follows when lattice is
                                  tilted */
          /* Initialize all booleans to false */
          bool _isNotInXnegtv = false;
          bool _isNotInXpostv = false;
          bool _isNotInYnegtv = false;
          bool _isNotInYpostv = false;
          bool _isNotInZnegtv = false;
          bool _isNotInZpostv = false;

          this->isNbrOut(l, m, n, &_isNotInXnegtv, &_isNotInXpostv,
                         &_isNotInYnegtv, &_isNotInYpostv, &_isNotInZnegtv,
                         &_isNotInZpostv);

          if (!_isNotInXnegtv && !_isNotInXpostv && !_isNotInYnegtv &&
              !_isNotInYpostv && !_isNotInZnegtv && !_isNotInZpostv)
            continue;

          /* if neighbor is out is negtv or postv direction and lattice is not
             periodic in this direction, no need to connect */
          if (((_isNotInXnegtv || _isNotInXpostv) && !this->isPeriodicX()) ||
              ((_isNotInYnegtv || _isNotInYpostv) && !this->isPeriodicY()) ||
              ((_isNotInZnegtv || _isNotInZpostv) && !this->isPeriodicZ())) {
            continue;
          }

          double delta = 1e-6;
          int thickNoY =
              int((this->lattice_y *
                       (this->cell_ymax - this->cell_ymin + 1 + this->dy) +
                   delta /*- this->adjThicknessY*/) /
                  this->lattice_y);
          int thickNoZ =
              int((this->lattice_z *
                       (this->cell_zmax - this->cell_zmin + 1 + this->dz) +
                   delta - this->adjThicknessZ) /
                  this->lattice_z);

          int noUnitCellsX =
              int((this->cell_xmax - this->cell_xmin + 1) / this->dx);

          if (_isNotInXnegtv) {
            in = i - (i_nbr - i) * (this->cell_xmax - this->cell_xmin) /
                         this->unitcelldivide;
            jn = j_nbr + noUnitCellsX * this->dy;
            kn = k_nbr + noUnitCellsX * this->dz;
          } else if (_isNotInXpostv) {
            in = i - (i_nbr - i) * (this->cell_xmax - this->cell_xmin) /
                         this->unitcelldivide;
            jn = j_nbr - noUnitCellsX * this->dy;
            kn = k_nbr - noUnitCellsX * this->dz;
          } else {
            in = i_nbr;
            jn = j_nbr;
            kn = k_nbr;
          }

          if (_isNotInYnegtv)
            jn = jn + thickNoY / this->unitcelldivide;
          else if (_isNotInYpostv)
            jn = jn - thickNoY / this->unitcelldivide;
          else
            jn = jn;

          if (_isNotInZnegtv)
            kn = kn + thickNoZ / this->unitcelldivide;
          else if (_isNotInZpostv)
            kn = kn - thickNoZ / this->unitcelldivide;
          else
            kn = kn;
        } /*if is_tilted */

#ifdef ELIMINATE_SSMAP
        if (!this->is_tilted)
          this->ijk__cell[i_nbr][j_nbr][k_nbr] =
              this->ijk__cell[in][jn][kn] - (this->N_Cell) - 1;
#else  /* ELIMINATE_SSMAP */
        this->ssmap[l][m][n] = this->ijk__cell[in][jn][kn];
// printf("ssmap[%d][%d][%d] = %d in=%d jn=%d
// kn=%d\n",l,m,n,this->ssmap[l][m][n],in,jn,kn);
#endif /* ELIMINATE_SSMAP */
      }
    }
  }
  nml_memory_report("SupercellReInit - end");
}

int geo_struct::get_cindx_tilted(int l, int m, int n) {
  int i = this->cell__ijk[l][0];
  int j = this->cell__ijk[l][1];
  int k = this->cell__ijk[l][2];
  int i_nbr = i + this->NbrCell(m, n, 0);
  int j_nbr = j + this->NbrCell(m, n, 1);
  int k_nbr = k + this->NbrCell(m, n, 2);

  int in, jn, kn;

  /* Initialize all booleans to false */
  bool _isNotInXnegtv = false;
  bool _isNotInXpostv = false;
  bool _isNotInYnegtv = false;
  bool _isNotInYpostv = false;
  bool _isNotInZnegtv = false;
  bool _isNotInZpostv = false;

  this->isNbrOut(l, m, n, &_isNotInXnegtv, &_isNotInXpostv, &_isNotInYnegtv,
                 &_isNotInYpostv, &_isNotInZnegtv, &_isNotInZpostv);

  if (!_isNotInXnegtv && !_isNotInXpostv && !_isNotInYnegtv &&
      !_isNotInYpostv && !_isNotInZnegtv && !_isNotInZpostv)
    return (this->ijk__cell[i_nbr][j_nbr][k_nbr]);

  /* if neighbor is out is negtv or postv direction and lattice is not
     periodic in this direction, no need to connect */
  if (((_isNotInXnegtv || _isNotInXpostv) && !this->isPeriodicX()) ||
      ((_isNotInYnegtv || _isNotInYpostv) && !this->isPeriodicY()) ||
      ((_isNotInZnegtv || _isNotInZpostv) && !this->isPeriodicZ())) {
    return (this->ijk__cell[i_nbr][j_nbr][k_nbr]);
  }

  double delta = 1e-6;
  int thickNoY = int(
      (this->lattice_y * (this->cell_ymax - this->cell_ymin + 1 + this->dy) +
       delta /*- this->adjThicknessY*/) /
      this->lattice_y);
  int thickNoZ = int(
      (this->lattice_z * (this->cell_zmax - this->cell_zmin + 1 + this->dz) +
       delta - this->adjThicknessZ) /
      this->lattice_z);

  int noUnitCellsX = int((this->cell_xmax - this->cell_xmin + 1) / this->dx);

  if (_isNotInXnegtv) {
    in = i - (i_nbr - i) * (this->cell_xmax - this->cell_xmin) /
                 this->unitcelldivide;
    jn = j_nbr + noUnitCellsX * this->dy;
    kn = k_nbr + noUnitCellsX * this->dz;
  } else if (_isNotInXpostv) {
    in = i - (i_nbr - i) * (this->cell_xmax - this->cell_xmin) /
                 this->unitcelldivide;
    jn = j_nbr - noUnitCellsX * this->dy;
    kn = k_nbr - noUnitCellsX * this->dz;
  } else {
    in = i_nbr;
    jn = j_nbr;
    kn = k_nbr;
  }

  if (_isNotInYnegtv)
    jn = jn + thickNoY / this->unitcelldivide;
  else if (_isNotInYpostv)
    jn = jn - thickNoY / this->unitcelldivide;
  else
    jn = jn;

  if (_isNotInZnegtv)
    kn = kn + thickNoZ / this->unitcelldivide;
  else if (_isNotInZpostv)
    kn = kn - thickNoZ / this->unitcelldivide;
  else
    kn = kn;

  return (this->ijk__cell[in][jn][kn]);

} /* geo_struct::get_cindx_tilted */

/* this method sets this->surfaceInfo[][] */
void geo_struct::FindSurfaceAtoms() {
  // Loop through all unit cells in the structure
  for (int l = 0; l < this->N_Cell; l++) {

    // These are real space indexes for the unit cell
    int i = this->cell__ijk[l][0];
    int j = this->cell__ijk[l][1];
    int k = this->cell__ijk[l][2];

    // Loop through all the atoms in the unit cell
    for (int m = 0; m < this->AtomsPerCellMax(); m++) {
      this->surfaceInfo[l][m] = 0; // initialize;
      if (!this->AtomType[l][m])
        continue;

      if (i == this->cell_xmin)
        this->surfaceInfo[l][m] += 4;
      if (i == this->cell_xmax)
        this->surfaceInfo[l][m] += 8;
      if (j == this->cell_ymin)
        this->surfaceInfo[l][m] += 16;
      if (j == this->cell_ymax)
        this->surfaceInfo[l][m] += 32;
      if (k == this->cell_zmin)
        this->surfaceInfo[l][m] += 64;
      if (k == this->cell_zmax)
        this->surfaceInfo[l][m] += 128;

      /* Added for tilt */
      //         if (j == this->cell_ymin+1)  this->surfaceInfo[l][m] += 16;
      //         if (j == this->cell_ymax-1)  this->surfaceInfo[l][m] += 32;
      if (k == this->cell_zmin + 1)
        this->surfaceInfo[l][m] += 64;
      if (k == this->cell_zmax - 1)
        this->surfaceInfo[l][m] += 128;
      /* End Added for tilt */

      int nbrcnt1 = 0, nbrcnt2 = 0;
      for (int n = 0; n < this->Neighbors(m); n++) {

        // Get relative position of unit cell where this neighbor lives.
        int i_nbr = i + this->NbrCell(m, n, 0);
        int j_nbr = j + this->NbrCell(m, n, 1);
        int k_nbr = k + this->NbrCell(m, n, 2);

        // Get the unit cell and atomic index for this neighbor
        int aindx = this->NbrCell(m, n, 3);

        int cindx = this->ijk__cell[i_nbr][j_nbr][k_nbr];

        if (cindx >= 0 && this->AtomType[cindx][aindx]) {
          nbrcnt1++;
          nbrcnt2++;
        } else {

#ifdef ELIMINATE_SSMAP
          if (hasPeriodicity() && ((cindx + this->N_Cell + 1) >= 0)) {
            cindx = cindx + this->N_Cell + 1;
            if (this->is_tilted)
              cindx = this->get_cindx_tilted(l, m, n);
#else  /* ELIMINATE_SSMAP */
          if (hasPeriodicity() && this->ssmap[l][m][n] != -1) {
            cindx = this->ssmap[l][m][n];
#endif /* ELIMINATE_SSMAP */

            if (cindx <= this->N_Cell) { /* Addition for tilt */
              if (this->AtomType[cindx][aindx])
                nbrcnt1++;
            } /* Addition for tilt */
          }
        }
      }

      if (nbrcnt1 != this->cellgeom.N_MaxNeighbor) {
        this->surfaceInfo[l][m] += 1;
      }
      if (nbrcnt2 != this->cellgeom.N_MaxNeighbor) {
        this->surfaceInfo[l][m] += 2;
      }
      /*
               if ( (!nbrcnt1 || !nbrcnt2) &&
                    (this->N_Cell>1 || this->AtomsPerCellMax()>1) )
                  die("unbonded atom encountered in FindSurfaceAtoms()\n");
      */
    }
  }
}

/* Fill the 2-d array geo.is_DiffNeighbor which is used in the Hamiltonian
   construction to check whether or not the parameter vector needs to be
   averaged over different kinds of neighbors.  */
void geo_struct::fill_is_DiffNeighbor() {
  // Loop through all unit cells in the structure
  for (int l = 0; l < this->N_Cell; l++) {
    // These are real space indexes for the unit cell
    int i = this->cell__ijk[l][0];
    int j = this->cell__ijk[l][1];
    int k = this->cell__ijk[l][2];

    // Loop through all the real atoms in the unit cell
    for (int m = 0; m < this->AtomsPerCellMax(); m++) {
      if (!this->AtomType[l][m])
        continue;

      /* Hamiltonian diagonal.  If anion send anion basis vector,
         otherwise use cation vector . */
      if (this->AtomType[l][m] == 12 || this->AtomType[l][m] == 16) {
        int atom_id_first = 0;
        // For anion we need a weighted avg of the orbital energies
        int nbrcnt = 0;
        int n = 0;
        // find the first neighbor
        while (n < this->Neighbors(m)) {
          /* Get relative position of unit cell where this
             neighbor lives. */
          int id = this->NbrCell(m, n, 0);
          int jd = this->NbrCell(m, n, 1);
          int kd = this->NbrCell(m, n, 2);
          // Get the unit cell and atomic index for this neighbor
          int cindx = this->ijk__cell[i + id][j + jd][k + kd];
          int aindx = this->NbrCell(m, n, 3);
          if (cindx >= 0) {
            if (this->AtomType[cindx][aindx]) {
              nbrcnt++;
              if (!atom_id_first) {
                atom_id_first = this->AtomType[cindx][aindx];
              }
              if (this->AtomType[cindx][aindx] != atom_id_first) {
                nbrcnt++;
                // NoNeed_isDiffNbr this->isDiffNbr[l][m]=1;
                break;
              }
            }
          }
          n++;
        }
        if (!nbrcnt && (this->N_Cell > 1 || this->AtomsPerCellMax() > 1)) {
          die("unbonded atom encountered in fill_is_DiffNeighbor()\n");
        }
      } else {
        /* For cation just use any neigbor for now...
           still thinking about this one */
        /* do not examine the neighbors at all for now */
      }
    }
  }
}

/* Function scans through the whole structure and removes atoms that
   have only a certain set of bonds.  if "bond_to_remove" is set to 0
   we will remove all the atoms that have no neighbors.   if "bond_to_remove"
   is set to 1 we will remove the atoms that have only a single neighbor.
   The function returns the number of removed atoms. */
int geo_struct::remove_few_bonded_atoms(int bond_to_remove) {
  int removed_atoms = 0;

  if (this->hasPeriodicity() && bond_to_remove > 0) {
    die("May not remove single-bonded atoms under periodic boundary conditions "
        "(remove_few_bonded_atoms)");
  }

  if (bond_to_remove != 0 && bond_to_remove != 1) {
    die("Invalid argument:  remove_few_bonded_atoms(d, bond_to_remove=%d)\n",
        bond_to_remove);
  }

  /* Check for unbonded atoms and remove them from the list */
  if (this->N_Cell > 1 || this->AtomsPerCellMax() > 1) {
    for (int i = 0; i < this->N_Cell; i++) {
      for (int j = 0; j < this->AtomsPerCellMax(); j++) {
        if (!this->AtomType[i][j])
          continue;

        int nb_flag = 0;
        for (int k = 0; k < this->Neighbors(j); k++) {
          int id = this->cell__ijk[i][0] + this->NbrCell(j, k, 0);
          int jd = this->cell__ijk[i][1] + this->NbrCell(j, k, 1);
          int kd = this->cell__ijk[i][2] + this->NbrCell(j, k, 2);
          int l = this->ijk__cell[id][jd][kd];
          if (l >= 0)
            if (this->AtomType[l][this->NbrCell(j, k, 3)] != 0)
              nb_flag++;
        }
        if (nb_flag == bond_to_remove) {
          this->AtomType[i][j] = 0;
          removed_atoms++;
        }
      }
    }
  }

  if (removed_atoms > 0)
    cout_master << "Removed " << removed_atoms << " atoms with "
                << bond_to_remove << " bond from list\n";

#undef _DO_REMOVE_RECURSE
#ifdef _DO_REMOVE_RECURSE
  /* if there were single bonded atoms that were removed in one pass,
     we must call this routine again to remove further stubs.
     This routine would ultimately eliminate a chain of atoms. */
  if (bond_to_remove == 1 && removed_atoms > 0)
    return (removed_atoms + remove_few_bonded_atoms(bond_to_remove));
#endif

  return removed_atoms;
}

int geo_struct::remove_user_defined_atoms(void) {
  int removed_atoms = 0;
  for (int Zc = 0; Zc < this->N_Cell; Zc++) {
    for (int Za = 0; Za < this->AtomsPerCellMax(); Za++) {
      if (Zc == 0 && Za == 1)
        continue;
      else if (Zc == 2 && Za == 0)
        continue;
      else if (Zc == 2 && Za == 1)
        continue;
      else if (Zc == 2 && Za == 2)
        continue;
      else if (Zc == 6 && Za == 1)
        continue;
      else {
        this->AtomType[Zc][Za] = 0;
        removed_atoms++;
      }
    }
  }
  return removed_atoms;
}

/* Function scans through the whole structure and removes
   a certain species of atoms from the surface.
   */
int geo_struct::remove_other_surface_species() {
  int removed_atoms = 0;

  if (this->hasPeriodicity()) {
    die("ERROR in remove_other_surface_species:\n"
        "Must not remove atoms from the structure\n"
        "if periodic boundary conditions are active\n");
  }

  cout_master << "homogenize the surface\n";
  /* Check for unbonded atoms and remove them from the list */
  if (this->N_Cell > 1 || this->AtomsPerCellMax() > 1) {
    for (int i = 0; i < this->N_Cell; i++) {
      for (int j = 0; j < this->AtomsPerCellMax(); j++) {
        if (!this->isOnSurface(i, j))
          continue;

        if (this->AtomType[i][j] == EL_As || this->AtomType[i][j] == EL_Sb ||
            this->AtomType[i][j] == EL_Se || this->AtomType[i][j] == EL_Ge) {
          // homo 2
          // if (this->AtomType[i][j]!= EL_As && this->AtomType[i][j]!=EL_Sb){
          this->AtomType[i][j] = 0;
          removed_atoms++;
        }
      }
    }
  }

  /* if there were  atoms removed in one pass, we must call this
     routine again to remove further stubs. */
  if (removed_atoms > 0) {
    cout_master << "Removed " << removed_atoms
                << " atoms with from the surface\n";
    return (removed_atoms + this->remove_few_bonded_atoms(0));
  }

  return removed_atoms;
}

bool geo_struct::isInvalidPeriodicity() {
  bool flag = false;

  for (int l = 0; l < this->N_Cell; l++) {
    /* These are real space indexes for the unit cell */
    int i = this->cell__ijk[l][0];
    int j = this->cell__ijk[l][1];
    int k = this->cell__ijk[l][2];

    for (int m = 0; m < this->AtomsPerCellMax(); m++) {
      if (!this->AtomType[l][m])
        continue;

      // Loop through nearest neighbors.
      for (int n = 0; n < this->Neighbors(m); n++) {
        // Get relative position of unit cell of virtual neighbor
        int i_nbr = i + this->NbrCell(m, n, 0);
        int j_nbr = j + this->NbrCell(m, n, 1);
        int k_nbr = k + this->NbrCell(m, n, 2);

        // Get the unit cell and atomic index for this neighbor

#ifdef ELIMINATE_SSMAP
        int l_nbr;
        if (this->ijk__cell[i_nbr][j_nbr][k_nbr] < 0 &&
            ((this->ijk__cell[i_nbr][j_nbr][k_nbr] + this->N_Cell + 1) >= 0)) {
          l_nbr = this->ijk__cell[i_nbr][j_nbr][k_nbr] + this->N_Cell + 1;
          if (this->is_tilted)
            l_nbr = this->get_cindx_tilted(l, m, n);
#else  /* ELIMINATE_SSMAP */
        if (this->ijk__cell[i_nbr][j_nbr][k_nbr] < 0 &&
            this->ssmap[l][m][n] != -1) {
          int l_nbr = this->ssmap[l][m][n];
#endif /* ELIMINATE_SSMAP */

          int m_nbr = this->NbrCell(m, n, 3);

#ifdef DEBUG_NEW_SHAPE
          printf("Checking %d %d %d:  ", l, m, n);
#endif
          if (!this->AtomType[l_nbr][m_nbr]) {
            flag = true;
#ifdef DEBUG_NEW_SHAPE
            printf("PROBLEM!!!");
#endif
          }
#ifdef DEBUG_NEW_SHAPE
          printf("\n");
#endif
        }
      }
    }
  }
  return flag;
}

// This may be broken -- fix
/* Re-adjusts ijk__cell to reflect periodic boundary conditions. */
void geo_struct::remap_periodic_bc() {
  int i, j, k, l, m, n, in, jn, kn;
  i3matrix buf;

  if (!this->hasPeriodicity())
    return; // no periodicity in any direction, so do nothing

  cout_master << "Accomodating periodicity.\n";

  /* Loop through all unit cells in the structure */
  for (l = 0; l < this->N_Cell; l++) {
    int ii, jj;

    /* create buff matrix and initialize its elements to -1 */
    buf = I3matrix(this->AtomsPerCellMax(), this->NeighborsMax(), 2);
    for (ii = 0; ii < this->AtomsPerCellMax(); ii++) {
      for (jj = 0; jj < this->NeighborsMax(); jj++) {
        buf[ii][jj][0] = buf[ii][jj][1] = -1;
      }
    }

    /* These are real space indexes for the unit cell */
    i = this->cell__ijk[l][0];
    j = this->cell__ijk[l][1];
    k = this->cell__ijk[l][2];

    /* Loop through all the atoms in the unit cell within simulation domain */
    for (m = 0; m < this->AtomsPerCellMax(); m++) {
      if (!this->AtomType[l][m])
        continue;

      /* Loop through nearest neighbors. */
      for (n = 0; n < this->Neighbors(m); n++) {
        short _isNotInX, _isNotInY, _isNotInZ;

        /* Get relative position of unit cell of virtual neighbor */
        int i_nbr = i + this->NbrCell(m, n, 0);
        int j_nbr = j + this->NbrCell(m, n, 1);
        int k_nbr = k + this->NbrCell(m, n, 2);

        /* Get the unit cell and atomic index for this neighbor */
        int cindx = this->ijk__cell[i_nbr][j_nbr][k_nbr];
        int aindx = this->NbrCell(m, n, 3);

        /* if virtual cell is in fact real, no need to look further */
        if (cindx >= 0)
          continue;

        _isNotInX = (i_nbr > this->cell_xmax || i_nbr < this->cell_xmin);
        _isNotInY = (j_nbr > this->cell_ymax || j_nbr < this->cell_ymin);
        _isNotInZ = (k_nbr > this->cell_zmax || k_nbr < this->cell_zmin);

        /* if neighboring cell doesn't exist and we are not periodic,
           throw out */
        if ((_isNotInX && !this->isPeriodicX()) ||
            (_isNotInY && !this->isPeriodicY()) ||
            (_isNotInZ && !this->isPeriodicZ())) {
          continue;
        }

        /* Cell (i_nbr, j_nbr, k_nbr) is outside of simulation domain.
           Therefore, we record the periodic boundary coupling information
           into the super cell map ssmap.  */
        in = _isNotInX ? i - (i_nbr - i) * (this->cell_xmax - this->cell_xmin) /
                                 this->unitcelldivide
                       : i_nbr;
        jn = _isNotInY ? j - (j_nbr - j) * (this->cell_ymax - this->cell_ymin) /
                                 this->unitcelldivide
                       : j_nbr;
        kn = _isNotInZ ? k - (k_nbr - k) * (this->cell_zmax - this->cell_zmin) /
                                 this->unitcelldivide
                       : k_nbr;

        buf[m][n][0] = this->ijk__cell[in][jn][kn];
        buf[m][n][1] = aindx;

#ifndef ELIMINATE_SSMAP
        this->ijk__cell[i_nbr][j_nbr][k_nbr] = this->ijk__cell[in][jn][kn];
#endif /* ELIMINATE_SSMAP */
      }
    }
  }
}

// useful if want to use periodic bc for an electronic calculation with
// a more sensible choice of lattice constant used for the (unperiodic)
// strain calculation
void geo_struct::changeLattConst(real ax_new, real ay_new, real az_new) {

#ifdef DISPL_SINGLE_PRECISION
  f3tensor dsp = F3tensor(0, N_Cell - 1, 0, this->AtomsPerCellMax() - 1, 0, 2);
#else  /* DISPL_SINGLE_PRECISION */
  r3tensor dsp = R3tensor(0, N_Cell - 1, 0, this->AtomsPerCellMax() - 1, 0, 2);
#endif /* DISPL_SINGLE_PRECISION */

  for (int Zcell = 0; Zcell < this->N_Cell; Zcell++) {
    for (int Zatom = 0; Zatom < this->AtomsPerCellMax(); Zatom++) {
      if (!this->AtomType[Zcell][Zatom])
        continue;

      int *Rcell = cell__ijk[Zcell];
      real *Ratom = PositionInCell(Zatom);
#ifdef DISPL_SINGLE_PRECISION
      nml_float *dR = l_m_2_atomdsp[Zcell][Zatom];
      nml_float *dR_new = dsp[Zcell][Zatom];
#else  /* DISPL_SINGLE_PRECISION */
      real *dR = l_m_2_atomdsp[Zcell][Zatom];
      real *dR_new = dsp[Zcell][Zatom];
#endif /* DISPL_SINGLE_PRECISION */

      (*dR_new++) = (lattice_x - ax_new) * ((*Rcell++) + (*Ratom++)) + (*dR++);
      (*dR_new++) = (lattice_y - ay_new) * ((*Rcell++) + (*Ratom++)) + (*dR++);
      (*dR_new++) = (lattice_z - az_new) * ((*Rcell) + (*Ratom)) + (*dR);
    }
  }

  this->lattice_x = ax_new;
  this->lattice_y = ay_new;
  this->lattice_z = az_new;

#ifdef DISPL_SINGLE_PRECISION
  rm_f3tensor(&this->l_m_2_atomdsp, 0);
#else  /* DISPL_SINGLE_PRECISION */
  rm_r3tensor(&this->l_m_2_atomdsp, 0);
#endif /* DISPL_SINGLE_PRECISION */

  this->l_m_2_atomdsp = dsp;
}

void geo_struct::updateOnlyLattConst(real ax_new, real ay_new, real az_new) {
  this->lattice_x = ax_new;
  this->lattice_y = ay_new;
  this->lattice_z = az_new;
}

/* temporary kludge to see if we can do away with ssmap */
void geo_struct::adjust_for_periodic_bc() {
  if (!this->hasPeriodicity())
    return;

  /* Loop through all unit cells in the structure */
  for (int l = 0; l < this->N_Cell; l++) {
    /* These are real space indexes for the unit cell */
    int i = this->cell__ijk[l][0];
    int j = this->cell__ijk[l][1];
    int k = this->cell__ijk[l][2];

    // Loop through all the atoms in the unit cell within simulation domain
    for (int m = 0; m < this->AtomsPerCellMax(); m++) {
      if (!this->AtomType[l][m])
        continue;

      /* Loop through nearest neighbors. */
      for (int n = 0; n < this->Neighbors(m); n++) {
        bool _isNotInX, _isNotInY, _isNotInZ;

        // Get relative position of unit cell of virtual neighbor
        int i_nbr = i + this->NbrCell(m, n, 0);
        int j_nbr = j + this->NbrCell(m, n, 1);
        int k_nbr = k + this->NbrCell(m, n, 2);

        // Get the unit cell and atomic index for this neighbor
        int cindx = this->ijk__cell[i_nbr][j_nbr][k_nbr];

        // if virtual cell is in fact real, no need to look further
        if (cindx >= 0)
          continue;

        _isNotInX = (i_nbr > this->cell_xmax || i_nbr < this->cell_xmin);
        _isNotInY = (j_nbr > this->cell_ymax || j_nbr < this->cell_ymin);
        _isNotInZ = (k_nbr > this->cell_zmax || k_nbr < this->cell_zmin);

        /* if neighboring cell doesn't exist and we are not periodic,
           throw out */
        if ((_isNotInX && !this->isPeriodicX()) ||
            (_isNotInY && !this->isPeriodicY()) ||
            (_isNotInZ && !this->isPeriodicZ())) {
          continue;
        }

        /* Cell (i_nbr, j_nbr, k_nbr) is outside of simulation domain.
           Therefore, we record the periodic boundary coupling information
           into the super cell map ssmap.  */
        int in = _isNotInX
                     ? i - (i_nbr - i) * (this->cell_xmax - this->cell_xmin) /
                               this->unitcelldivide
                     : i_nbr;
        int jn = _isNotInY
                     ? j - (j_nbr - j) * (this->cell_ymax - this->cell_ymin) /
                               this->unitcelldivide
                     : j_nbr;
        int kn = _isNotInZ
                     ? k - (k_nbr - k) * (this->cell_zmax - this->cell_zmin) /
                               this->unitcelldivide
                     : k_nbr;

#ifndef ELIMINATE_SSMAP
        this->ijk__cell[i_nbr][j_nbr][k_nbr] = this->ijk__cell[in][jn][kn];
#endif /* ELIMINATE_SSMAP */
      }
    }
  }
}

/* This parses a hyperchem file for atom locations, number
   of neighbors and total number of atoms */
#define DELIMITER_HC1 " \n"
void geo_struct::InitializeForHyperChem(const char *inputfile) {
  this->cellgeom.InitializeForHyperChem(inputfile);

  char line[100], item[50], element[3], ldesc[20], *tail, *token;
  int str_to_atomID(const char *);

  char *filename = NULL;
  filename = n3d_strdup_n(inputfile);
  n3d_FileTypeSet(&filename, "nd_hin", 1);

  /* Initialize arrays */
  this->N_Cell = 1;
  this->cell__ijk = I2tensor(0, 0, 0, 2);
#ifdef MEMORY_SAVE_NEWSTRUCT
  if (this->AtomType)
    rm_b2tensor(&this->AtomType, 0);
  this->AtomType = B2tensor(0, 0, 0, this->AtomsPerCellMax() - 1);
#else  /* MEMORY_SAVE_NEWSTRUCT  */
  if (this->AtomType)
    rm_i2tensor(&this->AtomType, 0);
  this->AtomType = I2tensor(0, 0, 0, this->AtomsPerCellMax() - 1);
#endif /* MEMORY_SAVE_NEWSTRUCT  */
  this->ijk__cell = I3tensor(0, 0, 0, 0, 0, 0);

  FILE *fhin = fopen(filename, "r");
  while (fgets(line, 100, fhin)) {
    /* Read line descriptor */
    token = strtok(line, DELIMITER_HC1);
    strcpy(ldesc, token);
    if (strcmp(ldesc, "atom") == 0) {

      /* Second column is atom number, save this one */
      token = strtok(NULL, DELIMITER_HC1);
      strcpy(item, token);
      int j = strtol(item, &tail, 0) - 1;

      /* 3rd column is junk read and forget */
      token = strtok(NULL, DELIMITER_HC1);

      /* 4th column is element, we need this */
      token = strtok(NULL, DELIMITER_HC1);
      strcpy(element, token);
      this->AtomType[0][j] = str_to_atomID(element);

      /* don't need rest of the columns */
      token = strtok(NULL, DELIMITER_HC1);
      token = strtok(NULL, DELIMITER_HC1);
      token = strtok(NULL, DELIMITER_HC1);
      token = strtok(NULL, DELIMITER_HC1);
      token = strtok(NULL, DELIMITER_HC1);
      token = strtok(NULL, DELIMITER_HC1);
      token = strtok(NULL, DELIMITER_HC1);

      /* Loop through neighbors and get their index */
      for (int i = 0; i < this->Neighbors(j); i++) {
        /* Get neighbor index */
        token = strtok(NULL, DELIMITER_HC1);
        /* Next column is bond type, forget it for now */
        token = strtok(NULL, DELIMITER_HC1);
      }
    }
  }
  str_free(filename);
}

ostream &operator<<(ostream &o, const geo_struct &g) {
  if (mpi_n3d_id)
    return o;

  const char LINE[] = "========================================================"
                      "================";

  int Nx = d_1h(g.ijk__cell) - d_1l(g.ijk__cell) + 1;
  int Ny = d_2h(g.ijk__cell) - d_2l(g.ijk__cell) + 1;
  int Nz = d_3h(g.ijk__cell) - d_3l(g.ijk__cell) + 1;

  o << LINE << "\n   GEO_STRUCT\n" << LINE << "\n";
  o << "Actual number of cells = " << g.N_Cell << "\n";
  o << "X = " << d_1l(g.ijk__cell) << " ... " << d_1h(g.ijk__cell) << "\n";
  o << "Y = " << d_2l(g.ijk__cell) << " ... " << d_2h(g.ijk__cell) << "\n";
  o << "Z = " << d_3l(g.ijk__cell) << " ... " << d_3h(g.ijk__cell) << "\n";
  o << "xmin = " << g.cell_xmin << "   xmax = " << g.cell_xmax << "\n";
  o << "ymin = " << g.cell_ymin << "   ymax = " << g.cell_ymax << "\n";
  o << "zmin = " << g.cell_zmin << "   zmax = " << g.cell_zmax << "\n";
  o << "lattice constant (x,y,z) = (" << g.lattice_x << "," << g.lattice_y
    << "," << g.lattice_z << ")\n";
  o << "unitcelldivide = " << g.unitcelldivide << "\n";
  o << "periodicity:  " << (g.isPeriodicX() ? "X " : " ")
    << (g.isPeriodicY() ? "Y " : " ") << (g.isPeriodicZ() ? "Z " : " ")
    << (g.StrainIsPeriodicX() ? "str.X " : " ")
    << (g.StrainIsPeriodicY() ? "str.Y " : " ")
    << (g.StrainIsPeriodicZ() ? "str.Z " : " ")
    << (g.ElIsPeriodicX() ? "el.X " : " ")
    << (g.ElIsPeriodicY() ? "el.Y " : " ")
    << (g.ElIsPeriodicZ() ? "el.Z " : " ") << "\n";

  o << g.cellgeom;

  o << LINE << "\n"
    << "cell(X,Y,Z):  atomtype(atom in cell)\n"
    << LINE << "\n";
  for (int l = 0; l < g.N_Cell; l++) {
    for (int m = 0; m < g.AtomsPerCellMax(); m++) {
      o << setw(4) << (g.atom_strain_only[l][m] ? 'X' : ' ') << setw(4)
        << g.AtomType[l][m]
        << (g.isOnSurface(l, m) ? 'S'
                                : (g.isOnSurfaceIgnoreBC(l, m) ? 's' : ' '));
    }
    o << "\n";
  }

  o.setf(ios::scientific, ios::floatfield);
  o.precision(6);

  o << LINE << "\n"
    << "  i j k m l_m_2_offset_ham_tot[l][m] l_m_2_atomdsp[l][m]\n"
    << LINE << "\n";
  for (int l = 0; l < g.N_Cell; l++) {
    int *Zr = g.cell__ijk[l];

    for (int m = 0; m < g.AtomsPerCellMax(); m++) {
      o << setw(6) << Zr[0] << "  " << setw(6) << Zr[1] << "  " << setw(6)
        << Zr[2] << "     " << setw(6) << m << "  " << setw(6)
        << g.l_m_2_offset_ham_tot[l][m] << "  " << setw(10)
        << g.l_m_2_atomdsp[l][m][0] << "  " << setw(10)
        << g.l_m_2_atomdsp[l][m][1] << "  " << setw(10)
        << g.l_m_2_atomdsp[l][m][2] << "\n";
    }
    o << "\n";
  }

  o.precision(0);
  o.setf(ios::floatfield);

  o << LINE << "\n"
    << "  i j k g.ijk__cell[i][j][k]\n"
    << LINE << "\n";
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        o << setw(3) << i << "  " << setw(3) << j << "  " << setw(3) << k
          << "  " << setw(6) << g.ijk__cell[i][j][k] << "\n";
      }
    }
  }

  o << LINE << "\n"
#ifdef ELIMINATE_SSMAP
    << "  (cell, atom, nbr) -> cell__ijk[i_nbr][j_nbr][k_nbr] + N_Cell + 1\n"
#else  /* ELIMINATE_SSMAP */
    << "  (cell, atom, nbr) -> ssmap[cell][atom][nbr]\n"
#endif /* ELIMINATE_SSMAP */
    << LINE << "\n";

#ifndef ELIMINATE_SSMAP
  for (int c = 0; c < g.N_Cell; c++) {
    for (int a = 0; a < g.cellgeom.AtomsPerCellMax(); a++) {
      for (int n = 0; n < g.cellgeom.NeighborsMax(); n++) {
        if (g.ssmap[c][a][n] == -1)
          continue;

        o << setw(3) << c << "  " << setw(3) << a << "  " << setw(3) << n
          << "  " << setw(3) << g.ssmap[c][a][n] << "\n";
      }
    }
  }
#endif /* ELIMINATE_SSMAP */

  return o;
}

void geo_struct::getPosition_unstrained(real *pos, real unstrnd_latt, int Zcell,
                                        int Zatom) {
  int *Rcell = cell__ijk[Zcell];
  real *Ratom = PositionInCell(Zatom);
  if (isCellWurtzite()) {
    pos[0] = unstrnd_latt * (0.5 * (Rcell[0] + Rcell[1]) + Ratom[0]);
    pos[1] =
        unstrnd_latt * (0.5 * sqrt(3.0) * (Rcell[0] - Rcell[1]) + Ratom[1]);
    pos[2] = unstrnd_latt * sqrt(8.0 / 3.0) * (Rcell[2] + Ratom[2]);
  }

  /* Added for ZB_110 */
  else if (isCellZB_110() || isCellZB_110_small() || isCellZB_111() ||
           isCellUserDef()) {
    pos[0] = get_a_latt_rot_norm(0) * unstrnd_latt * (Rcell[0] + Ratom[0]);
    pos[1] = get_a_latt_rot_norm(1) * unstrnd_latt * (Rcell[1] + Ratom[1]);
    pos[2] = get_a_latt_rot_norm(2) * unstrnd_latt * (Rcell[2] + Ratom[2]);
  }
  /* End Added for ZB_110 */

  else {
    pos[0] = unstrnd_latt * (Rcell[0] + Ratom[0]);
    pos[1] = unstrnd_latt * (Rcell[1] + Ratom[1]);
    pos[2] = unstrnd_latt * (Rcell[2] + Ratom[2]);
  }
}

void geo_struct::getPosition_strained(real *pos, int Zcell, int Zatom) {

  int *Rcell = cell__ijk[Zcell];
  real *Ratom = PositionInCell(Zatom);
#define MEMCHECK2
#ifdef MEMCHECK2
#ifdef DISPL_SINGLE_PRECISION
  nml_float *dR = NULL;
#else  /* DISPL_SINGLE_PRECISION */
  real *dR = NULL;
#endif /* DISPL_SINGLE_PRECISION */

  if (l_m_2_atomdsp) {
    dR = l_m_2_atomdsp[Zcell][Zatom];
    if (isCellWurtzite()) {
      pos[0] = lattice_x * (0.5 * (Rcell[0] + Rcell[1]) + Ratom[0]) + dR[0];
      pos[1] =
          lattice_y * (0.5 * sqrt(3.0) * (Rcell[0] - Rcell[1]) + Ratom[1]) +
          dR[1];
      pos[2] = lattice_z * (Rcell[2] + Ratom[2]) + dR[2];
    }

    /* Added for ZB_110 */
    else if (isCellZB_110() || isCellZB_110_small() || isCellZB_111() ||
             isCellUserDef()) {
      pos[0] = lattice_x * (Rcell[0] + Ratom[0]) + dR[0];
      pos[1] = lattice_y * (Rcell[1] + Ratom[1]) + dR[1];
      pos[2] = lattice_z * (Rcell[2] + Ratom[2]) + dR[2];
    }
    /* End Added for ZB_110 */

    else {
      pos[0] = lattice_x * (Rcell[0] + Ratom[0]) + dR[0];
      pos[1] = lattice_y * (Rcell[1] + Ratom[1]) + dR[1];
      pos[2] = lattice_z * (Rcell[2] + Ratom[2]) + dR[2];
    }
  } else {
    if (isCellWurtzite()) {
      pos[0] = lattice_x * (0.5 * (Rcell[0] + Rcell[1]) + Ratom[0]);
      pos[1] = lattice_y * (0.5 * sqrt(3.0) * (Rcell[0] - Rcell[1]) + Ratom[1]);
      pos[2] = lattice_z * (Rcell[2] + Ratom[2]);
    }

    /* Added for ZB_110 */
    else if (isCellZB_110() || isCellZB_110_small() || isCellZB_111() ||
             isCellUserDef()) {
      pos[0] = lattice_x * (Rcell[0] + Ratom[0]);
      pos[1] = lattice_y * (Rcell[1] + Ratom[1]);
      pos[2] = lattice_z * (Rcell[2] + Ratom[2]);
    }
    /* End Added for ZB_110 */

    else {
      pos[0] = lattice_x * (Rcell[0] + Ratom[0]);
      pos[1] = lattice_y * (Rcell[1] + Ratom[1]);
      pos[2] = lattice_z * (Rcell[2] + Ratom[2]);
    }
  }
#else
#ifdef DISPL_SINGLE_PRECISION
  nml_float *dR = l_m_2_atomdsp[Zcell][Zatom];
#else  /* DISPL_SINGLE_PRECISION */
  real *dR = l_m_2_atomdsp[Zcell][Zatom];
#endif /* DISPL_SINGLE_PRECISION */

  if (isCellWurtzite()) {
    pos[0] = lattice_x * (0.5 * (Rcell[0] + Rcell[1]) + Ratom[0]) + dR[0];
    pos[1] = lattice_y * (0.5 * sqrt(3.0) * (Rcell[0] - Rcell[1]) + Ratom[1]) +
             dR[1];
    pos[2] = lattice_z * (Rcell[2] + Ratom[2]) + dR[2];
  }

  /* Added for ZB_110 */
  else if (isCellZB_110() || isCellZB_110_small() || isCellZB_111() ||
           isCellUserDef()) {
    pos[0] = lattice_x * (Rcell[0] + Ratom[0]) + dR[0];
    pos[1] = lattice_y * (Rcell[1] + Ratom[1]) + dR[1];
    pos[2] = lattice_z * (Rcell[2] + Ratom[2]) + dR[2];
  }
  /* End Added for ZB_110 */

  else {
    pos[0] = lattice_x * (Rcell[0] + Ratom[0]) + dR[0];
    pos[1] = lattice_y * (Rcell[1] + Ratom[1]) + dR[1];
    pos[2] = lattice_z * (Rcell[2] + Ratom[2]) + dR[2];
  }
#endif
}

void geo_struct::getPosition_equilibrium(real *pos, int Zcell, int Zatom) {

  int *Rcell = cell__ijk[Zcell];
  real *Ratom = PositionInCell(Zatom);
  if (isCellWurtzite()) {
    pos[0] = lattice_x * (0.5 * (Rcell[0] + Rcell[1]) + Ratom[0]);
    pos[1] = lattice_y * (0.5 * sqrt(3.0) * (Rcell[0] - Rcell[1]) + Ratom[1]);
    pos[2] = lattice_z * (Rcell[2] + Ratom[2]);
  }

  /* Added for ZB_110 */
  else if (isCellZB_110() || isCellZB_110_small() || isCellZB_111() ||
           isCellUserDef()) {
    pos[0] = lattice_x * (Rcell[0] + Ratom[0]);
    pos[1] = lattice_y * (Rcell[1] + Ratom[1]);
    pos[2] = lattice_z * (Rcell[2] + Ratom[2]);
  }
  /* End Added for ZB_110 */

  else {
    pos[0] = lattice_x * (Rcell[0] + Ratom[0]);
    pos[1] = lattice_y * (Rcell[1] + Ratom[1]);
    pos[2] = lattice_z * (Rcell[2] + Ratom[2]);
  }
}

void geo_struct::getRelativePosVec(real *nnv0, real *nnv, real unstrnd_latt,
                                   int c1, int a1, int c2, int a2,
                                   int neighbor) {
  if (cellgeom.IsCellWurtzite) {
    real normalized[3];
    normalized[0] =
        0.5 * (NbrCell(a2, neighbor, 0) + NbrCell(a2, neighbor, 1)) +
        PositionInCell(a1, 0) - PositionInCell(a2, 0);
    normalized[1] = 0.5 * sqrt(3.0) *
                        (NbrCell(a2, neighbor, 0) - NbrCell(a2, neighbor, 1)) +
                    PositionInCell(a1, 1) - PositionInCell(a2, 1);
    normalized[2] = NbrCell(a2, neighbor, 2) + PositionInCell(a1, 2) -
                    PositionInCell(a2, 2);
    nnv0[0] = unstrnd_latt * normalized[0];
    nnv0[1] = unstrnd_latt * normalized[1];
    nnv0[2] = unstrnd_latt * sqrt(8.0 / 3.0) * normalized[2];

    nnv[0] = lattice_x * normalized[0] + l_m_2_atomdsp[c1][a1][0] -
             l_m_2_atomdsp[c2][a2][0];
    nnv[1] = lattice_y * normalized[1] + l_m_2_atomdsp[c1][a1][1] -
             l_m_2_atomdsp[c2][a2][1];
    nnv[2] = lattice_z * normalized[2] + l_m_2_atomdsp[c1][a1][2] -
             l_m_2_atomdsp[c2][a2][2];
    //      cout<<"nnv0, NN distance: "<<nnv0[0]<<" "<<nnv0[1]<<" "<<nnv0[2]
    //          <<"
    //          "<<sqrt(nnv0[0]*nnv0[0]+nnv0[1]*nnv0[1]+nnv0[2]*nnv0[2])<<endl<<endl;
    //      cout<<"nnv, NN distance: "<<nnv[0]<<" "<<nnv[1]<<" "<<nnv[2]
    //          <<"
    //          "<<sqrt(nnv[0]*nnv[0]+nnv[1]*nnv[1]+nnv[2]*nnv[2])<<endl<<endl;
  }

  /* New ZB_110 unit cell */
  else if (cellgeom.IsCellZB_110 || cellgeom.IsCellZB_110_small ||
           cellgeom.IsCellZB_111 ||
           cellgeom.IsCellUserDef) { /* Added for ZB_110 */
    for (int p = 0; p < 3; p++) {
      real normalized;
      normalized = NbrCell(a2, neighbor, p) + PositionInCell(a1, p) -
                   PositionInCell(a2, p);
      if (p == 0) {
        nnv0[p] = get_a_latt_rot_norm(0) * unstrnd_latt * normalized;
        nnv[p] = lattice_x * normalized + l_m_2_atomdsp[c1][a1][p] -
                 l_m_2_atomdsp[c2][a2][p];
      } else if (p == 1) {
        nnv0[p] = get_a_latt_rot_norm(1) * unstrnd_latt * normalized;
        nnv[p] = lattice_y * normalized + l_m_2_atomdsp[c1][a1][p] -
                 l_m_2_atomdsp[c2][a2][p];
      } else if (p == 2) {
        nnv0[p] = get_a_latt_rot_norm(2) * unstrnd_latt * normalized;
        nnv[p] = lattice_z * normalized + l_m_2_atomdsp[c1][a1][p] -
                 l_m_2_atomdsp[c2][a2][p];
      }
    }
  } /* if(cellgeom.IsCellZB_110) */
  /* End New ZB_110 unit cell */

  else {
    for (int p = 0; p < 3; p++) {
      real normalized;
      normalized = NbrCell(a2, neighbor, p) + PositionInCell(a1, p) -
                   PositionInCell(a2, p);
      nnv0[p] = unstrnd_latt * normalized;
      if (p == 0)
        nnv[p] = lattice_x * normalized + l_m_2_atomdsp[c1][a1][p] -
                 l_m_2_atomdsp[c2][a2][p];
      else if (p == 1)
        nnv[p] = lattice_y * normalized + l_m_2_atomdsp[c1][a1][p] -
                 l_m_2_atomdsp[c2][a2][p];
      else if (p == 2)
        nnv[p] = lattice_z * normalized + l_m_2_atomdsp[c1][a1][p] -
                 l_m_2_atomdsp[c2][a2][p];
    }
  }
}
