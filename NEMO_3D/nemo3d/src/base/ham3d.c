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

#include "ham3d.h"
#include "Cell.h"
#include "eigsys3d_par.h"
#include "ham3d_final.h"

ostream &operator<<(ostream &o, const MaterialHandle &mh) {
  o << "matid3d=" << mh.matid3d << "\nanion_id=";
  for (int i = 0; i < (signed int)mh.anion_id.size(); i++)
    o << mh.anion_id[i] << "  " << endl;
  o << "\ncation_id=";
  for (int i = 0; i < (signed int)mh.cation_id.size(); i++)
    o << mh.cation_id[i] << "  ";
  o << "\nanion=" << mh.anion << "\ncation=" << mh.cation
    << "\ncation_weight=" << mh.cation_weight
    << "\nanion_weight=" << mh.anion_weight << "\n"
    << mh.hndl->toString(0) << endl;

  return o;
}

static bool isNotInBox(real *posAtom, real *dimensions) {
  real *rMin = dimensions;
  real *rMax = dimensions + 3;

  if (*posAtom < *rMin || *posAtom > *rMax)
    return true;
  posAtom++;
  rMin++;
  rMax++;
  if (*posAtom < *rMin || *posAtom > *rMax)
    return true;
  posAtom++;
  rMin++;
  rMax++;
  if (*posAtom < *rMin || *posAtom > *rMax)
    return true;

  return false;
}

/* Defines boundaries of a tilted box */
static bool isNotInTiltedBox(real *posAtom, real *dimensions, int **delta_step,
                             int shape_rank) {
  real *rMin = dimensions;
  real *rMax = dimensions + 3;
  double delta = 1e-6;
  double ymin, ymax;
  ymin = (((posAtom[0] - dimensions[0]) / dimensions[12]) * dimensions[13]) *
             (dimensions[7] / dimensions[6]) +
         dimensions[1] + delta;
  /* exclude all atoms on lower boundary */
  ymax = ymin + dimensions[9];
  /* include all atoms on upper boundary */

  double zmin, zmax;
  zmin = (((posAtom[0] - dimensions[0]) / dimensions[12]) * dimensions[14]) *
             (dimensions[8] / dimensions[6]) +
         dimensions[2] + delta;
  /* exclude all atoms on lower boundary */

  if (shape_rank == 2) { /* If this line is uncommented roughness is introduced
                            in Si QW only */
    /* adjust step edges for roughness */
    double stepLen = (dimensions[12] / 4) / (dimensions[8] / dimensions[6]);
    int delta_indx =
        int(posAtom[1] / (dimensions[13] / 4)); /* row index in delta_step */
    double Delta;

    int stepNo = int((posAtom[0] - delta) / stepLen);
    int stepNo2 = int((posAtom[0] - delta) / (2 * stepLen));
    Delta = delta_step[delta_indx - 1][stepNo2] * dimensions[12] / 4;

    if (posAtom[0] > (double(stepNo2)) * 2 * stepLen &&
        posAtom[0] < (double(stepNo2) + 1) * 2 * stepLen) {
      if (Delta > 0 && stepNo == (2 * stepNo2 + 1)) {
        if (posAtom[0] > (double(stepNo)) * stepLen &&
            posAtom[0] < (double(stepNo)) * stepLen + Delta) {
          zmin = (stepNo)*dimensions[14] / 4 + dimensions[2] - delta;
        }
      }
      if (Delta < 0 && stepNo == 2 * stepNo2) {
        if (posAtom[0] > (stepNo + 1) * stepLen + Delta &&
            posAtom[0] < (stepNo + 1) * stepLen) {
          zmin = (stepNo + 2) * dimensions[14] / 4 + dimensions[2] - delta;
        }
      }
    }
  } /* if shape_rank */

  zmax = zmin + dimensions[10];
  /* include all atoms on upper boundary */

  if (posAtom[0] < dimensions[0] || posAtom[0] > dimensions[3]) {
    return true;
  }

  if (posAtom[1] < ymin) {
    return true;
  }
  if (posAtom[1] > ymax) {
    return true;
  }

  if (posAtom[2] < zmin) {
    return true;
  }
  if (posAtom[2] > zmax) {
    return true;
  }

  return false;
} /* isNotInTiltedBox */

static bool isNotInRhombohedron(real *posAtom, real *dimensions)
// This rhombohedron is made of the unit cell of wurtzite structure
// The primitive vectors of the unit cell are
//(a/2, a*sqrt(3)/2, 0), (a/2,-a*sqrt(3)/2, 0), (0, 0, c)
{
  real *rMin = dimensions;
  real *rMax = dimensions + 3;
  double transform_pos[3];

  transform_pos[0] = posAtom[0] + posAtom[1] / sqrt(3.0);
  transform_pos[1] = posAtom[0] - posAtom[1] / sqrt(3.0);
  transform_pos[2] = posAtom[2];
  /*
     cout<<"transform_pos "
         <<transform_pos[0]<<" "<<transform_pos[1]<<" "<<transform_pos[2]<<endl;
     cout<<"rMin, rMax "<<rMin[0]<<" "<<rMin[1]<<" "<<rMin[2]<<" "
          <<rMax[0]<<" "<<rMax[1]<<" "<<rMax[2]<<endl;
  */
  if (transform_pos[0] < rMin[0] || transform_pos[0] > rMax[0])
    return true;
  if (transform_pos[1] < rMin[1] || transform_pos[1] > rMax[1])
    return true;
  if (transform_pos[2] < rMin[2] || transform_pos[2] > rMax[2])
    return true;
  return false;
}

static bool isNotInCylinder_x(real *posAtom, real *dimensions) {
  real *rMin = dimensions;
  real *rMax = dimensions + 3;
  real *rCenter = dimensions + 6;
  real Radius2 = dimensions[8];

  if (posAtom[0] < rMin[0] || posAtom[0] > rMax[0])
    return true;

  real sum = 0.0;
  sum += pow((posAtom[2] - rCenter[0]), 2.0);
  sum += pow((posAtom[1] - rCenter[1]), 2.0);

  if (double(sum - Radius2) > -1e-6)
    return true;
  return false;
}

static bool isNotInCylinder(real *posAtom, real *dimensions) {
  real *rMin = dimensions;
  real *rMax = dimensions + 3;
  real *rCenter = dimensions + 6;
  real Radius2 = dimensions[8];

  if (posAtom[2] < rMin[2] || posAtom[2] > rMax[2])
    return true;

  real sum = 0.0;
  sum += pow((posAtom[0] - rCenter[0]), 2.0);
  sum += pow((posAtom[1] - rCenter[1]), 2.0);
  if (sum > Radius2)
    return true;

  return false;
}

static bool isNotInPyramid(real *posAtom, real *dimensions) {
  real *rMin = dimensions;
  real *rMax = dimensions + 3;
  const real &dzdx = dimensions[6];
  const real &dzdy = dimensions[7];

  if (posAtom[2] < rMin[2])
    return true;
  if (posAtom[2] > rMin[2] + dzdx * (posAtom[0] - rMin[0]))
    return true;
  if (posAtom[2] > rMin[2] - dzdx * (posAtom[0] - rMax[0]))
    return true;
  if (posAtom[2] > rMin[2] + dzdy * (posAtom[1] - rMin[1]))
    return true;
  if (posAtom[2] > rMin[2] - dzdy * (posAtom[1] - rMax[1]))
    return true;

  return false;
}

static bool isNotInDome(real *posAtom, real *dimensions) {
  real *rMin = dimensions;
  real *rMax = dimensions + 3;
  real *rCenter = dimensions + 6;
  real *Radius_Square = dimensions + 9;

  if (posAtom[0] < rMin[0] || posAtom[0] > rMax[0])
    return true;
  if (posAtom[1] < rMin[1] || posAtom[1] > rMax[1])
    return true;
  if (posAtom[2] < rMin[2] || posAtom[2] > rMax[2])
    return true;

  real sum = 0.0;
  sum += pow((*posAtom) - rCenter[0], 2.0) / Radius_Square[0];
  posAtom++;
  sum += pow((*posAtom) - rCenter[1], 2.0) / Radius_Square[1];
  posAtom++;
  sum += pow((*posAtom) - rCenter[2], 2.0) / Radius_Square[2];
  if (sum > 1.0)
    return true;
  else
    return false;
}

static bool isNotInEllipsoid(real *posAtom, real *dimensions) {
  real *rMin = dimensions;
  real *rMax = dimensions + 3;
  real *rCenter = dimensions + 6;
  real *Radius_Square = dimensions + 9;

  if (posAtom[0] < rMin[0] || posAtom[0] > rMax[0])
    return true;
  if (posAtom[1] < rMin[1] || posAtom[1] > rMax[1])
    return true;
  if (posAtom[2] < rMin[2] || posAtom[2] > rMax[2])
    return true;

  real sum = 0.0;
  sum += pow((*posAtom) - rCenter[0], 2.0) / Radius_Square[0];
  posAtom++;
  sum += pow((*posAtom) - rCenter[1], 2.0) / Radius_Square[1];
  posAtom++;
  sum += pow((*posAtom) - rCenter[2], 2.0) / Radius_Square[2];
  if (sum > 1.0)
    return true;
  else
    return false;
}

static void PrintShapeInfo(const Shape_struct &shape, const MaterialHandle &mh,
                           int *species_count, real *rMax, real *rMin,
                           int *cMin, int *cMax) {

  if (mpi_n3d_id != mpi_n3d_masterid)
    return;
  string shapeHdr;
  switch (shape.Shape) {
  case Shape_struct::Box:
    shapeHdr = "[BOX] -------------------------";
    break;
  case Shape_struct::TiltedBox:
    shapeHdr = "[TILTEDBOX] -------------------------";
    break;
  case Shape_struct::Cylinder:
    shapeHdr = "[Cylinder] -------------------------";
    break;
  case Shape_struct::Cylinder_x:
    shapeHdr = "[Cylinder_x] -------------------------";
    break;
  case Shape_struct::Pyramid:
    shapeHdr = "[PYRAMID] -------------------------";
    break;
  case Shape_struct::Dome:
    shapeHdr = "[DOME] -------------------------";
    break;
  case Shape_struct::Ellipsoid:
    shapeHdr = "[Ellipsoid] -------------------------";
    break;
  case Shape_struct::Rhombohedron:
    shapeHdr = "[RHOMBOHEDRON] -------------------------";
    break;
  default:
    shapeHdr = "   UNKNOWN SHAPE  ";
    break;
  }

  int NatomsInShape = 0;
  for (int i = 0; i < EL_TOTAL; i++) {
    if (species_count[i] == 0)
      continue;
    NatomsInShape += species_count[i];
  }

  cout << "\n   " << shapeHdr << endl;
  cout << "      material:  " << mh.cation << mh.anion << endl;
  if (shape.has_cell_granularity)
    cout << "      will have cell granularity " << endl;
  else
    cout << "      will not have cell granularity " << endl;
  cout << "      will use this shape for strain "
       << (shape.use_for_strain_and_electronic ? "and electronic calculation"
                                               : "calculation only")
       << endl;
  cout << "      requested dimensions:  " << shape.x << "   " << shape.y
       << "   " << shape.z << "   " << endl;

  if (NatomsInShape == 0) {
    cout << "      !!! THIS SHAPE AND FIRST SHAPE HAVE NO ATOMS IN COMMON !!!"
         << endl;
  } else {
    cout << "      adjusted dimensions:  " << rMax[0] - rMin[0] << "   "
         << rMax[1] - rMin[1] << "   " << rMax[2] - rMin[2] << "   " << endl;
    cout << "      cell index boundaries:  {"
         << "  [" << cMin[0] << "," << cMax[0] + 1 << ");"
         << "  [" << cMin[1] << "," << cMax[1] + 1 << ");"
         << "  [" << cMin[2] << "," << cMax[2] + 1 << ")" << " }" << endl;

    for (int i = 0; i < EL_TOTAL; i++) {
      if (species_count[i] == 0)
        continue;
      printf("         * number of %s atoms = %d\n", atomID_to_str(i),
             species_count[i]);
    }
    cout << "         * total number of atoms = " << NatomsInShape << endl;
  }
  if (shape.Shape == Shape_struct::Rhombohedron) {
    cout << "      important note: the dimensions of rhombohedron are defined "
            "as "
         << endl;
    cout << "                      the lengths along the wurtzite primitive "
            "vectors "
         << endl;
    cout << "                      A1=(a/2,  a*sqrt(3)/2, 0) " << endl;
    cout << "                      A2=(a/2, -a*sqrt(3)/2, 0) " << endl;
    cout << "                      A3=(  0,            0, c) " << endl;
  }
}

static int shape_rank = 0;
void geo_struct::FillShape(const Shape_struct &shape,
                           const vector<MaterialHandle> &sMatList,
                           bool has_cell_granularity) {

  /* The new data storage can crash on the dereferencing of the casted data */
  /* It appears to crash on the new gcc compilers                   */
  /* const int& NatomsPerCell = this->AtomsPerCellMax(); */
  int NatomsPerCell = this->AtomsPerCellMax();
  const int useOnlyForStrain = shape.use_for_strain_and_electronic ? 0 : 1;
  const MaterialHandle &mh = sMatList[shape.matid3d];
  bool isNotRandomAlloy = (mh.cation_id.size() == 1 && mh.anion_id.size() == 1);

  int species_count[EL_TOTAL];
  for (int i = 0; i < EL_TOTAL; i++)
    species_count[i] = 0;

  // nml_memory_report("FillShape beginning");
#ifdef MEMORY_SAVE_NEWSTRUCT
  bmatrix inSet = Bmatrix(this->N_Cell, NatomsPerCell);
#else  /* MEMORY_SAVE_NEWSTRUCT  */
  imatrix inSet = Imatrix(this->N_Cell, NatomsPerCell);
#endif /* MEMORY_SAVE_NEWSTRUCT  */

  // nml_memory_report("FillShape after allocation of inSet");

#undef MEMCHECK
#ifdef MEMCHECK
  double int3_memory = 0.0, tl_memory = 0.0;
  ;
  int3_memory = this->N_Cell * NatomsPerCell * 4 / 1.0e6;
  tl_memory += int3_memory;
  printf("\n FilleShape memory allocation in MB int3_memory=%g tl_memory=%g\n",
         int3_memory, tl_memory),
      fflush(stdout);
#endif

  bool (*isNotInShape)(real *, real *) = NULL;
  real shDim[15];
  shDim[0] = shape.xorigin;
  shDim[1] = shape.yorigin;
  shDim[2] = shape.zorigin;
  shDim[3] = shape.xorigin + shape.x;
  shDim[4] = shape.yorigin + shape.y;
  shDim[5] = shape.zorigin + shape.z;
  switch (shape.Shape) {
  case Shape_struct::Box:
    isNotInShape = isNotInBox;
    break;
  case Shape_struct::TiltedBox: /* shDim for tilted Box */
                                //      isNotInShape = isNotInTiltedBox;
    this->calc_tilt_params();
    /* Adjust origin and dimensions of the shape and make them integer multiples
     * of ax, ay, az */
    shDim[0] = ((int)(shape.xorigin / this->lattice_x)) * this->lattice_x;
    shDim[1] = ((int)(shape.yorigin / this->lattice_y)) * this->lattice_y;
    shDim[2] = ((int)(shape.zorigin / this->lattice_z)) * this->lattice_z;
    shDim[3] = shDim[0] + ((int)(shape.x / this->lattice_x)) * this->lattice_x;
    shDim[4] = shDim[1] + ((int)(shape.y / this->lattice_y)) * this->lattice_y;
    shDim[5] = shDim[2] + ((int)(shape.z / this->lattice_z)) * this->lattice_z;
    shDim[6] = this->dx;                                        // tiltdir[0];
    shDim[7] = this->dy;                                        // tiltdir[1];
    shDim[8] = this->dz;                                        // tiltdir[2];
    shDim[9] = (shDim[4] - shDim[1]) /*- this->adjThicknessY*/; /* thickInY */
    shDim[10] = (shDim[5] - shDim[2]) - this->adjThicknessZ;    /* thickInZ */
    //     double tiltedUnitCellX  = this->dx*this->lattice_x;
    shDim[11] =
        ((int)((shDim[3] - shDim[0]) / (this->dx * this->lattice_x))) *
        (this->dx * this->lattice_x); /* Adjusted x-dimension. It is made
                                         commensurate with tilt direction */
    shDim[12] = this->lattice_x;
    shDim[13] = this->lattice_y;
    shDim[14] = this->lattice_z;
    break;
  case Shape_struct::Rhombohedron:
    isNotInShape = isNotInRhombohedron;
    break;
  case Shape_struct::Cylinder_x:
    isNotInShape = isNotInCylinder_x;
    shDim[6] = (shDim[1] + shDim[4]) * 0.5;
    shDim[7] = (shDim[2] + shDim[5]) * 0.5;
    shDim[8] = shape.y * shape.y * 0.25;
    break;
  case Shape_struct::Cylinder:
    isNotInShape = isNotInCylinder;
    shDim[6] = (shDim[0] + shDim[3]) * 0.5;
    shDim[7] = (shDim[1] + shDim[4]) * 0.5;
    shDim[8] = shape.x * shape.x * 0.25;
    break;
  case Shape_struct::Pyramid:
    isNotInShape = isNotInPyramid;
    shDim[6] = 2.0 * shape.z / shape.x; //  dz/dx
    shDim[7] = 2.0 * shape.z / shape.y; //  dz/dy
    break;
  case Shape_struct::Dome:
    shDim[6] = shape.xorigin + shape.x / 2.0;
    shDim[7] = shape.yorigin + shape.y / 2.0;
    shDim[8] = shape.zorigin;
    shDim[9] = 0.25 * shape.x * shape.x;
    shDim[10] = 0.25 * shape.y * shape.y;
    shDim[11] = shape.z * shape.z;
    isNotInShape = isNotInDome;
    break;
  case Shape_struct::Ellipsoid:
    shDim[6] = shape.xorigin + shape.x / 2.0;
    shDim[7] = shape.yorigin + shape.y / 2.0;
    shDim[8] = shape.zorigin + shape.z / 2.0;
    shDim[9] = 0.25 * shape.x * shape.x;
    shDim[10] = 0.25 * shape.y * shape.y;
    shDim[11] = 0.25 * shape.z * shape.z;
    isNotInShape = isNotInEllipsoid;
    break;
  default:
    die("unimplemented shape");
    break;
    break;
  }

  // loop over atoms; if an atom is in the shape mark it.
  real rMin[3], rMax[3];
  int cMin[3], cMax[3];

  rMin[0] = rMin[1] = rMin[2] = 9999;
  rMax[0] = rMax[1] = rMax[2] = -9999;
  cMin[0] = cMin[1] = cMin[2] = 9999;
  cMax[0] = cMax[1] = cMax[2] = -9999;

  real pos[3];
  int NatomsInShape = 0;
  for (int Zc = 0; Zc < this->N_Cell; Zc++) {
    bool invalidate_this_cell = false;
    int NatomsInCell = 0;
    for (int Za = 0; Za < NatomsPerCell; Za++) {
      this->getPosition_strained(pos, Zc, Za);

      /* If lattice is tilted call isNotInTiltedBox function */
      if (this->is_tilted && shape.Shape == Shape_struct::TiltedBox) {
        if (isNotInTiltedBox(pos, shDim, this->delta_step, shape_rank)) {
          if (has_cell_granularity)
            invalidate_this_cell = true;
          continue;
        }
      } /* if is_tilted */
      else { // !is_tilted
        if (isNotInShape(pos, shDim)) {
          if (has_cell_granularity)
            invalidate_this_cell = true;
          continue;
        }
      } /* else is_tilted. */

      if (shape.Shape == Shape_struct::Rhombohedron) {
        double transform_pos[3];
        transform_pos[0] = pos[0] + pos[1] / sqrt(3.0);
        transform_pos[1] = pos[0] - pos[1] / sqrt(3.0);
        transform_pos[2] = pos[2];
        if (rMin[0] > transform_pos[0])
          rMin[0] = transform_pos[0];
        if (rMax[0] < transform_pos[0])
          rMax[0] = transform_pos[0];
        if (rMin[1] > transform_pos[1])
          rMin[1] = transform_pos[1];
        if (rMax[1] < transform_pos[1])
          rMax[1] = transform_pos[1];
        if (rMin[2] > transform_pos[2])
          rMin[2] = transform_pos[2];
        if (rMax[2] < transform_pos[2])
          rMax[2] = transform_pos[2];
      } else {
        if (rMin[0] > pos[0])
          rMin[0] = pos[0];
        if (rMax[0] < pos[0])
          rMax[0] = pos[0];
        if (rMin[1] > pos[1])
          rMin[1] = pos[1];
        if (rMax[1] < pos[1])
          rMax[1] = pos[1];
        if (rMin[2] > pos[2])
          rMin[2] = pos[2];
        if (rMax[2] < pos[2])
          rMax[2] = pos[2];
      }
      if (cMin[0] > this->cell__ijk[Zc][0])
        cMin[0] = this->cell__ijk[Zc][0];
      if (cMax[0] < this->cell__ijk[Zc][0])
        cMax[0] = this->cell__ijk[Zc][0];
      if (cMin[1] > this->cell__ijk[Zc][1])
        cMin[1] = this->cell__ijk[Zc][1];
      if (cMax[1] < this->cell__ijk[Zc][1])
        cMax[1] = this->cell__ijk[Zc][1];
      if (cMin[2] > this->cell__ijk[Zc][2])
        cMin[2] = this->cell__ijk[Zc][2];
      if (cMax[2] < this->cell__ijk[Zc][2])
        cMax[2] = this->cell__ijk[Zc][2];

      inSet[Zc][Za] = 1;
      NatomsInCell++;
    }
    if (invalidate_this_cell && !Shape_struct::Cylinder_x)
      for (int Za = 0; Za < NatomsPerCell; Za++)
        inSet[Zc][Za] = 0;
    else
      NatomsInShape += NatomsInCell;
    /* If lattice is tilted tilt then unit cells have to be cut */
    if (this->is_tilted)
      NatomsInShape += NatomsInCell;
    else { /* !is_tilted */
      if (invalidate_this_cell)
        for (int Za = 0; Za < NatomsPerCell; Za++)
          inSet[Zc][Za] = 0;
      else
        NatomsInShape += NatomsInCell;
    } /* else is_tilted */
  }

  for (int Zc = 0; Zc < this->N_Cell; Zc++) {
    for (int Za = 0; Za < NatomsPerCell; Za++) {
      if (inSet[Zc][Za])
        this->ShapeType[Zc][Za] = shape_rank;
    }
  }

  // now we populate the marked atoms with particular atomic speicies
  const vector<int> &cation = NatomsPerCell != 1 ? mh.cation_id : mh.anion_id;
  const vector<int> &anion = mh.anion_id;
  real prob_cation = (cation.size() == 2 ? mh.cation_weight : 1.0);
  real prob_anion = (anion.size() == 2 ? mh.anion_weight : 1.0);

  if (shape.AlloyFixConcentration) {
    if (shape.AlloyTreatment != Shape_struct::RandomInUnitCell ||
        isNotRandomAlloy)
      die("Cannot only use fixed alloy concentration with a random alloy and "
          "atomic granularity");

    real conc = prob_cation;
    int Nc = int(0.5 * conc * this->N_Cell * this->AtomsPerCellMax());
    ivectr Ga = Ivectr(Nc);

    printf("%d %e\n", Nc, prob_cation);

    /* construct array Ga containing atom tags and overwrite */
    for (int k = 0; k < Nc; k++) {
      real rnd = rand() / real(RAND_MAX);
      int Zrnd = (int)(this->N_Cell * this->AtomsPerCellMax() * rnd);
      int cell = Zrnd / this->AtomsPerCellMax();
      int m = Zrnd % this->AtomsPerCellMax();

      if (m > 3 && shape.AlloyClusterRadius % 2) {
        /* CASE: odd clusterlevel => select a cation (throw out anions) */
        k--;
        continue;
      } else if (m < 4 && shape.AlloyClusterRadius % 2 == 0) {
        /* CASE: even clusterlevel => select an anion (throw out cations) */
        k--;
        continue;
      }

      /* cluster size of a single atom */
      if (shape.AlloyClusterRadius == 1) {
        int j;
        /* search for duplicate */
        for (j = 0; j < k && Zrnd != Ga[j]; j++)
          ;

        if (j != k || cation.size() <= 1)
          k--; /* throw out duplicate */
        else {
          /* no duplicate */
          Ga[k] = Zrnd;

          species_count[cation[0]]++;
          species_count[cation[1]]--;

#ifdef MEMORY_SAVE_NEWSTRUCT
          this->AtomType[cell][m] = (unsigned char)cation[0];
#else  /* MEMORY_SAVE_NEWSTRUCT  */
          this->AtomType[cell][m] = cation[0];
#endif /* MEMORY_SAVE_NEWSTRUCT  */
        }
      }

      /* cluster size of a single primitive cell */
      if (shape.AlloyClusterRadius == 2) {
        int n;
        // int Zx = this->cell__ijk[cell][0];
        // int Zy = this->cell__ijk[cell][1];
        // int Zz = this->cell__ijk[cell][2];
        for (n = 0; n < 4; n++) {
          // int Zx_nbr = Zx + this->NbrCell(m,n,0);
          // int Zy_nbr = Zx + this->NbrCell(m,n,1);
          // int Zz_nbr = Zx + this->NbrCell(m,n,2);
          // int m_nbr = this->NbrCell(m,n,3);
        }
      }
    }
    rm_ivectr(&Ga);
  } else {
    //      reset_random_seed_with_time();
    for (int Zc = 0; Zc < this->N_Cell; Zc++) {
      real rnd = real(rand()) / real(RAND_MAX);
      int Za = 0;
      // unused   int count_atom=0;
      for (Za = 0; Za < NatomsPerCell / 2; Za++) {
        if (!inSet[Zc][Za])
          continue;

        if (shape.AlloyTreatment == Shape_struct::RandomInUnitCell)
          rnd = real(rand()) / real(RAND_MAX);

#ifdef MEMORY_SAVE_NEWSTRUCT
        this->AtomType[Zc][Za] =
            (unsigned char)cation[rnd <= prob_cation ? 0 : 1];
        this->atom_strain_only[Zc][Za] = (unsigned char)useOnlyForStrain;
        species_count[(int)this->AtomType[Zc][Za]]++;
#else  /* MEMORY_SAVE_NEWSTRUCT  */
        this->AtomType[Zc][Za] = cation[rnd <= prob_cation ? 0 : 1];
        this->atom_strain_only[Zc][Za] = useOnlyForStrain;
        species_count[this->AtomType[Zc][Za]]++;
#endif /* MEMORY_SAVE_NEWSTRUCT  */
      }

      for (; Za < NatomsPerCell; Za++) {
        if (!inSet[Zc][Za])
          continue;

        if (shape.AlloyTreatment == Shape_struct::RandomInUnitCell)
          rnd = real(rand()) / real(RAND_MAX);

#ifdef MEMORY_SAVE_NEWSTRUCT
        this->AtomType[Zc][Za] =
            (unsigned char)anion[rnd <= prob_anion ? 0 : 1];
        this->atom_strain_only[Zc][Za] = (unsigned char)useOnlyForStrain;
        species_count[(int)this->AtomType[Zc][Za]]++;
#else  /* MEMORY_SAVE_NEWSTRUCT  */
        this->AtomType[Zc][Za] = anion[rnd <= prob_anion ? 0 : 1];
        this->atom_strain_only[Zc][Za] = useOnlyForStrain;
        species_count[this->AtomType[Zc][Za]]++;
#endif /* MEMORY_SAVE_NEWSTRUCT  */
      }
    }
  }

  PrintShapeInfo(shape, mh, species_count, rMax, rMin, cMin, cMax);

  // cout << mh << endl;
  // rm_ivectr(&inSetAtom);
  // rm_ivectr(&inSetCell);
#ifdef MEMORY_SAVE_NEWSTRUCT
  rm_bmatrix(&inSet);
#else  /* MEMORY_SAVE_NEWSTRUCT  */
  rm_imatrix(&inSet);
#endif /* MEMORY_SAVE_NEWSTRUCT  */

  shape_rank++;
#undef MEMCHECK
#ifdef MEMCHECK
  {
    double int2_memory = 0.0, tl_memory = 0.0;
    ;
    int3_memory = this->N_Cell * NatomsPerCell * 4 / 1.0e6;
    tl_memory -= int3_memory;
    printf("\n FilleShape memory deallocationallocation in MB int3_memory=%g "
           "tl_memory=%g\n",
           int3_memory, tl_memory),
        fflush(stdout);
  }
#endif
  // nml_memory_report("FillShape at the end");
}

/* The first shape MUST be a box for cubic/zincblende structure
   or rhombohedron for wurtzite crystal structure
with minimal granularity on the cellular level */
void geo_struct::FillFirstShape(const Shape_struct &shape,
                                const top_struct &opt,
                                const vector<MaterialHandle> &sMatList) {
  int Zcell = 0;

  Cell cell;
  cell.Initialize(opt);

  // pad domain with a shell one cell thick
  int Zx0 = 1 + (int)(shape.xorigin / this->lattice_x);
  int Zy0 = 1 + (int)(shape.yorigin / this->lattice_y);
  int Zz0 = 1 + (int)(shape.zorigin / this->lattice_z);

  int Zxf = Zx0 + (int)(shape.x / this->lattice_x);
  int Zyf = Zy0 + (int)(shape.y / this->lattice_y);
  int Zzf = Zz0 + (int)(shape.z / this->lattice_z);

  int Nx = 2 + (int)(shape.x / this->lattice_x);
  int Ny = 2 + (int)(shape.y / this->lattice_y);
  int Nz = 2 + (int)(shape.z / this->lattice_z);

  if (Zxf <= Zx0 || Zyf <= Zy0 || Zzf <= Zz0)
    die("Domain shape has zero volume!");

  // Call allocate for the very first time!
  this->Allocate(Nx * Ny * Nz, Nx, Ny, Nz, 1, cell);

  if (opt.Dev.CrystalStruct == Dev_struct::Diamond_FCC ||
      opt.Dev.CrystalStruct == Dev_struct::Zincblende_FCC) {

    /* This branch is for the (currently broken) case of FCC symmetry
       (i.e a rhombohedron).  To build an FCC lattice we only define atoms
       at every other site of the cubic lattice in the z-direction.
       The phase of the on/off pattern is shifted by Pi for consecutive
       z-loops.  For BCC one can construct a similar algorithm. */

    const string bullet = "      * ";
    int sign, evenodd;
    /* Initialize Hamiltonian Data Structures */
    Zcell = 0;
    sign = -1;
    for (int i = 1; i < Nx - 1; i++) {
      // real x = i*this->lattice_x;
      for (int j = 1; j < Ny - 1; j++) {
        // real y = j*this->lattice_y;
        sign *= -1;
        for (int k = 1; k < Nz - 1; k++) {
          // real z = k*this->lattice_z;
          evenodd = (sign > 0 ? 1 : 0);
          /* Grid point inside domain */
          if (i >= Zx0 && i < Zxf && j >= Zy0 && j < Zyf && k >= Zz0 &&
              k < Zzf && k % 2 == evenodd) {
            this->cell__ijk[Zcell][0] = i;
            this->cell__ijk[Zcell][1] = j;
            this->cell__ijk[Zcell][2] = k;
            if (i < this->cell_xmin)
              this->cell_xmin = i;
            if (i > this->cell_xmax)
              this->cell_xmax = i;
            if (j < this->cell_ymin)
              this->cell_ymin = j;
            if (j > this->cell_ymax)
              this->cell_ymax = j;
            if (k < this->cell_zmin)
              this->cell_zmin = k;
            if (k > this->cell_zmax)
              this->cell_zmax = k;

            // cell_atom_fill(shape, di, Zcell, matinit);

            Zcell++;
          }
        }
      }
    }

    this->N_Cell = Zcell;

    if (mpi_n3d_id == mpi_n3d_masterid) {
      cout << "   [BOX (domain)] --------------------" << endl
           << bullet << "requested dimensions:  " << shape.x << "   " << shape.y
           << "   " << shape.z << "   " << endl
           << bullet
           << "adjusted dimensions:  " << (Zxf - Zx0) * this->lattice_x << "   "
           << (Zyf - Zy0) * this->lattice_y << "   "
           << (Zzf - Zz0) * this->lattice_z << "   " << endl
           << bullet << "Cell index boundaries:  {"
           << "  [" << Zx0 << "," << Zxf << ");"
           << "  [" << Zy0 << "," << Zyf << ");"
           << "  [" << Zz0 << "," << Zzf << ")" << " }" << endl;
    }
  } else {
    for (int i = Zx0; i < Zxf; i++) {
      for (int j = Zy0; j < Zyf; j++) {
        for (int k = Zz0; k < Zzf; k++) {
          this->cell__ijk[Zcell][0] = i;
          this->cell__ijk[Zcell][1] = j;
          this->cell__ijk[Zcell][2] = k;
          Zcell++;
        }
      }
    }

    this->cell_xmin = Zx0;
    this->cell_ymin = Zy0;
    this->cell_zmin = Zz0;
    this->cell_xmax = Zxf - 1;
    this->cell_ymax = Zyf - 1;
    this->cell_zmax = Zzf - 1;
    this->N_Cell = (Zxf - Zx0) * (Zyf - Zy0) * (Zzf - Zz0);
    this->FillShape(shape, sMatList, true);
  }
}
