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
$Header: /repo/nemo3d/src/base/geo_struct.h,v 1.36 2008/07/23 16:47:30 nkharche
Exp $
*****************************************************************************/

#ifndef GEO_STRUCT_H
#define GEO_STRUCT_H

#define MEMORY_SAVE_NEWSTRUCT
#define DISPL_SINGLE_PRECISION
#define ELIMINATE_SSMAP // Comment this line for tilting the crystal.

#include "b2tensor.h"
#include "constants_nemo.h"
#include "f3tensor.h"
#include "i2tensor.h"
#include "i3matrix.h"
#include "i3tensor.h"
#include "nml_global.h"
#include "r3tensor.h"
#include "s2tensor.h"
#include <ivector.h>
#include <rmatrix.h> /* Added for ZB_110 and ZB_111 */
#include <rvector.h> /* Added for ZB_110 and ZB_111 */

#include "Cell.h"
#include "MaterialHandle.h"
#include "mat_def.h"
#include "top_struct.H"

/*! \class geo_struct
 *  \brief describes geometry of device.

 *  contains all information necessary to describe physical properties
 *  of the device.
 */

class geo_struct {
  int periodicity;
  /*!< This is a bit (boolean) vector containing periodicity info.
   * What follows is a description of what each bit means:
   * - bit 0 determines whether the structure has periodicity along x
   * - bit 1 determines whether the structure has periodicity along y
   * - bit 2 determines whether the structure has periodicity along z
   * - bit 3 determines whether the strain calculation has periodicity along x
   * - bit 4 determines whether the strain calculation has periodicity along y
   * - bit 5 determines whether the strain calculation has periodicity along z
   * - bit 6 determines whether the elect. calculation has periodicity along x
   * - bit 7 determines whether the elect. calculation has periodicity along y
   * - bit 8 determines whether the elect. calculation has periodicity along z.
   */

  int unitcelldivide;
  /*!< Seems that # of so-called elementary crystal structure (i.e : fcc)
   * that's needed for construction of a single unit cell. The value = 2 when
   * the structure is zincblend or diamond fcc and 1 otherwise.*/

#ifdef MEMORY_SAVE_NEWSTRUCT
  b2tensor surfaceInfo;
  /*!< @see i2tensor surfaceInfo*/
#else  /* MEMORY_SAVE_NEWSTRUCT  */
  i2tensor surfaceInfo;
  /*!< This is a bit (boolean) vector containing surface info.
   * What follows is a description of what each bit means:
   *  - bit 0 determines whether or not an atom is on the surface (with bc)
   *  - bit 1 determines whether or not an atom is on the surface (sans bc)
   *  - bit 2 determines whether the x-position of the cell we're in is a min
   *  - bit 3 determines whether the x-position of the cell we're in is a max
   *  - bit 4 determines whether the y-position of the cell we're in is a min
   *  - bit 5 determines whether the y-position of the cell we're in is a max
   *  - bit 6 determines whether the z-position of the cell we're in is a min
   *  - bit 7 determines whether the z-position of the cell we're in is a max.
   *  - If MEMORY_SAVE_NEWSTRUCT is defined, boolean type is used.*/
#endif /* MEMORY_SAVE_NEWSTRUCT  */

  void set_geo_struct(
      int periodicity, int N_Cell, int cell_xmin, int cell_ymin, int cell_zmin,
      int cell_xmax, int cell_ymax, int cell_zmax,
      /* Additional variables to tilt the lattice */
      //                     ivectr tiltdir,
      int dx, int dy, int dz, int is_tilted,
      //                       int NaddALY,
      int NaddALZ,
      //                     double adjThicknessY,
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
#endif /* MEMORY_SAVE_NEWSTRUCT */
      i2tensor l_m_2_offset_ham_tot,
#ifdef DISPL_SINGLE_PRECISION
      f3tensor l_m_2_atomdsp
#else  /* DISPL_SINGLE_PRECISION */
      r3tensor l_m_2_atomdsp
#endif /* DISPL_SINGLE_PRECISION */
#ifndef ELIMINATE_SSMAP
      ,
      i3matrix ssmap /* Note that the comma before i3matrix is not a bug! */
#endif               /* ELIMINATE_SSMAP */
  );
  /*!< Private helper function for copy constructor and assignment operator */

public:
  int N_Cell; /*!< Number of cells in device. */
  /* MEMCHECKmods   */
  int Nx_Cell; /*!< Number of cells in x-direction.  */
  int Ny_Cell; /*!< Number of cells in y-direction.  */
  int Nz_Cell; /*!< Number of cells in z-direction.  */
  /* MEMCHECKmods   */

  int cell_xmin; /*!< Minimum value of index in x direction */
  int cell_ymin; /*!< Minimum value of index in y direction */
  int cell_zmin; /*!< Minimum value of index in z direction */
  int cell_xmax; /*!< Maximum value of index in x direction */
  int cell_ymax; /*!< Maximum value of index in y direction */
  int cell_zmax; /*!< Maximum value of index in z direction */
                 /* Additional variables to tilt the lattice */
                 //                       ivectr tiltdir;
  int dx;
  int dy;
  int dz;
  int is_tilted;
  //                       int NaddALY;
  int NaddALZ;
  //                       double adjThicknessY;
  double adjThicknessZ;
  i2tensor delta_step;
  //                       double adjustedLenX;
  /* End: Additional variables to tilt the lattice */
  real lattice_x; /*!< Lattice constant in x direction for a cubic cell */
  real lattice_y; /*!< Lattice constant in y direction for a cubic cell */
  real lattice_z; /*!< Lattice constant in z direction for a cubic cell */

  Cell cellgeom;
  /*!< Cell information.
   *  @see Cell::Cell()
   */
  i3tensor ijk__cell;
  /*!< Integer 3D array. ijk_cell[i][j][k] returns the index of cell whose
   * relative location of unitcell is (i,j,k) */
  i2tensor cell__ijk;
  /*!< Integer 2D array. If index of a cell is indx, cell_ijk[indx][0] returns
   * i, cell[indx][1] returns j and cell[indx][2] returns k. */

#ifdef MEMORY_SAVE_NEWSTRUCT
  b2tensor AtomType;         /*!< see i2tensor AtomType */
  s2tensor ShapeType;        /*!< see i2tensor ShapeType */
  b2tensor isDiffNbr;        /*!< see i2tensor isDiffNbr */
  b2tensor atom_strain_only; /*!< @see i2tensor atom_strain_only */
#else                        /* MEMORY_SAVE_NEWSTRUCT  */
  i2tensor AtomType;
  /*!< 2D integer array. Atoptype[Cell index][Atom index] = 1 if atom is filled
   * in. See mat_def.h for more information. If MEMORY_SAVE_NEWSTRUCT is
   * defined, boolean type is used. */
  i2tensor ShapeType;
  /*!< No more use in Nemo3D. ShapeType[Cell index][Atom index] = 1 if atom is
   * filled in This was originally used to determine shape of structure.
   * Shapetype start from box shape. Then, according to the defined shape,
   * NEMO3D fills 1 in part of this array where should be filled by atom. If
   * MEMORY_SAVE_NEWSTRUCT is defined, boolean type is used.
   */

  i2tensor isDiffNbr;
  /*!< 2D interger array.
   * - isDiffNbr[cindx][aindx] = 1
   *   - the atom is surrounded by atoms of different material. Parameters for
   * Hamiltonian setup should be averaged.
   * - isDiffNbr[cindx][aindx] = 0
   *   - the atom is surrounded by atoms of same material. Parameters for
   * Hamiltonian doens't need to be averaged here.
   * - If MEMORY_SAVE_NEWSTRUCT is defined, boolean type is used. */

  i2tensor atom_strain_only;
  /*!< 2D interger array.
   * - atom_strain_only[cindx][aindx] = 1
   *   - the atom is in the region of so called 'Strain calculation only'
   * region.
   * - atom_strain_only[cindx][aindx] = 0
   *   - atom is not only in 'Strain calculation' region, but also in
   * 'Electronic calculation' region.
   * - If MEMORY_SAVE_NEWSTRUCT is defined, boolean type is used. */

#endif /* MEMORY_SAVE_NEWSTRUCT  */
  i2tensor l_m_2_offset_ham_tot;
  /*!< Offset index for atom whose (cell index = l, atomic index = m).
   *   Wave function array in NEMO3D usually has extremely large size.
   *   (i.e : number of basis * number of atoms in unit cell * number of
   * unitcell) Therefore, when we want to access wave function array using a
   * single index (l,m), it could be extremely difficult to calculate
   * correspondent index. l_m_2_offset_ham_tot[cindx][aindx] has offset index of
   * the atom with cell index of cindx and atomic index of aindex, therefore we
   * can access wave function array by using simple index such as
   * (l_m_2_offset_ham_tot[l][m] + i), where i is number of i-th basis. */

#ifdef DISPL_SINGLE_PRECISION
  f3tensor l_m_2_atomdsp;
  /*!< 3D float array.
   * displacement info of atomic positions caused from strain calculation.
   * - l_m_2_atomdsp[cidx][aidx][d]
   *   - displacement in x-axis when d=1
   *   - displacement in y-axis when d=2
   *   - displacement in z-axis when d=3
   * - If DISPL_SINGLE_PRESISTION isn't defined, real type array will be used.
   */

  f3tensor l_m_2_atomdsp_onedomain;
  /*!< Displacement info in one domain. Will not be used if
   * DISPL_SINGLE_PRESISTION isn't defined. */
#else  /* DISPL_SINGLE_PRECISION */
  r3tensor l_m_2_atomdsp;
  /*!< See f3tensor l_m_2_atomdsp */
#endif /* DISPL_SINGLE_PRECISION */

#ifndef ELIMINATE_SSMAP
  i3matrix ssmap;
  /*!< 3D integer array.
   * ssmap[c][a][n] have index of a unit cell where n-th neighbor atom of a-th
   * atom in this cell belongs. */

#endif /* ELIMINATE_SSMAP */

  geo_struct()
      : periodicity(0), unitcelldivide(0), surfaceInfo(NULL), N_Cell(0),
        cell_xmin(0), cell_ymin(0), cell_zmin(0), cell_xmax(0), cell_ymax(0),
        cell_zmax(0),
        /* Additional variables to tilt the lattice */
        /*tiltdir(NULL),*/ dx(0), dy(0), dz(0), is_tilted(0),
        /*NaddALY(0),*/ NaddALZ(0), /* adjThicknessY(0),*/ adjThicknessZ(0),
        delta_step(NULL), /*adjustedLenX(0),*/
        /* End: Additional variables to tilt the lattice */
        lattice_x(0), lattice_y(0), lattice_z(0), ijk__cell(NULL),
        cell__ijk(NULL), AtomType(NULL), ShapeType(NULL), isDiffNbr(NULL),
        atom_strain_only(NULL), l_m_2_offset_ham_tot(NULL), l_m_2_atomdsp(NULL)
#ifndef ELIMINATE_SSMAP
        ,
        ssmap(NULL) /* Note that the comma before ssmap is NOT a bug! */
#endif              /* ELIMINATE_SSMAP */
  {
    ;
  }
  /*!< default constructor; all data get initialized to 0 or NULL. */

  geo_struct(const geo_struct &g) {
    set_geo_struct(
        g.periodicity, g.N_Cell, g.cell_xmin, g.cell_ymin, g.cell_zmin,
        g.cell_xmax, g.cell_ymax, g.cell_zmax,
        /* Additional variables to tilt the lattice */
        /*g.tiltdir,*/ g.dx, g.dy, g.dz, g.is_tilted, /* g.NaddALY,*/ g.NaddALZ,
        /* g.adjThicknessY,*/ g.adjThicknessZ, g.delta_step, /*g.adjustedLenX,*/
        /* End: Additional variables to tilt the lattice */
        g.lattice_x, g.lattice_y, g.lattice_z, g.unitcelldivide, g.cellgeom,
        g.ijk__cell, g.cell__ijk, g.AtomType, g.ShapeType, g.isDiffNbr,
        g.surfaceInfo, g.atom_strain_only, g.l_m_2_offset_ham_tot,
        g.l_m_2_atomdsp
#ifndef ELIMINATE_SSMAP
        ,
        g.ssmap /* the comma is NOT a bug! */
#endif          /* ELIMINATE_SSMAP  */
    );
  }
  /*!< overrided constructor; all data get initialized to parameter values. */

  geo_struct &operator=(const geo_struct &g) {
    set_geo_struct(g.periodicity, g.N_Cell, g.cell_xmin, g.cell_ymin,
                   g.cell_zmin, g.cell_xmax, g.cell_ymax, g.cell_zmax,
                   /* Additional variables to tilt the lattice */
                   /*g.tiltdir,*/ g.dx, g.dy, g.dz, g.is_tilted,
                   /*g.NaddALY,*/ g.NaddALZ, /*adjThicknessY,*/ g.adjThicknessZ,
                   g.delta_step, /*g.adjustedLenX,*/
                   /* End: Additional variables to tilt the lattice */
                   g.lattice_x, g.lattice_y, g.lattice_z, g.unitcelldivide,
                   g.cellgeom, g.ijk__cell, g.cell__ijk, g.AtomType,
                   g.ShapeType, g.isDiffNbr, g.surfaceInfo, g.atom_strain_only,
                   g.l_m_2_offset_ham_tot, g.l_m_2_atomdsp
#ifndef ELIMINATE_SSMAP
                   ,
                   g.ssmap /* the comma is NOT a bug! */
#endif                     /* ELIMINATE_SSMAP  */
    );

    return *this;
  }

  void Allocate(int, int, int, int, int, const Cell &);
  /*!< Allocate memories and Set initial values to the structure. */

  /* Additional functions to tilt the lattice */
  //   void readTiltParams(int,int); /* tilt_fix */
  void readStepRoughness(int, int); /* tilt_fix */
  void calc_tilt_params();
  /* End: Additional functions to tilt the lattice */

  void Allocate_NOSTRAIN(int, int, int, int, int, const Cell &);
  /*!< Allocate memories and Set initial values to the structure. But variables
   * for strain calculation, (i.e: 3D array l_m_2_atomdsp) is not initialized
   * here  */

  void InitializeFromInputFile(const top_struct &opt);
  /*!< Information such as lattice constant and periodicity should be directly
   * read from input xml code. This information comes to top_struct and is saved
   * in dev_struct, which is a member class of top_struct. This function
   * allocate such information to geo_struct.
   * @see top_struct::top_struct()
   * @see dev_struct::dev_struct() */

  void Deallocate();
  /*!< Deallocate all the memory. Will be used in destructor. */

  void zero();
  /*!< Set all member variables to zero. */
  ~geo_struct() { Deallocate(); }
  /*!< Destructor */

  bool InvalidateStrainOnly();
  /*!< Exclude atoms in so-called strain-only region for electonic caculation */

  void Reduce(const geo_struct *const di);
  /*!< Allocate space for geo_struct and copies information from parameter. It
   * copies geo_struct::ijk__cell, geo_struct::Atoptype,
   * geo_struct::atom_strain_only and geo_struct::cell_(x,y,z)(min,max).
   * @param di object of class geo_struct which has geometric information to be
   * used here.
   * @see geo_struct::geo_struct() */

  void
  isNbrOut(int l, int m, int n, bool *_isNotInXnegtv, bool *_isNotInXpostv,
           bool *_isNotInYnegtv, bool *_isNotInYpostv, bool *_isNotInZnegtv,
           bool *_isNotInZpostv); /* Additional function to tilt the lattice */

  void SupercellReInit();
  /*!< Initialize the super cell coupling map, which contains three indices
   * ssmap[i][j][k]. i&j index an atome's unit cell and atom repectively. k is a
   * nearest neighbor's index for the atom indexed by i&j. The mapping is to a
   * unit cell. */

  int get_cindx_tilted(int l, int m, int n);
  /*!< Determine the unit cell number of neighbor n of atom (l,m) when lattice
   * is tilted */

  void fill_is_DiffNeighbor();
  /*!< Fill the 2D array geo_struct::isDiffNbr which is used in the hamiltonian
   * construction to check whether or not the parameter vector needs to be
   * averaged over different neighbors.
   * @see geo_struct::isDiffNbr */

  void FindSurfaceAtoms();
  /*!< scan all atoms in structure and chedck if an atom is on the surface of
   * the structure. It fills geo_struct::surfaceInfo.
   * @see geo_struct::surfaceInfo */

  bool isInvalidPeriodicity();
  /*!< checks if periodicity is valid in the structure. If structure was formed
   * such that periodicity is valid, will return true. Otherwise return false.
   * @return true if structure is periodic. */

  int remove_few_bonded_atoms(int bond_to_remove);
  /*!< Function scans through the whole structure and removes atoms that have
   * only a certain set of bonds.  if "bond_to_remove" is set to 0, we will
   * remove all the atoms that have no neighbors. if "bond_to_remove" is set to
   * 1 we will remove the atoms that have only a single neighbor. The function
   * returns the number of removed atoms.
   * @param bond_to_remove # of bonds which the atoms that will be deleted has
   * @return # of removed atoms. */

  int remove_other_surface_species();
  /*!< Function scans through the whole structure and removes a certain species
   * of atoms from the surface. If the structure's periodic condition is open
   * and the type of atom on the surface is one of EL_As, EL_Sb, EL_Se and
   * EL_Ge, the atom will be removed from the surface.
   * @return # of deleted atoms from the surface of structure.
   * @see mat_def.h for more information about material types */

  int remove_user_defined_atoms(void);
  /*!< remove user-defined atoms. Criteria is defined with use of set of [cell
   * index][atom index].
   * @return # of removed atoms */

  void adjust_for_periodic_bc();
  /*!< This function records the periodic boundary coupling information into the
   * super cell map. if neighbor cell does't exist in real, we create virtual
   * cell and set the index such that periodic boundary condition is virtually
   * satisfied. */

  void remap_periodic_bc();
  /*!< readjusts ijk_cell to reflect periodic BC. */

  void FillFirstShape(const Shape_struct &shape, const top_struct &opt,
                      const vector<MaterialHandle> &);
  /*!< this function constructs geometric information such as cell_ijk,
   * cell_(x,y,z)(min,max) for this structure. This function is operated on each
   * shape defined in the input xml file. After excution, it calls
   * geo_struct::FillShape. This, should be seen as sorts of preparation for
   * 'real' filling of structure.
   * @param shape This will give information of shape of each structure defined
   * as "ShapeName" in input deck.
   * @param opt This has information of crystal structure type
   * @param mat vector of objects of MaterialHandle class. It has various
   * information of materials.
   * @see class Shape_struct
   * @see class top_struct
   * @see class MaterialHandle
   * @warning This function is defined in nemo3d/base/src/ham3d.c
   */

  void FillShape(const Shape_struct &, const vector<MaterialHandle> &, bool);
  /*!< Loops over unitcells & atoms. assigns atom types if atom is in the shape.
   * To check this, IsNotInShape function is used.
   * @param shape This will give information of shape of each structure defined
   * as "ShapeName" in input deck.
   * @param mat vector of objects of MaterialHandle class. It has various
   * information of materials.
   * @param granularity If 1, the shape has granularity.
   * @see class Shape_struct
   * @see class MaterialHandle
   * @warning This function is defined in nemo3d/base/src/ham3d.c
   */

  void InitializeForHyperChem(const char *inputfile);
  /*!< @See Cell::IntializeForHyperChem(const char* inputfile) */

  void changeLattConst(real ax, real ay, real az);
  /*!< Update lattice constants and calculate corresponding displacement.
   * This function could be quite useful if you want to use periodic BC for
   * an electronic calculation with a more sensible choice of lattice constant
   * used for the strain calculation.
   * @param ax lattice constant in x-direction
   * @param ay lattice constant in y-direction
   * @param az lattice constant in z-direction */

  void updateOnlyLattConst(real ax, real ay, real az);
  /*!< Update lattice constants only with given parameters.
   * @param ax lattice constant in x-direction
   * @param ay lattice constant in y-direction
   * @param az lattice constant in z-direction
   */

  // short functions that will be inlined
  int AtomsPerCellMax() const { return cellgeom.N_MaxAtomCell; }
  /*!< @return Maximum # of atoms a cell can have. */
  int NeighborsMax() const { return cellgeom.N_MaxNeighbor; }
  /*!< @return Maximum # of bonds that an atom can have. */
  int Neighbors(int m) const { return cellgeom.N_Neighbor[m]; }
  /*!< @param m atomic index
   * @return actual # of neighbors this atom has. */
  real *PositionInCell(int atom) { return cellgeom.CellAtom_2_ijk[atom]; }
  /*!< @param atom atomic index
   * @return position vector of this atom in the cell cellgeom.
   * @see Cell:CellAtom_2_ijk
   */
  real PositionInCell(int atom, int dir) {
    return cellgeom.CellAtom_2_ijk[atom][dir];
  }
  /*!< @param atom atomic index
   * @param dir direction. 0 for x, 1 for y and 2 for z.
   * @return position of selected direction of this atom in the cell cellgeom.
   * @see Cell:CellAtom_2_ijk
   */

  int *NbrCell(int m, int n) { return cellgeom.nnmap[m][n]; }
  /*!< @param m atomic index
   * @param n neighbor index
   * @return returns position vector of cell where n-th neighbor of m-th atom in
   * this cell is belonged.
   * @see Cell::nnmap
   */

  int NbrCell(int m, int n, int d) const { return cellgeom.nnmap[m][n][d]; }
  /*!< @param m atomic index
   * @param n neighbor index
   * @param d direction. 0 for x, 1 for y, 2 for z and 3 for the actual index of
   * atom.
   * @return position in the selected direction. Nbr(m,n,3) returns atomic
   * index. */

  double SP3_hybrid(int Za, int i, int j) const {
    return cellgeom.SP3_hybrid[Za][i][j];
  }
  /*!< @param Za atomic index
   * @param i row # of matrix.
   * @param j colomn # of matrix.
   * @return component at (i,j) of a matrix which transforms initial basis
   * (S,Px,Py,Px) to the set of hybridized orbital SP3. (i.e: weighted
   * combination of initial basis.) */

  /* Rotation matrix from new unit cell to standard Zincblende basis. */
  double get_T_rot(int i, int j) const { return cellgeom.T_rot[i][j]; }
  /* End Rotation matrix from new unit cell to standard Zincblende basis. */

  /* Normalized lattice constant for rotated unit cells. */
  double get_a_latt_rot_norm(int i) const {
    return cellgeom.a_latt_rot_norm[i];
  }
  /* End Normalized lattice constant for rotated unit cells. */

  bool isCellZB_110() const {
    return cellgeom.IsCellZB_110;
  } /* Added for ZB_110 */
  /*!< @return true if the crystal structure is ZB_110. */
  bool isCellZB_110_small() const {
    return cellgeom.IsCellZB_110_small;
  } /* Added for ZB_110_small */
  /*!< @return true if the crystal structure is ZB_110_small. */
  bool isCellZB_111() const {
    return cellgeom.IsCellZB_111;
  } /* Added for ZB_111 */
  /*!< @return true if the crystal structure is ZB_111. */
  bool isCellUserDef() const {
    return cellgeom.IsCellUserDef;
  } /* Added for CellUserDef */
  /*!< @return true if the crystal structure is CellUserDef. */
  bool isCellWurtzite() const { return cellgeom.IsCellWurtzite; }
  /*!< @return true if the crystal structure is Wurtzite. */
  bool isPeriodicX() const { return periodicity & 1; }
  /*!< @return 1 if structure is periodic in x-direction. */
  bool isPeriodicY() const { return periodicity & 2; }
  /*!< @return 1 if structure is periodic in y-direction. */
  bool isPeriodicZ() const { return periodicity & 4; }
  /*!< @return 1 if structure is periodic in z-direction. */
  bool isPeriodicXYZ() const { return (periodicity & 7) == 7; }
  /*!< @return 1 if the structure periodic in xyz-direction.  */
  bool hasPeriodicity() const { return periodicity & 7; }
  /*!< @return 1 if the structure is periodic in any direction. */
  bool StrainIsPeriodicX() const { return periodicity & 8; }
  /*!< @return 1 if strain calculation is periodic in x-direction. */
  bool StrainIsPeriodicY() const { return periodicity & 16; }
  /*!< @return 1 if strain calculation is periodic in y-direction. */
  bool StrainIsPeriodicZ() const { return periodicity & 32; }
  /*!< @return 1 if strain calculation is periodic in z-direction. */
  bool StrainIsPeriodicXYZ() const { return (periodicity & 56) == 56; }
  /*!< @return 1 if strain calculation is periodic in xyz-direction. */
  bool StrainHasPeriodicity() const { return periodicity & 56; }
  /*!< @return 1 if strain calculation is periodic in any direction. */
  bool ElIsPeriodicX() const { return periodicity & 64; }
  /*!< @return 1 if electronic calculation is periodic in x-direction. */
  bool ElIsPeriodicY() const { return periodicity & 128; }
  /*!< @return 1 if electronic calculation is periodic in y-direction. */
  bool ElIsPeriodicZ() const { return periodicity & 256; }
  /*!< @return 1 if electronic calculation is periodic in z-direction. */
  bool ElIsPeriodicXYZ() const { return (periodicity & 448) == 448; }
  /*!< @return 1 if electronic calculation is periodic in xyz-direction. */
  bool ElHasPeriodicity() const { return periodicity & 448; }
  /*!< @return 1 if electronic calculation is periodic in any direction.*/
  void zeroPeriodicity() { periodicity = (periodicity >> 3) << 3; }
  /*!< set the first 3 bits of periodicity to zero. This invalidate structural
   * periodicity.
   * @see geo_struct::periodicity */
  void setPeriodicityToEl() {
    zeroPeriodicity();
    periodicity += (periodicity & 448) >> 6;
  }
  /*!< Make periodicity of structure same as that of electric calculation. */
  void setPeriodicityToStrain() {
    zeroPeriodicity();
    periodicity += (periodicity & 56) >> 3;
  }
  /*!< Make periodicity of structure same as that of strain calculation. */
  bool isCation(int l, int m) const {
    /* The new data storage can crash on the dereferencing of the casted data */
    /* It appears to crash on the new gcc compilers                   */
    /* const int& AtomType_lm = AtomType[l][m]; */
    int AtomType_lm = AtomType[l][m];
    if (AtomType_lm == EL_As || AtomType_lm == EL_Sb || AtomType_lm == EL_Se ||
        AtomType_lm == EL_N || AtomType_lm == EL_O || AtomType_lm == EL_P ||
        AtomType_lm == EL_S)
      return false;
    else if (AtomType_lm == EL_Si || AtomType_lm == EL_Ge)
      return (2 * m < this->cellgeom.AtomsPerCellMax()) ? true : false;
    else
      return true;
  }
  /*!< Given an atom, returns boolean indicating whether atom is a cation.
   * As, Sb are automatically treated as anions. Else, if atom index is in
   * the upper half of allowed indices treat as anion. Consider everything
   * else to be a cation.
   * @param l cell index.
   * @param m atom index within cell.
   * @see geo_struct()
   * @return TRUE/FALSE if atom is a cation/anion
   * @warning The determination of anion/cation is a bit of a kludge right
   * now.  See implementation details.
   */
  bool isOnSurface(int l, int m) const { return surfaceInfo[l][m] & 1; }
  /*!< @param l cell index
   * @param m atomic index
   * @return 1 if atom(l,m) is on surface of the structure.
   * @see geo_struct::surfaceInfo */
  bool isOnSurfaceIgnoreBC(int l, int m) const { return surfaceInfo[l][m] & 2; }
  /*!< @param l cell index
   * @param m atomic index
   * @return 1 if atom(l,m) is on surface of the structure with closed BC.
   * @see geo_struct::surfaceInfo */
  bool isInX_min(int l, int m) const { return surfaceInfo[l][m] & 4; }
  /*!< @param l cell index
   * @param m atomic index
   * @return 1 if atom(l,m) is on the position of minimum x of this structure.
   * @see geo_struct::surfaceInfo */
  bool isInX_max(int l, int m) const { return surfaceInfo[l][m] & 8; }
  /*!< @param l cell index
   * @param m atomic index
   * @return 1 if atom(l,m) is on the position of maximum x of this structure.
   * @see geo_struct::surfaceInfo */
  bool isInY_min(int l, int m) const { return surfaceInfo[l][m] & 16; }
  /*!< @param l cell index
   * @param m atomic index
   * @return 1 if atom(l,m) is on the position of minimum y of this structure.
   * @see geo_struct::surfaceInfo */
  bool isInY_max(int l, int m) const { return surfaceInfo[l][m] & 32; }
  /*!< @param l cell index
   * @param m atomic index
   * @return 1 if atom(l,m) is on the position of maximum y of this structure.
   * @see geo_struct::surfaceInfo */
  bool isInZ_min(int l, int m) const { return surfaceInfo[l][m] & 64; }
  /*!< @param l cell index
   * @param m atomic index
   * @return 1 if atom(l,m) is on the position of minimum z of this structure.
   * @see geo_struct::surfaceInfo */
  bool isInZ_max(int l, int m) const { return surfaceInfo[l][m] & 128; }
  /*!< @param l cell index
   * @param m atomic index
   * @return 1 if atom(l,m) is on the position of maximum z of this structure.
   * @see geo_struct::surfaceInfo */
  int phaseInfoForNemo1D(int l, int m, int n) {
    const int *dCell = NbrCell(m, n);
    int phaseInfo = 0;

    if (isInX_min(l, m) && dCell[0] == -1)
      phaseInfo |= (ElIsPeriodicX() ? 01 : 0100);
    if (isInX_max(l, m) && dCell[0] == 1)
      phaseInfo |= (ElIsPeriodicX() ? 02 : 0200);
    if (isInY_min(l, m) && dCell[1] == -1)
      phaseInfo |= (ElIsPeriodicY() ? 04 : 0400);
    if (isInY_max(l, m) && dCell[1] == 1)
      phaseInfo |= (ElIsPeriodicY() ? 010 : 01000);
    if (isInZ_min(l, m) && dCell[2] == -1)
      phaseInfo |= (ElIsPeriodicZ() ? 020 : 02000);
    if (isInZ_max(l, m) && dCell[2] == 1)
      phaseInfo |= (ElIsPeriodicZ() ? 040 : 04000);

    return phaseInfo;
  }
  /*!< Calculate and allocate value of PhaseInfo. It determines phaseInfo as
   * follows.
   *  - Relative x location of neighbor atom is -1 and current atom is on xmin
   * (cell) of structure
   *   - Electronic calculation is periodic in x direction : PhaseInfo =
   * (PhaseInfo) | 01
   *   - Electronic calculation is not periodic in x direction : PhaseInfo =
   * (PhaseInfo) | 0100
   *  - Relative x location of neighbor atom is 1 and current atom is on xmax
   * (cell) of structure
   *   - Electronic calculation is periodic in x direction : PhaseInfo =
   * (PhaseInfo) | 02
   *   - Electronic calculation is not periodic in x direction : PhaseInfo =
   * (PhaseInfo) | 0200
   *  - Relative y location of neighbor atom is -1 and current atom is on ymin
   * (cell) of structure
   *   - Electronic calculation is periodic in y direction : PhaseInfo =
   * (PhaseInfo) | 04
   *   - Electronic calculation is not periodic in y direction : PhaseInfo =
   * (PhaseInfo) | 04
   *  - Relative y location of neighbor atom is 1 and current atom is on ymax
   * (cell) of structure
   *   - Electronic calculation is periodic in y direction : PhaseInfo =
   * (PhaseInfo) | 010
   *   - Electronic calculation is not periodic in y direction : PhaseInfo =
   * (PhaseInfo) | 01000
   *  - Relative z location of neighbor atom is -1 and current atom is on zmin
   * (cell) of structure
   *   - Electronic calculation is periodic in z direction : PhaseInfo =
   * (PhaseInfo) | 020
   *   - Electronic calculation is not periodic in z direction : PhaseInfo =
   * (PhaseInfo) | 02000
   *  - Relative z location of neighbor atom is 1 and current atom is on zmax
   * (cell) of structure
   *   - Electronic calculation is periodic in z direction : PhaseInfo =
   * (PhaseInfo) | 040
   *   - Electronic calculation is not periodic in z direction : PhaseInfo =
   * (PhaseInfo) | 04000
   *
   * @param l cell index
   * @param m atomic index
   * @param n neighbor index
   * @return phase information.
   * @brief functional description only. Don't understand physical meaning quite
   * well. */

  real phaseBloch(real kxL, real kyL, real kzL, int l, int m, int n) {
    real phase = 0.0;
    if (!is_tilted) { /* Calculate Bloch phase as follows when lattice is NOT
                         tilted */
      const int *dCell = NbrCell(m, n);
      if (isInX_min(l, m) && dCell[0] == -1)
        phase -= kxL;
      if (isInX_max(l, m) && dCell[0] == 1)
        phase += kxL;
      if (isInY_min(l, m) && dCell[1] == -1)
        phase -= kyL;
      if (isInY_max(l, m) && dCell[1] == 1)
        phase += kyL;
      if (isInZ_min(l, m) && dCell[2] == -1)
        phase -= kzL;
      if (isInZ_max(l, m) && dCell[2] == 1)
        phase += kzL;
    } /* if !is_tilted */

    if (is_tilted) { /* Calculate Bloch phase as follows when lattice is tilted
                      */
                     /* Initialize all booleans to false */
      bool _isNotInXnegtv = false;
      bool _isNotInXpostv = false;
      bool _isNotInYnegtv = false;
      bool _isNotInYpostv = false;
      bool _isNotInZnegtv = false;
      bool _isNotInZpostv = false;

      isNbrOut(l, m, n, &_isNotInXnegtv, &_isNotInXpostv, &_isNotInYnegtv,
               &_isNotInYpostv, &_isNotInZnegtv, &_isNotInZpostv);

      if (_isNotInXnegtv) {
        phase -= kxL;
        //           if(isPeriodicY()) phase -= ((double)dy/(double)dx)*kxL;
        if (isPeriodicZ())
          phase -= (((double)dz * lattice_z) / ((double)dx * lattice_x)) *
                   kzL; /* Modified for ZB_111: scaled by lattice constants */
      }
      if (_isNotInXpostv) {
        phase += kxL;
        //           if(isPeriodicY()) phase += ((double)dy/(double)dx)*kxL;
        if (isPeriodicZ())
          phase += (((double)dz * lattice_z) / ((double)dx * lattice_x)) *
                   kzL; /* Modified for ZB_111: scaled by lattice constants */
      }
      if (_isNotInYnegtv) {
        phase -= kyL;
      }
      if (_isNotInYpostv) {
        phase += kyL;
      }
      if (_isNotInZnegtv) {
        phase -= kzL;
      }
      if (_isNotInZpostv) {
        phase += kzL;
      }
    } /* if is_tilted */

    return -Pi * (phase - 2 * int(0.5 * phase));
  }
  /*!< Seems that this function shifts phase by some constant value if neighbor
   * cell is outside real structure. According to the bloch theorm, kxL,kyL,kzL
   * should be integer multiple of 2PI. But this function multiply PI before
   * returning its final phase value. Therefore kxL,kyL,kzL here have a little
   * different value from what we understand from books. Need more study to
   * understand what 'real value' of these are.
   *
   * @param kxL not sure exactly. Seems that phase shift in x-direction between
   * two adjacent unit-cell.
   * @param kyL not sure exactly. Seems that phase shift in x-direction between
   * two adjacent unit-cell.
   * @param kzL not sure exactly. Seems that phase shift in x-direction between
   * two adjacent unit-cell.
   *
   * @param l cell index
   * @param m atomic index
   * @param n neighbor index
   * @return amount of phase shift. can't understand equation for conversion
   * exactly.
   * @brief Don't understand physical meaning quite well. */

  real phasePeriodic(real kxL, real kyL, real kzL, int i_nbr, int j_nbr,
                     int k_nbr) {
    real phase = 0.0;
    if (i_nbr < this->cell_xmin)
      phase -= kxL;
    if (i_nbr > this->cell_xmax)
      phase += kxL;
    if (j_nbr < this->cell_ymin)
      phase -= kyL;
    if (j_nbr > this->cell_ymax)
      phase += kyL;
    if (k_nbr < this->cell_zmin)
      phase -= kzL;
    if (k_nbr > this->cell_zmax)
      phase += kzL;
    return -Pi * (phase - 2 * int(0.5 * phase));
  }
  /*!< What this function does is exactly same (functionally) as that
   * geo_struct::phaseBloch does. But while phasBloch function shifts phase with
   * use of (cindx, aindx, nindx), this does so with use of (i,j,k) of neighbor
   * cell.
   *
   * @param kxL not sure exactly. Seems that phase shift in x-direction between
   * two adjacent unit-cell.
   * @param kyL not sure exactly. Seems that phase shift in x-direction between
   * two adjacent unit-cell.
   * @param kzL not sure exactly. Seems that phase shift in x-direction between
   * two adjacent unit-cell.
   *
   * @param i_nbr i of neighber cell
   * @param j_nbr j of neighbor cell
   * @param k_nbr k of neighbor cell
   * @return amount of phase shift. can't understand equation for conversion
   * exactly.
   * @brief Don't understand physical meaning quite well. */

  void getPosition_unstrained(real *pos, real unstrnd_latt, int Zcell,
                              int Zatom);
  /*!< Calculate atomic position using unstrained lattice constant.
   * @param pos array that will have current atomic position info.
   * @param unstrnd_latt lattice constant without strain effect.
   * @param Zcell cell index
   * @param Zatom atomic index
   * @return position vector. pos[0],pos[1],pos[2] will have current atomic
   * position in x,y,z respectively. */
  void getPosition_strained(real *pos, int Zcell, int Zatom);
  /*!< Calculate strain effect and get new position of current atom.
   * @param pos array that will have current atomic position info.
   * @param Zcell cell index
   * @param Zatom atomic index
   * @return position vector. pos[0],pos[1],pos[2] will have current atomic
   * position in x,y,z respectively. */
  void getPosition_equilibrium(real *pos, int Zcell, int Zatom);
  /*!< calculate atomic position in equillibrium
   * @see geo_struct::getPosition_strained for description of parameters and
   * returned value. */
  void getRelativePosVec(real *nnv0, real *nnv, real unstrnd_latt, int c1,
                         int a1, int c2, int a2, int neighbor);
  /*!< get the relative position vector between two atoms.
   * - Function first will calculate relative atomic position of a1,a2 in a unit
   * cell. If neigbor atom of a2 and atom a2 are in the same unit cell, result
   * is simply the relative posion of atom a1 to atom a2 in a unit cell.
   * Otherwise, function will add +1 or -1 to the result, which becomes then
   * relative position of atom a1 in a unitcell to atom a2 in the neighbor
   * unitcell. Normalized result will be converted in a unit of unstraind
   * lattice constant and stored in nnv0.
   * - In unstrained environment, all the relative positions of atom X's are
   * same in any unitcells. But once strain starts to affect the structure, the
   * displacements of atoms X due to strain become different with different
   * indice of unitcell. Therefore, we need to get displacement with
   * consideration of cell index of each unit cell. The result in a
   * strained-structure will be saved in nnv.
   *
   *
   * @param nnv0 will have the result in an unstrained structure.
   * @param nnv will have the result in a strain effect.
   * @param unstrnd_latt lattice constant in an unstarined structure.
   * @param c1 cell index 1
   * @param a1 atomic index 1
   * @param c2 cell index 2
   * @param a2 atomic index 2
   * @param neighbor neighbor index of atom a2 (in unitcell c2)
   */
  friend ostream &operator<<(ostream &o, const geo_struct &g);
};

#endif
