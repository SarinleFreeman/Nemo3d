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
$Header: /repo/nemo3d/src/base/Cell.h,v 1.10 2007/06/14 20:08:20 nkharche Exp $
*****************************************************************************/

#ifndef CELL_H
#define CELL_H

#include <iostream>

#include "realtype.h"
#include <ivector.h>
#include <rvector.h> /* Added for ZB_111 */
#include <rmatrix.h>
#include "i3tensor.h"
#include "r3tensor.h"

#include "top_struct.H"
#include "nml_global.h"

// need to declare this class, since it is declared to be a friend
class geo_struct;


/*! \class Cell
 *  \brief Represents a single unit cell.

 *  This class encapsulates data relevant to a single cell, including its topology, neighbor information, and equilibrium (unstrained) atomic positions within the cell. */

class Cell 
{
   int N_MaxAtomCell; 
   /*!< Maximum # of atoms that current unit cell can have. */
   int N_MaxNeighbor;
   /*!< Maximum # of neighbor atoms that an atom in this unitcell can have. */
   bool IsCellZB_110; /* Added for ZB_110 */
   /*!< Returns 1 if the structure of current unit cell is ZB_110. */
   bool IsCellZB_110_small; /* Added for ZB_110_small */
   /*!< Returns 1 if the structure of current unit cell is ZB_110_small. */
   bool IsCellZB_111; /* Added for ZB_111 */
   /*!< Returns 1 if the structure of current unit cell is ZB_111. */
   bool IsCellUserDef; /* Added by for CellUserDef */
   /*!< Returns 1 if the structure of current unit cell is IsCellUserDef. */
   bool IsCellWurtzite;
   /*!< Returns 1 if the structure of current unit cell is Wurtzite. */
   ivectr N_Neighbor;
   /*!< 1D integer array. N_Neighbor[m] has the # of bonds of the atom with index m. */ 
   rmatrix CellAtom_2_ijk;
   /*!< 2D real array. More exactly, (atom # in an unit cell) by 3 matrix. CellAtom_2_ijk[m][i] has position data
    * of m-th atom in this unitcell. i(0,1,2) corresponds to x,y and z, respectively. */
   i3tensor nnmap;
   /*!< 3D integer arrray. nnmap[m][n][d] has (i,j,k) information of the unicell where n-th neighbor atom of m-th
    * atom in current unitcell belongs. (i,j,k) info's are mapped onto i=(0,1,2). nnmap[m][n][3] has tje atomic index
    * of the n-th neighbor atom. */
   r3tensor SP3_hybrid;
   /*!< 2D real array. SP3_hybrid[Za][i][j] is a matrix that transfroms orbital (s,px,py,pz) to SP3 hybridized alighend 
    * along bonds for atom Za. This matrix will be used to construct new Hamiltonian based on the new basis resulted from 
    * hybridization. Za indexs atomic sequence, i indexs initial set of basis (s,px,py,pz), and j indexes the hybridized 
    * orbital along a bond ranked by the nearest location. */
   rmatrix T_rot; /* ZB_111: Rotation matrix from new unit cell to standard Zincblende basis. */
   rvectr a_latt_rot_norm; /* ZB_111: Normalized lattice constant for rotated unit cells. */
   double a_latt_unstr; /* ZB_111: Lattice constant of the underlying Zincblende crystal. */
   void set_Cell(
                 int N_MaxAtomCell, 
                 int N_MaxNeighbor, 
                 bool IsCellZB_110, /* Added for ZB_110 */
                 bool IsCellZB_110_small, /* Added for ZB_110_small */
                 bool IsCellZB_111, /* Added for ZB_111 */
                 bool IsCellUserDef, /* Added for CellUserDef */
                 bool IsCellWurtzite,
                 ivectr N_Neighbor,
                 rmatrix CellAtom_2_ijk,
                 i3tensor nnmap,
		 r3tensor SP3_hybrid,
                 rmatrix T_rot, /* ZB_111: Rotation matrix from new unit cell to standard Zincblende basis. */
                 rvectr a_latt_rot_norm, /* ZB_111: Normalized lattice constant for rotated unit cells. */
                 double a_latt_unstr /* ZB_111: Lattice constant of the underlying Zincblende crystal. */
                 );
   /*!< Private helper function for copy constructor and assignment operator */
   
   friend class geo_struct;
 public:
   // default constructor; make sure all data get initialized to something
   Cell() : N_MaxAtomCell(0), N_MaxNeighbor(0), 
            IsCellZB_110(false), /* Added for ZB_110 */
            IsCellZB_110_small(false), /* Added for ZB_110_small */
            IsCellZB_111(false), /* Added for ZB_111 */
            IsCellUserDef(false), /* Added for CellUserDef */
            IsCellWurtzite(false), 
            N_Neighbor(NULL), CellAtom_2_ijk(NULL), 
            nnmap(NULL), SP3_hybrid(NULL),
            T_rot(NULL), /* ZB_111: Rotation matrix from new unit cell to standard Zincblende basis. */
            a_latt_rot_norm(NULL), /* ZB_111: Normalized lattice constant for rotated unit cells. */
            a_latt_unstr(0) /* ZB_111: Lattice constant of the underlying Zincblende crystal. */
{;}
   /*!< default constructor */
   Cell(const Cell &c) { 
      set_Cell(c.N_MaxAtomCell, 
               c.N_MaxNeighbor, 
               c.IsCellZB_110, /* Added for ZB_110 */
               c.IsCellZB_110_small, /* Added for ZB_110_small */
               c.IsCellZB_111, /* Added for ZB_111 */
               c.IsCellUserDef, /* Added for CellUserDef */
               c.IsCellWurtzite,
               c.N_Neighbor,
               c.CellAtom_2_ijk,
               c.nnmap,
               c.SP3_hybrid,
               c.T_rot, /* ZB_111: Rotation matrix from new unit cell to standard Zincblende basis. */
               c.a_latt_rot_norm, /* ZB_111: Normalized lattice constant for rotated unit cells. */
               c.a_latt_unstr /* ZB_111: Lattice constant of the underlying Zincblende crystal. */
               );
   }
   /*!< overrided constructor */
   Cell& operator =(const Cell &c) { 
      set_Cell(c.N_MaxAtomCell, 
               c.N_MaxNeighbor, 
               c.IsCellZB_110, /* Added for ZB_110 */
               c.IsCellZB_110_small, /* Added for ZB_110_small */
               c.IsCellZB_111, /* Added for ZB_111 */
               c.IsCellUserDef, /* Added for CellUserDef */
               c.IsCellWurtzite,
               c.N_Neighbor,
               c.CellAtom_2_ijk,
               c.nnmap,
               c.SP3_hybrid,
               c.T_rot, /* ZB_111: Rotation matrix from new unit cell to standard Zincblende basis. */
               c.a_latt_rot_norm, /* ZB_111: Normalized lattice constant for rotated unit cells. */
               c.a_latt_unstr /* ZB_111: Lattice constant of the underlying Zincblende crystal. */
               );
      return *this;
   }

   void Initialize(const top_struct& opt);
   /*!< Initialize nnmap, N_MaxAtomCell, N_Maxneighbor, N_Neighbor, CellAtom_2_ijk and SP3_hybrid 
    * (in the case of Wurtzite).   
    * @param opt it offers information of type of crystal structure.
    * @see top_struct::top_struct() */
   void InitializeForHyperChem(const char* inputfile);
    /*!< This function parses a hyperchem file for atomic locations, # of neighbers 
    * and total # of atomes in current cell.
    * @param inputfile name of input file */
   
   int AtomsPerCellMax() const { return N_MaxAtomCell;}
   /*!< @return Max # of atoms that current cell can have */
   int NeighborsMax() const { return N_MaxNeighbor;}
   /*!< @return Max # of neighbor atoms that an atom in this cell can have */
   int Neighbors(int m) const { return N_Neighbor[m];}
   /*!< @param m index of current atom
    * @return Actual # of neighbor atoms that m-th atom of this cell has */ 
   real* PositionInCell(int atom) { return CellAtom_2_ijk[atom]; }
   /*!< @param atom index of current atom 
    * @return position vector of m-th atom of this cell has. */
   real  PositionInCell(int atom, int d) { return CellAtom_2_ijk[atom][d]; }
   /*!< @param atom index of current atom
    * @param d data index.
    * @return each component of Cell_2_ijk[atom]. (d=0,1,2) correspond to (x,y,z) respectively. */
   int* NbrCell(int m, int n) { return nnmap[m][n]; }
   /*!< @param m index of current atom
    * @param n index of neighbor atom
    * @return position vector of the cell where n-th neighbor atom of m-th atom in current cell belongs. */
   int  NbrCell(int m, int n, int d) { return nnmap[m][n][d]; }
   /*!< @param m index of current atom
    * @param n index of neighbor atom
    * @param d data index
    * @return each element of nnmap[d]. (d=0,1,2,3) correspond to (x,y,z,atomic index of the neighbor). */
   friend ostream& operator<<(ostream& o, const Cell& c);

   void Deallocate();
   /*!< Cleanup function that will be used in destructor */
   ~Cell() { Deallocate(); }
   /*!< destructor */
};

#endif
