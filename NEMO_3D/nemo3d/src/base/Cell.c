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
$Header: /repo/nemo3d/src/base/Cell.c,v 1.13 2007/06/14 20:11:28 nkharche Exp $
*****************************************************************************/

#include "Cell.h"
#include "io_utils.h"
#include "util_stringops.h"

void Cell::set_Cell(
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
              rmatrix T_rot, /* Rotation matrix from new unit cell to standard Zincblende basis. */
              rvectr a_latt_rot_norm, /* Normalized lattice constant for rotated unit cells. */
              double a_latt_unstr /* Lattice constant of the underlying Zincblende crystal. */
              )
{
   this->N_MaxAtomCell = N_MaxAtomCell;
   this->N_MaxNeighbor = N_MaxNeighbor;
   this->IsCellZB_110 = IsCellZB_110; /* Added for ZB_110 */
   this->IsCellZB_110_small = IsCellZB_110_small; /* Added for ZB_110_small */
   this->IsCellZB_111 = IsCellZB_111; /* Added for ZB_111 */
   this->IsCellUserDef = IsCellUserDef; /* Added for CellUserDef */
   this->IsCellWurtzite = IsCellWurtzite;
   this->N_Neighbor = Ivectr(N_MaxAtomCell);
   for (int Za=0; Za < N_MaxAtomCell; Za++)
      this->N_Neighbor[Za] = N_Neighbor[Za];
   
   this->CellAtom_2_ijk = Rmatrix(N_MaxAtomCell, 3);
   for ( int Za=0; Za < N_MaxAtomCell; Za++ ) {
      this->CellAtom_2_ijk[Za][0] = CellAtom_2_ijk[Za][0];
      this->CellAtom_2_ijk[Za][1] = CellAtom_2_ijk[Za][1];
      this->CellAtom_2_ijk[Za][2] = CellAtom_2_ijk[Za][2];
   }

   this->nnmap = I3tensor(0, N_MaxAtomCell-1, 0, N_MaxNeighbor-1, 0, 3);
   for ( int Za=0; Za < N_MaxAtomCell; Za++ ) {
   for ( int Zn=0; Zn < N_MaxNeighbor; Zn++ ) {
      this->nnmap[Za][Zn][0] = nnmap[Za][Zn][0];
      this->nnmap[Za][Zn][1] = nnmap[Za][Zn][1];
      this->nnmap[Za][Zn][2] = nnmap[Za][Zn][2];
      this->nnmap[Za][Zn][3] = nnmap[Za][Zn][3];
   }}

   if(IsCellZB_110 || IsCellZB_110_small || IsCellZB_111 || IsCellUserDef) { /* Added for ZB_ and CellUserDef */
   this->T_rot = Rmatrix(3, 3);
   this->a_latt_rot_norm = Rvectr(3);
      for(int i=0; i<3; i++){
          for(int j=0; j<3; j++){
              this->T_rot[i][j] = T_rot[i][j];
              this->a_latt_rot_norm[j] = a_latt_rot_norm[j];
             }
         }
     }

   if(IsCellUserDef) { /* Added for CellUserDef */
   this->a_latt_unstr = a_latt_unstr; /* Lattice constant of the underlying Zincblende crystal. */
   }

   if(IsCellWurtzite) {
     this->SP3_hybrid = R3tensor(0, N_MaxAtomCell-1, 0, 3, 0, 3);
     for (int Za=0; Za< N_MaxAtomCell; Za++){
     for (int i=0; i<4; i++){
     for(int j=0; j<4; j++){
       this->SP3_hybrid[Za][i][j] = SP3_hybrid[Za][i][j];
     }}}
   }
}


/* Initialize the unit cell structure.  Arrays are defined:
   N_MaxAtomCell = # of atoms in unit cell,
   N_Neighbor[i] is the number of bonds associated ith atom the unit cell.
   CellAtom_2_ijk[i][j] contains the relative locations of the atoms in the
   unit cell - length units are the lattice constant d->a.  i references
   atom #, and j label the cartesian axes x=0,y=1,z=2. */

/* Initialize the atomic connectivity within and among unit cells
   nnmap[i][j][k] contains all information. i indexes
   atom # and j indexes atoms i's neigbors.  k=0, k=1, and
   k=2 indexes the relative unit cell location for atom i's
   jth neighbor.  k=3 contains the atom number for the jth
   neighbor. Obvious isn't it? */

/* SP3_hybrid[Za][i][j] is a matrix that transforms orbital s, px, py, pz to sp3
   hybridized orbitals aligned along bonds for atom Za. 
   i indexes s, px, py, pz and j indexes the hybridized orbital along 
   a bond ranked by the nearest neighbor location. */  

void Cell::Initialize(const top_struct& opt)
{
    int i, j;

    // now do neighbor topology
   if (opt.Dev.StructFile != "NULL") {
      cout_master << "Reading first shape from file is not fully implemented.\n";
      die("");
   }

    switch(opt.Dev.CrystalStruct) {
    case Dev_struct::Cubic:
       this->IsCellZB_110 = false; /* Added for ZB_110 */
       this->IsCellZB_110_small = false; /* Added for ZB_110_small */
       this->IsCellZB_111 = false; /* Added for ZB_111 */
       this->IsCellUserDef = false; /* Added for CellUserDef */
       this->IsCellWurtzite = false;
       this->N_MaxAtomCell = 1;
       this->N_MaxNeighbor = 6;
       this->N_Neighbor = Ivectr(this->N_MaxAtomCell);
       this->N_Neighbor[0] = 6;
       this->CellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);
       this->CellAtom_2_ijk[0][0] = 
          this->CellAtom_2_ijk[0][1] = this->CellAtom_2_ijk[0][2] = 0;
       break;
    case Dev_struct::Diamond:
       this->IsCellZB_110 = false; /* Added for ZB_110 */
       this->IsCellZB_110_small = false; /* Added for ZB_110_small */
       this->IsCellZB_111 = false; /* Added for ZB_111 */
       this->IsCellUserDef = false; /* Added for CellUserDef */
       this->IsCellWurtzite = false;
       this->N_MaxAtomCell = 8;
       this->N_MaxNeighbor = 4;
       this->N_Neighbor = Ivectr(this->N_MaxAtomCell);
       for ( i=0; i < this->N_MaxAtomCell; i++ )
          this->N_Neighbor[i] = this->N_MaxNeighbor;
       this->CellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);
       this->CellAtom_2_ijk[0][0] = this->CellAtom_2_ijk[0][1] = 
          this->CellAtom_2_ijk[0][2] = 0;
       
       this->CellAtom_2_ijk[1][0] = 0;
       this->CellAtom_2_ijk[1][1] = this->CellAtom_2_ijk[1][2] = -0.5;
       
       this->CellAtom_2_ijk[2][0] = this->CellAtom_2_ijk[2][1] = -0.5;
       this->CellAtom_2_ijk[2][2] = 0;
       
       this->CellAtom_2_ijk[3][0] = this->CellAtom_2_ijk[3][2] = -0.5;
       this->CellAtom_2_ijk[3][1] = 0;
       
       for ( i=0; i < 4; i++ )
          for ( j=0; j < 3; j++ )
             this->CellAtom_2_ijk[i+4][j] = this->CellAtom_2_ijk[i][j] - 0.25;
       break;	    
    case Dev_struct::Zincblende:
       this->IsCellZB_110 = false; /* Added for ZB_110 */
       this->IsCellZB_110_small = false; /* Added for ZB_110_small */
       this->IsCellZB_111 = false; /* Added for ZB_111 */
       this->IsCellUserDef = false; /* Added for CellUserDef */
       this->IsCellWurtzite = false;
       this->N_MaxAtomCell = 8;
       this->N_MaxNeighbor = 4;
       this->N_Neighbor = Ivectr(this->N_MaxAtomCell);
       for ( i=0; i < this->N_MaxAtomCell; i++ ) 
          this->N_Neighbor[i] = this->N_MaxNeighbor;
       this->CellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);
       this->CellAtom_2_ijk[0][0] = 
          this->CellAtom_2_ijk[0][1] = this->CellAtom_2_ijk[0][2] = 0;
       
       this->CellAtom_2_ijk[1][0] = 0;
       this->CellAtom_2_ijk[1][1] = this->CellAtom_2_ijk[1][2] = -0.5;
       
       this->CellAtom_2_ijk[2][0] = this->CellAtom_2_ijk[2][1] = -0.5;
       this->CellAtom_2_ijk[2][2] = 0;
       
       this->CellAtom_2_ijk[3][0] = this->CellAtom_2_ijk[3][2] = -0.5;
       this->CellAtom_2_ijk[3][1] = 0;
       
       for ( i=0; i < 4; i++ )
          for ( j=0; j < 3; j++ )
             this->CellAtom_2_ijk[i+4][j] = this->CellAtom_2_ijk[i][j] - 0.25;
       break;
    case Dev_struct::Diamond_FCC:
       this->IsCellZB_110 = false; /* Added for ZB_110 */
       this->IsCellZB_110_small = false; /* Added for ZB_110_small */
       this->IsCellZB_111 = false; /* Added for ZB_111 */
       this->IsCellUserDef = false; /* Added for CellUserDef */
       this->IsCellWurtzite = false;
       this->N_MaxAtomCell = 2;
       this->N_MaxNeighbor = 4;
       this->N_Neighbor = Ivectr(this->N_MaxAtomCell);
       for ( i=0; i < this->N_MaxAtomCell; i++ )
          this->N_Neighbor[i] = this->N_MaxNeighbor;
       this->CellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);
       this->CellAtom_2_ijk[0][0] = this->CellAtom_2_ijk[0][1] = 
          this->CellAtom_2_ijk[0][2] = 0;
       this->CellAtom_2_ijk[1][0] = this->CellAtom_2_ijk[1][1] = 
          this->CellAtom_2_ijk[1][2] = -0.5;	    
       break;	    
    case Dev_struct::Wurtzite:
       this->IsCellZB_110 = false; /* Added for ZB_110 */
       this->IsCellZB_110_small = false; /* Added for ZB_110_small */
       this->IsCellZB_111 = false; /* Added for ZB_111 */
       this->IsCellUserDef = false; /* Added for CellUserDef */
       this->IsCellWurtzite = true;
       this->N_MaxAtomCell = 4;
       this->N_MaxNeighbor = 4;
       this->N_Neighbor = Ivectr(this->N_MaxAtomCell);		
       for ( i=0; i < this->N_MaxAtomCell; i++ )
          this->N_Neighbor[i] = this->N_MaxNeighbor;
       this->CellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);

       this->CellAtom_2_ijk[0][0] = -0.5; 
       this->CellAtom_2_ijk[0][1] = 0.5/sqrt(3.0);
       this->CellAtom_2_ijk[0][2] = -1.0;

       this->CellAtom_2_ijk[1][0] = -0.5; 
       this->CellAtom_2_ijk[1][1] = -0.5/sqrt(3.0);
       this->CellAtom_2_ijk[1][2] = -0.5;

       this->CellAtom_2_ijk[2][0] = -0.5; 
       this->CellAtom_2_ijk[2][1] = 0.5/sqrt(3.0);
       this->CellAtom_2_ijk[2][2] = -0.625;

       this->CellAtom_2_ijk[3][0] = -0.5; 
       this->CellAtom_2_ijk[3][1] = -0.5/sqrt(3.0);
       this->CellAtom_2_ijk[3][2] = -0.125;
       break;

    /* -------- New unit cell (ZB_110)  --------
       This unit cell is useful to simulate [110] directed wires
       and quantum wells */
    case Dev_struct::ZB_110:
       this->IsCellZB_110 = true; /* Added for ZB_110 */
       this->IsCellZB_110_small = false; /* Added for ZB_110_small */
       this->IsCellZB_111 = false; /* Added for ZB_111 */
       this->IsCellUserDef = false; /* Added for CellUserDef */
       this->IsCellWurtzite = false;
       this->N_MaxAtomCell = 16;
       this->N_MaxNeighbor = 4;
       this->N_Neighbor = Ivectr(this->N_MaxAtomCell);
       for ( i=0; i < this->N_MaxAtomCell; i++ )
          this->N_Neighbor[i] = this->N_MaxNeighbor;
       this->CellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);
       this->CellAtom_2_ijk[0][0] =    0;  this->CellAtom_2_ijk[0][1] =    0;  this->CellAtom_2_ijk[0][2] =    0;
       this->CellAtom_2_ijk[1][0] = -0.5;  this->CellAtom_2_ijk[1][1] = -0.5;  this->CellAtom_2_ijk[1][2] =    0;
       this->CellAtom_2_ijk[2][0] = -0.5;  this->CellAtom_2_ijk[2][1] =    0;  this->CellAtom_2_ijk[2][2] =    0;
       this->CellAtom_2_ijk[3][0] =    0;  this->CellAtom_2_ijk[3][1] = -0.5;  this->CellAtom_2_ijk[3][2] =    0;
       this->CellAtom_2_ijk[4][0] = -0.75; this->CellAtom_2_ijk[4][1] = -0.75; this->CellAtom_2_ijk[4][2] = -0.5;
       this->CellAtom_2_ijk[5][0] = -0.25; this->CellAtom_2_ijk[5][1] = -0.25; this->CellAtom_2_ijk[5][2] = -0.5;
       this->CellAtom_2_ijk[6][0] = -0.75; this->CellAtom_2_ijk[6][1] = -0.25; this->CellAtom_2_ijk[6][2] = -0.5;
       this->CellAtom_2_ijk[7][0] = -0.25; this->CellAtom_2_ijk[7][1] = -0.75; this->CellAtom_2_ijk[7][2] = -0.5;

       for ( i=0; i < (this->N_MaxAtomCell)/2; i++ )
           {
//          for ( j=0; j < 3; j++ )
//             this->CellAtom_2_ijk[i+4][j] = this->CellAtom_2_ijk[i][j] - 0.25;
            this->CellAtom_2_ijk[i+(this->N_MaxAtomCell)/2/*8*/][0] = this->CellAtom_2_ijk[i][0] - 0.25;
            this->CellAtom_2_ijk[i+(this->N_MaxAtomCell)/2/*8*/][1] = this->CellAtom_2_ijk[i][1] - 0.00;
            this->CellAtom_2_ijk[i+(this->N_MaxAtomCell)/2/*8*/][2] = this->CellAtom_2_ijk[i][2] - 0.25;
           }

       /* Test print */
       for ( i=0; i < (this->N_MaxAtomCell); i++ )
           printf("\n%f %f %f",this->CellAtom_2_ijk[i][0],this->CellAtom_2_ijk[i][1],this->CellAtom_2_ijk[i][2]);

       break;
    /* -------- End New unit cell (ZB_110) -------- */


    /* -------- New unit cell (ZB_110_small)  --------
       This unit cell is useful to simulate [110] directed wires
       and quantum wells */
    case Dev_struct::ZB_110_small:
       this->IsCellZB_110 = false; /* Added for ZB_110 */
       this->IsCellZB_110_small = true; /* Added for ZB_110_small */
       this->IsCellZB_111 = false; /* Added for ZB_111 */
       this->IsCellUserDef = false; /* Added for CellUserDef */
       this->IsCellWurtzite = false;
       this->N_MaxAtomCell = 4;
       this->N_MaxNeighbor = 4;
       this->N_Neighbor = Ivectr(this->N_MaxAtomCell);
       for ( i=0; i < this->N_MaxAtomCell; i++ )
          this->N_Neighbor[i] = this->N_MaxNeighbor;
       this->CellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);
       this->CellAtom_2_ijk[0][0] =    0;  this->CellAtom_2_ijk[0][1] =    0;  this->CellAtom_2_ijk[0][2] =    0;
       this->CellAtom_2_ijk[1][0] = -0.5;  this->CellAtom_2_ijk[1][1] = -0.5;  this->CellAtom_2_ijk[1][2] = -0.5;

       for ( i=0; i < (this->N_MaxAtomCell)/2; i++ )
           {
//          for ( j=0; j < 3; j++ )
//             this->CellAtom_2_ijk[i+4][j] = this->CellAtom_2_ijk[i][j] - 0.25;
            this->CellAtom_2_ijk[i+(this->N_MaxAtomCell)/2/*8*/][0] = this->CellAtom_2_ijk[i][0] - 0.50;
            this->CellAtom_2_ijk[i+(this->N_MaxAtomCell)/2/*8*/][1] = this->CellAtom_2_ijk[i][1] - 0.00;
            this->CellAtom_2_ijk[i+(this->N_MaxAtomCell)/2/*8*/][2] = this->CellAtom_2_ijk[i][2] - 0.25;
           }

       /* Test print */
       for ( i=0; i < (this->N_MaxAtomCell); i++ )
           printf("\n%f %f %f",this->CellAtom_2_ijk[i][0],this->CellAtom_2_ijk[i][1],this->CellAtom_2_ijk[i][2]);

       break;
    /* -------- End New unit cell (ZB_110_small) -------- */


    /* -------- New unit cell (ZB_111) --------
       This unit cell is useful to simulate [111] directed wires
       and quantum wells */
    case Dev_struct::ZB_111:
       this->IsCellZB_110 = false; /* Added for ZB_110 */
       this->IsCellZB_110_small = false; /* Added for ZB_110_small */
       this->IsCellZB_111 = true; /* Added for ZB_111 */
       this->IsCellUserDef = false; /* Added for CellUserDef */
       this->IsCellWurtzite = false;
       this->N_MaxAtomCell = 12;
       this->N_MaxNeighbor = 4;
       this->N_Neighbor = Ivectr(this->N_MaxAtomCell);
       for ( i=0; i < this->N_MaxAtomCell; i++ )
          this->N_Neighbor[i] = this->N_MaxNeighbor;
       this->CellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);
       this->CellAtom_2_ijk[0][0] =      0;  this->CellAtom_2_ijk[0][1] =      0;  this->CellAtom_2_ijk[0][2] =      0;
       this->CellAtom_2_ijk[1][0] = -1.0/2;  this->CellAtom_2_ijk[1][1] = -1.0/2;  this->CellAtom_2_ijk[1][2] =      0;
       this->CellAtom_2_ijk[2][0] = -1.0/3;  this->CellAtom_2_ijk[2][1] =      0;  this->CellAtom_2_ijk[2][2] = -1.0/3;
       this->CellAtom_2_ijk[3][0] = -5.0/6;  this->CellAtom_2_ijk[3][1] = -1.0/2;  this->CellAtom_2_ijk[3][2] = -1.0/3;
       this->CellAtom_2_ijk[4][0] = -2.0/3;  this->CellAtom_2_ijk[4][1] =      0;  this->CellAtom_2_ijk[4][2] = -2.0/3;
       this->CellAtom_2_ijk[5][0] = -1.0/6;  this->CellAtom_2_ijk[5][1] = -1.0/2;  this->CellAtom_2_ijk[5][2] = -2.0/3;

       for ( i=0; i < (this->N_MaxAtomCell)/2; i++ )
           {
//          for ( j=0; j < 3; j++ )
//             this->CellAtom_2_ijk[i+4][j] = this->CellAtom_2_ijk[i][j] - 0.25;
            this->CellAtom_2_ijk[i+(this->N_MaxAtomCell)/2][0] = this->CellAtom_2_ijk[i][0] -     0;
            this->CellAtom_2_ijk[i+(this->N_MaxAtomCell)/2][1] = this->CellAtom_2_ijk[i][1] -     0;
            this->CellAtom_2_ijk[i+(this->N_MaxAtomCell)/2][2] = this->CellAtom_2_ijk[i][2] - 1.0/4;
           }

       /* Test print */
       for ( i=0; i < (this->N_MaxAtomCell); i++ )
           printf("\n%f %f %f",this->CellAtom_2_ijk[i][0],this->CellAtom_2_ijk[i][1],this->CellAtom_2_ijk[i][2]);

       break;
    /* -------- End New unit cell (ZB_111) -------- */

    /* -------- New unit cell (CellUserDef) --------
       This unit cell is useful to simulate arbitrary directed wires
       and quantum wells */
    case Dev_struct::CellUserDef:
       this->IsCellZB_110 = false; /* Added for ZB_110 */
       this->IsCellZB_110_small = false; /* Added for ZB_110 */
       this->IsCellZB_111 = false; /* Added for ZB_111 */
       this->IsCellUserDef = true; /* Added for CellUserDef */
       this->IsCellWurtzite = false;

       int m,n;
       double temp_coord[3];
       FILE *fp_pdb_in;
       char s[1000];
       /* Get pdb input file */
//       printf("Enter name of pdb file containing unit cell atomic positions> ");
//       scanf("%s",ifile);

       fp_pdb_in = fopen("pdb_in.txt","r");
       if (fp_pdb_in == NULL)
           die("\nERROR: No pdb_in.txt file in the current directory.\n");

       fgets(s,1000,fp_pdb_in);
       sscanf(s,"%d",&m);
       fgets(s,1000,fp_pdb_in);
       sscanf(s,"%d",&n);
       this->N_MaxAtomCell = m;
       this->N_MaxNeighbor = n;
       /* Test print */
       printf("\nN_MaxAtomCell = %d N_MaxNeighbor = %d\n",this->N_MaxAtomCell,this->N_MaxNeighbor);

       this->N_Neighbor = Ivectr(this->N_MaxAtomCell);
       for ( i=0; i < this->N_MaxAtomCell; i++ )
          this->N_Neighbor[i] = this->N_MaxNeighbor;

       this->CellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);

       this->T_rot = Rmatrix(3,3);
       this->a_latt_rot_norm = Rvectr(3);

       /* Discard next 2 lines in pdb_in.txt */
       for(int i=0;i<2;i++)
           fgets(s,1000,fp_pdb_in);

       /* Read rotation matrix from pdb_in.txt */
       for(int i=0;i<3;i++) {
           fgets(s,1000,fp_pdb_in);
           sscanf(s,"%lf%lf%lf",&temp_coord[0],&temp_coord[1],&temp_coord[2]);
           for(int j=0;j<3;j++)
               this->T_rot[i][j] = temp_coord[j];
          }

       /* Test print */
       printf("\nT_rot read from pdb file:\n");
       for(int i=0;i<3; i++ )
           printf("\n%12.10le %12.10le %12.10le",this->T_rot[i][0],this->T_rot[i][1],this->T_rot[i][2]);

       /* Discard next 2 lines in pdb_in.txt */
       for(int i=0;i<2;i++)
           fgets(s,1000,fp_pdb_in);

       /* Read  lenghts of unit cell in xyz directions (lx, ly, lz) from pdb_in.txt */
       fgets(s,1000,fp_pdb_in);
       sscanf(s,"%lf%lf%lf",&temp_coord[0],&temp_coord[1],&temp_coord[2]);
       for(int j=0;j<3;j++)
           this->a_latt_rot_norm[j] = temp_coord[j];

       /* Test print */
       printf("\na_latt_rot_norm read from pdb file:\n");
       printf("\n%12.10le %12.10le %12.10le",this->a_latt_rot_norm[0],this->a_latt_rot_norm[1],this->a_latt_rot_norm[2]);

       /* Discard next 2 lines in pdb_in.txt */
       for(int i=0;i<2;i++)
           fgets(s,1000,fp_pdb_in);

       /* Read unstrained lattice constant of the underlying Zincblende crystal */
       fgets(s,1000,fp_pdb_in);
       sscanf(s,"%lf",&temp_coord[0]);
       this->a_latt_unstr = temp_coord[0];

       /* Discard next 2 lines in pdb_in.txt */
       for(int i=0;i<2;i++)
           fgets(s,1000,fp_pdb_in);

       /* Read atomic positions */
       for(int i=0;i<this->N_MaxAtomCell;i++) {
           fgets(s,1000,fp_pdb_in);
               sscanf(s,"%lf%lf%lf",&temp_coord[0],&temp_coord[1],&temp_coord[2]);
           for(int j=0;j<3;j++)
               this->CellAtom_2_ijk[i][j] = temp_coord[j];
          }

       /* Test print */
       printf("\nCellAtom_2_ijk read from pdb file:\n");
       for(int i=0;i<this->N_MaxAtomCell; i++ )
           printf("\n%12.10le %12.10le %12.10le",this->CellAtom_2_ijk[i][0],this->CellAtom_2_ijk[i][1],this->CellAtom_2_ijk[i][2]);

       printf("\n\nCellAtom_2_ijk normalized and shifted back by one unit cell \n");
       for(int i=0;i<this->N_MaxAtomCell;i++) {
           printf("%d  ",i);
           for(int j=0;j<3;j++) {
               this->CellAtom_2_ijk[i][j] = this->CellAtom_2_ijk[i][j]/this->a_latt_rot_norm[j] - 1.0;
               printf("%12.10le  ",this->CellAtom_2_ijk[i][j]);
           }
       printf("\n");
       }
       printf("\n\n\n");


       break;
    /* -------- End New unit cell (ZB_111) -------- */

    default:
       break;
    }

    switch(opt.Dev.CrystalStruct) {
    case(Dev_struct::Cubic):
       this->nnmap = I3tensor(0, N_MaxAtomCell-1, 0, N_MaxNeighbor-1, 0, 3);
       
       this->nnmap[0][0][0] = 1;
       this->nnmap[0][0][1] = this->nnmap[0][0][2] = 0;
       this->nnmap[0][0][3] = 0;
       
       this->nnmap[0][1][0] = -1;
       this->nnmap[0][1][1] = this->nnmap[0][0][2] = 0;
       this->nnmap[0][1][3] = 0;
       
       this->nnmap[0][2][1] = 1;
       this->nnmap[0][2][0] = this->nnmap[0][0][2] = 0;
       this->nnmap[0][2][3] = 0;
       
       this->nnmap[0][3][1] = -1;
       this->nnmap[0][3][0] = this->nnmap[0][0][2] = 0;
       this->nnmap[0][3][3] = 0;
       
       this->nnmap[0][4][2] = 1;
       this->nnmap[0][4][0] = this->nnmap[0][0][1] = 0;
       this->nnmap[0][4][3] = 0;
       
       this->nnmap[0][5][2] = -1;
       this->nnmap[0][5][0] = this->nnmap[0][0][1] = 0;
       this->nnmap[0][5][3] = 0;
       break;
    case(Dev_struct::Zincblende_FCC):
    case(Dev_struct::Diamond_FCC):
       this->nnmap = I3tensor(0, N_MaxAtomCell-1, 0, N_MaxNeighbor-1, 0, 3);
       this->nnmap[0][0][0] = this->nnmap[0][0][1] = this->nnmap[0][0][2] = 0;
       this->nnmap[0][0][3] = 1;
       
       this->nnmap[0][1][0] = 0;
       this->nnmap[0][1][1] = this->nnmap[0][1][2] = 1;
       this->nnmap[0][1][3] = 1;
       
       this->nnmap[0][2][0] = this->nnmap[0][2][2] = 1;
       this->nnmap[0][2][1] = 0;
       this->nnmap[0][2][3] = 1;
       
       this->nnmap[0][3][0] = this->nnmap[0][3][1] = 1;
       this->nnmap[0][3][2] = 0;
       this->nnmap[0][3][3] = 1;
       
       this->nnmap[1][0][0] = this->nnmap[1][0][1] = this->nnmap[1][0][2] = 0;
       this->nnmap[1][0][3] = 0;
       
       this->nnmap[1][1][0] = 0;
       this->nnmap[1][1][1] = this->nnmap[1][1][2] = -1;
       this->nnmap[1][1][3] = 0;
       
       this->nnmap[1][2][0] = this->nnmap[1][2][2] = -1;
       this->nnmap[1][2][1] = 0;
       this->nnmap[1][2][3] = 0;	    
       
       this->nnmap[1][3][0] = this->nnmap[1][3][1] = -1;
       this->nnmap[1][3][2] = 0;
       this->nnmap[1][3][3] = 0;	    
       break;
    case(Dev_struct::Zincblende):
    case(Dev_struct::Diamond):
       this->nnmap = I3tensor(0, N_MaxAtomCell-1, 0, N_MaxNeighbor-1, 0, 3); 
       this->nnmap[0][0][0] = this->nnmap[0][0][1] = this->nnmap[0][0][2] = 0;
       this->nnmap[0][0][3] = 4;
       
       this->nnmap[0][1][0] = this->nnmap[0][1][2] = 1;
       this->nnmap[0][1][1] = 0;
       this->nnmap[0][1][3] = 7;
       
       this->nnmap[0][2][0] = this->nnmap[0][2][1] = 1;
       this->nnmap[0][2][2] = 0;
       this->nnmap[0][2][3] = 6;
       
       this->nnmap[0][3][0] = 0;
       this->nnmap[0][3][1] = this->nnmap[0][3][2] = 1;
       this->nnmap[0][3][3] = 5;
       
       this->nnmap[1][0][0] = this->nnmap[1][0][1] = this->nnmap[1][0][2] = 0;
       this->nnmap[1][0][3] = 5;
       
       this->nnmap[1][1][0] = this->nnmap[1][1][1] = this->nnmap[1][1][2] = 0;
       this->nnmap[1][1][3] = 4;
       
       this->nnmap[1][2][0] = 1;
       this->nnmap[1][2][1] = this->nnmap[1][2][2] = 0;
       this->nnmap[1][2][3] = 6;
       
       this->nnmap[1][3][0] = 1;
       this->nnmap[1][3][1] = this->nnmap[1][3][2] = 0;
       this->nnmap[1][3][3] = 7;
       
       this->nnmap[2][0][0] = this->nnmap[2][0][1] = this->nnmap[2][0][2] = 0;
       this->nnmap[2][0][3] = 4;
       
       this->nnmap[2][1][0] = this->nnmap[2][1][1] = this->nnmap[2][1][2] = 0;
       this->nnmap[2][1][3] = 6;
       
       this->nnmap[2][2][0] = this->nnmap[2][2][1] = 0;
       this->nnmap[2][2][2] = 1;
       this->nnmap[2][2][3] = 5;
       
       this->nnmap[2][3][0] = this->nnmap[2][3][1] = 0;
       this->nnmap[2][3][2] = 1;
       this->nnmap[2][3][3] = 7;
       
       this->nnmap[3][0][0] = this->nnmap[3][0][1] = this->nnmap[3][0][2] = 0;
       this->nnmap[3][0][3] = 4;
       
       this->nnmap[3][1][0] = this->nnmap[3][1][1] = this->nnmap[3][1][2] = 0;
       this->nnmap[3][1][3] = 7;
       
       this->nnmap[3][2][0] = this->nnmap[3][2][2] = 0;
       this->nnmap[3][2][1] = 1;
       this->nnmap[3][2][3] = 5;
       
       this->nnmap[3][3][0] = this->nnmap[3][3][2] = 0;
       this->nnmap[3][3][1] = 1;
       this->nnmap[3][3][3] = 6;
       
       this->nnmap[4][0][0] = this->nnmap[4][0][1] = this->nnmap[4][0][2] = 0;
       this->nnmap[4][0][3] = 0;
       
       this->nnmap[4][1][0] = this->nnmap[4][1][1] = this->nnmap[4][1][2] = 0;
       this->nnmap[4][1][3] = 1;
       
       this->nnmap[4][2][0] = this->nnmap[4][2][1] = this->nnmap[4][2][2] = 0;
       this->nnmap[4][2][3] = 2;
       
       this->nnmap[4][3][0] = this->nnmap[4][3][1] = this->nnmap[4][3][2] = 0;
       this->nnmap[4][3][3] = 3;
       
       this->nnmap[5][0][0] = this->nnmap[5][0][1] = this->nnmap[5][0][2] = 0;
       this->nnmap[5][0][3] = 1;
       
       this->nnmap[5][1][0] = this->nnmap[5][1][2] = 0;
       this->nnmap[5][1][1] = -1;
       this->nnmap[5][1][3] = 3;
       
       this->nnmap[5][2][0] = this->nnmap[5][2][1] = 0;
       this->nnmap[5][2][2] = -1;
       this->nnmap[5][2][3] = 2;
       
       this->nnmap[5][3][0] = 0;
       this->nnmap[5][3][1] = this->nnmap[5][3][2] = -1;
       this->nnmap[5][3][3] = 0;
       
       this->nnmap[6][0][0] = this->nnmap[6][0][1] = this->nnmap[6][0][2] = 0;
       this->nnmap[6][0][3] = 2;
       
       this->nnmap[6][1][0] = this->nnmap[6][1][2] = 0;
       this->nnmap[6][1][1] = -1;
       this->nnmap[6][1][3] = 3;
       
       this->nnmap[6][2][0] = -1;
       this->nnmap[6][2][1] = this->nnmap[6][2][2] = 0;
       this->nnmap[6][2][3] = 1;
       
       this->nnmap[6][3][0] = this->nnmap[6][3][1] = -1;
       this->nnmap[6][3][2] = 0;
       this->nnmap[6][3][3] = 0;
       
       this->nnmap[7][0][0] = this->nnmap[7][0][1] = this->nnmap[7][0][2] = 0;
       this->nnmap[7][0][3] = 3;
       
       this->nnmap[7][1][0] = this->nnmap[7][1][1] = 0;
       this->nnmap[7][1][2] = -1;
       this->nnmap[7][1][3] = 2;
       
       this->nnmap[7][2][0] = -1;
       this->nnmap[7][2][1] = this->nnmap[7][2][2] = 0;
       this->nnmap[7][2][3] = 1;
       
       this->nnmap[7][3][0] = this->nnmap[7][3][2] = -1;
       this->nnmap[7][3][1] = 0;
       this->nnmap[7][3][3] = 0;
       break;
    case(Dev_struct::Wurtzite):
       this->nnmap = I3tensor(0, N_MaxAtomCell-1, 0, N_MaxNeighbor-1, 0, 3);

       this->nnmap[0][0][0] = 0; this->nnmap[0][0][1] = 0; 
       this->nnmap[0][0][2] = 0; this->nnmap[0][0][3] = 2;

       this->nnmap[0][1][0] = 0; this->nnmap[0][1][1] = 0;
       this->nnmap[0][1][2] = -1; this->nnmap[0][1][3] = 3;

       this->nnmap[0][2][0] = 0; this->nnmap[0][2][1] = -1;
       this->nnmap[0][2][2] = -1; this->nnmap[0][2][3] = 3;

       this->nnmap[0][3][0] = 1; this->nnmap[0][3][1] = 0;
       this->nnmap[0][3][2] = -1; this->nnmap[0][3][3] = 3;

       this->nnmap[1][0][0] = 0; this->nnmap[1][0][1] = 0; 
       this->nnmap[1][0][2] = 0; this->nnmap[1][0][3] = 3;

       this->nnmap[1][1][0] = 0; this->nnmap[1][1][1] = 0;
       this->nnmap[1][1][2] = 0; this->nnmap[1][1][3] = 2;

       this->nnmap[1][2][0] = 0; this->nnmap[1][2][1] = 1;
       this->nnmap[1][2][2] = 0; this->nnmap[1][2][3] = 2;

       this->nnmap[1][3][0] = -1; this->nnmap[1][3][1] = 0;
       this->nnmap[1][3][2] = 0; this->nnmap[1][3][3] = 2;

       this->nnmap[2][0][0] = 0; this->nnmap[2][0][1] = 0; 
       this->nnmap[2][0][2] = 0; this->nnmap[2][0][3] = 0;

       this->nnmap[2][1][0] = 0; this->nnmap[2][1][1] = 0;
       this->nnmap[2][1][2] = 0; this->nnmap[2][1][3] = 1;

       this->nnmap[2][2][0] = 0; this->nnmap[2][2][1] = -1;
       this->nnmap[2][2][2] = 0; this->nnmap[2][2][3] = 1;

       this->nnmap[2][3][0] = 1; this->nnmap[2][3][1] = 0;
       this->nnmap[2][3][2] = 0; this->nnmap[2][3][3] = 1;

       this->nnmap[3][0][0] = 0; this->nnmap[3][0][1] = 0; 
       this->nnmap[3][0][2] = 0; this->nnmap[3][0][3] = 1;

       this->nnmap[3][1][0] = 0; this->nnmap[3][1][1] = 0;
       this->nnmap[3][1][2] = 1; this->nnmap[3][1][3] = 0;

       this->nnmap[3][2][0] = 0; this->nnmap[3][2][1] = 1;
       this->nnmap[3][2][2] = 1; this->nnmap[3][2][3] = 0;

       this->nnmap[3][3][0] = -1; this->nnmap[3][3][1] = 0;
       this->nnmap[3][3][2] = 1; this->nnmap[3][3][3] = 0;

       this->SP3_hybrid = R3tensor(0, N_MaxAtomCell-1, 0, 3, 0, 3);
   
       this->SP3_hybrid[0][0][0] = 0.5; 
       this->SP3_hybrid[0][0][1] = 0.0; 
       this->SP3_hybrid[0][0][2] = 0.0; 
       this->SP3_hybrid[0][0][3] = sqrt(3.0)*0.5;
       this->SP3_hybrid[0][1][0] = 0.5; 
       this->SP3_hybrid[0][1][1] = 0.0; 
       this->SP3_hybrid[0][1][2] = -sqrt(6.0)/3.0; 
       this->SP3_hybrid[0][1][3] = -sqrt(3.0)/6.0;
       this->SP3_hybrid[0][2][0] = 0.5; 
       this->SP3_hybrid[0][2][1] = -1.0/sqrt(2.0); 
       this->SP3_hybrid[0][2][2] = 1.0/sqrt(6.0); 
       this->SP3_hybrid[0][2][3] = -sqrt(3.0)/6.0;
       this->SP3_hybrid[0][3][0] = 0.5; 
       this->SP3_hybrid[0][3][1] = 1.0/sqrt(2.0); 
       this->SP3_hybrid[0][3][2] = 1.0/sqrt(6.0); 
       this->SP3_hybrid[0][3][3] = -sqrt(3.0)/6.0;

       this->SP3_hybrid[1][0][0] = 0.5; 
       this->SP3_hybrid[1][0][1] = 0.0; 
       this->SP3_hybrid[1][0][2] = 0.0; 
       this->SP3_hybrid[1][0][3] = sqrt(3.0)*0.5;
       this->SP3_hybrid[1][1][0] = 0.5; 
       this->SP3_hybrid[1][1][1] = 0.0; 
       this->SP3_hybrid[1][1][2] = sqrt(6.0)/3.0; 
       this->SP3_hybrid[1][1][3] = -sqrt(3.0)/6.0;
       this->SP3_hybrid[1][2][0] = 0.5; 
       this->SP3_hybrid[1][2][1] = 1.0/sqrt(2.0); 
       this->SP3_hybrid[1][2][2] = -1.0/sqrt(6.0); 
       this->SP3_hybrid[1][2][3] = -sqrt(3.0)/6.0;
       this->SP3_hybrid[1][3][0] = 0.5; 
       this->SP3_hybrid[1][3][1] = -1.0/sqrt(2.0); 
       this->SP3_hybrid[1][3][2] = -1.0/sqrt(6.0); 
       this->SP3_hybrid[1][3][3] = -sqrt(3.0)/6.0;

       this->SP3_hybrid[2][0][0] = 0.5; 
       this->SP3_hybrid[2][0][1] = 0.0; 
       this->SP3_hybrid[2][0][2] = 0.0; 
       this->SP3_hybrid[2][0][3] = -sqrt(3.0)*0.5;
       this->SP3_hybrid[2][1][0] = 0.5; 
       this->SP3_hybrid[2][1][1] = 0.0; 
       this->SP3_hybrid[2][1][2] = -sqrt(6.0)/3.0; 
       this->SP3_hybrid[2][1][3] = sqrt(3.0)/6.0;
       this->SP3_hybrid[2][2][0] = 0.5; 
       this->SP3_hybrid[2][2][1] = -1.0/sqrt(2.0); 
       this->SP3_hybrid[2][2][2] = 1.0/sqrt(6.0); 
       this->SP3_hybrid[2][2][3] = sqrt(3.0)/6.0;
       this->SP3_hybrid[2][3][0] = 0.5; 
       this->SP3_hybrid[2][3][1] = 1.0/sqrt(2.0); 
       this->SP3_hybrid[2][3][2] = 1.0/sqrt(6.0); 
       this->SP3_hybrid[2][3][3] = sqrt(3.0)/6.0;

       this->SP3_hybrid[3][0][0] = 0.5; 
       this->SP3_hybrid[3][0][1] = 0.0; 
       this->SP3_hybrid[3][0][2] = 0.0; 
       this->SP3_hybrid[3][0][3] = -sqrt(3.0)*0.5;
       this->SP3_hybrid[3][1][0] = 0.5; 
       this->SP3_hybrid[3][1][1] = 0.0; 
       this->SP3_hybrid[3][1][2] = sqrt(6.0)/3.0; 
       this->SP3_hybrid[3][1][3] = sqrt(3.0)/6.0;
       this->SP3_hybrid[3][2][0] = 0.5; 
       this->SP3_hybrid[3][2][1] = 1.0/sqrt(2.0); 
       this->SP3_hybrid[3][2][2] = -1.0/sqrt(6.0); 
       this->SP3_hybrid[3][2][3] = sqrt(3.0)/6.0;
       this->SP3_hybrid[3][3][0] = 0.5; 
       this->SP3_hybrid[3][3][1] = -1.0/sqrt(2.0); 
       this->SP3_hybrid[3][3][2] = -1.0/sqrt(6.0); 
       this->SP3_hybrid[3][3][3] = sqrt(3.0)/6.0;

       break;

    /* ---------- New unit cell (ZB_110) ----------- */
    case(Dev_struct::ZB_110):
     this->nnmap = I3tensor(0, N_MaxAtomCell-1, 0, N_MaxNeighbor-1, 0, 3);
     /* Atom 0 */
     this->nnmap[0][0][0] = 0; this->nnmap[0][0][1] = 0; this->nnmap[0][0][2] = 0; this->nnmap[0][0][3] =  8;//9;
     this->nnmap[0][1][0] = 1; this->nnmap[0][1][1] = 0; this->nnmap[0][1][2] = 1; this->nnmap[0][1][3] = 14;//15;
     this->nnmap[0][2][0] = 1; this->nnmap[0][2][1] = 0; this->nnmap[0][2][2] = 0; this->nnmap[0][2][3] = 10;//11;
     this->nnmap[0][3][0] = 1; this->nnmap[0][3][1] = 1; this->nnmap[0][3][2] = 1; this->nnmap[0][3][3] = 12;//13;

     /* Atom 1 */
     this->nnmap[1][0][0] = 0; this->nnmap[1][0][1] = 0; this->nnmap[1][0][2] = 0; this->nnmap[1][0][3] = 9;//10;
     this->nnmap[1][1][0] = 0; this->nnmap[1][1][1] = 0; this->nnmap[1][1][2] = 1; this->nnmap[1][1][3] = 15;//16;
     this->nnmap[1][2][0] = 0; this->nnmap[1][2][1] = 0; this->nnmap[1][2][2] = 0; this->nnmap[1][2][3] = 11;//12;
     this->nnmap[1][3][0] = 0; this->nnmap[1][3][1] = 0; this->nnmap[1][3][2] = 1; this->nnmap[1][3][3] = 13;//14;

     /* Atom 2 */
     this->nnmap[2][0][0] = 0; this->nnmap[2][0][1] = 0; this->nnmap[2][0][2] = 0; this->nnmap[2][0][3] = 10;//11;
     this->nnmap[2][1][0] = 0; this->nnmap[2][1][1] = 0; this->nnmap[2][1][2] = 1; this->nnmap[2][1][3] = 13;//14;
     this->nnmap[2][2][0] = 0; this->nnmap[2][2][1] = 0; this->nnmap[2][2][2] = 0; this->nnmap[2][2][3] =  8;//9;
     this->nnmap[2][3][0] = 0; this->nnmap[2][3][1] = 1; this->nnmap[2][3][2] = 1; this->nnmap[2][3][3] = 15;//16;

     /* Atom 3 */
     this->nnmap[3][0][0] = 0; this->nnmap[3][0][1] = 0; this->nnmap[3][0][2] = 0; this->nnmap[3][0][3] = 11;//12;
     this->nnmap[3][1][0] = 1; this->nnmap[3][1][1] = 0; this->nnmap[3][1][2] = 1; this->nnmap[3][1][3] = 12;//13;
     this->nnmap[3][2][0] = 1; this->nnmap[3][2][1] = 0; this->nnmap[3][2][2] = 0; this->nnmap[3][2][3] = 9;//10;
     this->nnmap[3][3][0] = 1; this->nnmap[3][3][1] = 0; this->nnmap[3][3][2] = 1; this->nnmap[3][3][3] = 14;//15;

     /* Atom 4 */
     this->nnmap[4][0][0] = 0; this->nnmap[4][0][1] =  0; this->nnmap[4][0][2] = 0; this->nnmap[4][0][3] = 12;//13;
     this->nnmap[4][1][0] = 0; this->nnmap[4][1][1] = -1; this->nnmap[4][1][2] = 0; this->nnmap[4][1][3] = 10;//11;
     this->nnmap[4][2][0] = 0; this->nnmap[4][2][1] =  0; this->nnmap[4][2][2] = 0; this->nnmap[4][2][3] = 15;//16;
     this->nnmap[4][3][0] = 0; this->nnmap[4][3][1] =  0; this->nnmap[4][3][2] = 0; this->nnmap[4][3][3] = 9;//10;

     /* Atom 5 */
     this->nnmap[5][0][0] = 0; this->nnmap[5][0][1] = 0; this->nnmap[5][0][2] = 0; this->nnmap[5][0][3] = 13;//14;
     this->nnmap[5][1][0] = 0; this->nnmap[5][1][1] = 0; this->nnmap[5][1][2] = 0; this->nnmap[5][1][3] = 11;//12;
     this->nnmap[5][2][0] = 1; this->nnmap[5][2][1] = 0; this->nnmap[5][2][2] = 0; this->nnmap[5][2][3] = 14;//15;
     this->nnmap[5][3][0] = 0; this->nnmap[5][3][1] = 0; this->nnmap[5][3][2] = 0; this->nnmap[5][3][3] =  8;//9;

     /* Atom 6 */
     this->nnmap[6][0][0] = 0; this->nnmap[6][0][1] = 0; this->nnmap[6][0][2] = 0; this->nnmap[6][0][3] = 14;//15;
     this->nnmap[6][1][0] = 0; this->nnmap[6][1][1] = 0; this->nnmap[6][1][2] = 0; this->nnmap[6][1][3] = 9;//10;
     this->nnmap[6][2][0] = 0; this->nnmap[6][2][1] = 0; this->nnmap[6][2][2] = 0; this->nnmap[6][2][3] = 13;//14;
     this->nnmap[6][3][0] = 0; this->nnmap[6][3][1] = 0; this->nnmap[6][3][2] = 0; this->nnmap[6][3][3] = 10;//11;

     /* Atom 7 */
     this->nnmap[7][0][0] = 0; this->nnmap[7][0][1] =  0; this->nnmap[7][0][2] = 0; this->nnmap[7][0][3] = 15;//16;
     this->nnmap[7][1][0] = 0; this->nnmap[7][1][1] = -1; this->nnmap[7][1][2] = 0; this->nnmap[7][1][3] =  8;//9;
     this->nnmap[7][2][0] = 1; this->nnmap[7][2][1] =  0; this->nnmap[7][2][2] = 0; this->nnmap[7][2][3] = 12;//13;
     this->nnmap[7][3][0] = 0; this->nnmap[7][3][1] =  0; this->nnmap[7][3][2] = 0; this->nnmap[7][3][3] = 11;//12;

     /* Atom 8 */
     this->nnmap[8][0][0] = 0; this->nnmap[8][0][1] = 0; this->nnmap[8][0][2] = 0; this->nnmap[8][0][3] = 0;//1;
     this->nnmap[8][1][0] = 0; this->nnmap[8][1][1] = 1; this->nnmap[8][1][2] = 0; this->nnmap[8][1][3] = 7;//8;
     this->nnmap[8][2][0] = 0; this->nnmap[8][2][1] = 0; this->nnmap[8][2][2] = 0; this->nnmap[8][2][3] = 2;//3;
     this->nnmap[8][3][0] = 0; this->nnmap[8][3][1] = 0; this->nnmap[8][3][2] = 0; this->nnmap[8][3][3] = 5;//6;

     /* Atom 9 */
     this->nnmap[9][0][0] =  0; this->nnmap[9][0][1] = 0; this->nnmap[9][0][2] = 0; this->nnmap[9][0][3] = 1;//2;
     this->nnmap[9][1][0] =  0; this->nnmap[9][1][1] = 0; this->nnmap[9][1][2] = 0; this->nnmap[9][1][3] = 6;//7;
     this->nnmap[9][2][0] = -1; this->nnmap[9][2][1] = 0; this->nnmap[9][2][2] = 0; this->nnmap[9][2][3] = 3;//4;
     this->nnmap[9][3][0] =  0; this->nnmap[9][3][1] = 0; this->nnmap[9][3][2] = 0; this->nnmap[9][3][3] = 4;//5;

     /* Atom 10 */
     this->nnmap[10][0][0] =  0; this->nnmap[10][0][1] = 0; this->nnmap[10][0][2] = 0; this->nnmap[10][0][3] = 2;//3;
     this->nnmap[10][1][0] =  0; this->nnmap[10][1][1] = 1; this->nnmap[10][1][2] = 0; this->nnmap[10][1][3] = 4;//5;
     this->nnmap[10][2][0] = -1; this->nnmap[10][2][1] = 0; this->nnmap[10][2][2] = 0; this->nnmap[10][2][3] = 0;//1;
     this->nnmap[10][3][0] =  0; this->nnmap[10][3][1] = 0; this->nnmap[10][3][2] = 0; this->nnmap[10][3][3] = 6;//7;

     /* Atom 11 */
     this->nnmap[11][0][0] = 0; this->nnmap[11][0][1] = 0; this->nnmap[11][0][2] = 0; this->nnmap[11][0][3] = 3;//4;
     this->nnmap[11][1][0] = 0; this->nnmap[11][1][1] = 0; this->nnmap[11][1][2] = 0; this->nnmap[11][1][3] = 5;//6;
     this->nnmap[11][2][0] = 0; this->nnmap[11][2][1] = 0; this->nnmap[11][2][2] = 0; this->nnmap[11][2][3] = 1;//2;
     this->nnmap[11][3][0] = 0; this->nnmap[11][3][1] = 0; this->nnmap[11][3][2] = 0; this->nnmap[11][3][3] = 7;//8;

     /* Atom 12 */
     this->nnmap[12][0][0] =  0; this->nnmap[12][0][1] =  0; this->nnmap[12][0][2] =  0; this->nnmap[12][0][3] = 4;//5;
     this->nnmap[12][1][0] = -1; this->nnmap[12][1][1] =  0; this->nnmap[12][1][2] = -1; this->nnmap[12][1][3] = 3;//4;
     this->nnmap[12][2][0] = -1; this->nnmap[12][2][1] =  0; this->nnmap[12][2][2] =  0; this->nnmap[12][2][3] = 7;//8;
     this->nnmap[12][3][0] = -1; this->nnmap[12][3][1] = -1; this->nnmap[12][3][2] = -1; this->nnmap[12][3][3] = 0;//1;

     /* Atom 13 */
     this->nnmap[13][0][0] = 0; this->nnmap[13][0][1] = 0; this->nnmap[13][0][2] =  0; this->nnmap[13][0][3] = 5;//6;
     this->nnmap[13][1][0] = 0; this->nnmap[13][1][1] = 0; this->nnmap[13][1][2] = -1; this->nnmap[13][1][3] = 2;//3;
     this->nnmap[13][2][0] = 0; this->nnmap[13][2][1] = 0; this->nnmap[13][2][2] =  0; this->nnmap[13][2][3] = 6;//7;
     this->nnmap[13][3][0] = 0; this->nnmap[13][3][1] = 0; this->nnmap[13][3][2] = -1; this->nnmap[13][3][3] = 1;//2;

     /* Atom 14 */
     this->nnmap[14][0][0] =  0; this->nnmap[14][0][1] = 0; this->nnmap[14][0][2] =  0; this->nnmap[14][0][3] = 6;//7;
     this->nnmap[14][1][0] = -1; this->nnmap[14][1][1] = 0; this->nnmap[14][1][2] = -1; this->nnmap[14][1][3] = 0;//1;
     this->nnmap[14][2][0] = -1; this->nnmap[14][2][1] = 0; this->nnmap[14][2][2] =  0; this->nnmap[14][2][3] = 5;//6;
     this->nnmap[14][3][0] = -1; this->nnmap[14][3][1] = 0; this->nnmap[14][3][2] = -1; this->nnmap[14][3][3] = 3;//4;

     /* Atom 15 */
     this->nnmap[15][0][0] = 0; this->nnmap[15][0][1] =  0; this->nnmap[15][0][2] =  0; this->nnmap[15][0][3] = 7;//8;
     this->nnmap[15][1][0] = 0; this->nnmap[15][1][1] =  0; this->nnmap[15][1][2] = -1; this->nnmap[15][1][3] = 1;//2;
     this->nnmap[15][2][0] = 0; this->nnmap[15][2][1] =  0; this->nnmap[15][2][2] =  0; this->nnmap[15][2][3] = 4;//5;
     this->nnmap[15][3][0] = 0; this->nnmap[15][3][1] = -1; this->nnmap[15][3][2] = -1; this->nnmap[15][3][3] = 2;//3;

/* Rotation matrix from new unit cell to standard Zincblende basis. */
    this->T_rot = Rmatrix(3,3);
    this->T_rot[0][0]=  sqrt(1.0/2.0); this->T_rot[0][1]=  sqrt(1.0/2.0); this->T_rot[0][2]= 0.0;
    this->T_rot[1][0]= -sqrt(1.0/2.0); this->T_rot[1][1]=  sqrt(1.0/2.0); this->T_rot[1][2]= 0.0;
    this->T_rot[2][0]=            0.0; this->T_rot[2][1]=            0.0; this->T_rot[2][2]= 1.0;

    this->a_latt_rot_norm = Rvectr(3);
    this->a_latt_rot_norm[0] = sqrt(2.0);
    this->a_latt_rot_norm[1] = sqrt(2.0);
    this->a_latt_rot_norm[2] = 1.0;
/* End Rotation matrix from new unit cell to standard Zincblende basis. */

       break;
/* ---------- End New unit cell (ZB_110) ----------- */


    /* ---------- New unit cell (ZB_110_small) ----------- */
    case(Dev_struct::ZB_110_small):
     this->nnmap = I3tensor(0, N_MaxAtomCell-1, 0, N_MaxNeighbor-1, 0, 3);
     /* Atom 0 */
     this->nnmap[0][0][0] = 1; this->nnmap[0][0][1] = 1; this->nnmap[0][0][2] = 1; this->nnmap[0][0][3] = 3;
     this->nnmap[0][1][0] = 0; this->nnmap[0][1][1] = 0; this->nnmap[0][1][2] = 0; this->nnmap[0][1][3] = 2;
     this->nnmap[0][2][0] = 1; this->nnmap[0][2][1] = 0; this->nnmap[0][2][2] = 1; this->nnmap[0][2][3] = 3;
     this->nnmap[0][3][0] = 1; this->nnmap[0][3][1] = 0; this->nnmap[0][3][2] = 0; this->nnmap[0][3][3] = 2;

     /* Atom 1 */
     this->nnmap[1][0][0] = 0; this->nnmap[1][0][1] = 0; this->nnmap[1][0][2] = 0; this->nnmap[1][0][3] = 2; 
     this->nnmap[1][1][0] = 0; this->nnmap[1][1][1] = 0; this->nnmap[1][1][2] = 0; this->nnmap[1][1][3] = 3;
     this->nnmap[1][2][0] = 0; this->nnmap[1][2][1] =-1; this->nnmap[1][2][2] = 0; this->nnmap[1][2][3] = 2;
     this->nnmap[1][3][0] = 1; this->nnmap[1][3][1] = 0; this->nnmap[1][3][2] = 0; this->nnmap[1][3][3] = 3;

     /* Atom 2 */
     this->nnmap[2][0][0] = 0; this->nnmap[2][0][1] = 0; this->nnmap[2][0][2] = 0; this->nnmap[2][0][3] = 1;
     this->nnmap[2][1][0] = 0; this->nnmap[2][1][1] = 0; this->nnmap[2][1][2] = 0; this->nnmap[2][1][3] = 0;
     this->nnmap[2][2][0] = 0; this->nnmap[2][2][1] = 1; this->nnmap[2][2][2] = 0; this->nnmap[2][2][3] = 1;
     this->nnmap[2][3][0] =-1; this->nnmap[2][3][1] = 0; this->nnmap[2][3][2] = 0; this->nnmap[2][3][3] = 0;

     /* Atom 3 */
     this->nnmap[3][0][0] =-1; this->nnmap[3][0][1] =-1; this->nnmap[3][0][2] =-1; this->nnmap[3][0][3] = 0;
     this->nnmap[3][1][0] = 0; this->nnmap[3][1][1] = 0; this->nnmap[3][1][2] = 0; this->nnmap[3][1][3] = 1;
     this->nnmap[3][2][0] =-1; this->nnmap[3][2][1] = 0; this->nnmap[3][2][2] =-1; this->nnmap[3][2][3] = 0;
     this->nnmap[3][3][0] =-1; this->nnmap[3][3][1] = 0; this->nnmap[3][3][2] = 0; this->nnmap[3][3][3] = 1;

/* Rotation matrix from new unit cell to standard Zincblende basis. */
    this->T_rot = Rmatrix(3,3);
    this->T_rot[0][0]=  sqrt(1.0/2.0); this->T_rot[0][1]=  sqrt(1.0/2.0); this->T_rot[0][2]= 0.0;
    this->T_rot[1][0]= -sqrt(1.0/2.0); this->T_rot[1][1]=  sqrt(1.0/2.0); this->T_rot[1][2]= 0.0;
    this->T_rot[2][0]=            0.0; this->T_rot[2][1]=            0.0; this->T_rot[2][2]= 1.0;

    this->a_latt_rot_norm = Rvectr(3);
    this->a_latt_rot_norm[0] = 1/sqrt(2.0);
    this->a_latt_rot_norm[1] = 1/sqrt(2.0);
    this->a_latt_rot_norm[2] = 1.0;
/* End Rotation matrix from new unit cell to standard Zincblende basis. */

       break;
/* ---------- End New unit cell (ZB_110_small) ----------- */

    /* ---------- New unit cell (ZB_111) ----------- */
    case(Dev_struct::ZB_111):
     this->nnmap = I3tensor(0, N_MaxAtomCell-1, 0, N_MaxNeighbor-1, 0, 3);
     /* Atom 0 */
     this->nnmap[0][0][0] =  0; this->nnmap[0][0][1] =  0; this->nnmap[0][0][2] =  0; this->nnmap[0][0][3] =  6;
     this->nnmap[0][1][0] =  0; this->nnmap[0][1][1] =  0; this->nnmap[0][1][2] =  1; this->nnmap[0][1][3] = 11;
     this->nnmap[0][2][0] =  1; this->nnmap[0][2][1] =  0; this->nnmap[0][2][2] =  1; this->nnmap[0][2][3] = 10;
     this->nnmap[0][3][0] =  0; this->nnmap[0][3][1] =  1; this->nnmap[0][3][2] =  1; this->nnmap[0][3][3] = 11;

     /* Atom 1 */
     this->nnmap[1][0][0] =  0; this->nnmap[1][0][1] =  0; this->nnmap[1][0][2] =  0; this->nnmap[1][0][3] =  7;
     this->nnmap[1][1][0] =  0; this->nnmap[1][1][1] = -1; this->nnmap[1][1][2] =  1; this->nnmap[1][1][3] = 10;
     this->nnmap[1][2][0] =  0; this->nnmap[1][2][1] =  0; this->nnmap[1][2][2] =  1; this->nnmap[1][2][3] = 11;
     this->nnmap[1][3][0] =  0; this->nnmap[1][3][1] =  0; this->nnmap[1][3][2] =  1; this->nnmap[1][3][3] = 10;

     /* Atom 2 */
     this->nnmap[2][0][0] =  0; this->nnmap[2][0][1] =  0; this->nnmap[2][0][2] =  0; this->nnmap[2][0][3] =  8;
     this->nnmap[2][1][0] =  0; this->nnmap[2][1][1] =  0; this->nnmap[2][1][2] =  0; this->nnmap[2][1][3] =  7;
     this->nnmap[2][2][0] =  0; this->nnmap[2][2][1] =  0; this->nnmap[2][2][2] =  0; this->nnmap[2][2][3] =  6;
     this->nnmap[2][3][0] =  0; this->nnmap[2][3][1] =  1; this->nnmap[2][3][2] =  0; this->nnmap[2][3][3] =  7;

     /* Atom 3 */
     this->nnmap[3][0][0] =  0; this->nnmap[3][0][1] =  0; this->nnmap[3][0][2] =  0; this->nnmap[3][0][3] =  9;
     this->nnmap[3][1][0] = -1; this->nnmap[3][1][1] = -1; this->nnmap[3][1][2] =  0; this->nnmap[3][1][3] =  6;
     this->nnmap[3][2][0] =  0; this->nnmap[3][2][1] =  0; this->nnmap[3][2][2] =  0; this->nnmap[3][2][3] =  7;
     this->nnmap[3][3][0] = -1; this->nnmap[3][3][1] =  0; this->nnmap[3][3][2] =  0; this->nnmap[3][3][3] =  6;

     /* Atom 4 */
     this->nnmap[4][0][0] =  0; this->nnmap[4][0][1] =  0; this->nnmap[4][0][2] =  0; this->nnmap[4][0][3] = 10;
     this->nnmap[4][1][0] =  0; this->nnmap[4][1][1] =  0; this->nnmap[4][1][2] =  0; this->nnmap[4][1][3] = 9;
     this->nnmap[4][2][0] =  0; this->nnmap[4][2][1] =  0; this->nnmap[4][2][2] =  0; this->nnmap[4][2][3] = 8;
     this->nnmap[4][3][0] =  0; this->nnmap[4][3][1] =  1; this->nnmap[4][3][2] =  0; this->nnmap[4][3][3] = 9;

     /* Atom 5 */
     this->nnmap[5][0][0] =  0; this->nnmap[5][0][1] =  0; this->nnmap[5][0][2] =  0; this->nnmap[5][0][3] = 11;
     this->nnmap[5][1][0] =  0; this->nnmap[5][1][1] = -1; this->nnmap[5][1][2] =  0; this->nnmap[5][1][3] = 8;
     this->nnmap[5][2][0] =  1; this->nnmap[5][2][1] =  0; this->nnmap[5][2][2] =  0; this->nnmap[5][2][3] = 9;
     this->nnmap[5][3][0] =  0; this->nnmap[5][3][1] =  0; this->nnmap[5][3][2] =  0; this->nnmap[5][3][3] = 8;

     /* Atom 6 */
     this->nnmap[6][0][0] =  0; this->nnmap[6][0][1] =  0; this->nnmap[6][0][2] =  0; this->nnmap[6][0][3] = 0;
     this->nnmap[6][1][0] =  1; this->nnmap[6][1][1] =  1; this->nnmap[6][1][2] =  0; this->nnmap[6][1][3] = 3;
     this->nnmap[6][2][0] =  0; this->nnmap[6][2][1] =  0; this->nnmap[6][2][2] =  0; this->nnmap[6][2][3] = 2;
     this->nnmap[6][3][0] =  1; this->nnmap[6][3][1] =  0; this->nnmap[6][3][2] =  0; this->nnmap[6][3][3] = 3;
     /* Atom 7 */
     this->nnmap[7][0][0] =  0; this->nnmap[7][0][1] =  0; this->nnmap[7][0][2] =  0; this->nnmap[7][0][3] = 1;
     this->nnmap[7][1][0] =  0; this->nnmap[7][1][1] =  0; this->nnmap[7][1][2] =  0; this->nnmap[7][1][3] = 2;
     this->nnmap[7][2][0] =  0; this->nnmap[7][2][1] =  0; this->nnmap[7][2][2] =  0; this->nnmap[7][2][3] = 3;
     this->nnmap[7][3][0] =  0; this->nnmap[7][3][1] = -1; this->nnmap[7][3][2] =  0; this->nnmap[7][3][3] = 2;

     /* Atom 8 */
     this->nnmap[8][0][0] =  0; this->nnmap[8][0][1] =  0; this->nnmap[8][0][2] =  0; this->nnmap[8][0][3] = 2;
     this->nnmap[8][1][0] =  0; this->nnmap[8][1][1] =  1; this->nnmap[8][1][2] =  0; this->nnmap[8][1][3] = 5;
     this->nnmap[8][2][0] =  0; this->nnmap[8][2][1] =  0; this->nnmap[8][2][2] =  0; this->nnmap[8][2][3] = 4;
     this->nnmap[8][3][0] =  0; this->nnmap[8][3][1] =  0; this->nnmap[8][3][2] =  0; this->nnmap[8][3][3] = 5;

     /* Atom 9 */
     this->nnmap[9][0][0] =  0; this->nnmap[9][0][1] =  0; this->nnmap[9][0][2] =  0; this->nnmap[9][0][3] = 3;
     this->nnmap[9][1][0] =  0; this->nnmap[9][1][1] =  0; this->nnmap[9][1][2] =  0; this->nnmap[9][1][3] = 4;
     this->nnmap[9][2][0] = -1; this->nnmap[9][2][1] =  0; this->nnmap[9][2][2] =  0; this->nnmap[9][2][3] = 5;
     this->nnmap[9][3][0] =  0; this->nnmap[9][3][1] = -1; this->nnmap[9][3][2] =  0; this->nnmap[9][3][3] = 4;

     /* Atom 10 */
     this->nnmap[10][0][0] =  0; this->nnmap[10][0][1] =  0; this->nnmap[10][0][2] =  0; this->nnmap[10][0][3] = 4;
     this->nnmap[10][1][0] =  0; this->nnmap[10][1][1] =  1; this->nnmap[10][1][2] = -1; this->nnmap[10][1][3] = 1;
     this->nnmap[10][2][0] = -1; this->nnmap[10][2][1] =  0; this->nnmap[10][2][2] = -1; this->nnmap[10][2][3] = 0;
     this->nnmap[10][3][0] =  0; this->nnmap[10][3][1] =  0; this->nnmap[10][3][2] = -1; this->nnmap[10][3][3] = 1;
     /* Atom 11 */
     this->nnmap[11][0][0] =  0; this->nnmap[11][0][1] =  0; this->nnmap[11][0][2] =  0; this->nnmap[11][0][3] = 5;
     this->nnmap[11][1][0] =  0; this->nnmap[11][1][1] =  0; this->nnmap[11][1][2] = -1; this->nnmap[11][1][3] = 0;
     this->nnmap[11][2][0] =  0; this->nnmap[11][2][1] =  0; this->nnmap[11][2][2] = -1; this->nnmap[11][2][3] = 1;
     this->nnmap[11][3][0] =  0; this->nnmap[11][3][1] = -1; this->nnmap[11][3][2] = -1; this->nnmap[11][3][3] = 0;

/* Rotation matrix from new unit cell to standard Zincblende basis. */
    this->T_rot = Rmatrix(3,3);
    this->T_rot[0][0]=  sqrt(1.0/6.0); this->T_rot[0][1]=  sqrt(1.0/6.0); this->T_rot[0][2]= -sqrt(2.0/3.0);
    this->T_rot[1][0]= -sqrt(1.0/2.0); this->T_rot[1][1]=  sqrt(1.0/2.0); this->T_rot[1][2]=            0.0;
    this->T_rot[2][0]=  sqrt(1.0/3.0); this->T_rot[2][1]=  sqrt(1.0/3.0); this->T_rot[2][2]=  sqrt(1.0/3.0);
    this->a_latt_rot_norm = Rvectr(3);
    this->a_latt_rot_norm[0] = sqrt(3.0/2.0);
    this->a_latt_rot_norm[1] = sqrt(1.0/2.0);
    this->a_latt_rot_norm[2] = sqrt(3.0);
/* End Rotation matrix from new unit cell to standard Zincblende basis. */

       break;
/* ---------- End New unit cell (ZB_111) ----------- */

    /* ---------- New unit cell (CellUserDef) ----------- */
    case(Dev_struct::CellUserDef):
     this->nnmap = I3tensor(0, N_MaxAtomCell-1, 0, N_MaxNeighbor-1, 0, 3);

    /* Following code calculates nnmap for this unit cell */
    FILE *fp_nnmap_out;
    fp_nnmap_out = fopen("nnmap_test.txt","w");

    double **nnmapCellAtom_2_ijk;
    nnmapCellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);

    printf("\n\nnnmapCellAtom_2_ijk = CellAtom_2_ijk + 1.0\n");
    for(int i=0;i<this->N_MaxAtomCell;i++) {
        printf("%d  ",i);
        for(int j=0;j<3;j++) {
            nnmapCellAtom_2_ijk[i][j] = this->CellAtom_2_ijk[i][j] + 1.0;///this->a_latt_rot_norm[j];
            printf("%12.10le  ",nnmapCellAtom_2_ijk[i][j]);
        }
    printf("\n");
    }
    printf("\n\n\n");

    /* relative vectors of neighbors in original co-ordinate system */
    double **nbr_a_2_c,**nbr_c_2_a, **nbr;
    nbr_a_2_c = Rmatrix(this->N_MaxNeighbor,3);
    nbr_c_2_a = Rmatrix(this->N_MaxNeighbor,3);
    nbr       = Rmatrix(this->N_MaxNeighbor,3);
    /* Hardcode lattice constants for now */
    nbr_a_2_c[0][0] = -1.0/4.0; nbr_a_2_c[0][1] = -1.0/4.0; nbr_a_2_c[0][2] = -1.0/4.0;
    nbr_a_2_c[1][0] =  1.0/4.0; nbr_a_2_c[1][1] = -1.0/4.0; nbr_a_2_c[1][2] =  1.0/4.0;
    nbr_a_2_c[2][0] =  1.0/4.0; nbr_a_2_c[2][1] =  1.0/4.0; nbr_a_2_c[2][2] = -1.0/4.0;
    nbr_a_2_c[3][0] = -1.0/4.0; nbr_a_2_c[3][1] =  1.0/4.0; nbr_a_2_c[3][2] =  1.0/4.0;

    /* Rotate relative neighbor vectors nbr_a_2_c' = T*nbr_a_2_c' */
    double **nbr_temp;
    nbr_temp = Rmatrix(this->N_MaxNeighbor,3);
    for(int i=0; i<this->N_MaxNeighbor; i++){
      for(int j=0; j<3; j++){
        nbr_temp[i][j]=0.0;
         for(int k=0; k<3; k++){
           nbr_temp[i][j]+=this->T_rot[j][k]*nbr_a_2_c[i][k];
         }
      }
    }

    for(int i=0;i<this->N_MaxNeighbor;i++)
        for(int j=0;j<3;j++)
            nbr_a_2_c[i][j] = nbr_temp[i][j];

    /* Scale by unit cell vectors */
    for(int i=0;i<this->N_MaxNeighbor;i++)
        for(int j=0;j<3;j++)
               nbr_a_2_c[i][j] = nbr_a_2_c[i][j]/(this->a_latt_rot_norm[j]/this->a_latt_unstr);

    for(int i=0;i<this->N_MaxNeighbor;i++)
        for(int j=0;j<3;j++)
            nbr_c_2_a[i][j] = -nbr_a_2_c[i][j];


    fp_nnmap_out = fopen("nnmap_test.txt","w");

    double *r0, *rv, *rv_temp, *r_cell_nbr;
    r0 = Rvectr(3);
    rv = Rvectr(3);
    rv_temp = Rvectr(3);
    r_cell_nbr = Rvectr(3);

    for(int m=0;m<this->N_MaxAtomCell;m++)
       {
        for(int j=0;j<3;j++)
            r0[j] = nnmapCellAtom_2_ijk[m][j];
//        fprintf(fp_nnmap_out,"%0.2d  %6.2f %6.2f %6.2f\n",m,r0[0],r0[1],r0[2]);
        fprintf(fp_nnmap_out,"Atom %0.2d: \n",m);
        if(m<this->N_MaxAtomCell/2)
          {
           for(int i=0;i<this->N_MaxNeighbor;i++)
               for(int j=0;j<3;j++)
                   nbr[i][j] = nbr_a_2_c[i][j];
          }
        else
          {
           for(int i=0;i<this->N_MaxNeighbor;i++)
               for(int j=0;j<3;j++)
                   nbr[i][j] = nbr_c_2_a[i][j];
          }
        for(int n=0;n<this->N_MaxNeighbor;n++)
           {
            for(int j=0;j<3;j++) {
                rv[j] = r0[j] + nbr[n][j];
                if(fabs(rv[j] - int(rv[j]))<1e-5) rv[j] = int(rv[j]);
               }
//            fprintf(fp_nnmap_out,"%4.2d  %12.10le %12.10le %12.10le   ",n,rv[0],rv[1],rv[2]);
            for(int j=0;j<3;j++)
               {
                r_cell_nbr[j] = ceil(rv[j]) - 1.0;
                this->nnmap[m][n][j] = (int)r_cell_nbr[j];
               }
            int indx = (this->N_MaxAtomCell/2)*m + n;
//            fprintf(fp_nnmap_out,"%4.2d  %6.2f %6.2f %6.2f   ",indx,nnmap[n][0][m],nnmap[n][1][m],nnmap[n][2][m]);

            for(int j=0;j<3;j++)
                r_cell_nbr[j] = rv[j] - r_cell_nbr[j];

            int nbr_number = -1;
            for(int l=0;l<this->N_MaxAtomCell;l++)
               {
                if(fabs(nnmapCellAtom_2_ijk[l][0]-r_cell_nbr[0])<1e-3 &&
                   fabs(nnmapCellAtom_2_ijk[l][1]-r_cell_nbr[1])<1e-3 &&
                   fabs(nnmapCellAtom_2_ijk[l][2]-r_cell_nbr[2])<1e-3)
                  { nbr_number = l;
                    break;
                  }
               }
            this->nnmap[m][n][3] = nbr_number;
//          fprintf(fp_nnmap_out,"%6.2f %6.2f %6.2f  %d\n",r_cell_nbr[0],r_cell_nbr[1],r_cell_nbr[2],nbr_number);
            fprintf(fp_nnmap_out,"%4.1d  %4.1d %4.1d %4.1d %4.1d  \n",indx,this->nnmap[m][n][0],this->nnmap[m][n][1],this->nnmap[m][n][2],this->nnmap[m][n][3]);
           }
           fprintf(fp_nnmap_out,"\n");
       }

    fclose(fp_nnmap_out);

    rm_rvectr(&r0);
    rm_rvectr(&rv);
    rm_rvectr(&rv_temp);
    rm_rvectr(&r_cell_nbr);

    rm_rmatrix(&nnmapCellAtom_2_ijk);
    rm_rmatrix(&nbr_a_2_c);
    rm_rmatrix(&nbr_c_2_a);

    /* Normalize  lenghts of unit cell in xyz directions (lx, ly, lz) from pdb_in.txt */
    for(int j=0;j<3;j++)
        this->a_latt_rot_norm[j] = this->a_latt_rot_norm[j]/this->a_latt_unstr;


       break;
/* ---------- End New unit cell (CellUserDef) ----------- */

    default:
       break;
    }
}


/* This parses a hyperchem file for atom locations, number
   of neighbors and total number of atoms */
#define DELIMITER_HC2 " \n"
void Cell::InitializeForHyperChem(const char* inputfile)
{
   int i, j, nnbr, nbr_indx, nbr_rm;
   char line[100], item[50], ldesc[20], *tail, *token;

   FILE *fhin;
   char *filename=NULL;
   
   filename = n3d_strdup_n(inputfile);
   n3d_FileTypeSet(&filename, "nd_hin", 1);
   fhin = fopen(filename,"r");
   
   this->N_MaxAtomCell=0;
   this->N_MaxNeighbor = 0;

   while(fgets(line,100,fhin)) {
      /* Read line descriptor */
      token = strtok(line,DELIMITER_HC2);     
      strcpy(ldesc,token);                  
      if ( strcmp(ldesc,"atom") == 0 ) {
	 
	 /* Second column is atom number, save this one */      
	 token = strtok(NULL,DELIMITER_HC2);
	 strcpy(item,token);
         this->N_MaxAtomCell = strtol(item,&tail,0);

	 /* Skip over columns 3-10  */
	 for ( i=3; i <= 10; i++ )
	    token = strtok(NULL,DELIMITER_HC2);
	 
	 /* Get number of neighbors for this atom */
	 token = strtok(NULL,DELIMITER_HC2);
	 strcpy(item,token);
	 nnbr = strtol(item,&tail,0);
	 
	 /* Update maximum number of neighbors in structure */
	 if ( nnbr > this->N_MaxNeighbor )
            this->N_MaxNeighbor = nnbr;
      }
   }
   fclose(fhin);
   
   /* Initialize arrays */
   this->N_Neighbor = Ivectr(this->N_MaxAtomCell);
   this->CellAtom_2_ijk = Rmatrix(this->N_MaxAtomCell,3);
   this->nnmap = I3tensor(0,this->N_MaxAtomCell-1,0,this->N_MaxNeighbor-1,0,3);
   
   fhin = fopen(filename,"r");
   while(fgets(line,100,fhin)) {
      /* Read line descriptor */
      token = strtok(line,DELIMITER_HC2);     
      strcpy(ldesc,token);                  
      if ( strcmp(ldesc,"atom") == 0 ) {
	 
	 /* Second column is atom number, save this one */      
	 token = strtok(NULL,DELIMITER_HC2);
	 strcpy(item,token);
	 j = strtol(item,&tail,0)-1;   
	 
	 /* 3rd column is junk read and forget */
	 token = strtok(NULL,DELIMITER_HC2);
	 
	 /* 4th column is used by geo_struct */
	 token = strtok(NULL,DELIMITER_HC2);
	 
	 /* 5th, 6th, and 7th columns are junk read and forget */
	 token = strtok(NULL,DELIMITER_HC2);
	 token = strtok(NULL,DELIMITER_HC2);
	 token = strtok(NULL,DELIMITER_HC2);
	 
	 /* 8th column is the x-coordinate in A, convert to nm */
	 token = strtok(NULL,DELIMITER_HC2);
	 strcpy(item,token);
	 this->CellAtom_2_ijk[j][0] = 0.1*strtod(item,&tail);
	 
	 /* 9th column is the y-coordinate in A, convert to nm */
	 token = strtok(NULL,DELIMITER_HC2);
	 strcpy(item,token);
	 this->CellAtom_2_ijk[j][1] = 0.1*strtod(item,&tail);
	 
	 /* 10th column is the z-coordinate in A, convert to nm */
	 token = strtok(NULL,DELIMITER_HC2);
	 strcpy(item,token);
	 this->CellAtom_2_ijk[j][2] = 0.1*strtod(item,&tail);    
	 
	 /* 11th column is the number of neighbors */
	 token = strtok(NULL,DELIMITER_HC2);
	 strcpy(item,token);
	 this->N_Neighbor[j] = strtol(item,&tail,0);
	 
	 /* Loop through neighbors and get their index */
	 nbr_rm=0;
	 for ( i=0; i < this->N_Neighbor[j]; i++ ) {
	    
	    /* Get neighbor index */
	    token = strtok(NULL,DELIMITER_HC2);
	    strcpy(item,token);
	    nbr_indx = strtol(item,&tail,0)-1;
	    /* If the neighbor exists add it to the list.
	       Else remove it. */
	    if (nbr_indx < this->N_MaxAtomCell)
	       this->nnmap[j][i-nbr_rm][3] = nbr_indx;
	    else
	       nbr_rm++;
	    
	    /* Next column is bond type, forget it for now */
	    token = strtok(NULL,DELIMITER_HC2);
	 }
      }
   }
   str_free(filename);
}     


void Cell::Deallocate() 
{
   if (this->N_Neighbor){
      rm_ivectr(&this->N_Neighbor);
      this->N_Neighbor=NULL;
   }

   if (this->CellAtom_2_ijk){
      rm_rmatrix(&this->CellAtom_2_ijk);
      this->CellAtom_2_ijk=NULL;
   }
 
   if ( this->nnmap){
      rm_i3tensor(& this->nnmap,0);
      this->nnmap=NULL;
   }

   if ( this->SP3_hybrid){
      rm_r3tensor(& this->SP3_hybrid,0);
      this->SP3_hybrid=NULL;
   }
}


ostream& operator << (ostream& o, const Cell& c)
{
   if (mpi_n3d_id) return o;

   o << "=================================================\n"
     << "Max number of atoms per cell = " << c.N_MaxAtomCell << "\n"
     << "Max number of neighbors =      " << c.N_MaxNeighbor << "\n"
     << "=================================================\n";

   o << "Max number of neighbors (atom):  \n";
   for (int Za=0; Za < c.N_MaxAtomCell; Za++)
      o << "   atom=" << Za << "  " << c.N_Neighbor[Za] << "\n";
   
   o << "Equilibrium atomic positions within cell:\n";
   o.setf(ios::fixed, ios::floatfield);
   o.precision(4);
   for (int Za=0; Za < c.N_MaxAtomCell; Za++)
      o << "   atom=" << Za << "  "
        << setw(9) << c.CellAtom_2_ijk[Za][0]
        << setw(9) << c.CellAtom_2_ijk[Za][1]
        << setw(9) << c.CellAtom_2_ijk[Za][2]
        << "\n";

   // set back to default
   o.precision(0);
   o.setf(ios::floatfield);

   o << "Information regarding neighbor's cell:\n";
   for (int Za=0; Za < c.N_MaxAtomCell; Za++)
   for (int Zn=0; Zn < c.N_MaxNeighbor; Zn++)
      o << "   atom=" << Za 
        << "   nbr_indx=" << Zn
        << "  " << setw(2) << c.nnmap[Za][Zn][0]
        << "  " << setw(2) << c.nnmap[Za][Zn][1]
        << "  " << setw(2) << c.nnmap[Za][Zn][2]
        << "  " << setw(2) << c.nnmap[Za][Zn][3]
        << "\n";
   
   o << flush;

   return o;
}
