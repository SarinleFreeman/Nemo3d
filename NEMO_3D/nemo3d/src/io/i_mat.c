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
$Header: /repo/nemo3d/src/io/i_mat.c,v 1.16 2006/01/17 02:59:40 gekco Exp $ 
*****************************************************************************/

/*
 * T_HEADER_FILE{i_mat3d_comp.c}
 * This file contains functions that deal with the completeness
 * checking and computation of material parameters.
 */

#include "i_mat.h"

#ifdef CYGWIN
/* 
   gekco:  Jan 16, 2006
   somehow the CYGWIN installation chokes on the use of strcasecmp
   so I am using the simple strcmp instead.  If an input deck does 
   use some strange identifications of atoms on this system we will have
   trouble.  Really this cygwin installation is for development of a few 
   people only anyhow!
*/
#define strcasecmp strcmp
#endif

static int new_atom_count=0;
static char **new_atom_list=NULL;

/* Return an integer value corresponding to a particular atom.
   If the atom is not in the list associate the character string with 
   a particularly assigned number */
int str_to_atomID(const char * atomname)
{
   int result=-1;
  
   if ( (strcasecmp(atomname,"NULL") == 0) || 
        (strcasecmp(atomname,"XX") == 0) )
      result = EL_None;
   else if (!strcasecmp(atomname,"B"))
      result = EL_B;
   else if (!strcasecmp(atomname,"C"))
      result = EL_C;
   else if (!strcasecmp(atomname,"N"))
      result = EL_N;
   else if (!strcasecmp(atomname,"O"))
      result = EL_O;
   else if (!strcasecmp(atomname,"F"))
      result = EL_F;
   else if (!strcasecmp(atomname,"Al"))
      result = EL_Al;
   else if (!strcasecmp(atomname,"Si"))
      result = EL_Si;
   else if (!strcasecmp(atomname,"P"))
      result = EL_P;
   else if (!strcasecmp(atomname,"Zn"))
      result = EL_Zn;
   else if (!strcasecmp(atomname,"Ga"))
      result = EL_Ga;
   else if (!strcasecmp(atomname,"Ge"))
      result = EL_Ge;
   else if (!strcasecmp(atomname,"As"))
      result = EL_As;
   else if (!strcasecmp(atomname,"Cd"))
      result = EL_Cd;
   else if (!strcasecmp(atomname,"In"))
      result = EL_In;
   else if (!strcasecmp(atomname,"Sn"))
      result = EL_Sn;
   else if (!strcasecmp(atomname,"Sb"))
      result = EL_Sb;
   else if (!strcasecmp(atomname,"H"))
      result = EL_H;
   else if (!strcasecmp(atomname,"S"))
      result = EL_S;
   else if (!strcasecmp(atomname,"Se"))
      result = EL_Se;
   else{
      new_atom_count = n3d_StringListNum(new_atom_list);
      for (int i=0; i<new_atom_count; i++) {
         if (!strcasecmp(atomname,new_atom_list[i]))
            return EL_EXPLICIT+i;
      }
      if (EL_EXPLICIT+new_atom_count<EL_TOTAL){
         /* Can add another user defined material */
         new_atom_list = n3d_StringListAdd(&new_atom_list,atomname);
         result = EL_EXPLICIT+new_atom_count;
         new_atom_count++;
         return result;
      } else {
         die("Ran out of space to create the user defined material \"%s\"\n",
             atomname);
         return -1;
      }
   }
   
   return result;
}

/* Return a pointer to a character string that expresses 
   the atom corresponding to the input index. */
char * atomID_to_str(int mat_index)
{
   if ( mat_index == EL_None )
      return "NULL";
   else if (mat_index == EL_B )
      return "B";
   else if (mat_index == EL_C)
      return "C";
   else if (mat_index == EL_N)
      return "N";
   else if (mat_index == EL_O)
      return "O";
   else if (mat_index == EL_F)
      return "F";
   else if (mat_index == EL_Al)
      return "Al";
   else if (mat_index == EL_Si)
      return "Si";
   else if (mat_index == EL_P)
      return "P";
   else if (mat_index == EL_Zn)
      return "Zn";
   else if (mat_index == EL_Ga)
      return "Ga";
   else if (mat_index == EL_Ge)
      return "Ge";
   else if (mat_index == EL_As)
      return "As";
   else if (mat_index == EL_Cd)
      return "Cd";
   else if (mat_index == EL_In)
      return "In";
   else if (mat_index == EL_Sn)
      return "Sn";
   else if (mat_index == EL_Sb)
      return "Sb";
   else if (mat_index == EL_H)
      return "H";
   else if (mat_index == EL_S)
      return "S";
   else if (mat_index == EL_Se)
      return "Se";
   else if (mat_index-EL_EXPLICIT < new_atom_count)
      return new_atom_list[mat_index-EL_EXPLICIT];
   else{
      die("Unidentifyable material index: \"%d\"\n",mat_index);
   }
   
   return NULL;
}


// given a vector of material structs, find the element that
// satisfies the properties of the subsequent arguements
static Material_struct* find_Material_struct(
                             vector<Material_struct>& m, 
                             string material, 
                             string band_model,
                             string comp1_mat, real comp1_val,
                             string comp2_mat, real comp2_val,
                             real tkel, bool get_strain ) 
{
   for (vector<Material_struct>::iterator mPtr = m.begin(); 
        mPtr != m.end(); mPtr++) {
      
//#define _PRINT_DEBUG
#ifdef _PRINT_DEBUG
      cout << "--------------------------" << endl;
      cout << material << "  " << mPtr->material << endl;
      cout << tkel << "  " << mPtr->tkel << endl;
      cout << comp1_mat  << "  " <<  mPtr->comp1_mat << endl;
      cout << comp1_val << "  " << mPtr->xval << endl;
      cout << comp2_mat  << "  " <<  mPtr->comp2_mat << endl;
      cout << comp2_val << "  " << mPtr->yval << endl;
      cout << band_model  << "  " <<  mPtr->band_model << endl;
#endif

      if (material!=mPtr->material) continue;
      if (tkel!=mPtr->tkel && tkel>=0)	continue;
      if (comp1_mat != mPtr->comp1_mat) continue;
      if ((comp1_val!=mPtr->xval) && comp1_val>=0) continue;
      if (comp2_mat != mPtr->comp2_mat) continue;
      if ((comp2_val!=mPtr->yval) && comp2_val>=0) continue;
      if (band_model != mPtr->band_model) continue;
      
      return &*mPtr;
   }

   // could not find material in list; return NULL
   return NULL;
   // access_database();
}


static void printMaterialList(vector<MaterialHandle> &sml) 
{
   if (mpi_n3d_id) return;  // only the master checks

   printf("\nSummarizing Material List:\n");
   for (vector<MaterialHandle>::iterator p=sml.begin(); p!=sml.end(); p++) {
      printf("   Material #%d:  %s + %s", 
             p->matid3d, p->cation.c_str(), p->anion.c_str());
      printf("  cations={ ");
      for (vector<int>::iterator c=p->cation_id.begin(); 
           c!=p->cation_id.end(); c++)
         printf("%d ", *c);
      printf("}  anions={ ");
      for (vector<int>::iterator a=p->anion_id.begin(); 
           a!=p->anion_id.end(); a++)
         printf("%d ", *a);
      printf("}\n");
   }
   fflush(stdout);
}


void QD_struct::set_sMatList()
{
   real tkel = 300.0;
   real comp1_val=0.0, comp2_val=0.0;
   
   string bandModel;
   if      (this->opt.Dev.band_model==Dev_struct::Bands_20_sp3d5ss_spin)
      bandModel = string("Bands_20_sp3d5ss_spin");
   else if (this->opt.Dev.band_model==Dev_struct::Bands_10_sp3ss_spin)
      bandModel = string("Bands_10_sp3ss_spin");
   else if (this->opt.Dev.band_model==Dev_struct::Bands_10_sp3d5ss_nospin)
      bandModel = string("Bands_10_sp3d5ss_nospin");
   else if (this->opt.Dev.band_model==Dev_struct::Bands_1_s_nospin)
      bandModel = string("Bands_1_s_nospin");
   else
      die("ERROR in QD_struct::set_sMatList: Found an unmatched BandModel\n");

   cout_master << "\nMatching shape materials with material database:\n";
   
   /* Loop through the general shapes and check for material 
      parameter list completeness. */
   int Zmat_valid = 0;
   for (vector<Shape_struct>::iterator forme=this->opt.Dev.ShapeList.begin();
        forme != this->opt.Dev.ShapeList.end(); forme++) {
      
      cout_master << "   Discovered a shape:\n";
      
      // count the number of polar molecules
      int Npolar=0;
      for (vector<int>::iterator c=forme->pol.begin(); 
           c!=forme->pol.end(); c++) {
         if (*c != 0) Npolar++;
      }

      bool is_monoatomic = (forme->mat.size()==1);
      bool is_alloy = (Npolar!=2);
      bool is_vca = forme->AlloyTreatment==Shape_struct::VCA;
      
      // case:  monoatomic
      if (is_monoatomic) {
         cout_master << "      Shape material is monoatomic.\n";
         
         Material_struct* m = 
            find_Material_struct(this->opt.material_list,
                                 forme->mat[0],
                                 bandModel,
                                 string("NULL"), 
                                 0.0, 
                                 string("NULL"), 
                                 0.0,
                                 tkel,
                                 this->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian
                                 );
         
         if (!m)  die("Could not find requested material %s with band model %s.", forme->mat[0].c_str(), bandModel.c_str());
         
         // search sMatList for a match with this Material_struct m
         vector<MaterialHandle>::iterator h;
         for (h = this->sMatList.begin(); 
              h != this->sMatList.end() && h->hndl != m; h++);
         
         if (h == this->sMatList.end()) {
            MaterialHandle mh(m);
            
            mh.anion_id = vector<int>(1);
            mh.anion = forme->mat[0];
            mh.anion_id[0]  = str_to_atomID(forme->mat[0].c_str());
            mh.anion_weight = 1.0;

            mh.cation_id = vector<int>(1);
            mh.cation = forme->mat[0];
            mh.cation_id[0] = str_to_atomID(forme->mat[0].c_str());
            mh.cation_weight = 1.0;
            
            mh.matid3d = Zmat_valid;
            forme->matid3d = Zmat_valid;
            Zmat_valid++;
            
            this->sMatList.push_back(mh);
         }
         else {
            forme->matid3d = h->matid3d;
         }
      } 
      //   case:  binary material or alloy
      else {
         int alloy_indx = -1; // pointer to a MaterialHandle
         
         if (is_alloy) {
            // shape material is a not a simple binary but an alloy
            char mat_val_str[20];
            string material, anion_composite, cation_composite;
            string comp1_mat, comp2_mat;
            vector<int> cation_list, anion_list;

            // set comp1_mat, comp1_val, 
            // cation_composite, material,cation_list
            for (unsigned int cation=0; cation<forme->pol.size(); cation++) {
               if (forme->pol[cation] < 0) continue;
               
               if (cation_list.size()==0) {
                  comp1_mat = forme->mat[cation];
                  comp1_val = forme->mol[cation];
               } 
               else if (cation_list.size()==1) {
               //else if (cation_list.size()>=1)
                  sprintf(mat_val_str,"%1.3f",comp1_val);
                  cation_composite += string(mat_val_str);
               }
               else {
                  die("Material cannot have more than two cations!\n");
               }
               cation_composite += forme->mat[cation];
               material += forme->mat[cation];
               cation_list.push_back(str_to_atomID(forme->mat[cation].c_str()));
            }
            if (cation_list.size()==1){
               comp1_mat=string("NULL"); comp1_val=0.0; 
            }
            
            // set comp2_mat, comp2_val, 
            // anion_composite, material, anion_list
            for (unsigned int anion=0; anion<forme->pol.size(); anion++) {
               if (forme->pol[anion] > 0) continue;
               
               if (anion_list.size()==0){
                  comp2_mat = forme->mat[anion];
                  comp2_val = forme->mol[anion];
               } 
               else if (anion_list.size()==1) {
                 //else if (anion_list.size()>=1)
                  sprintf(mat_val_str,"%1.3f",comp2_val);
                  anion_composite += string(mat_val_str);
               }
               else {
                  die("Material cannot have more than two anions!\n");
               }
               anion_composite += forme->mat[anion];
               material += forme->mat[anion];
               anion_list.push_back(str_to_atomID(forme->mat[anion].c_str()));
            }
            if (anion_list.size()==1){
               comp2_mat=string("NULL"); comp2_val=0.0; 
            }

            // look for a match among material list; 
            // if not found, m is NULL
            Material_struct* m = 
               find_Material_struct(this->opt.material_list,
                                    material,
                                    bandModel,
                                    comp1_mat,
                                    comp1_val, 
                                    comp2_mat,
                                    comp2_val, 
                                    tkel,
                                    this->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian
                                    );
            
            // search sMatList for a match with this Material_struct
            vector<MaterialHandle>::iterator h;
            for (h = this->sMatList.begin(); 
                 h != this->sMatList.end() && h->hndl != m; h++);
               
            
            if (!m || h==this->sMatList.end()) {
               // no match was found => append to MaterialHandle list
               MaterialHandle mh(m);
               
               mh.anion = anion_composite;
               mh.cation = cation_composite;
               mh.cation_weight = comp1_val;
               mh.anion_weight = comp2_val;
               
               if (forme->AlloyTreatment==Shape_struct::VCA) {
                  mh.anion_id     = vector<int>(1);
                  mh.cation_id    = vector<int>(1);
                  mh.anion_id[0]  = str_to_atomID(mh.anion.c_str());
                  mh.cation_id[0] = str_to_atomID(mh.cation.c_str());
               } else {
                  mh.anion_id     = vector<int>(anion_list.size());
                  mh.cation_id    = vector<int>(cation_list.size());
                  for (unsigned int i=0;i<anion_list.size();i++)   
                     mh.anion_id[i]=-1;
                  for (unsigned int i=0;i<cation_list.size();i++)  
                     mh.cation_id[i]=-1;

                  for (unsigned int i=0;i<anion_list.size();i++)   
                     mh.anion_id[i]=anion_list[i];
                  for (unsigned int i=0;i<cation_list.size();i++)  
                     mh.cation_id[i]=cation_list[i];
               }
               
               alloy_indx = Zmat_valid;
               
               mh.matid3d = Zmat_valid;
               forme->matid3d = Zmat_valid;
               
               this->sMatList.push_back(mh);
               Zmat_valid++;
               
               cout_master << "      Shape material is an alloy (newly found):  "
                           << mh.cation << " + " << mh.anion
                           << "  (matid3d=" << mh.matid3d << ")\n";
            }
            else {  // we have a match
               forme->matid3d = h->matid3d;
               alloy_indx = h->matid3d;

               cout_master<< "      Shape material is an alloy (previously found):  "
                          << h->cation << " + " << h->anion
                          << "   (matid3d=" << h->matid3d << ")\n";
            }
         }
         
         /* Make up the binary materials that are in this shape if 
            they are needed for an atomistic representation of an 
            alloy, or if this shape only consists of a binary. */
         if ( !(is_vca && is_alloy) ) {
          
            for (unsigned int cation=0; cation<forme->pol.size(); cation++) {
               if (forme->pol[cation] < 0) continue;  // ignore anions
               
               // Find a corresponding anion.
               for (unsigned int anion=0; anion<forme->pol.size(); anion++) {
                  if (forme->pol[anion] > 0) continue;  // ignore cations
                  
                  string material = (forme->mat[cation]!=forme->mat[anion] 
                                     ? forme->mat[cation] + forme->mat[anion]
                                     : forme->mat[cation]);
                  
                  /* leave out any VCA treatment and look for the 
                    simple pairs only. */
                  Material_struct* m = 
                     find_Material_struct(this->opt.material_list,
                                          material,
                                          bandModel,
                                          string("NULL"), 
                                          0.0, 
                                          string("NULL"), 
                                          0.0,
                                          tkel,
                                          this->opt.ExecParam.ElCalc.include_strain_in_Hamiltonian
                                          );
                  
                  if (!m)  
                     die("Could not find requested material %s with band model %s.", material.c_str(), bandModel.c_str());
                  
                  // search sMatList for a match with this Material_struct m
                  vector<MaterialHandle>::iterator h;
                  for (h = this->sMatList.begin(); 
                       h != this->sMatList.end() && h->hndl != m; h++);
               
                  if (h == this->sMatList.end()) { // no match; must be new

                     MaterialHandle mh(m);
                     
                     mh.anion = forme->mat[anion];
                     mh.cation = forme->mat[cation];

                     
                     if (alloy_indx < 0){
                        // shape material is a new binary
                        forme->matid3d=Zmat_valid;
                        cout_master 
                           << "      Shape material is a binary (newly found):  "
                           << forme->mat[cation]  << " + " << forme->mat[anion]
                           << "  (matid3d=" << forme->matid3d << ")\n";
                     }
                     
                     mh.anion_id  = vector<int>(1);
                     mh.anion_id[0]  = str_to_atomID(mh.anion.c_str());
                     mh.cation_id = vector<int>(1);
                     mh.cation_id[0] = str_to_atomID(mh.cation.c_str());
                     mh.cation_weight = mh.anion_weight = 1.0;
                     
                     mh.matid3d=Zmat_valid;
                     Zmat_valid++;
                     this->sMatList.push_back(mh);

                     cout_master
                        << "      Adding a new binary material:  "
                        << m->material
                        << "  (matid3d=" << mh.matid3d << ")\n";
                  }
                  // we have a match with a previously found binary
                  else {
                     if (alloy_indx < 0){
                        // shape material is a previously found binary
                        forme->matid3d=h->matid3d;
                        cout_master 
                           <<"      Shape material is a binary (previously found):  "
                           << forme->mat[cation]  << " + " << forme->mat[anion]
                           << "  (matid3d=" << forme->matid3d << ")\n";
                     }
                  }
               }
            }
         }
      }
   }

   printMaterialList(this->sMatList);
}


