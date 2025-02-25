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
$Header: /repo/nemo3d/src/base/matsite_init.c,v 1.8 2004/12/02 12:09:56 marek
Exp $
*****************************************************************************/

#include "matsite_init.h"

void rm_matsite_init_struct(matsite_init_struct *mat_ptr) {
  matsite_init_struct mat = NULL;

  if (!mat_ptr)
    return;

  mat = *mat_ptr;

  rm_ivectrvectr_shellonly(&mat->cellatom);
  rm_rvectr(&mat->randweight);
  rm_ivectr(&mat->species_count);

  nml_free(mat);
  *mat_ptr = NULL;

  return;
}

matsite_init_struct make_matsite_init_struct(qd_struct d, const Cell &cell,
                                             int matid3d) {
  const int atoms_per_cell_max = cell.AtomsPerCellMax();

  matsite_init_struct mat =
      (matsite_init_struct)nml_calloc(1, sizeof(struct MatSite_Init_struct));

  mat->imat = matid3d;

  mat->cation = Ivectr(d->sMatList[matid3d].cation_id.size());
  for (int i = 0; i < nml_iv_extent(mat->cation); i++)
    mat->cation[i] = d->sMatList[matid3d].cation_id[i];

  mat->anion = Ivectr(d->sMatList[matid3d].anion_id.size());
  for (int i = 0; i < nml_iv_extent(mat->anion); i++)
    mat->anion[i] = d->sMatList[matid3d].anion_id[i];

  mat->cation_weight = d->sMatList[matid3d].cation_weight;
  mat->anion_weight = d->sMatList[matid3d].anion_weight;

  mat->cellatom = Ivectrvectr(atoms_per_cell_max);
  mat->randweight = Rvectr(atoms_per_cell_max);
  mat->species_count = Ivectr(EL_TOTAL);

  if (nml_iv_extent(mat->cation) > 1 || nml_iv_extent(mat->anion) > 1) {
    mat->random_alloy = 1;
  }

  if (atoms_per_cell_max == 8) {
    for (int i = 0; i < 4; i++) {
      add_elem_ivectrvectr(mat->cation, i, mat->cellatom);
      mat->randweight[i] = mat->cation_weight;
    }
    for (int i = 4; i < 8; i++) {
      add_elem_ivectrvectr(mat->anion, i, mat->cellatom);
      mat->randweight[i] = mat->anion_weight;
    }
  } else if (atoms_per_cell_max == 2) {
    add_elem_ivectrvectr(mat->cation, 0, mat->cellatom);
    mat->randweight[0] = mat->cation_weight;
    add_elem_ivectrvectr(mat->anion, 1, mat->cellatom);
    mat->randweight[1] = mat->anion_weight;
  } else
    add_elem_ivectrvectr(mat->anion, 0, mat->cellatom);

  return mat;
}

void print_species_count(qd_struct d, matsite_init_struct mat) {
  if (mpi_n3d_id != mpi_n3d_masterid)
    return;

  FILE *fp;
  char *filename = n3d_strdup_n(d->inputfile);
  n3d_FileTypeSet(&filename, "nd_species_count", TRUE);
  fp = fopen(filename, "w");
  str_free(filename);

  int species = 0;
  for (int i = 0; i < EL_TOTAL; i++) {
    if (mat->species_count[i] <= 0)
      continue;

    species += mat->species_count[i];
    printf("      * number of %s atoms = %d\n", atomID_to_str(i),
           mat->species_count[i]);
    fprintf(fp, "%s  %d\n", atomID_to_str(i), mat->species_count[i]);
  }

  printf("      * total number of atoms = %d\n", species);
  fclose(fp);
}
