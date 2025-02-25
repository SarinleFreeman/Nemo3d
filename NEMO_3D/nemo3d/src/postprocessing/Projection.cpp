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

Written by:  Seungwon Lee

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/postprocessing/Projection.cpp,v 1.2 2004/10/27
19:43:10 gekco Exp $
*****************************************************************************/

#include "Projection.h"
#include "SimpleParser.h"

void Projection::get_data(complex *coef_e, complex *coef_h, int *atom_id,
                          int *shape_rank, int *nn, double *lat) {
  _cf_e = coef_e;
  _cf_h = coef_h;
  _atom = atom_id;
  _shape = shape_rank;
  _nn = nn;
  _lattice = lat;
}

void Projection::compute_weight_on_shape(int n1) {

  double weight[6];
  weight[0] = 0.0;
  weight[1] = 0.0;
  weight[2] = 0.0;
  weight[3] = 0.0;
  weight[4] = 0.0;
  weight[5] = 0.0;

  for (int g = g_min; g < g_max; g++) {
    double tmp_weight = 0.0;
    for (int spin = 0; spin < 2; spin++) {
      for (int orbital = 0; orbital < o_max; orbital++) {
        complex c1 = cf_e(n1, g, spin, orbital);
        tmp_weight += c1.r * c1.r + c1.i * c1.i;
      }
    }
    if (_shape[g] == 0)
      weight[0] += tmp_weight;
    else if (_shape[g] == 1)
      weight[1] += tmp_weight;
    else if (_shape[g] == 2)
      weight[2] += tmp_weight;
    else if (_shape[g] == 3)
      weight[3] += tmp_weight;
    else if (_shape[g] == 4)
      weight[4] += tmp_weight;
    else if (_shape[g] == 5)
      weight[5] += tmp_weight;
  }

  double Weight[6];
  MPI_Allreduce(weight, Weight, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  cout_master << n1 << " Wavefunction,  Weight on each shape 0, 1, 2, 3, ... "
              << Weight[0] << " " << Weight[1] << " " << Weight[2] << " "
              << Weight[3] << " " << Weight[4] << " " << Weight[5] << " "
              << endl;
}

void Projection::compute_electron_hole_alignment(int n1, int n2) {
  double X, Y, Z, x, y, z;

  x = 0.0;
  y = 0.0;
  z = 0.0;
  for (int g = g_min; g < g_max; g++) {
    double tmp_weight = 0.0;
    for (int spin = 0; spin < 2; spin++) {
      for (int orbital = 0; orbital < o_max; orbital++) {
        complex c1 = cf_e(n1, g, spin, orbital);
        tmp_weight += c1.r * c1.r + c1.i * c1.i;
      }
    }
    x += tmp_weight * lattice(g, 0);
    y += tmp_weight * lattice(g, 1);
    z += tmp_weight * lattice(g, 2);
  }
  MPI_Allreduce(&x, &X, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&y, &Y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&z, &Z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  cout_master << "Electron state, center of mass X, Y, Z: " << n1 << " " << X
              << " " << Y << " " << Z << endl;

  double x_div = 0.0, y_div = 0.0, z_div = 0.0, r_div = 0.0;
  for (int g = g_min; g < g_max; g++) {
    double tmp_weight = 0.0;
    for (int spin = 0; spin < 2; spin++) {
      for (int orbital = 0; orbital < o_max; orbital++) {
        complex c1 = cf_e(n1, g, spin, orbital);
        tmp_weight += c1.r * c1.r + c1.i * c1.i;
      }
    }
    x_div += tmp_weight * (lattice(g, 0) - X) * (lattice(g, 0) - X);
    y_div += tmp_weight * (lattice(g, 1) - Y) * (lattice(g, 1) - Y);
    z_div += tmp_weight * (lattice(g, 2) - Z) * (lattice(g, 2) - Z);
    r_div += tmp_weight * ((lattice(g, 0) - X) * (lattice(g, 0) - X) +
                           (lattice(g, 1) - Y) * (lattice(g, 1) - Y) +
                           (lattice(g, 2) - Z) * (lattice(g, 2) - Z));
  }

  double X_div, Y_div, Z_div, R_div;
  MPI_Allreduce(&x_div, &X_div, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&y_div, &Y_div, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&z_div, &Z_div, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&r_div, &R_div, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  cout_master << "Electron state, diviation R, X, Y, Z " << n1 << " "
              << sqrt(R_div) << " " << sqrt(X_div) << " " << sqrt(Y_div) << " "
              << sqrt(Z_div) << endl;

  x = 0.0;
  y = 0.0;
  z = 0.0;
  for (int g = g_min; g < g_max; g++) {
    double tmp_weight = 0.0;
    for (int spin = 0; spin < 2; spin++) {
      for (int orbital = 0; orbital < o_max; orbital++) {
        complex c2 = cf_h(n2, g, spin, orbital);
        tmp_weight += c2.r * c2.r + c2.i * c2.i;
      }
    }
    x += tmp_weight * lattice(g, 0);
    y += tmp_weight * lattice(g, 1);
    z += tmp_weight * lattice(g, 2);
  }
  MPI_Allreduce(&x, &X, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&y, &Y, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&z, &Z, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  cout_master << "Hole state, center of mass X, Y, Z: " << n1 << " " << X << " "
              << Y << " " << Z << endl;

  x_div = 0.0, y_div = 0.0, z_div = 0.0, r_div = 0.0;
  for (int g = g_min; g < g_max; g++) {
    double tmp_weight = 0.0;
    for (int spin = 0; spin < 2; spin++) {
      for (int orbital = 0; orbital < o_max; orbital++) {
        complex c2 = cf_h(n2, g, spin, orbital);
        tmp_weight += c2.r * c2.r + c2.i * c2.i;
      }
    }
    x_div += tmp_weight * (lattice(g, 0) - X) * (lattice(g, 0) - X);
    y_div += tmp_weight * (lattice(g, 1) - Y) * (lattice(g, 1) - Y);
    z_div += tmp_weight * (lattice(g, 2) - Z) * (lattice(g, 2) - Z);
    r_div += tmp_weight * ((lattice(g, 0) - X) * (lattice(g, 0) - X) +
                           (lattice(g, 1) - Y) * (lattice(g, 1) - Y) +
                           (lattice(g, 2) - Z) * (lattice(g, 2) - Z));
  }
  MPI_Allreduce(&x_div, &X_div, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&y_div, &Y_div, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&z_div, &Z_div, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&r_div, &R_div, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  cout_master << "Hole state, diviation R, X, Y, Z " << n2 << " " << sqrt(R_div)
              << " " << sqrt(X_div) << " " << sqrt(Y_div) << " " << sqrt(Z_div)
              << endl;
}

void Projection::compute_qd_wl_extent(int n1) {
  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  double *weight = new double[g_max - g_min];
  double *max_weight = new double[num_procs];
  int *max_atom = new int[num_procs];
  max_weight[my_id] = 0.0;
  max_atom[my_id] = 0;

  for (int g = g_min; g < g_max; g++) {
    weight[g - g_min] = 0.0;
    for (int spin = 0; spin < 2; spin++) {
      for (int orbital = 0; orbital < o_max; orbital++) {
        complex c1 = cf_e(n1, g, spin, orbital);
        weight[g - g_min] += c1.r * c1.r + c1.i * c1.i;
      }
    }
    if (weight[g - g_min] > max_weight[my_id]) {
      max_weight[my_id] = weight[g - g_min];
      max_atom[my_id] = g;
    }
  }

  cout << "My id, maximum atom, local maximum wf " << my_id << " "
       << max_atom[my_id] << " " << max_weight[my_id] << endl;

  for (int p = 0; p < num_procs; p++) {
    MPI_Bcast(&max_weight[p], 1, MPI_DOUBLE, p, MPI_COMM_WORLD);
    MPI_Bcast(&max_atom[p], 1, MPI_INT, p, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  double global_max_weight = -1000.0;
  int global_max_atom = 0;

  for (int p = 0; p < num_procs; p++) {
    if (max_weight[p] > global_max_weight) {
      global_max_weight = max_weight[p];
      global_max_atom = max_atom[p];
    }
  }

  cout_master << "Center atom location: " << lattice(global_max_atom, 0) << ", "
              << lattice(global_max_atom, 1) << ", "
              << lattice(global_max_atom, 2) << endl;
  double r_extent = 0.0;
  double xy_extent = 0.0;
  double z_extent = 0.0;
  double xy_WL_extent = 0.0;
  double xy_QD_extent = 0.0;
  double WL = 0.0;
  double QD = 0.0;
  double atom_WL = 0;
  double atom_QD = 0;
  for (int g = g_min; g < g_max; g++) {
    double distance = pow(lattice(g, 0) - lattice(global_max_atom, 0), 2.0);
    distance += pow(lattice(g, 1) - lattice(global_max_atom, 1), 2.0);
    double xy_distance = sqrt(distance);
    double z_distance = fabs(lattice(g, 2) - lattice(global_max_atom, 2));
    distance += pow(z_distance, 2.0);
    distance = sqrt(distance);
    r_extent += weight[g - g_min] * distance;
    xy_extent += weight[g - g_min] * xy_distance;
    z_extent += weight[g - g_min] * z_distance;
    if (_shape[g] == 3) {
      WL += weight[g - g_min];
      xy_WL_extent += weight[g - g_min] * xy_distance;
      atom_WL++;
    }
    if (_shape[g] == 4) {
      QD += weight[g - g_min];
      xy_QD_extent += weight[g - g_min] * xy_distance;
      atom_QD++;
    }
  }
  double R_extent, XY_extent, Z_extent, XY_WL_extent, XY_QD_extent;
  double Total_WL, Total_atom_WL, Total_QD, Total_atom_QD;
  MPI_Allreduce(&r_extent, &R_extent, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&xy_extent, &XY_extent, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&xy_WL_extent, &XY_WL_extent, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&xy_QD_extent, &XY_QD_extent, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&z_extent, &Z_extent, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&WL, &Total_WL, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&QD, &Total_QD, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&atom_WL, &Total_atom_WL, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&atom_QD, &Total_atom_QD, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  cout_master << "R, XY, Z: " << R_extent << " " << XY_extent << " " << Z_extent
              << endl;
  cout_master << "charge density in wetting layer, Number of Atom in WL: "
              << Total_WL << " " << Total_atom_WL << endl;
  cout_master << "charge density in quantum dot, Number of Atom in QD: "
              << Total_QD << " " << Total_atom_QD << endl;
  cout_master << "XY in WL: " << XY_WL_extent / Total_WL << endl;
  cout_master << "XY in QD: " << XY_QD_extent / Total_QD << endl;
}

void Projection::compute_extent(int n1) {
  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  double *weight = new double[g_max - g_min];
  double *max_weight = new double[num_procs];
  int *max_atom = new int[num_procs];
  max_weight[my_id] = 0.0;
  max_atom[my_id] = 0;

  for (int g = g_min; g < g_max; g++) {
    weight[g - g_min] = 0.0;
    for (int spin = 0; spin < 2; spin++) {
      for (int orbital = 0; orbital < o_max; orbital++) {
        complex c1 = cf_e(n1, g, spin, orbital);
        weight[g - g_min] += c1.r * c1.r + c1.i * c1.i;
      }
    }
    if (weight[g - g_min] > max_weight[my_id]) {
      max_weight[my_id] = weight[g - g_min];
      max_atom[my_id] = g;
    }
  }

  cout << "My id, maximum atom, local maximum wf " << my_id << " "
       << max_atom[my_id] << " " << max_weight[my_id] << endl;

  for (int p = 0; p < num_procs; p++) {
    MPI_Bcast(&max_weight[p], 1, MPI_DOUBLE, p, MPI_COMM_WORLD);
    MPI_Bcast(&max_atom[p], 1, MPI_INT, p, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  double global_max_weight = -1000.0;
  int global_max_atom = 0;

  for (int p = 0; p < num_procs; p++) {
    if (max_weight[p] > global_max_weight) {
      global_max_weight = max_weight[p];
      global_max_atom = max_atom[p];
    }
  }

  cout_master << "Center atom location: " << lattice(global_max_atom, 0) << ", "
              << lattice(global_max_atom, 1) << ", "
              << lattice(global_max_atom, 2) << endl;

  double r_extent = 0.0;
  double xy_extent = 0.0;
  double z_extent = 0.0;
  for (int g = g_min; g < g_max; g++) {
    double distance = pow(lattice(g, 0) - lattice(global_max_atom, 0), 2.0);
    distance += pow(lattice(g, 1) - lattice(global_max_atom, 1), 2.0);
    double xy_distance = sqrt(distance);
    double z_distance = fabs(lattice(g, 2) - lattice(global_max_atom, 2));
    distance += pow(z_distance, 2.0);
    distance = sqrt(distance);
    r_extent += weight[g - g_min] * distance;
    xy_extent += weight[g - g_min] * xy_distance;
    z_extent += weight[g - g_min] * z_distance;
  }
  double R_extent, XY_extent, Z_extent;
  MPI_Allreduce(&r_extent, &R_extent, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&xy_extent, &XY_extent, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&z_extent, &Z_extent, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  cout_master << "R, XY, Z: " << R_extent << " " << XY_extent << " " << Z_extent
              << endl;
}

void Projection::compute_overlap(int n1, int n2) {
  double overlap_r, overlap_i;
  overlap_r = 0.0;
  overlap_i = 0.0;
  for (int g = g_min; g < g_max; g++) {
    for (int spin = 0; spin < 2; spin++) {
      for (int orbital = 0; orbital < o_max; orbital++) {
        complex c1 = cf_e(n1, g, spin, orbital);
        complex c2 = cf_e(n2, g, spin, orbital);
        overlap_r += c1.r * c2.r + c1.i * c2.i;
        overlap_i += c1.r * c2.i - c1.i * c2.r;
      }
    }
  }
  double Overlap_r, Overlap_i;
  MPI_Allreduce(&overlap_r, &Overlap_r, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&overlap_i, &Overlap_i, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  cout << "Overlap between wave functions " << n1 << " and " << n2 << ": "
       << "( " << Overlap_r << ", " << Overlap_i << " )" << endl;
}

void Projection::compute_overlap_without_spin(int n1, int n2) {
  double overlap = 0.0;
  double fluctuation = 0.0;
  for (int g = g_min; g < g_max; g++) {
    for (int orbital = 0; orbital < o_max; orbital++) {
      complex c1_up = cf_e(n1, g, 0, orbital);
      complex c2_up = cf_e(n2, g, 0, orbital);
      complex c1_down = cf_e(n1, g, 1, orbital);
      complex c2_down = cf_e(n2, g, 1, orbital);
      double density1 = c1_up.r * c1_up.r + c1_up.i * c1_up.i;
      density1 += c1_down.r * c1_down.r + c1_down.i * c1_down.i;
      double density2 = c2_up.r * c2_up.r + c2_up.i * c2_up.i;
      density2 += c2_down.r * c2_down.r + c2_down.i * c2_down.i;
      overlap += sqrt(density1 * density2);
      fluctuation += fabs(density1 - density2);
    }
  }
  double Overlap;
  double Fluctuation;
  MPI_Allreduce(&overlap, &Overlap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&fluctuation, &Fluctuation, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  cout << "Overlap between wave functions without spin " << n1 << " and " << n2
       << ": " << Overlap << endl;
  cout << "Fluctuation between wave functions without spin " << n1 << " and "
       << n2 << ": " << Fluctuation / (numAtoms * o_max) << endl;
}

void Projection::print_amplitude(int n) {
  double *density = new double[g_max - g_min];
  for (int g = g_min; g < g_max; g++) {
    density[g - g_min] = 0.0;
    for (int orbital = 0; orbital < o_max; orbital++) {
      complex c_up = cf_e(n, g, 0, orbital);
      complex c_down = cf_e(n, g, 1, orbital);
      density[g - g_min] += c_up.r * c_up.r + c_up.i * c_up.i;
      density[g - g_min] += c_down.r * c_down.r + c_down.i * c_down.i;
    }
  }

  char *filename1 = "projection";
  char label[100];
  sprintf(label, "projection_along_z_%d", my_id);
  filename1 = n3d_strdup_n(filename1);
  n3d_FileTypeSet(&filename1, label, true);
  ofstream DensityOut(filename1);
  str_free(filename1);
  DensityOut.setf(ios::scientific);
  for (int g = g_min; g < g_max; g++) {
    if (lattice(g, 0) > 17.5 && lattice(g, 0) < 18.0 && lattice(g, 1) > 17.5 &&
        lattice(g, 1) < 18.0)
      DensityOut << lattice(g, 1) << "\t" << density[g - g_min] << "\t"
                 << _atom[g] << endl;
  }

  delete density;
}

void Projection::compute_projection(void) {
  double one_over_sqrt_three = 1.0 / sqrt(3.0);
  double sqrt_two_over_sqrt_three = sqrt(2.0) / sqrt(3.0);
  double one_over_sqrt_six = 1.0 / sqrt(6.0);

  for (int n = 0; n < e_max + h_max; n++) {

    double Tmp0 = 0.0, Tmp1 = 0.0, Tmp2 = 0.0, Tmp3 = 0.0, Tmp4 = 0.0;
    double Tmp5 = 0.0, Tmp6 = 0.0, Tmp7 = 0.0;
    double Tmp8 = 0.0, Tmp9 = 0.0, Tmp10 = 0.0;

    for (int g = g_min; g < g_max; g++) {

      for (int s = 0; s < 2; s++) {
        for (int orbital = 0; orbital < o_max; orbital++) {
          if (n < e_max)
            coef(s, orbital) = cf_e(n, g, s, orbital);
          else
            coef(s, orbital) = cf_h(n - e_max, g, s, orbital);
        }
      }

      Tmp0 += coef(0, 0).r * coef(0, 0).r + coef(0, 0).i * coef(0, 0).i;
      Tmp0 +=
          coef(1, 0).r * coef(1, 0).r + coef(1, 0).i * coef(1, 0).i; // s*-band

      Tmp1 += coef(0, 1).r * coef(0, 1).r + coef(0, 1).i * coef(0, 1).i;
      Tmp1 +=
          coef(1, 1).r * coef(1, 1).r + coef(1, 1).i * coef(1, 1).i; // s-band

      complex hh_u, hh_d, lh_u, lh_d, ss_u, ss_d;
      hh_u.r = coef(0, 2).r + coef(0, 3).i; // heavy-hole J=3/2, Jz=3/2
      hh_u.i = coef(0, 2).i - coef(0, 3).r; // heavy-hole J=3/2, Jz=3/2
      hh_d.r = coef(1, 2).r - coef(1, 3).i; // heavy-hole J=3/2, Jz=-3/2
      hh_d.i = coef(1, 2).i + coef(1, 3).r; // heavy-hole J=3/2, Jz=-3/2
      lh_u.r = sqrt_two_over_sqrt_three * coef(0, 4).r -
               one_over_sqrt_six *
                   (coef(1, 2).r + coef(1, 3).i); // light-hole J=3/2, Jz=1/2
      lh_u.i = sqrt_two_over_sqrt_three * coef(0, 4).i -
               one_over_sqrt_six *
                   (coef(1, 2).i - coef(1, 3).r); // light-hole J=3/2, Jz=1/2
      lh_d.r = sqrt_two_over_sqrt_three * coef(1, 4).r +
               one_over_sqrt_six *
                   (coef(0, 2).r - coef(0, 3).i); // light-hole J=3/2, Jz=-1/2
      lh_d.i = sqrt_two_over_sqrt_three * coef(1, 4).i +
               one_over_sqrt_six *
                   (coef(0, 2).i + coef(0, 3).r); // light-hole J=3/2, Jz=-1/2
      ss_u.r = one_over_sqrt_three * (coef(0, 4).r + coef(1, 2).r +
                                      coef(1, 3).i); // split-off J=1/2, Jz=1/2
      ss_u.i = one_over_sqrt_three * (coef(0, 4).i + coef(1, 2).i -
                                      coef(1, 3).r); // split-off J=1/2, Jz=1/2
      ss_d.r = one_over_sqrt_three * (coef(0, 2).r - coef(0, 3).i -
                                      coef(1, 4).r); // split-off J=1/2, Jz=-1/2
      ss_d.i = one_over_sqrt_three * (coef(0, 2).i + coef(0, 3).r -
                                      coef(1, 4).i); // split-off J=1/2, Jz=-1/2

      Tmp2 +=
          hh_u.r * hh_u.r + hh_u.i * hh_u.i + hh_d.r * hh_d.r + hh_d.i * hh_d.i;
      Tmp3 +=
          lh_u.r * lh_u.r + lh_u.i * lh_u.i + lh_d.r * lh_d.r + lh_d.i * lh_d.i;
      Tmp4 +=
          ss_u.r * ss_u.r + ss_u.i * ss_u.i + ss_d.r * ss_d.r + ss_d.i * ss_d.i;

      Tmp8 += coef(0, 2).r * coef(0, 2).r + coef(0, 2).i * coef(0, 2).i;
      Tmp8 +=
          coef(1, 2).r * coef(1, 2).r + coef(1, 2).i * coef(1, 2).i; //  px band
      Tmp9 += coef(0, 3).r * coef(0, 3).r + coef(0, 3).i * coef(0, 3).i;
      Tmp9 +=
          coef(1, 3).r * coef(1, 3).r + coef(1, 3).i * coef(1, 3).i; //  py band
      Tmp10 += coef(0, 4).r * coef(0, 4).r + coef(0, 4).i * coef(0, 4).i;
      Tmp10 +=
          coef(1, 4).r * coef(1, 4).r + coef(1, 4).i * coef(1, 4).i; // pz band

      if (o_max == 10) {
        for (int o1 = 5; o1 < 10; o1++) {
          Tmp5 += coef(0, o1).r * coef(0, o1).r + coef(0, o1).i * coef(0, o1).i;
          Tmp5 += coef(1, o1).r * coef(1, o1).r + coef(1, o1).i * coef(1, o1).i;
        }
      }

      double local_amplitude_up = 0.0;
      double local_amplitude_down = 0.0;
      for (int orb = 0; orb < o_max; orb++) {
        double amplitude_up = coef(0, orb).r * coef(0, orb).r +
                              coef(0, orb).i * coef(0, orb).i; // spin up
        double amplitude_down = coef(1, orb).r * coef(1, orb).r +
                                coef(1, orb).i * coef(1, orb).i; // spin down
        Tmp6 += amplitude_up;
        Tmp7 += amplitude_down;
        local_amplitude_up += amplitude_up;
        local_amplitude_down += amplitude_down;
      }
    }

    double projection[11];

    projection[0] = Tmp0; // total amplitude of s* orbital in envelope functions
    projection[1] = Tmp1; // total amplitude of s  orbital in envelope functions
    projection[2] =
        Tmp2 *
        0.5; // total amplitude of heavy-hole Bloch function (J=3/2, Jz=+/-3/2)
    projection[3] =
        Tmp3; // total amplitude of light-hole Bloch function (J=3/2, Jz=+/-1/2)
    projection[4] =
        Tmp4; // total amplitude of spin-split Bloch function (J=1/2, Jz=+/-1/2)
    projection[5] = Tmp5; // total amplitude of d orbital  in envelope functions
    projection[6] = Tmp6; // total amplitude of spin up  in envelope functions
    projection[7] = Tmp7; // total amplitude of spin down  in envelope functions
    projection[8] = Tmp8; // total amplitude of  px orbital
    projection[9] = Tmp9; // total amplitude of  py orbital
    projection[10] = Tmp10; // total amplitude of  pz orbital

    MPI_Barrier(MPI_COMM_WORLD);

    double total_projection[11];
    MPI_Allreduce(projection, total_projection, 11, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    if (mpi_n3d_id == 0) {
      if (n < e_max)
        cout << "Electron Level: " << setw(5) << n + 1 << endl;
      else
        cout << "Hole Level: " << setw(5) << n - e_max + 1 << endl;
      cout << "Projection onto s* band:          " << setw(10)
           << total_projection[0] << endl;
      cout << "Projection onto s  band:          " << setw(10)
           << total_projection[1] << endl;
      cout << "Projection onto heavy-hole band:  " << setw(10)
           << total_projection[2] << endl;
      cout << "Projection onto light-hole band:  " << setw(10)
           << total_projection[3] << endl;
      cout << "Projection onto split-off  band:  " << setw(10)
           << total_projection[4] << endl;
      cout << "Projection onto d  band:          " << setw(10)
           << total_projection[5] << endl;
      cout << "Projection onto spin up in z:     " << setw(10)
           << total_projection[6] << endl;
      cout << "Projection onto spin down in z:   " << setw(10)
           << total_projection[7] << endl;
      cout << "Projection onto px band:          " << setw(10)
           << total_projection[8] << endl;
      cout << "Projection onto py band:          " << setw(10)
           << total_projection[9] << endl;
      cout << "Projection onto pz band:          " << setw(10)
           << total_projection[10] << endl;
      cout << endl;
    }
  }
}
