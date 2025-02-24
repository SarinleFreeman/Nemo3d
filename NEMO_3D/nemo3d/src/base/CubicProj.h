#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "r2tensor.h"
#include "rvector.h"
#include "ivector.h"

#define system_error printf
#define EPS_UNITARY     1.0e-10 /* effective 0 for unitary test */
#define EPS             1.0e-10
#define A_SM_SC 1.0

int check_unitarity(double **Ur, double **Ui, int n, double eps);

void define_small_cell_of_nanostructure(int *n_nemo, int *n_nanostruct, int *NBASIS, int *periodic_y, int *periodic_z, int primitive_p_cell);

void gen_Rvecs_cubic(double **R_a, int Nrect[]);

void gen_SmCell_Gvecs(double **Gvecs,int *Nrect,double bz_max);

void gen_SmCell_kpGvecs(double **Gvecs,double *k,double **kpGvecs,
                           double bz_max, int Nc);

void projection_matrix(double **Gv, double **Rv, double **Ur, double **Ui, int Nc);

void project_probs(double **Ekmat, double **Ur, double **Ui, 
                            double vr[], double vi[], double wvr[], 
                            double wvi[], int NBASIS, int Norb, int Nc_nano, 
                            int Negy, int eix, double k[], int *n_nanostruct);

