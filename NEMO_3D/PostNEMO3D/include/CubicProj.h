#ifndef CUBICPROJ
#define CUBICPROJ

void get_nemo3d_energies_cubic(FILE *fp, double emin, double emax, double *all_energies,
                                int max_states, int *negy, int *min_ix, int *max_ix,
                                double *k);   

void get_nemo3d_vector_cubic(FILE *fp, double vr[], double vi[]);              

void gen_Rvecs_cubic(double **R_a, int Nrect[]);

void gen_SmCell_kpGvecs_cubic(double **kpGvecs, double *k, int *Nrect, double bz_max);

void projection_matrix(double **Gv, double **Rv, double **Ur, double **Ui, int Nc);

void U_dagger_v(double **Ur, double **Ui, double vr[], double vi[],
                                double uvr[], double uvi[], int start, int nth, int nn);

int check_unitarity(double **Ur, double **Ui, int n, double eps);

void project_probs_cubic(double **Ekmat, double **Ur, double **Ui, double vr[], double vi[],
				double wvr[], double wvi[], int NBASIS, int Norb, int Nc_nano, int Negy, int eix);

void define_small_cell_of_nanostructure(int *n_nemo, int *n_nanostruct, int *NBASIS, int *periodic_y, int *periodic_z);

int cubicproj(void);

#endif

