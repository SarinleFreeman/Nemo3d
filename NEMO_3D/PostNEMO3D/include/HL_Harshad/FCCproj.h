#ifndef FCCPROJ
#define FCCPROJ

void pos_itoa(int n, char s[]);

void full_filename(char template1[], int num, char suffix[], char fname[]);

void hoover(FILE *fp);

void clean_to_eq(FILE *fp);

void get_header(FILE *fp, int *Norb, int *Natom);

int maximum_states(FILE *fp);

void get_nemo3d_energies(FILE *fp, double emin, double emax, double *all_energies,
                                int max_states, int *negy, int *min_ix, int *max_ix);

void get_nemo3d_vector(FILE *fp, double vr[], double vi[], int Ncube, int Norb);

long get_index(int g,int s,int orb, int Norb);

void compute_tr_state(double *vr,double *vi,double *vr_tr,double *vi_tr,int Norb,int nrec);

void real_space_coords(double **av, int n[], double **Rv, int row);

void gen_Rvecs(double **R_a, int Nrect[]);

void gen_SmCell_k(double **Gv, int Nrect[]);

void projection_matrix(double **Gv, double **Rv, double **Ur, double **Ui, int Nc);

void U_dagger_v(double **Ur, double **Ui, double vr[], double vi[],
                                double uvr[], double uvi[], int start, int nth, int nn);

int check_unitarity(double **Ur, double **Ui, int n, double eps);

void project_probs(double **Ekmat, double **Ur, double **Ui, double vr[], double vi[],
			double wvr[], double wvi[], int Ncube, int Norb, int Nc, 
                        int Negy, int eix);

int fccproj(void);

#endif
