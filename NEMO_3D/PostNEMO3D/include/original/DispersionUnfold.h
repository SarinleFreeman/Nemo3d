#ifndef DISPERSIONUNFOLD
#define DISPERSIONUNFOLD

void piksrtix(int n,double arr[],int index[]);
void normalize_all(double **Wr, double **Wi, int n);
void dir_cos(double r0[], double r1[], double dc[]);
double mx_elt_nn(int orb0, double r0[], int orb1, double r1[], double param[]);
void ham_mtx_elt(int orb0, int orb1, int r0[], double k[], int Nc[], double a_sm, double param[], double **Hr, double **Hi);
void make_hermitian(double **Hr, double **Hi, int n);
void init_ham_sc_sp3(double **Hr, double **Hi, double param[], double k[], double a_sm, int Nc[]);
void rect_to_polar_mat(double **Wr, double **Wi, double **Wm, double **Wp, int n);
void init_exp_mat(double **Ur, double **Ui, int Nc[]);
void coeffs_super2bulk_nm(double **Wr, double **Wi, double **Ur, double **Ui, double **Bcr, double **Bci, double k[], double **norm_mat, int Nc[]);
int unfold_dispersion(void);

#endif

