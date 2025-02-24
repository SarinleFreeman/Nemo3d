#ifndef CUBICGUESS
#define CUBICGUESS

int comp_dir(double k[], double k_unit[], double eps);

void band_guess(double energies[], double **Ekmat, double control[], int ns,
                double **begymat, double **bsdevmat, double **bprobmat,
                int nb[], int kix);

int cubicguess(void);

#endif

