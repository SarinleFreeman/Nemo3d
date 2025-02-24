#ifndef FCCGUESS
#define FCCGUESS 

int comp_dirFCC(double k[], double k_unit[], double eps);

void band_guess(double energies[], double **Ekmat, double control[], int ns,
                double **begymat, double **bsdevmat, double **bprobmat,
                int nb[], int kix);

int fccguess(void);

#endif
