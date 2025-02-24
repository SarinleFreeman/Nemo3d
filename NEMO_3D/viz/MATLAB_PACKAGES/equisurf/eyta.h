#undef NULL
#define NULL (void *)0
#define malloc(size) mxCalloc(size, 1)
#define calloc mxCalloc
#define free mxFree

extern char konf[128][4][6];

struct point {
    double x,y,z;
};

struct polygon {
    int sidor;
    struct polygon *next;
    struct point horn[6];
};

void *allok(int size, char *first_free);

struct polygon **listcontour3d(double *d, double u, int in, int jn, int kn,
                              int i1, int i2, int j1, int j2, int k1, int k2,
			      struct polygon **poly);

