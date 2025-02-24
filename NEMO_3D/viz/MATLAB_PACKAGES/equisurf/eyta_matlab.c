#include "LICENSE.h"
#include "mex.h"
#include "eyta.h"

void sort_poly(struct polygon *polylista, mxArray *plhs[])
{
    int antal[4] = {0,0,0,0};
    struct polygon *ptr = polylista;
    double *ptr_x[4], *ptr_y[4], *ptr_z[4], *tabs[4];
    int i, index;

    while (ptr) {
	antal[ptr -> sidor - 3]++;
	ptr = ptr->next;
    }
    for (i = 0; i < 4; i++) {
	plhs[i] = mxCreateDoubleMatrix(i + 3, antal[i] * 3, mxREAL);
	tabs[i] = mxGetPr(plhs[i]);
    }
    if ((antal[0] && !tabs[0]) || (antal[1] && !tabs[1]) || 
	(antal[2] && !tabs[2]) || (antal[3] && !tabs[3]))
	mexErrMsgTxt("Out of memory when allocating output matrices");

    for (i = 0; i < 4; i++) {
	ptr_x[i] = tabs[i];
	ptr_y[i] = tabs[i] + (i + 3) * antal[i];
	ptr_z[i] = tabs[i] + (i + 3) * 2 *antal[i];
    }
    ptr = polylista;
    while (ptr) {
	index = ptr -> sidor - 3;
	for (i = 0; i < index + 3; i++) {
	    *ptr_x[index]++ = ptr -> horn[i].x;
	    *ptr_y[index]++ = ptr -> horn[i].y;
	    *ptr_z[index]++ = ptr -> horn[i].z;
	}
	ptr = ptr->next;
    }
}

#define Mat prhs[0]
#define Lim prhs[1]
#define Const prhs[2]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i1, i2, in, j1, j2, jn, k1, k2, kn, i, ant_konst;
    double *gr;
    int ndim;
    const int *dim_array;
    struct polygon *polylista = NULL, **next_ptr = &polylista;

    if (nrhs != 3 || nlhs != 4)
	mexErrMsgTxt("Usage: [p3,p4,p5,p6] = polygons(3Dmatrix, "
	    "limits, constant)");
    ndim = mxGetNumberOfDimensions(Mat);
    if (ndim != 3)
	mexErrMsgTxt("The first argument should be a 3D matrix");
    dim_array = mxGetDimensions(Mat);
/* meshgrid swaps the order of the two first indices by convention */
    in = dim_array[1];
    jn = dim_array[0];
    kn = dim_array[2];

    if (mxGetM(Lim) == 0) {
	i1 = 0;
	i2 = in - 1;
	j1 = 0;
	j2 = jn - 1;
	k1 = 0;
	k2 = kn - 1;
    } else if (mxGetM(Lim) * mxGetN(Lim) == 6) {
	gr = mxGetPr(Lim);
	i1 = gr[0] - 1;
	i2 = gr[1] - 1;
	j1 = gr[2] - 1;
	j2 = gr[3] - 1;
	k1 = gr[4] - 1;
	k2 = gr[5] - 1;
    } else
	mexErrMsgTxt("The limits vector should have six elements or be empty");

    if (mxGetM(Const) == 0)
	mexErrMsgTxt("No constant supplied");
    if (i1 < 0 || i1 >= i2 || i2 >= in ||
	j1 < 0 || j1 >= j2 || j2 >= jn ||
	k1 < 0 || k1 >= k2 || k2 >= kn)
	mexErrMsgTxt("The limits are wrong");

    ant_konst = mxGetM(Const) * mxGetN(Const);

    for (i = 0; i < ant_konst; i++)
	next_ptr = listcontour3d(mxGetPr(Mat), mxGetPr(Const)[i],
	    in, jn, kn, i1, i2, j1, j2, k1, k2, next_ptr);

    if (!polylista)
	mexErrMsgTxt("No output");
    sort_poly(polylista, plhs);
    allok(-1, NULL);	
}


