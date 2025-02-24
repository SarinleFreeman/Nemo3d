#include "LICENSE.h"
#include <stdlib.h>
#include <math.h>
#include "eyta.h"
#include "mex.h"

const double eps = 0.000001;

struct key {
    unsigned int index : 24;
    unsigned int horn  :  4; 
    unsigned int list  :  4;
};

struct properties {
    double ambient;
    double diffuse;
    double specular;
    double spread;
};

static double *horn[4][3];

int compare(const struct key *k1, const struct key *k2)
{
    double diff;
    double x1, y1, z1, x2, y2, z2;

    x1 = horn[k1->list][0][k1->index];
    y1 = horn[k1->list][1][k1->index];
    z1 = horn[k1->list][2][k1->index];
    x2 = horn[k2->list][0][k2->index];
    y2 = horn[k2->list][1][k2->index];
    z2 = horn[k2->list][2][k2->index];

    diff = x1 - x2;
    if (diff < -eps)
	return -1;
    else if (diff > eps)
	return 1;
    else {
	diff = y1 - y2;
	if (diff < -eps)
	    return -1;
	else if (diff > eps)
	    return 1;
	else {
	    diff = z1 - z2;
	    if (diff < -eps)
		return -1;
	    else if (diff > eps)
		return 1;
	    else
		return 0;
	}
    }
}

double calc_colour(struct point *normal, struct point *view_dir,
    struct point *lightning_dir, struct properties *params, double vlangle)
{
    double vangle, langle, diffuse_light, specular_light;

    langle = normal->x * lightning_dir->x + normal->y * lightning_dir->y + 
	normal->z * lightning_dir->z;
    if (langle > 0) {		/* Framsidan belyst */
	diffuse_light = params->diffuse * langle;
	vangle = normal->x * view_dir->x + normal->y * view_dir->y +
        normal->z * view_dir->z;
	specular_light = params->specular * (2*vangle*langle - vlangle);
	if (specular_light < 0.0)
	    specular_light = 0.0;
    } else
	diffuse_light = specular_light = 0.0;
    return params->ambient + diffuse_light +pow(specular_light, params->spread);
} 

void average_normals(int first, int last, struct key keylist[], 
    double *horn[4][3], struct point *view_dir, struct point *normal)
{
    int i, list, current, previous, corner, next;
    struct point v1, v2, n;
    double norm;

    normal->x = normal->y = normal->z = 0.0;
    for (i = first; i <= last; i++) {
	list = keylist[i].list;
	current = keylist[i].index;
	corner = keylist[i].horn;
	if (corner == 0)
	    previous = current + list + 2;
	else
	    previous = current - 1;
	if (corner == list + 2)
	    next = current - list - 2;
	else
	    next = current + 1;
	v1.x = horn[list][0][current] - horn[list][0][previous];
        v2.x = horn[list][0][current] - horn[list][0][next];
        v1.y = horn[list][1][current] - horn[list][1][previous];
        v2.y = horn[list][1][current] - horn[list][1][next];
        v1.z = horn[list][2][current] - horn[list][2][previous];
        v2.z = horn[list][2][current] - horn[list][2][next];

	n.x = v1.y * v2.z - v1.z * v2.y;
	n.y = v1.z * v2.x - v1.x * v2.z;
	n.z = v1.x * v2.y - v1.y * v2.x;

	if (n.x * view_dir->x + n.y * view_dir->y + n.z * view_dir->z < 0) {
	    normal->x -= n.x;
	    normal->y -= n.y;
	    normal->z -= n.z;
	} else {
	    normal->x += n.x;
	    normal->y += n.y;
	    normal->z += n.z;
	}
    }
    norm = sqrt(normal->x * normal->x + normal->y * normal->y +
	normal->z * normal->z);
    if (norm > 0) {
	normal->x /= norm;
	normal->y /= norm;
	normal->z /= norm;
    }
}	
 
void shade(double *corner[4], int len[4], struct point *view_dir,
    struct point *lightning_dir, struct properties *params, mxArray *plhs[4])
{
    double *colour[4], corner_colour, vlangle;
    int corners = 0;
    struct key *keylist, *keyptr;
    struct point normal;
    int i, j, k, index, first, keyindex;

    for (i = 0; i < 4; i++) {
	plhs[i] = mxCreateDoubleMatrix(i+3, len[i], mxREAL);
	if (!plhs[i])
	    mexErrMsgTxt("Out of memory");
	colour[i] = mxGetPr(plhs[i]);
	corners += len[i] * (i+3);
	horn[i][0] = corner[i];
	horn[i][1] = corner[i] + len[i] * (i+3);
	horn[i][2] = corner[i] + 2 * len[i] * (i+3);
    }
    if (corners == 0)
	return;
    keylist = mxCalloc(corners, sizeof (struct key));
    if (!keylist)
	mexErrMsgTxt("Out of memory");

    keyptr = keylist;
    for (i = 0; i < 4; i++)
	for (j = 0, index = 0; j < len[i]; j++)
	    for (k = 0; k < i + 3; k++, index++) {
		keyptr->list = i;
		keyptr->horn = k;
		keyptr++ ->index = index;
	    }

    qsort(keylist, corners, sizeof(struct key), (int (*)(const void*, const void*))
 	compare);

    vlangle = view_dir->x * lightning_dir->x + view_dir->y * lightning_dir->y +
	view_dir->z * lightning_dir->z;
    for (keyindex = 0; keyindex < corners; keyindex++) {
	first = keyindex;
	while (keyindex < corners - 1 && compare(&keylist[first],
	    &keylist[keyindex + 1]) == 0)
	    keyindex++;			/* The loop variable is modified */
	average_normals(first, keyindex, keylist, horn, view_dir, &normal);
	corner_colour = calc_colour(&normal, view_dir, lightning_dir, params,
	    vlangle);
	for (i = first; i <= keyindex; i++)
	    colour[keylist[i].list][keylist[i].index] = corner_colour;
    }
    mxFree(keylist);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, len[4];
    double *corner[4];

    if (nrhs != 7 || nlhs != 4)
	mexErrMsgTxt("Usage: [c3,c4,c5,c6] = interp(p3, p4, p5, p6,"
	    " view_dir, source_dir, params)");
    if (mxGetM(prhs[4])*mxGetN(prhs[4]) != 3)
	mexErrMsgTxt("The view direction should be a 3D unit vector");
    if (mxGetM(prhs[5])*mxGetN(prhs[5]) != 3)
        mexErrMsgTxt("The source direction should be a 3D unit vector");
    if (mxGetM(prhs[6])*mxGetN(prhs[6]) != 4)
	mexErrMsgTxt("The parameter vector should have four elements");

    for (i = 0; i < 4; i++) {
	corner[i] = mxGetPr(prhs[i]);
	len[i] = mxGetN(prhs[i]) / 3;
    }

    shade(corner, len, (struct point *)mxGetPr(prhs[4]),
	(struct point *)mxGetPr(prhs[5]), (struct properties *)mxGetPr(prhs[6]),
	plhs);
} 
