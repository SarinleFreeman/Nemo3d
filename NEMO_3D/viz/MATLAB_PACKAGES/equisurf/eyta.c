#include "LICENSE.h"
#include "mex.h"
#include "eyta.h"

#define ANT_AREOR 1000
#define AREA_SIZE 16384
#define INTERP(i, f1, f2) (a = (f1), i + (a - u)/(a - (f2)))

extern char konf[128][4][6];

/* Custom memory allocation routine.  */
/* Makes it unnessary to free each allocated item (a single call with size equal
   -1 frees everything.) Good for mex file use to avoid memory leaks.
*/
void *allok(int size, char *first_free)
{
    static int akt_area;
    static char *area[ANT_AREOR];
    int i;

    if (size < 0) {
	for (i = 0; i < ANT_AREOR; i++)
	    if (area[i])
		free(area[i]);
	return NULL;
    }
    if (first_free == NULL) {
	for (i = 0; i < ANT_AREOR; i++)
	    area[i] = NULL;
        akt_area = 0;
	area[akt_area] = malloc(AREA_SIZE);
	if (!area[akt_area])
	    return NULL;
	first_free = area[akt_area];
    }
    if (AREA_SIZE - (first_free - area[akt_area]) < size) {
	akt_area++;
	if (akt_area >= ANT_AREOR || !(area[akt_area] = malloc(AREA_SIZE)))
	    return NULL;
	first_free = area[akt_area];
    }
    return first_free;
}

struct polygon **listcontour3d(double *d, double u, int in, int jn, int kn,
			      int i1, int i2, int j1, int j2, int k1, int k2,
			      struct polygon **poly)
{
    int i, j, k;
    int left = -(j2 - j1 + 2);
    int down = -(j2 - j1 + 2)*(i2 - i1 + 2);
    int up = jn*in;
    unsigned char *maparea, *map;
    double *kol, a;
/*    struct polygon *head = NULL, **poly = &head;	*/
    struct point *horn = NULL;

    maparea = calloc(-down * (k2 - k1 + 2), sizeof(char));
    if (!maparea)
	return NULL;
    map = maparea-down-left+1;	/* The beginning is not used */
    for (k = k1; k <= k2; k++, map -= left)
	for (i = i1; i <= i2; i++, map++) {
	    kol = d + k*in*jn + i*jn;
	    for (j = j1; j <= j2; j++, map++)
		if (kol[j] >= u) {
		    map[0] |= 128;
		    map[-1] |= 32;
		    map[left] |= 8;
		    map[left-1] |= 2;
		    map[down] |= 64;
		    map[down-1] |= 16;
		    map[down+left] |= 4;
		    map[down+left-1] |= 1;
		}
	}
    map = maparea-down-left+1;
    for (k = k1; k < k2; k++, map -= 2*left)
        for (i = i1; i < i2; i++, map += 2) {
            kol = d + k*in*jn + i*jn;
            for (j = j1; j < j2; j++, map++) {
		int p, sida;
		unsigned char index = *map;
		if (index >= 128)
		    index ^= 255;
		if (!index)
		    continue;
		for (p = 0; p < 4 && konf[index][p][1]; p++) {
		    *poly = allok(sizeof(struct polygon), (char *)horn);
		    if (!*poly)
			return NULL;
		    horn = (*poly)->horn;
		    for (sida = 0; sida < 6 && konf[index][p][sida]; sida++) {
			switch (konf[index][p][sida]) {
			    case 1:
				horn -> x = INTERP(i, kol[j], kol[j+jn]);
				horn -> y = j;
				horn -> z = k;
				break;
			    case 2:
				horn -> x = INTERP(i, kol[j+up], kol[j+up+jn]);
				horn -> y = j;
				horn -> z = k + 1;
				break;
			    case 3:
			        horn -> x = INTERP(i, kol[j+1], kol[j+1+jn]);
				horn -> y = j + 1;
				horn -> z = k;
				break;
			    case 4:
			        horn -> x = INTERP(i, kol[j+1+up],
				    kol[j+1+jn+up]);
				horn -> y = j + 1;
				horn -> z = k + 1;
				break;
			    case 5:
			    	horn -> x = i;
				horn -> y = INTERP(j, kol[j], kol[j+1]);
				horn -> z = k;
				break;
			    case 6:
			    	horn -> x = i;
				horn -> y = INTERP(j, kol[j+up], kol[j+1+up]);
				horn -> z = k + 1;
				break;
			    case 7:
			    	horn -> x = i + 1;
				horn -> y = INTERP(j, kol[j+jn], kol[j+1+jn]);
				horn -> z = k;
				break;
			    case 8:
			        horn -> x = i + 1;
				horn -> y = INTERP(j, kol[j+jn+up],
				    kol[j+1+jn+up]);
				horn -> z = k + 1;
				break;
			    case 9:
			        horn -> x = i;
				horn -> y = j;
				horn -> z = INTERP(k, kol[j], kol[j+up]);
				break;
			    case 10:
			        horn -> x = i;
				horn -> y = j + 1;
				horn -> z = INTERP(k, kol[j+1], kol[j+1+up]);
				break;
			    case 11:
			        horn -> x = i + 1;
				horn -> y = j;
				horn -> z = INTERP(k, kol[j+jn], kol[j+jn+up]);
				break;
			    case 12:
			        horn -> x = i + 1;
				horn -> y = j + 1;
				horn -> z = INTERP(k, kol[j+1+jn],
				    kol[j+1+jn+up]);
				break;
			}
			horn++;
		    }
		    (*poly)->sidor = sida;
		    poly = &(*poly)->next;
		}
	    }
	}
    *poly = NULL;
    free(maparea);
    return poly;	/* Return the address of the last element's 'next' */
}			/* pointer */
