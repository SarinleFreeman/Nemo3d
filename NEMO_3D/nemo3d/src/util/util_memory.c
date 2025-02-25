/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This application is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)

Written by:  Chris Bowen
             Gerhard Klimeck
             Fabiano Oyafuso
             Seungwon Lee
             Olga Lazarenkova
             Hook Hua

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/util/util_memory.c,v 1.3 2004/10/27 19:43:10 gekco Exp
$
*****************************************************************************/

#if 0
#include "util_memory.h"

#define NEW_PSIZE

void *n3d_malloc(size_t size) {
    void *p = NULL;
    if ((size_t) NULL != size) {
        if (NULL == (p = (void*) malloc(size)) ) {
            printf("\n\nIn function n3d_malloc(size_t):: Error - n3d_malloc() could not allocate memory.\n");
            exit(1);
        }
    }
    return p;
}

void *n3d_calloc(int n, size_t size) {
    void *p = NULL;
    if ( (size_t) NULL != size) {
        if (NULL == (p = (void*) calloc(n, size)) ) {
            printf("\n\nIn function n3d_calloc(size_t):: Error - n3d_calloc() could not allocate memory.\n");
            exit(1);
        }
    }
    return p;
}

void n3d_free(void *p) {
    if (p != NULL) {
        /* free storage */
        free(p);
    }
    return;
}

#endif
