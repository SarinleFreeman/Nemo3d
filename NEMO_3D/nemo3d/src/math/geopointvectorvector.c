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
$Header: /repo/nemo3d/src/math/geopointvectorvector.c,v 1.3 2004/10/27 19:43:10 gekco Exp $ 
*****************************************************************************/

#include "geopointvectorvector.h"

/* Create memory to a vector of Geopointvector */
geopointvectrvectr Geopointvectrvectr(int n){
  int i;
  geopointvectrvectr I = (geopointvectr *) nml_calloc ((n+1), sizeof (geopointvectr));
  if (I == NULL){
    die("unable to allocate memory in Geopointvectrvectr.\n");
  }
  I++;
  vdim(I) = n;

  for (i=0;i<n;i++){
    I[i] = NULL;
  }

  return I;
}


/* Remove the memory associated with the vector to Geopointvectr. */
void rm_geopointvectrvectr (geopointvectrvectr *I_ptr)
{
  geopointvectrvectr I;
  int i;

  if (!I_ptr || !(I=*I_ptr))
    return;

  for (i=0;i<vdim(I);i++)
    rm_geopointvectr(&I[i]);

  nml_free(--I);

  *I_ptr=NULL;

  return;
}

/* Remove the memory associated with the vector to Geopointvectr. 
   Only remove the memory associated with the 
   holding vector; not the datavectors. */
void rm_geopointvectrvectr_shellonly (geopointvectrvectr *I_ptr)
{
  geopointvectrvectr I;

  if (!I_ptr || !(I=*I_ptr))
    return;

  nml_free(--I);

  *I_ptr=NULL;

  return;
}

/* Attach a data vector to the holding structure.  Do not create memory for d
   the data vector itself.  */
int add_elem_geopointvectrvectr(geopointvectr data_vectr, int pos, geopointvectrvectr hold_vectr)
{
  if (pos>vdim(hold_vectr)){
    die("Error in add_elem_geopointvectrvectr\n"
	  "Requested Position in the vector of vectors\n"
	  "exceeds holding capacity.\n"
	  "pos=%d dim=%d\n",pos,vdim(hold_vectr)); 
    return -1;
  }

  if (hold_vectr[pos]!=NULL){
    rm_geopointvectr(&hold_vectr[pos]);
  }

  hold_vectr[pos] = data_vectr;
  
  return 0;
}

geopoint min_gp_geopointvectrvectr(geopointvectrvectr gpvv)
{
    int i,j;
    geopoint gp;
    gp = gpvv[0][0];

    for (i=0;i<vdim(gpvv);i++){
	for (j=0;j<vdim(gpvv[i]);j++){
	    if (gpvv[i][j].x < gp.x &&
		gpvv[i][j].x < gp.x &&
		gpvv[i][j].x < gp.x  ){
		gp = gpvv[i][j];
	    }
	}
    }

    return gp;
}

geopoint max_gp_geopointvectrvectr(geopointvectrvectr gpvv)
{
    int i,j;
    geopoint gp;
    gp = gpvv[0][0];

    for (i=0;i<vdim(gpvv);i++){
	for (j=0;j<vdim(gpvv[i]);j++){
	    if (gpvv[i][j].x > gp.x &&
		gpvv[i][j].y > gp.y &&
		gpvv[i][j].z > gp.z  ){
		gp = gpvv[i][j];
	    }
	}
    }

    return gp;
}

geopoint max2_gp_geopointvectrvectr(geopointvectrvectr gpvv)
{
    int i,j;
    real x = gpvv[0][0].x;
    real y = gpvv[0][0].y;
    real z = gpvv[0][0].z;
    geopoint gp;

    for (i=0;i<vdim(gpvv);i++){
	for (j=0;j<vdim(gpvv[i]);j++){
	    if (gpvv[i][j].x > x ){
		x = gpvv[i][j].x;
	    }
	    if (gpvv[i][j].y > y ){
		y = gpvv[i][j].y;
	    }
	    if (gpvv[i][j].z > z ){
		z = gpvv[i][j].z;
	    }
	}
    }

    gp.x=x;
    gp.y=y;
    gp.z=z;

    return gp;
}

geopoint min2_gp_geopointvectrvectr(geopointvectrvectr gpvv)
{
    int i,j;
    real x = gpvv[0][0].x;
    real y = gpvv[0][0].y;
    real z = gpvv[0][0].z;
    geopoint gp;

    for (i=0;i<vdim(gpvv);i++){
	for (j=0;j<vdim(gpvv[i]);j++){
	    if (gpvv[i][j].x < x ){
		x = gpvv[i][j].x;
	    }
	    if (gpvv[i][j].y < y ){
		y = gpvv[i][j].y;
	    }
	    if (gpvv[i][j].z < z ){
		z = gpvv[i][j].z;
	    }
	}
    }

    gp.x=x;
    gp.y=y;
    gp.z=z;

    return gp;
}


