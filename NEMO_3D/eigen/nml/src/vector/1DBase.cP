dnl This is a class prototype input file for .../<ns>/bin/genclass.
pushdef(`index', `ifelse($#, 0, ``index'', ``index($*)'')')dnl
pushdef(`shift', `ifelse($#, 0, ``shift'', ``shift($*)'')')dnl

/*****************************************************************************
The <namespace>.
Copyright (C) 2002 California Institute of Technology (Caltech)

This file is part of
The NanoElectronic MOdeling (NEMO) Math Library.

This library is free software which you can redistribute and/or modify
under the terms of the GNU Library General Public License
as published by the Free Software Foundation;
either version 2, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

Written by E. Robert Tisdale

*****************************************************************************
$Header: /repo/nml/src/vector/1DBase.cP,v 1.4 2004/11/26 22:56:42 marek Exp $
*****************************************************************************/


#undef <NS>INLINE
#undef <NS>MACROS
#include<<ns><p><s>1DBase.h>

/* Functions								*/

#ifdef EXT_MEM_HANDLING
<ns>index		(<ns><p><s>1D_lower)(const <ns><p><s>1DBase* pv) {
   return(  
            ( (<ns><p><s>1D_darray *)
                     (    ((char *)pv) - <ns><p><s>1D_offsetLength     )
            ) -> lower
         ) ;
  }
#define <ns><p><s>1D_lower(pv) (((<ns><p><s>1D_darray *) \
     (((char *)pv) - <ns><p><s>1D_offsetLength )  ) -> lower  )

<ns>index		(<ns><p><s>1D_upper)(const <ns><p><s>1DBase* pv) {
   return(  
            ( (<ns><p><s>1D_darray *)
                     (    ((char *)pv) - <ns><p><s>1D_offsetLength     )
            ) -> upper
         ) ;
  }
#define <ns><p><s>1D_upper(pv) (((<ns><p><s>1D_darray *) \
     (((char *)pv) - <ns><p><s>1D_offsetLength )  ) -> upper  )

#else /*   EXT_MEM_HANDLING   */

<ns>index		(<ns><p><s>1D_lower)(const <ns><p><s>1DBase* pv) {
  return *((const <ns>index*)(pv - 1));
  }
#define <ns><p><s>1D_lower(pv) (*((const <ns>index*)((pv) - 1)))
<ns>index		(<ns><p><s>1D_upper)(const <ns><p><s>1DBase* pv) {
  return *((const <ns>index*)(pv - 2));
  }
#define <ns><p><s>1D_upper(pv) (*((const <ns>index*)((pv) - 2)))
#endif  /*   EXT_MEM_HANDLING   */

<ns>extent		(<ns><p><s>1D_extent)(const <ns><p><s>1DBase* pv) {
  return 1+<ns><p><s>1D_upper(pv)-<ns><p><s>1D_lower(pv);
  }
#define <ns><p><s>1D_extent(pv) (1+<ns><p><s>1D_upper(pv)-<ns><p><s>1D_lower(pv))

int
(<ns><p><s>1D_fprintf)(FILE* stream, const <ns><p><s>1DBase* pv, 
  int width, unsigned int precision, <ns>fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  <ns>extent	n = <ns><p><s>1D_extent(pv);
  <ns>offset	j = 0;
ifelse(system, complex,
`  columns = (0 < columns)? columns: 2;
',
`  columns = (0 < columns)? columns: 4;
')dnl
  for (j = 0; j < n; ++j) {
    if (0 < j) {
      characters
	= (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = <ns><p><s>fprintf(
      stream, pv[j], width, precision, flags);
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  characters = fprintf(stream, "\n");
  if (0 <= characters) { total += characters; } else { return -1; }
  return total;
  }



/* Constructor								*/
#ifdef EXT_MEM_HANDLING
<ns><p><s>1DBase*	(<ns><p><s>1D_new)(<ns>index lb, <ns>index ub) {
   <ns><p><s>1D_darray * p = NULL ;
   <ns>index * pointerCopy ;
   <ns>index		n = 1 + ub - lb;
  if (n < 0) {
    <ns>message("In function <ns><p><s>1D_new(<ns>index. <ns>index):\n"
      "extent is negative.");
    }
  else {		/* (0 <= n)		*/
    p = ( <ns><p><s>1D_darray * ) <ns>calloc(<ns><p><s>1D_offsetLength +
                                       n * sizeof(<ns><p><s>scalar) , 1 ) ;
    if (NULL == p) {
      <ns>message("In function <ns><p><s>1D_new(<ns>index. <ns>index):\n"
	"unable to allocate memory.");
      }
    else {		/* (NULL != p)		*/
      pointerCopy = ( <ns>index * ) p ;
      *pointerCopy = ub ; 
      ++pointerCopy ;
      *pointerCopy = lb ;
      }
    }
  return( (<ns><p><s>scalar*)  &p->start ) ;
  }
#else /*   EXT_MEM_HANDLING   */

<ns><p><s>1DBase*	(<ns><p><s>1D_new)(<ns>index lb, <ns>index ub) {
  <ns><p><s>scalar*	p = NULL;
  <ns>index		n = 1 + ub - lb;
  if (n < 0) {
    <ns>message("In function <ns><p><s>1D_new(<ns>index, <ns>index):\n"
      "extent is negative.");
    }
  else {		/* (0 <= n)		*/
    p = (<ns><p><s>scalar*)<ns>calloc(2+n, sizeof(<ns><p><s>scalar));
    if (NULL == p) {
      <ns>message("In function <ns><p><s>1D_new(<ns>index, <ns>index):\n"
	"unable to allocate memory.");
      }
    else {		/* (NULL != p)		*/
      *((<ns>index*)p)	= ub; ++p;
      *((<ns>index*)p)	= lb; ++p;
      }
    }
  return p;
  }
#endif /*   EXT_MEM_HANDLING   */

/* Destructor								*/

#ifdef EXT_MEM_HANDLING
void			(<ns><p><s>1D_delete)(<ns><p><s>1DBase* pv) {
   <ns><p><s>1D_darray * pointerToDarray ;
	if (NULL != pv) {
   pointerToDarray = 
   ( <ns><p><s>1D_darray * ) ( (char*) pv - <ns><p><s>1D_offsetLength ) ;
   <ns>free( pointerToDarray ) ;
	}
  }
#else /*   EXT_MEM_HANDLING   */
void			(<ns><p><s>1D_delete)(<ns><p><s>1DBase* pv) {
  <ns>free(pv - 2);
  }
#endif  /*   EXT_MEM_HANDLING   */


/* Reconstructor							*/
<ns><p><s>1DBase*	(<ns><p><s>1D_resize)(<ns><p><s>1DBase* *ppv,
    <ns>index lb, <ns>index ub) {
  if (NULL != ppv) {
    (<ns><p><s>1D_delete)(*ppv);
    *ppv = (<ns><p><s>1D_new)(lb, ub);
    return *ppv;
    }
  else {		/* (NULL == ppv)	*/
    <ns>message(
"In function <ns><p><s>1D_resize(<ns><p><s>1DBase**, <ns>index, <ns>index):\n"
"invalid pointer to <ns><p><s>1DBase*.");
    return NULL;
    }
  }

