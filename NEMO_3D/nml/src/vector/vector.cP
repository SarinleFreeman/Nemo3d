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
Modified by Marek J. Korkusinski, November 2004
*****************************************************************************
$Header: /repo/nml/src/vector/vector.cP,v 1.4 2004/11/26 22:56:42 marek Exp $
*****************************************************************************/


define(_point_,
ifelse(_precision_, long double, floating,
ifelse(_precision_,      double, floating,
ifelse(_precision_,       float, floating, fixed))))dnl

#undef <NS>INLINE
#undef <NS>MACROS
#include<<ns><p><s>vector.h>

/* Functions								*/
#ifdef EXT_MEM_HANDLING
<ns>extent		(<ns><p><s>v_extent)(const <ns><p><s>vector* pv) {
   return(  
            ( (<ns><p><s>darray *)
                     (    ((char *)pv) - <ns><p><s>offsetLength     )
            ) -> noElements
         ) ;
  }
#define <ns><p><s>v_extent(pv) (((<ns><p><s>darray *) \
   (((char *)pv) - <ns><p><s>offsetLength)) -> noElements )
#else   /* EXT_MEM_HANDLING */
<ns>extent		(<ns><p><s>v_extent)(const <ns><p><s>vector* pv) {
  return *((const <ns>extent*)(pv - 1));
  }
#define <ns><p><s>v_extent(pv) (*((const <ns>extent*)((pv) - 1)))
#endif   /* EXT_MEM_HANDLING */

<ns><p><s>vector*	(<ns><p><s>v_fill)(<ns><p><s>vector* pv,
    <ns><p><s>scalar x) {
  <ns>extent		n = <ns><p><s>v_extent(pv);
  <ns>offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = x;
  return pv;
  }

<ns><p><s>vector*	(<ns><p><s>v_smul)(<ns><p><s>vector* pv,
    <ns><p><s>scalar x) {
  <ns>extent		n = <ns><p><s>v_extent(pv);
  <ns>offset		j = 0;
  for (j = 0; j < n; ++j)
    <NS><P><S>MUL(pv[j], pv[j], x);
  return pv;
  }

<ns><p><s>vector*	(<ns><p><s>v_sdiv)(<ns><p><s>vector* pv,
    <ns><p><s>scalar x) {
  <ns>extent		n = <ns><p><s>v_extent(pv);
  <ns>offset		j = 0;
  for (j = 0; j < n; ++j)
    <NS><P><S>DIV(pv[j], pv[j], x);
  return pv;
  }

<ns><p><s>vector*	(<ns><p><s>v_sadd)(<ns><p><s>vector* pv,
    <ns><p><s>scalar x) {
  <ns>extent		n = <ns><p><s>v_extent(pv);
  <ns>offset		j = 0;
  for (j = 0; j < n; ++j)
    <NS><P><S>SUB(pv[j], pv[j], x);
  return pv;
  }

<ns><p><s>vector*	(<ns><p><s>v_ssub)(<ns><p><s>vector* pv,
    <ns><p><s>scalar x) {
  <ns>extent		n = <ns><p><s>v_extent(pv);
  <ns>offset		j = 0;
  for (j = 0; j < n; ++j)
    <NS><P><S>SUB(pv[j], pv[j], x);
  return pv;
  }

<ns><p>scalar		(<ns><p><s>v_sminmag)(const <ns><p><s>vector* pv) {
  <ns>extent		n = <ns><p><s>v_extent(pv);
  if (0 < n) {
    <ns><p>scalar	y = <NS><P><S>ABS(pv[0]);
    <ns>offset		j = 0;
    for (j = 0; j < n; ++j) {
      const
      <ns><p>scalar	t = <NS><P><S>ABS(pv[j]);
      if (t < y) y = t;
      }
    return y;
    }
#ifdef  <NS>ERRANT
  <ns>message("In function <ns><p><s>v_sminmag(const <ns><p><s>vector*),\n"
      "zero vector extent.");
#endif/*<NS>ERRANT              */
  return (<ns><p>scalar)0;
  }

<ns><p>scalar		(<ns><p><s>v_smaxmag)(const <ns><p><s>vector* pv) {
  <ns>extent		n = <ns><p><s>v_extent(pv);
  if (0 < n) {
    <ns><p>scalar	y = <NS><P><S>ABS(pv[0]);
    <ns>offset		j = 0;
    for (j = 1; j < n; ++j) {
      const
      <ns><p>scalar	t = <NS><P><S>ABS(pv[j]);
      if (y < t) y = t;
      }
    return y;
    }
#ifdef  <NS>ERRANT
  <ns>message("In function <ns><p><s>v_smaxmag(const <ns><p><s>vector*),\n"
      "zero vector extent.");
#endif/*<NS>ERRANT              */
  return (<ns><p>scalar)0;
  }

<ns><p>scalar		(<ns><p><s>v_ssqr)(const <ns><p><s>vector* pv) {
  <ns><p>scalar		t = 0;
  <ns>extent		n = <ns><p><s>v_extent(pv);
  <ns>offset		j = 0;
  for (j = 0; j < n; ++j)
    t += <NS><P><S>SQR(pv[j]);
  return t;
  }
ifelse(floating, _point_, `
<ns><p>scalar		(<ns><p><s>v_sabs)(const <ns><p><s>vector* pv) {
  return sqrt(<ns><p><s>v_ssqr(pv));
  }
')
<ns><p><s>vector*	(<ns><p><s>v_vcpy)(<ns><p><s>vector* pv,
    const <ns><p><s>vector* pw) {
  <ns>extent		n = <ns><p><s>v_extent(pv);
  <ns>offset		j = 0;
  for (j = 0; j < n; ++j)
    pv[j] = pw[j];
  return pv;
  }


int
(<ns><p><s>v_fprintf)(FILE* stream, const <ns><p><s>vector* pv, 
  int width, unsigned int precision, <ns>fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  <ns>extent	n = <ns><p><s>v_extent(pv);
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
<ns><p><s>vector*	(<ns><p><s>v_new)(<ns>extent n) {
   <ns><p><s>darray * p ;
   <ns>extent * pointerCopy ;

   p = ( <ns><p><s>darray * ) <ns>calloc(<ns><p><s>offsetLength +
                                       n * sizeof(<ns><p><s>scalar) , 1 ) ;
  if (NULL == p) {
    <ns>message("In function <ns><p><s>v_new(<ns>extent):\n"
      "unable to allocate memory.");
       }
  else {
      pointerCopy = (<ns>extent * ) p ;
      *pointerCopy = n ;
       } 
   return( (<ns><p><s>vector*) &p->start  ) ;
  }
#else  /*   EXT_MEM_HANDLING   */
<ns><p><s>vector*	(<ns><p><s>v_new)(<ns>extent n) {
     <ns><p><s>scalar*	p
       = (<ns><p><s>scalar*)<ns>calloc(1+n, sizeof(<ns><p><s>scalar));
  if (NULL == p) {
    <ns>message("In function <ns><p><s>v_new(<ns>extent):\n"
      "unable to allocate memory.");
    }
  else {
       *((<ns>extent*)p)	= n; ++p;
    }
  return p;
  }
#endif  /*   EXT_MEM_HANDLING   */

<ns><p><s>vector*	(<ns><p><s>v_clone)(const <ns><p><s>vector* pv) {
  <ns>extent		n = <ns><p><s>v_extent(pv);
  return <ns><p><s>v_vcpy(<ns><p><s>v_new(n), pv);
  }

/* Destructor								*/
#ifdef EXT_MEM_HANDLING
void			(<ns><p><s>v_delete)(<ns><p><s>vector* pv) {
   <ns><p><s>darray * pointerToDarray ;
	if (NULL != pv ) {
   pointerToDarray = 
   ( <ns><p><s>darray * ) ( (char*) pv - <ns><p><s>offsetLength ) ;
   <ns>free( pointerToDarray ) ;
	}
  }
#else  /*   EXT_MEM_HANDLING   */
void			(<ns><p><s>v_delete)(<ns><p><s>vector* pv) {
  if (NULL != pv)
  <ns>free(pv - 1);
  }
#endif  /*   EXT_MEM_HANDLING   */



/* Reconstructor							*/
<ns><p><s>vector*	(<ns><p><s>v_resize)(<ns><p><s>vector* *ppv, <ns>extent n) {
  if (NULL != ppv) {
    (<ns><p><s>v_delete)(*ppv);
    *ppv = (<ns><p><s>v_new)(n);
    return *ppv;
    }
  else {
    <ns>message(
"In function <ns><p><s>v_resize(<ns><p><s>vector**, <ns>extent):\n"
"invalid pointer to <ns><p><s>vector*.");
    return NULL;
    }
  }

