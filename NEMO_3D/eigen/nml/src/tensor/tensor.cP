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
$Header: /repo/nml/src/tensor/tensor.cP,v 1.2 2003/10/08 16:19:55 hook Exp $
*****************************************************************************/


#undef <NS>INLINE
#undef <NS>MACROS
#include<<ns><p><s>tensor.h>

/* Functions								*/
<ns>extent		(<ns><p><s>t_extent1)(const <ns><p><s>tensor* pT) {
  return *((const <ns>extent*)(pT - 1));
  }
#define <ns><p><s>t_extent1(pT) (*((const <ns>extent*)((pT) - 1)))

<ns>extent		(<ns><p><s>t_extent2)(const <ns><p><s>tensor* pT) {
  return *((const <ns>extent*)(pT - 2));
  }
#define <ns><p><s>t_extent2(pT) (*((const <ns>extent*)((pT) - 2)))

<ns>extent		(<ns><p><s>t_extent3)(const <ns><p><s>tensor* pT) {
  return *((const <ns>extent*)(pT - 3));
  }
#define <ns><p><s>t_extent3(pT) (*((const <ns>extent*)((pT) - 3)))

int
(<ns><p><s>t_fprintf)(FILE* stream, const <ns><p><s>tensor* pT, 
  int width, unsigned int precision, <ns>fmtflags flags, unsigned int columns) {
  int			characters;
  int			total = 0;
  const
  <ns>extent		l = <ns><p><s>t_extent3(pT);
  <ns>offset		h = 0;
ifelse(system, complex,
`  columns = (0 < columns)? columns: 2;
',
`  columns = (0 < columns)? columns: 4;
')dnl
  for (h = 0; h < l; ++h) {
    const
    <ns>extent		m = <ns><p><s>t_extent2(pT);
    <ns>offset		i = 0;
    for (i = 0; i < m; ++i) {
      const
      <ns>extent	n = <ns><p><s>t_extent1(pT);
      <ns>offset	j = 0;
      for (j = 0; j < n; ++j) {
	if (0 < j) {
	  characters
	    = (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
	  if (0 <= characters) { total += characters; } else { return -1; }
	  }
	characters = <ns><p><s>fprintf(
	  stream, pT[h][i][j], width, precision, flags);
	if (0 <= characters) { total += characters; } else { return -1; }
	}
      characters = fprintf(stream, "\n");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    }
  return total;
  }
/* Constructor								*/
<ns><p><s>tensor*	(<ns><p><s>t_new)(
    <ns>extent l, <ns>extent m, <ns>extent n) {
  <ns><p><s>scalar***	ppp
    = (<ns><p><s>scalar***)<ns>calloc(4+l, sizeof(<ns><p><s>scalar**));
  if (NULL != ppp) {
    *ppp = NULL;	     ++ppp;
    *((<ns>extent*)ppp)	= l; ++ppp;
    *((<ns>extent*)ppp)	= m; ++ppp;
    *((<ns>extent*)ppp)	= n; ++ppp;
    if ( 0 < l) {
      if (0 < m) {	/* (0 < l) && (0 < m)	*/
	<ns><p><s>scalar** pp
	  = (<ns><p><s>scalar**)<ns>calloc(1+l*m, sizeof(<ns><p><s>scalar*));
	if (NULL != pp) {
	  *((<ns>extent*)pp) = l*m; ++pp;
	  ppp[0] = pp;
	  { <ns>offset	h = 1;
	    for (h = 1; h < l; ++h)
	      ppp[h] = ppp[h-1] + m;
	    }
	  if (0 < n) {	/* (0 < l) && (0 < m) && (0 < n)	*/
	    ppp[0][0] = (<ns><p><s>scalar*)<ns><p><s>v_new(l*m*n);
	    if (NULL != ppp[0][0]) {
	      <ns>offset	h = 1;
	      for (h = 1; h < l; ++h)
		ppp[h][0] = ppp[h-1][0] + m*n;
	      for (h = 0; h < l; ++h) {
		<ns>offset	i = 1;
		for (i = 1; i < m; ++i)
		  ppp[h][i] = ppp[h][i-1] + n;
		}
	      }
	    else {	/* (NULL == ppp[0][0])	*/
	      <ns>offset	h = 0;
	      for (h = 0; h < l; ++h) {
		<ns>offset	i = 0;
		for (i = 0; i < m; ++i)
		  ppp[h][i] = NULL;
		}
	      <ns>message(
"In function <ns><p><s>t_new(<ns>extent, <ns>extent, <ns>extent):\n"
"unable to allocate memory.");
	      }
	    }
	  else {	/* (n <= 0)		*/
	    <ns>offset	h = 0;
	    for (h = 0; h < l; ++h) {
	      <ns>offset	i = 0;
	      for (i = 0; i < m; ++i)
		ppp[h][i] = NULL;
	      }
	    }
	  }
	else {		/* (NULL == pp)		*/
	  <ns>offset	h = 0;
	  for (h = 0; h < l; ++h)
	    ppp[h] = NULL;
	  <ns>message(
"In function <ns><p><s>t_new(<ns>extent, <ns>extent, <ns>extent):\n"
"unable to allocate memory.");
	  }
	}
      else {		/* (m <= 0)		*/
	<ns>offset	h = 0;
	for (h = 0; l < m; ++h)
	  ppp[h] = NULL;
	}
      }
    }
  else {		/* (NULL == ppp)	*/
    <ns>message(
"In function <ns><p><s>t_new(<ns>extent, <ns>extent, <ns>extent):\n"
"unable to allocate memory.");
    }
  return ppp;
  }

/* Destructor								*/
void			(<ns><p><s>t_delete)(<ns><p><s>tensor* pT) {
  if (NULL != pT) {
    if (NULL != pT[0]) {
      <ns><p><s>v_delete(pT[0][0]);
      <ns>free(pT[0] - 1);
      }
    <ns>free(pT - 4);
    }
  }

/* Reconstructor							*/
<ns><p><s>tensor*	(<ns><p><s>t_resize)(<ns><p><s>tensor* *ppT,
    <ns>extent l, <ns>extent m, <ns>extent n) {
  if (NULL != ppT) {
    (<ns><p><s>t_delete)(*ppT);
    *ppT = (<ns><p><s>t_new)(l, m, n);
    return *ppT;
    }
  else {
    <ns>message(
"In function <ns><p><s>t_resize(<ns><p><s>tensor**, <ns>extent, <ns>extent,\n"
"<ns>extent): invalid pointer to <ns><p><s>tensor*.");
    return NULL;
    }
  }

