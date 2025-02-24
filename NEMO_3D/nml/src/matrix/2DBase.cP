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
$Header: /repo/nml/src/matrix/2DBase.cP,v 1.2 2003/10/08 16:19:30 hook Exp $
*****************************************************************************/

#undef <NS>INLINE
#undef <NS>MACROS
#include<<ns><p><s>2DBase.h>

/* Functions								*/
<ns>index		(<ns><p><s>2D_lower2)(const <ns><p><s>2DBase* pM) {
  return *((const <ns>index*)(pM - 1));
  }
#define <ns><p><s>2D_lower2(pM) (*((const <ns>index*)((pM) - 1)))

<ns>index		(<ns><p><s>2D_upper2)(const <ns><p><s>2DBase* pM) {
  return *((const <ns>index*)(pM - 2));
  }
#define <ns><p><s>2D_upper2(pM) (*((const <ns>index*)((pM) - 2)))

<ns>extent		(<ns><p><s>2D_extent2)(const <ns><p><s>2DBase* pM) {
  return 1+<ns><p><s>2D_upper2(pM)-<ns><p><s>2D_lower2(pM);
  }
#define <ns><p><s>2D_extent2(pM) (1+<ns><p><s>2D_upper2(pM)-<ns><p><s>2D_lower2(pM))

<ns>index		(<ns><p><s>2D_lower1)(const <ns><p><s>2DBase* pM) {
  return *((const <ns>index*)(pM - 3));
  }
#define <ns><p><s>2D_lower1(pM) (*((const <ns>index*)((pM) - 3)))

<ns>index		(<ns><p><s>2D_upper1)(const <ns><p><s>2DBase* pM) {
  return *((const <ns>index*)(pM - 4));
  }
#define <ns><p><s>2D_upper1(pM) (*((const <ns>index*)((pM) - 4)))

<ns>extent		(<ns><p><s>2D_extent1)(const <ns><p><s>2DBase* pM) {
  return 1+<ns><p><s>2D_upper1(pM)-<ns><p><s>2D_lower1(pM);
  }
#define <ns><p><s>2D_extent1(pM) (1+<ns><p><s>2D_upper1(pM)-<ns><p><s>2D_lower1(pM))

int
(<ns><p><s>2D_fprintf)(FILE* stream, const <ns><p><s>2DBase* pM, 
  int width, unsigned int precision, <ns>fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  <ns>extent	m = <ns><p><s>2D_extent2(pM);
  <ns>offset	i = 0;
ifelse(system, complex,
`  columns = (0 < columns)? columns: 2;
',
`  columns = (0 < columns)? columns: 4;
')dnl
  for (i = 0; i < m; ++i) {
    const
    <ns>extent	n = <ns><p><s>2D_extent1(pM);
    <ns>offset	j = 0;
    for (j = 0; j < n; ++j) {
      <ns>index k = (<ns>index)j + <ns><p><s>2D_lower1(pM);
      if (0 < j) {
	characters
	  = (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
	if (0 <= characters) { total += characters; } else { return -1; }
	}
      characters = <ns><p><s>fprintf(
	stream, pM[i][k], width, precision, flags);
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = fprintf(stream, "\n");
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  return total;
  }

/* Constructor								*/
<ns><p><s>2DBase*	(<ns><p><s>2D_new)(
    <ns>index lb2, <ns>index ub2, <ns>index lb1, <ns>index ub1) {
  <ns><p><s>scalar**	pp = NULL;
  <ns>index		m = 1 + ub2 - lb2;
  <ns>index		n = 1 + ub1 - lb1;
  if (m < 0) {
    <ns>message(
"In function <ns><p><s>2D_new(<ns>index, <ns>index, <ns>index, <ns>index):\n"
"extent2 is negative.");
    }
  else
  if (n < 0) {
    <ns>message(
"In function <ns><p><s>2D_new(<ns>index, <ns>index, <ns>index, <ns>index):\n"
"extent1 is negative.");
    }
  else {		/* (0 <= m) && (0 <= n)	*/
    pp = (<ns><p><s>scalar**)<ns>calloc(4+m, sizeof(<ns><p><s>scalar*));
    if (NULL != pp) {
      *((<ns>index*)pp) = ub1; ++pp;
      *((<ns>index*)pp) = lb1; ++pp;
      *((<ns>index*)pp) = ub2; ++pp;
      *((<ns>index*)pp) = lb2; ++pp;
      if (0 < m) {
	if (0 < n) {
	  <ns><p><s>scalar*	p
	    = (<ns><p><s>scalar*)<ns>calloc(m*n, sizeof(<ns><p><s>scalar));
	  if (NULL != p) {
	    pp[0] = p - lb1;
	    { <ns>index	i = 1;
	      for (i = 1; i < m; ++i)
		pp[i] = pp[i-1] + n;
	      }
	    }
	  else {	/* (NULL == p)		*/
	    int	i = 0;
            for (i = 0; i < m; ++i)
	      pp[i] = NULL;
	    <ns>message(
"In function <ns><p><s>2D_new(<ns>index, <ns>index, <ns>index, <ns>index):\n"
"unable to allocate memory.");
	    }
	  }
	else {		/* (n <= 0)		*/
	  int	i = 0;
          for (i = 0; i < m; ++i)
	    pp[i] = NULL;
	  }
	}
      }
    else {		/* (NULL == pp)		*/
      <ns>message(
"In function <ns><p><s>2D_new(<ns>index, <ns>index, <ns>index, <ns>index):\n"
"unable to allocate memory.");
      }
    }
  return pp;
  }

/* Destructor								*/
void			(<ns><p><s>2D_delete)(<ns><p><s>2DBase* pM) {
  if (NULL != pM) {
    if (NULL != pM[0]) {
      <ns>index   lb1 = <ns><p><s>2D_lower1(pM);
      <ns>free(pM[0] + lb1);
      }
    <ns>free(pM - 4);
    }
  }

/* Reconstructor							*/
<ns><p><s>2DBase*	(<ns><p><s>2D_resize)(<ns><p><s>2DBase* *ppM,
    <ns>index lb2, <ns>index ub2, <ns>index lb1, <ns>index ub1) {
  if (NULL != ppM) {
    (<ns><p><s>2D_delete)(*ppM);
    *ppM = (<ns><p><s>2D_new)(lb2, ub2, lb1, ub1);
    return *ppM;
    }
  else {
    <ns>message(
"In function <ns><p><s>2D_resize(<ns><p><s>2DBase**, <ns>index, <ns>index,\n"
"  <ns>index, <ns>index): invalid pointer to <ns><p><s>2DBase*.");
    return NULL;
    }
  }

