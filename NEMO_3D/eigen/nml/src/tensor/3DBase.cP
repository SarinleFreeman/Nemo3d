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
$Header: /repo/nml/src/tensor/3DBase.cP,v 1.2 2003/10/08 16:19:53 hook Exp $
*****************************************************************************/


#undef <NS>INLINE
#undef <NS>MACROS
#include<<ns><p><s>3DBase.h>

/* Functions								*/
<ns>index		(<ns><p><s>3D_lower3)(const <ns><p><s>3DBase* pT) {
  return *((const <ns>index*)(pT - 1));
  }
#define <ns><p><s>3D_lower3(pT) (*((const <ns>index*)((pT) - 1)))

<ns>index		(<ns><p><s>3D_upper3)(const <ns><p><s>3DBase* pT) {
  return *((const <ns>index*)(pT - 2));
  }
#define <ns><p><s>3D_upper3(pT) (*((const <ns>index*)((pT) - 2)))

<ns>extent		(<ns><p><s>3D_extent3)(const <ns><p><s>3DBase* pT) {
  return 1+<ns><p><s>3D_upper3(pT)-<ns><p><s>3D_lower3(pT);
  }
#define <ns><p><s>3D_extent3(pT) (1+<ns><p><s>3D_upper3(pT)-<ns><p><s>3D_lower3(pT))

<ns>index		(<ns><p><s>3D_lower2)(const <ns><p><s>3DBase* pT) {
  return *((const <ns>index*)(pT - 3));
  }
#define <ns><p><s>3D_lower2(pT) (*((const <ns>index*)((pT) - 3)))

<ns>index		(<ns><p><s>3D_upper2)(const <ns><p><s>3DBase* pT) {
  return *((const <ns>index*)(pT - 4));
  }
#define <ns><p><s>3D_upper2(pT) (*((const <ns>index*)((pT) - 4)))

<ns>extent		(<ns><p><s>3D_extent2)(const <ns><p><s>3DBase* pT) {
  return 1+<ns><p><s>3D_upper2(pT)-<ns><p><s>3D_lower2(pT);
  }
#define <ns><p><s>3D_extent2(pT) (1+<ns><p><s>3D_upper2(pT)-<ns><p><s>3D_lower2(pT))

<ns>index		(<ns><p><s>3D_lower1)(const <ns><p><s>3DBase* pT) {
  return *((const <ns>index*)(pT - 5));
  }
#define <ns><p><s>3D_lower1(pT) (*((const <ns>index*)((pT) - 5)))

<ns>index		(<ns><p><s>3D_upper1)(const <ns><p><s>3DBase* pT) {
  return *((const <ns>index*)(pT - 6));
  }
#define <ns><p><s>3D_upper1(pT) (*((const <ns>index*)((pT) - 6)))

<ns>extent		(<ns><p><s>3D_extent1)(const <ns><p><s>3DBase* pT) {
  return 1+<ns><p><s>3D_upper1(pT)-<ns><p><s>3D_lower1(pT);
  }
#define <ns><p><s>3D_extent1(pT) (1+<ns><p><s>3D_upper1(pT)-<ns><p><s>3D_lower1(pT))
int
(<ns><p><s>3D_fprintf)(FILE* stream, const <ns><p><s>3DBase* pT, 
  int width, unsigned int precision, <ns>fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  <ns>extent	l = <ns><p><s>3D_extent3(pT);
  <ns>offset	h = 0;
ifelse(system, complex,
`  columns = (0 < columns)? columns: 2;
',
`  columns = (0 < columns)? columns: 4;
')dnl
  for (h = 0; h < l; ++h) {
    const
    <ns>extent	m = <ns><p><s>3D_extent2(pT);
    <ns>offset	i = 0;
    for (i = 0; i < m; ++i) {
      <ns>index k2 = (<ns>index)i + <ns><p><s>3D_lower2(pT);
      const
      <ns>extent	n = <ns><p><s>3D_extent1(pT);
      <ns>offset	j = 0;
      for (j = 0; j < n; ++j) {
	<ns>index k1 = (<ns>index)j + <ns><p><s>3D_lower1(pT);
	if (0 < j) {
	  characters
	    = (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
	  if (0 <= characters) { total += characters; } else { return -1; }
	  }
	characters = <ns><p><s>fprintf(
	  stream, pT[h][k2][k1], width, precision, flags);
	if (0 <= characters) { total += characters; } else { return -1; }
	}
      characters = fprintf(stream, "\n");
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    }
  return total;
  }

/* Constructor								*/
<ns><p><s>3DBase*	(<ns><p><s>3D_new)(<ns>index lb3, <ns>index ub3,
    <ns>index lb2, <ns>index ub2, <ns>index lb1, <ns>index ub1) {
  <ns><p><s>scalar***	ppp = NULL;
  <ns>index		l = 1 + ub3 - lb3;
  <ns>index		m = 1 + ub2 - lb2;
  <ns>index		n = 1 + ub1 - lb1;
  if (l < 0) {
    <ns>message(
"In function <ns><p><s>3D_new(<ns>index, <ns>index, <ns>index, <ns>index,\n"
"  <ns>index, <ns>index): extent3 is negative.");
    }
  else
  if (m < 0) {
    <ns>message(
"In function <ns><p><s>3D_new(<ns>index, <ns>index, <ns>index, <ns>index,\n"
"  <ns>index, <ns>index): extent2 is negative.");
    }
  else
  if (n < 0) {
    <ns>message(
"In function <ns><p><s>3D_new(<ns>index, <ns>index, <ns>index, <ns>index,\n"
"  <ns>index, <ns>index): extent1 is negative.");
    }
  else {		/* (0 <= l) && (0 <= m) && (0 <= n)	*/
    ppp = (<ns><p><s>scalar***)<ns>calloc(6 + l, sizeof(<ns><p><s>scalar**));
    if (NULL != ppp) {
      *((<ns>index*)ppp) = ub1; ++ppp;
      *((<ns>index*)ppp) = lb1; ++ppp;
      *((<ns>index*)ppp) = ub2; ++ppp;
      *((<ns>index*)ppp) = lb2; ++ppp;
      *((<ns>index*)ppp) = ub3; ++ppp;
      *((<ns>index*)ppp) = lb3; ++ppp;
      if (0 < l) {
	if (0 < m) {	/* (0 < l) && (0 < m)	*/
	  <ns><p><s>scalar** pp
	    = (<ns><p><s>scalar**)<ns>calloc(l*m, sizeof(<ns><p><s>scalar*));
	  if (NULL != pp) {
	    ppp[0] = pp - lb2;
	    { <ns>index	h = 1;
	      for (h = 1; h < l; ++h)
		ppp[h] = ppp[h-1] + m;
	      }
	    if (0 < n) {
	      <ns><p><s>scalar*	p
		= (<ns><p><s>scalar*)<ns>calloc(l*m*n, sizeof(<ns><p><s>scalar));
	      if (NULL != p) {
		ppp[0][lb2] = p - lb1;
		{ <ns>index	h = 1;
		  for (h = 1; h < l; ++h)
		    ppp[h][lb2] = ppp[h-1][lb2] + m*n;
		  for (h = 0; h < l; ++h) {
		    <ns>index	i = lb2 + 1;
		    for (i = lb2 + 1; i <= ub2; ++i)
		      ppp[h][i] = ppp[h][i-1] + n;
		    }
		  }
		}

	      else {	/* (NULL == p)		*/
		<ns>index	h = 0;
		for (h = 0; h < l; ++h) {
		  <ns>index	i = lb2;
		  for (i = lb2; i <= ub2; ++i)
		    ppp[h][i] = NULL;
		  }
		<ns>message(
"In function <ns><p><s>3D_new(<ns>index, <ns>index, <ns>index, <ns>index,\n"
"  <ns>index, <ns>index): unable to allocate memory.");
		}
	      }
	    else {	/* (n <= 0)		*/
	      <ns>index	h = 0;
	      for (h = 0; h < l; ++h) {
		<ns>index	i = lb2;
		for (i = lb2; i <= ub2; ++i)
		  ppp[h][i] = NULL;
		}
	      }
	    }
	  else {	/* (NULL == pp)		*/
	    <ns>index	h = 0;
	    for (h = 0; h < l; ++h)
	      ppp[h] = NULL;
	    <ns>message(
"In function <ns><p><s>3D_new(<ns>index, <ns>index, <ns>index, <ns>index,\n"
"  <ns>index, <ns>index): unable to allocate memory.");
	    }
	  }
	else {		/* (m <= 0)		*/
	  <ns>index	h = 0;
	  for (h = 0; h < l; ++h)
	    ppp[h] = NULL;
	  }
	}
      }
    else {		/* (NULL == ppp)	*/
      <ns>message(
"In function <ns><p><s>3D_new(<ns>index, <ns>index, <ns>index, <ns>index,\n"
"  <ns>index, <ns>index): unable to allocate memory.");
      }
    }
  return ppp;
  }

/* Destructor								*/
void			(<ns><p><s>3D_delete)(<ns><p><s>3DBase* pT) {
  if (NULL != pT) {
    if (NULL != pT[0]) {
      <ns>index	lb2 = <ns><p><s>3D_lower2(pT);
      if (NULL != pT[0][lb2]) {
	<ns>index	lb1 = <ns><p><s>3D_lower1(pT);
	<ns>free(pT[0][lb2] + lb1);
	}
      <ns>free(pT[0] + lb2);
      }
    <ns>free(pT - 6);
    }
  }

/* Reconstructor							*/
<ns><p><s>3DBase*	(<ns><p><s>3D_resize)(<ns><p><s>3DBase* *ppT,
    <ns>index lb3, <ns>index ub3, <ns>index lb2, <ns>index ub2,
    <ns>index lb1, <ns>index ub1) {
  if (NULL != ppT) {
    (<ns><p><s>3D_delete)(*ppT);
    *ppT = (<ns><p><s>3D_new)(lb3, ub3, lb2, ub2, lb1, ub1);
    return *ppT;
    }
  else {
    <ns>message(
"In function <ns><p><s>3D_resize(<ns><p><s>3DBase**,\n"
"  <ns>index, <ns>index, <ns>index, <ns>index, <ns>index, <ns>index):\n"
"invalid pointer to <ns><p><s>3DBase*.");
    return NULL;
    }
  }

