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
$Header: /repo/nml/src/matrix/matrix.cP,v 1.2 2003/10/08 16:19:34 hook Exp $
*****************************************************************************/

#undef <NS>INLINE
#undef <NS>MACROS
#include<<ns><p><s>matrix.h>

/* Functions								*/
<ns>extent		(<ns><p><s>m_extent1)(const <ns><p><s>matrix* pM) {
  return *((const <ns>extent*)(pM - 1));
  }
#define <ns><p><s>m_extent1(pM) (*((const <ns>extent*)((pM) - 1)))

<ns>extent		(<ns><p><s>m_extent2)(const <ns><p><s>matrix* pM) {
  return *((const <ns>extent*)(pM - 2));
  }
#define <ns><p><s>m_extent2(pM) (*((const <ns>extent*)((pM) - 2)))

int
(<ns><p><s>m_fprintf)(FILE* stream, const <ns><p><s>matrix* pM, 
  int width, unsigned int precision, <ns>fmtflags flags, unsigned int columns) {
  int		characters;
  int		total = 0;
  const
  <ns>extent	m = <ns><p><s>m_extent2(pM);
  <ns>offset	i = 0;
ifelse(system, complex,
`  columns = (0 < columns)? columns: 2;
',
`  columns = (0 < columns)? columns: 4;
')dnl
  for (i = 0; i < m; ++i) {
    const
    <ns>extent	n = <ns><p><s>m_extent1(pM);
    <ns>offset	j = 0;
    for (j = 0; j < n; ++j) {
      if (0 < j) {
	characters
	  = (0 == j%columns)? fprintf(stream, "\n"): fprintf(stream, " ");
	if (0 <= characters) { total += characters; } else { return -1; }
	}
      characters = <ns><p><s>fprintf(
	stream, pM[i][j], width, precision, flags);
      if (0 <= characters) { total += characters; } else { return -1; }
      }
    characters = fprintf(stream, "\n");
    if (0 <= characters) { total += characters; } else { return -1; }
    }
  return total;
  }

/* Constructor								*/
<ns><p><s>matrix*	(<ns><p><s>m_new)(<ns>extent m, <ns>extent n) {
  <ns><p><s>scalar**	pp
    = (<ns><p><s>scalar**)<ns>calloc(3+m, sizeof(<ns><p><s>scalar*));
  if (NULL == pp) {
    <ns>message("In function <ns><p><s>m_new(<ns>extent, <ns>extent):\n"
      "unable to allocate memory.");
    }
  else {		/* (NULL != pp)		*/
    *pp = NULL;		     ++pp;
    *((<ns>extent*)pp)	= m; ++pp;
    *((<ns>extent*)pp)	= n; ++pp;
    if (0 < m) {
      if (0 < n) {
	pp[0] = (<ns><p><s>scalar*)<ns><p><s>v_new(n*m);
	if (NULL == pp[0]) {
	  <ns>offset	i = 1;
          for (i = 1; i < m; ++i)
	    pp[i] = NULL;
	  <ns>message("In function <ns><p><s>m_new(<ns>extent, <ns>extent):\n"
	    "unable to allocate memory.");
	  }
	else {		/* (NULL != pp[0])	*/
	  <ns>offset	i = 1;
	  for (i = 1; i < m; ++i)
	    pp[i] = pp[i-1] + n;
	  }
	}
      else {		/* (n <= 0)		*/
	<ns>offset	i = 0;
        for (i = 0; i < m; ++i)
	  pp[i] = NULL;
	}
      }
    }
  return pp;
  }

/* Destructor								*/
void			(<ns><p><s>m_delete)(<ns><p><s>matrix* pM) {
  if (NULL != pM) {
    <ns><p><s>v_delete(pM[0]);
    <ns>free(pM - 3);
    }
  }

/* Reconstructor							*/
<ns><p><s>matrix*	(<ns><p><s>m_resize)(<ns><p><s>matrix* *ppM,
    <ns>extent m, <ns>extent n) {
  if (NULL != ppM) {
    (<ns><p><s>m_delete)(*ppM);
    *ppM = (<ns><p><s>m_new)(m, n);
    return *ppM;
    }
  else {
    <ns>message(
"In function <ns><p><s>m_resize(<ns><p><s>matrix**, <ns>extent, <ns>extent):\n"
"invalid pointer to <ns><p><s>matrix*.");
    return NULL;
    }
  }

