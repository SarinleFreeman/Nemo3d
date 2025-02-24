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
$Header: /repo/nml/src/matrix/3bands.cP,v 1.2 2003/10/08 16:19:32 hook Exp $
*****************************************************************************/

#undef	<NS>INLINE
#undef	<NS>MACROS
#include<<ns><p><s>3bands.h>

/* Functions								*/
<ns>extent		(<ns><p><s>3_extent)(const <ns><p><s>3bands* pT) {
  return *((const <ns>extent*)(pT - 2));
  }
#define <ns><p><s>3_extent(pT) (*((const <ns>extent*)((pT) - 2)))

<ns><p><s>3bands*	(<ns><p><s>3_fill)(
    <ns><p><s>3bands* pT, const <ns><p><s>scalar x) {
  <ns><p><s>v_fill(pT[+1], x);
  <ns><p><s>v_fill(pT[ 0], x);
  <ns><p><s>v_fill(pT[-1], x);
  return pT;
  }

<ns><p><s>3bands*	(<ns><p><s>3_smul)(
    <ns><p><s>3bands* pT, const <ns><p><s>scalar x) {
  <ns><p><s>v_smul(pT[+1], x);
  <ns><p><s>v_smul(pT[ 0], x);
  <ns><p><s>v_smul(pT[-1], x);
  return pT;
  }

<ns><p><s>3bands*	(<ns><p><s>3_sdiv)(
    <ns><p><s>3bands* pT, const <ns><p><s>scalar x) {
  <ns><p><s>v_sdiv(pT[+1], x);
  <ns><p><s>v_sdiv(pT[ 0], x);
  <ns><p><s>v_sdiv(pT[-1], x);
  return pT;
  }

<ns><p><s>3bands*	(<ns><p><s>3_sadd)(
    <ns><p><s>3bands* pT, const <ns><p><s>scalar x) {
  <ns><p><s>v_sadd(pT[+1], x);
  <ns><p><s>v_sadd(pT[ 0], x);
  <ns><p><s>v_sadd(pT[-1], x);
  return pT;
  }

<ns><p><s>3bands*	(<ns><p><s>3_ssub)(
    <ns><p><s>3bands* pT, const <ns><p><s>scalar x) {
  <ns><p><s>v_ssub(pT[+1], x);
  <ns><p><s>v_ssub(pT[ 0], x);
  <ns><p><s>v_ssub(pT[-1], x);
  return pT;
  }

<ns><p><s>3bands*	(<ns><p><s>3_3cpy)(
    <ns><p><s>3bands* pT, const <ns><p><s>3bands* pU) {
#ifdef	<NS>ERRANT
  if (<ns><p><s>3_extent(pT) != <ns><p><s>3_extent(pU))
    <ns>message(
"In function <ns><p><s>3_3cpy(<ns><p><s>3bands*, const <ns><p><s>3bands*),\n"
"unequal tridiagonal matrix extents.\n");
#endif/*<NS>ERRANT			*/
  <ns><p><s>v_vcpy(pT[+1], pU[+1]);
  <ns><p><s>v_vcpy(pT[ 0], pU[ 0]);
  <ns><p><s>v_vcpy(pT[-1], pU[-1]);
  return pT;
  }

<ns><p><s>vector* <ns><p><s>v_v3dot(
    <ns><p><s>vector*	py,	/*   out product     vector		*/
    const
    <ns><p><s>vector*	px,	/* in    multiplier  vector		*/
    const
    <ns><p><s>3bands*	pT	/* in    tridiagonal matrix		*/
    ) {				/* y <-- xT^{T}				*/
  const <ns>extent	m = <ns><p><s>v_extent(py);
#ifdef	<NS>ERRANT
  if (<ns><p><s>v_extent(px) != m)
    <ns>message(
	"In function <ns><p><s>v_3vdot(<ns><p><s>vector*,\n"
	"  const <ns><p><s>3bands*, const <ns><p><s>vector*),\n"
        "unequal product vector and multiplier vector extents.\n");
  if (<ns><p><s>3_extent(pT) != m)
    <ns>message(
	"In function <ns><p><s>v_3vdot(<ns><p><s>vector*,\n"
	"  const <ns><p><s>3bands*, const <ns><p><s>vector*),\n"
        "unequal product vector and tridiagonal matrix extents.\n");
#endif/*<NS>ERRANT			*/
  if (0 < m)
    if (1 < m) {
      <ns><p><s>scalar	p, q, r;		/* intermediate results	*/
      /* y[0] = T[0][0]*x[0] + T[1][0]*x[1];				*/
      <NS><P><S>MUL(r, pT[1][0], px[1]);	/* r = T[1][0]*x[1]	*/
      <NS><P><S>MUL(q, pT[0][0], px[0]);	/* q = T[0][0]*x[0]	*/
      <NS><P><S>ADD(py[0], q, r);		/* y[0] = q + r		*/
      { <ns>offset	j = 0;
	for (j = 1; j < m-1; ++j) {
	  /* y[j] = T[-1][j]*x[j-1] + T[0][j]*x[j] + T[1][j]*x[j+1];	*/
	  <NS><P><S>MUL(r, pT[+1][j], px[j+1]);	/* r = T[+1][j]*x[j+1]	*/
	  <NS><P><S>MUL(q, pT[ 0][j], px[j]);	/* q = T[ 0][j]*x[j]	*/
	  <NS><P><S>MUL(p, pT[-1][j], px[j-1]);	/* p = T[-1][j]*x[j-1]	*/
	  <NS><P><S>ADD(q, q, r);		/* q += r		*/
	  <NS><P><S>ADD(py[j], p, q);		/* y[j] = p + q		*/
	  }
	/* y[m-1] = T[-1][m-1]*x[m-2] + T[0][m-1]*x[m-1];		*/
	<NS><P><S>MUL(q, pT[ 0][m-1], px[m-1]);	/* q = T[ 0][m-1]*x[m-1]*/
	<NS><P><S>MUL(p, pT[-1][m-1], px[m-2]);	/* p = T[-1][m-1]*x[m-1]*/
	<NS><P><S>ADD(py[m-1], p, q);		/* y[m-1] = p + q	*/
	}
      }
    else {					/* (1 == m)		*/
      <NS><P><S>MUL(py[0], pT[0][0], px[0]);	/* y[0] = T[0][0]*x[0];	*/
      }
  return py;
  }

<ns><p><s>vector* <ns><p><s>v_solveTridiagonal(
    <ns><p><s>vector*	px,	/*   out    solution vector		*/
    const
    <ns><p><s>3bands*	pT,	/* in    tridiagonal matrix		*/
    const
    <ns><p><s>vector*	py	/* in        product vector		*/
    ) {	/* Solve xT^{T} = y  for x where T is a tridiagonal matrix.	*/
  const <ns>extent	m = <ns><p><s>v_extent(px);
#ifdef	<NS>ERRANT
  if (<ns><p><s>3_extent(pT) != m)
    <ns>message(
	"In function <ns><p><s>v_solveTridiagonal(<ns><p><s>vector*,\n"
	"  const <ns><p><s>3bands*, const <ns><p><s>vector*),\n"
        "unequal solution vector and tridiagonal matrix extents.\n");
  if (<ns><p><s>v_extent(py) != m)
    <ns>message(
	"In function <ns><p><s>v_solveTridiagonal(<ns><p><s>vector*,\n"
	"  const <ns><p><s>3bands*, const <ns><p><s>vector*),\n"
        "unequal solution vector and product vector extents.\n");
#endif/*<NS>ERRANT			*/
  if (<NS><P><S>NE(<NS><P><S>ZERO, pT[0][0])) {	/* (0 != pT[0][0])	*/
    <ns><p><s>scalar	d = pT[0][0];
    <ns><p><s>vector*	pg = <ns><p><s>v_new(m);
    <ns>offset		j = 0;
    <NS><P><S>DIV(px[0], py[0], d);		/* x[0] = y[0]/d;	*/
    for (j = 1; j < m; ++j) {
      /* d = T[0][j] - T[-1][j]*g[j];					*/
      <ns><p><s>scalar	t;			/* intermediate result	*/
      <NS><P><S>DIV(pg[j], pT[1][j-1], d);	/* g[j] = T[1][j-1]/d;	*/
      <NS><P><S>MUL(t, pT[-1][j], pg[j]);	/* t = T[-1][j]*g[j]	*/
      <NS><P><S>SUB(d, pT[0][j], d);		/* d = T[0][j] - t	*/
      if (<NS><P><S>NE(<NS><P><S>ZERO, d)) {	/* 0 == d		*/
	nml_message(
	"In function <ns><p><s>v_solveTridiagonal(<ns><p><s>vector*,\n"
	"  const <ns><p><s>3bands*, const <ns><p><s>vector*),\n"
        "singular tridiagonal matrix.\n");
	return px;
	}
      /* px[j] = (py[j] - pT[-1][j]*px[j-1])/d;				*/
      <NS><P><S>MUL(t, pT[-1][j], px[j-1]);	/* t = pT[-1][j]*px[j-1]*/
      <NS><P><S>SUB(t, py[j], t);		/* t = py[j] - t	*/
      <NS><P><S>DIV(px[j], t, d);		/* px[j] = t/d		*/ 
      }
    for (j = m-1; j > 0; --j) {
      /* x[j-1] -= g[j]*x[j];						*/
      <ns><p><s>scalar	t;			/* intermediate result	*/
      <NS><P><S>MUL(t, pg[j], px[j]);		/* t = g[j]*x[j]	*/
      <NS><P><S>SUB(px[j-1], px[j-1], t);	/* x[j-1] -= t		*/
      }

    <ns><p><s>v_delete(pg);
    }
  else {					/* (0 == pT[0][0])	*/
    <ns>message(
	"In function <ns><p><s>v_solveTridiagonal(<ns><p><s>vector*,\n"
	"  const <ns><p><s>3bands*, const <ns><p><s>vector*),\n"
        "singular tridiagonal matrix.\n");
    }
  return px;
  }


/* Constructor								*/
<ns><p><s>3bands*	(<ns><p><s>3_new)(<ns>extent n) {
  <ns><p><s>scalar**	pT
    = (<ns><p><s>scalar**)<ns>malloc(4*sizeof(<ns><p><s>scalar*));
  if (NULL == pT) {
    <ns>message("In function <ns><p><s>3_new(<ns>extent):\n"
      "unable to allocate memory.");
    }
  else {		/* (NULL != pT)		*/
    *((<ns>extent*)pT) = n;	++pT;
    *pT = <ns><p><s>v_new(n);	++pT;	/* T[-1]			*/
    *pT = <ns><p><s>v_new(n);	++pT;	/* T[ 0]			*/
    *pT = <ns><p><s>v_new(n);	--pT;	/* T[+1]			*/
    }
  return pT;
  }

<ns><p><s>3bands*	(<ns><p><s>3_clone)(const <ns><p><s>3bands* pU) {
  return <ns><p><s>3_3cpy(<ns><p><s>3_new(<ns><p><s>3_extent(pU)), pU);
  }
#define <ns><p><s>3_clone(pT) \
	<ns><p><s>3_3cpy(<ns><p><s>3_new(<ns><p><s>extent(pT)), (pT))

/* Destructor								*/
void			(<ns><p><s>3_delete)(<ns><p><s>3bands* pT) {
  if (NULL != pT) {
    ++pT; <ns><p><s>v_delete(*pT);	/* T[+1]			*/
    --pT; <ns><p><s>v_delete(*pT);	/* T[ 0]			*/
    --pT; <ns><p><s>v_delete(*pT);	/* T[-1]			*/
    --pT; <ns>free(pT);
    }
  }

/* Reconstructor							*/
<ns><p><s>3bands*	(<ns><p><s>3_resize)(<ns><p><s>3bands* *ppT,
    <ns>extent n) {
  if (NULL != ppT) {
    (<ns><p><s>3_delete)(*ppT);
    *ppT = (<ns><p><s>3_new)(n);
    return *ppT;
    }
  else {
    <ns>message(
	"In function <ns><p><s>3_resize(<ns><p><s>3bands**, <ns>extent):\n"
	"invalid pointer to <ns><p><s>3bands*.");
    return NULL;
    }
  }

