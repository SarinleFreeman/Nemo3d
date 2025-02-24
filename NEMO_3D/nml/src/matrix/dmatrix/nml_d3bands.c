
/*****************************************************************************
The NEMO Math Library.
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

#undef	NML_INLINE
#undef	NML_MACROS
#include<nml_d3bands.h>

/* Functions								*/
nml_extent		(nml_d3_extent)(const nml_d3bands* pT) {
  return *((const nml_extent*)(pT - 2));
  }
#define nml_d3_extent(pT) (*((const nml_extent*)((pT) - 2)))

nml_d3bands*	(nml_d3_fill)(
    nml_d3bands* pT, const nml_dscalar x) {
  nml_dv_fill(pT[+1], x);
  nml_dv_fill(pT[ 0], x);
  nml_dv_fill(pT[-1], x);
  return pT;
  }

nml_d3bands*	(nml_d3_smul)(
    nml_d3bands* pT, const nml_dscalar x) {
  nml_dv_smul(pT[+1], x);
  nml_dv_smul(pT[ 0], x);
  nml_dv_smul(pT[-1], x);
  return pT;
  }

nml_d3bands*	(nml_d3_sdiv)(
    nml_d3bands* pT, const nml_dscalar x) {
  nml_dv_sdiv(pT[+1], x);
  nml_dv_sdiv(pT[ 0], x);
  nml_dv_sdiv(pT[-1], x);
  return pT;
  }

nml_d3bands*	(nml_d3_sadd)(
    nml_d3bands* pT, const nml_dscalar x) {
  nml_dv_sadd(pT[+1], x);
  nml_dv_sadd(pT[ 0], x);
  nml_dv_sadd(pT[-1], x);
  return pT;
  }

nml_d3bands*	(nml_d3_ssub)(
    nml_d3bands* pT, const nml_dscalar x) {
  nml_dv_ssub(pT[+1], x);
  nml_dv_ssub(pT[ 0], x);
  nml_dv_ssub(pT[-1], x);
  return pT;
  }

nml_d3bands*	(nml_d3_3cpy)(
    nml_d3bands* pT, const nml_d3bands* pU) {
#ifdef	NML_ERRANT
  if (nml_d3_extent(pT) != nml_d3_extent(pU))
    nml_message(
"In function nml_d3_3cpy(nml_d3bands*, const nml_d3bands*),\n"
"unequal tridiagonal matrix extents.\n");
#endif/*NML_ERRANT			*/
  nml_dv_vcpy(pT[+1], pU[+1]);
  nml_dv_vcpy(pT[ 0], pU[ 0]);
  nml_dv_vcpy(pT[-1], pU[-1]);
  return pT;
  }

nml_dvector* nml_dv_v3dot(
    nml_dvector*	py,	/*   out product     vector		*/
    const
    nml_dvector*	px,	/* in    multiplier  vector		*/
    const
    nml_d3bands*	pT	/* in    tridiagonal matrix		*/
    ) {				/* y <-- xT^{T}				*/
  const nml_extent	m = nml_dv_extent(py);
#ifdef	NML_ERRANT
  if (nml_dv_extent(px) != m)
    nml_message(
	"In function nml_dv_3vdot(nml_dvector*,\n"
	"  const nml_d3bands*, const nml_dvector*),\n"
        "unequal product vector and multiplier vector extents.\n");
  if (nml_d3_extent(pT) != m)
    nml_message(
	"In function nml_dv_3vdot(nml_dvector*,\n"
	"  const nml_d3bands*, const nml_dvector*),\n"
        "unequal product vector and tridiagonal matrix extents.\n");
#endif/*NML_ERRANT			*/
  if (0 < m)
    if (1 < m) {
      nml_dscalar	p, q, r;		/* intermediate results	*/
      /* y[0] = T[0][0]*x[0] + T[1][0]*x[1];				*/
      NML_DMUL(r, pT[1][0], px[1]);	/* r = T[1][0]*x[1]	*/
      NML_DMUL(q, pT[0][0], px[0]);	/* q = T[0][0]*x[0]	*/
      NML_DADD(py[0], q, r);		/* y[0] = q + r		*/
      { nml_offset	j = 0;
	for (j = 1; j < m-1; ++j) {
	  /* y[j] = T[-1][j]*x[j-1] + T[0][j]*x[j] + T[1][j]*x[j+1];	*/
	  NML_DMUL(r, pT[+1][j], px[j+1]);	/* r = T[+1][j]*x[j+1]	*/
	  NML_DMUL(q, pT[ 0][j], px[j]);	/* q = T[ 0][j]*x[j]	*/
	  NML_DMUL(p, pT[-1][j], px[j-1]);	/* p = T[-1][j]*x[j-1]	*/
	  NML_DADD(q, q, r);		/* q += r		*/
	  NML_DADD(py[j], p, q);		/* y[j] = p + q		*/
	  }
	/* y[m-1] = T[-1][m-1]*x[m-2] + T[0][m-1]*x[m-1];		*/
	NML_DMUL(q, pT[ 0][m-1], px[m-1]);	/* q = T[ 0][m-1]*x[m-1]*/
	NML_DMUL(p, pT[-1][m-1], px[m-2]);	/* p = T[-1][m-1]*x[m-1]*/
	NML_DADD(py[m-1], p, q);		/* y[m-1] = p + q	*/
	}
      }
    else {					/* (1 == m)		*/
      NML_DMUL(py[0], pT[0][0], px[0]);	/* y[0] = T[0][0]*x[0];	*/
      }
  return py;
  }

nml_dvector* nml_dv_solveTridiagonal(
    nml_dvector*	px,	/*   out    solution vector		*/
    const
    nml_d3bands*	pT,	/* in    tridiagonal matrix		*/
    const
    nml_dvector*	py	/* in        product vector		*/
    ) {	/* Solve xT^{T} = y  for x where T is a tridiagonal matrix.	*/
  const nml_extent	m = nml_dv_extent(px);
#ifdef	NML_ERRANT
  if (nml_d3_extent(pT) != m)
    nml_message(
	"In function nml_dv_solveTridiagonal(nml_dvector*,\n"
	"  const nml_d3bands*, const nml_dvector*),\n"
        "unequal solution vector and tridiagonal matrix extents.\n");
  if (nml_dv_extent(py) != m)
    nml_message(
	"In function nml_dv_solveTridiagonal(nml_dvector*,\n"
	"  const nml_d3bands*, const nml_dvector*),\n"
        "unequal solution vector and product vector extents.\n");
#endif/*NML_ERRANT			*/
  if (NML_DNE(NML_DZERO, pT[0][0])) {	/* (0 != pT[0][0])	*/
    nml_dscalar	d = pT[0][0];
    nml_dvector*	pg = nml_dv_new(m);
    nml_offset		j = 0;
    NML_DDIV(px[0], py[0], d);		/* x[0] = y[0]/d;	*/
    for (j = 1; j < m; ++j) {
      /* d = T[0][j] - T[-1][j]*g[j];					*/
      nml_dscalar	t;			/* intermediate result	*/
      NML_DDIV(pg[j], pT[1][j-1], d);	/* g[j] = T[1][j-1]/d;	*/
      NML_DMUL(t, pT[-1][j], pg[j]);	/* t = T[-1][j]*g[j]	*/
      NML_DSUB(d, pT[0][j], d);		/* d = T[0][j] - t	*/
      if (NML_DNE(NML_DZERO, d)) {	/* 0 == d		*/
	nml_message(
	"In function nml_dv_solveTridiagonal(nml_dvector*,\n"
	"  const nml_d3bands*, const nml_dvector*),\n"
        "singular tridiagonal matrix.\n");
	return px;
	}
      /* px[j] = (py[j] - pT[-1][j]*px[j-1])/d;				*/
      NML_DMUL(t, pT[-1][j], px[j-1]);	/* t = pT[-1][j]*px[j-1]*/
      NML_DSUB(t, py[j], t);		/* t = py[j] - t	*/
      NML_DDIV(px[j], t, d);		/* px[j] = t/d		*/ 
      }
    for (j = m-1; j > 0; --j) {
      /* x[j-1] -= g[j]*x[j];						*/
      nml_dscalar	t;			/* intermediate result	*/
      NML_DMUL(t, pg[j], px[j]);		/* t = g[j]*x[j]	*/
      NML_DSUB(px[j-1], px[j-1], t);	/* x[j-1] -= t		*/
      }

    nml_dv_delete(pg);
    }
  else {					/* (0 == pT[0][0])	*/
    nml_message(
	"In function nml_dv_solveTridiagonal(nml_dvector*,\n"
	"  const nml_d3bands*, const nml_dvector*),\n"
        "singular tridiagonal matrix.\n");
    }
  return px;
  }


/* Constructor								*/
nml_d3bands*	(nml_d3_new)(nml_extent n) {
  nml_dscalar**	pT
    = (nml_dscalar**)nml_malloc(4*sizeof(nml_dscalar*));
  if (NULL == pT) {
    nml_message("In function nml_d3_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {		/* (NULL != pT)		*/
    *((nml_extent*)pT) = n;	++pT;
    *pT = nml_dv_new(n);	++pT;	/* T[-1]			*/
    *pT = nml_dv_new(n);	++pT;	/* T[ 0]			*/
    *pT = nml_dv_new(n);	--pT;	/* T[+1]			*/
    }
  return pT;
  }

nml_d3bands*	(nml_d3_clone)(const nml_d3bands* pU) {
  return nml_d3_3cpy(nml_d3_new(nml_d3_extent(pU)), pU);
  }
#define nml_d3_clone(pT) \
	nml_d3_3cpy(nml_d3_new(nml_dextent(pT)), (pT))

/* Destructor								*/
void			(nml_d3_delete)(nml_d3bands* pT) {
  if (NULL != pT) {
    ++pT; nml_dv_delete(*pT);	/* T[+1]			*/
    --pT; nml_dv_delete(*pT);	/* T[ 0]			*/
    --pT; nml_dv_delete(*pT);	/* T[-1]			*/
    --pT; nml_free(pT);
    }
  }

/* Reconstructor							*/
nml_d3bands*	(nml_d3_resize)(nml_d3bands* *ppT,
    nml_extent n) {
  if (NULL != ppT) {
    (nml_d3_delete)(*ppT);
    *ppT = (nml_d3_new)(n);
    return *ppT;
    }
  else {
    nml_message(
	"In function nml_d3_resize(nml_d3bands**, nml_extent):\n"
	"invalid pointer to nml_d3bands*.");
    return NULL;
    }
  }

