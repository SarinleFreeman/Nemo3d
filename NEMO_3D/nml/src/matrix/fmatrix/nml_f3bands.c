
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
#include<nml_f3bands.h>

/* Functions								*/
nml_extent		(nml_f3_extent)(const nml_f3bands* pT) {
  return *((const nml_extent*)(pT - 2));
  }
#define nml_f3_extent(pT) (*((const nml_extent*)((pT) - 2)))

nml_f3bands*	(nml_f3_fill)(
    nml_f3bands* pT, const nml_fscalar x) {
  nml_fv_fill(pT[+1], x);
  nml_fv_fill(pT[ 0], x);
  nml_fv_fill(pT[-1], x);
  return pT;
  }

nml_f3bands*	(nml_f3_smul)(
    nml_f3bands* pT, const nml_fscalar x) {
  nml_fv_smul(pT[+1], x);
  nml_fv_smul(pT[ 0], x);
  nml_fv_smul(pT[-1], x);
  return pT;
  }

nml_f3bands*	(nml_f3_sdiv)(
    nml_f3bands* pT, const nml_fscalar x) {
  nml_fv_sdiv(pT[+1], x);
  nml_fv_sdiv(pT[ 0], x);
  nml_fv_sdiv(pT[-1], x);
  return pT;
  }

nml_f3bands*	(nml_f3_sadd)(
    nml_f3bands* pT, const nml_fscalar x) {
  nml_fv_sadd(pT[+1], x);
  nml_fv_sadd(pT[ 0], x);
  nml_fv_sadd(pT[-1], x);
  return pT;
  }

nml_f3bands*	(nml_f3_ssub)(
    nml_f3bands* pT, const nml_fscalar x) {
  nml_fv_ssub(pT[+1], x);
  nml_fv_ssub(pT[ 0], x);
  nml_fv_ssub(pT[-1], x);
  return pT;
  }

nml_f3bands*	(nml_f3_3cpy)(
    nml_f3bands* pT, const nml_f3bands* pU) {
#ifdef	NML_ERRANT
  if (nml_f3_extent(pT) != nml_f3_extent(pU))
    nml_message(
"In function nml_f3_3cpy(nml_f3bands*, const nml_f3bands*),\n"
"unequal tridiagonal matrix extents.\n");
#endif/*NML_ERRANT			*/
  nml_fv_vcpy(pT[+1], pU[+1]);
  nml_fv_vcpy(pT[ 0], pU[ 0]);
  nml_fv_vcpy(pT[-1], pU[-1]);
  return pT;
  }

nml_fvector* nml_fv_v3dot(
    nml_fvector*	py,	/*   out product     vector		*/
    const
    nml_fvector*	px,	/* in    multiplier  vector		*/
    const
    nml_f3bands*	pT	/* in    tridiagonal matrix		*/
    ) {				/* y <-- xT^{T}				*/
  const nml_extent	m = nml_fv_extent(py);
#ifdef	NML_ERRANT
  if (nml_fv_extent(px) != m)
    nml_message(
	"In function nml_fv_3vdot(nml_fvector*,\n"
	"  const nml_f3bands*, const nml_fvector*),\n"
        "unequal product vector and multiplier vector extents.\n");
  if (nml_f3_extent(pT) != m)
    nml_message(
	"In function nml_fv_3vdot(nml_fvector*,\n"
	"  const nml_f3bands*, const nml_fvector*),\n"
        "unequal product vector and tridiagonal matrix extents.\n");
#endif/*NML_ERRANT			*/
  if (0 < m)
    if (1 < m) {
      nml_fscalar	p, q, r;		/* intermediate results	*/
      /* y[0] = T[0][0]*x[0] + T[1][0]*x[1];				*/
      NML_FMUL(r, pT[1][0], px[1]);	/* r = T[1][0]*x[1]	*/
      NML_FMUL(q, pT[0][0], px[0]);	/* q = T[0][0]*x[0]	*/
      NML_FADD(py[0], q, r);		/* y[0] = q + r		*/
      { nml_offset	j = 0;
	for (j = 1; j < m-1; ++j) {
	  /* y[j] = T[-1][j]*x[j-1] + T[0][j]*x[j] + T[1][j]*x[j+1];	*/
	  NML_FMUL(r, pT[+1][j], px[j+1]);	/* r = T[+1][j]*x[j+1]	*/
	  NML_FMUL(q, pT[ 0][j], px[j]);	/* q = T[ 0][j]*x[j]	*/
	  NML_FMUL(p, pT[-1][j], px[j-1]);	/* p = T[-1][j]*x[j-1]	*/
	  NML_FADD(q, q, r);		/* q += r		*/
	  NML_FADD(py[j], p, q);		/* y[j] = p + q		*/
	  }
	/* y[m-1] = T[-1][m-1]*x[m-2] + T[0][m-1]*x[m-1];		*/
	NML_FMUL(q, pT[ 0][m-1], px[m-1]);	/* q = T[ 0][m-1]*x[m-1]*/
	NML_FMUL(p, pT[-1][m-1], px[m-2]);	/* p = T[-1][m-1]*x[m-1]*/
	NML_FADD(py[m-1], p, q);		/* y[m-1] = p + q	*/
	}
      }
    else {					/* (1 == m)		*/
      NML_FMUL(py[0], pT[0][0], px[0]);	/* y[0] = T[0][0]*x[0];	*/
      }
  return py;
  }

nml_fvector* nml_fv_solveTridiagonal(
    nml_fvector*	px,	/*   out    solution vector		*/
    const
    nml_f3bands*	pT,	/* in    tridiagonal matrix		*/
    const
    nml_fvector*	py	/* in        product vector		*/
    ) {	/* Solve xT^{T} = y  for x where T is a tridiagonal matrix.	*/
  const nml_extent	m = nml_fv_extent(px);
#ifdef	NML_ERRANT
  if (nml_f3_extent(pT) != m)
    nml_message(
	"In function nml_fv_solveTridiagonal(nml_fvector*,\n"
	"  const nml_f3bands*, const nml_fvector*),\n"
        "unequal solution vector and tridiagonal matrix extents.\n");
  if (nml_fv_extent(py) != m)
    nml_message(
	"In function nml_fv_solveTridiagonal(nml_fvector*,\n"
	"  const nml_f3bands*, const nml_fvector*),\n"
        "unequal solution vector and product vector extents.\n");
#endif/*NML_ERRANT			*/
  if (NML_FNE(NML_FZERO, pT[0][0])) {	/* (0 != pT[0][0])	*/
    nml_fscalar	d = pT[0][0];
    nml_fvector*	pg = nml_fv_new(m);
    nml_offset		j = 0;
    NML_FDIV(px[0], py[0], d);		/* x[0] = y[0]/d;	*/
    for (j = 1; j < m; ++j) {
      /* d = T[0][j] - T[-1][j]*g[j];					*/
      nml_fscalar	t;			/* intermediate result	*/
      NML_FDIV(pg[j], pT[1][j-1], d);	/* g[j] = T[1][j-1]/d;	*/
      NML_FMUL(t, pT[-1][j], pg[j]);	/* t = T[-1][j]*g[j]	*/
      NML_FSUB(d, pT[0][j], d);		/* d = T[0][j] - t	*/
      if (NML_FNE(NML_FZERO, d)) {	/* 0 == d		*/
	nml_message(
	"In function nml_fv_solveTridiagonal(nml_fvector*,\n"
	"  const nml_f3bands*, const nml_fvector*),\n"
        "singular tridiagonal matrix.\n");
	return px;
	}
      /* px[j] = (py[j] - pT[-1][j]*px[j-1])/d;				*/
      NML_FMUL(t, pT[-1][j], px[j-1]);	/* t = pT[-1][j]*px[j-1]*/
      NML_FSUB(t, py[j], t);		/* t = py[j] - t	*/
      NML_FDIV(px[j], t, d);		/* px[j] = t/d		*/ 
      }
    for (j = m-1; j > 0; --j) {
      /* x[j-1] -= g[j]*x[j];						*/
      nml_fscalar	t;			/* intermediate result	*/
      NML_FMUL(t, pg[j], px[j]);		/* t = g[j]*x[j]	*/
      NML_FSUB(px[j-1], px[j-1], t);	/* x[j-1] -= t		*/
      }

    nml_fv_delete(pg);
    }
  else {					/* (0 == pT[0][0])	*/
    nml_message(
	"In function nml_fv_solveTridiagonal(nml_fvector*,\n"
	"  const nml_f3bands*, const nml_fvector*),\n"
        "singular tridiagonal matrix.\n");
    }
  return px;
  }


/* Constructor								*/
nml_f3bands*	(nml_f3_new)(nml_extent n) {
  nml_fscalar**	pT
    = (nml_fscalar**)nml_malloc(4*sizeof(nml_fscalar*));
  if (NULL == pT) {
    nml_message("In function nml_f3_new(nml_extent):\n"
      "unable to allocate memory.");
    }
  else {		/* (NULL != pT)		*/
    *((nml_extent*)pT) = n;	++pT;
    *pT = nml_fv_new(n);	++pT;	/* T[-1]			*/
    *pT = nml_fv_new(n);	++pT;	/* T[ 0]			*/
    *pT = nml_fv_new(n);	--pT;	/* T[+1]			*/
    }
  return pT;
  }

nml_f3bands*	(nml_f3_clone)(const nml_f3bands* pU) {
  return nml_f3_3cpy(nml_f3_new(nml_f3_extent(pU)), pU);
  }
#define nml_f3_clone(pT) \
	nml_f3_3cpy(nml_f3_new(nml_fextent(pT)), (pT))

/* Destructor								*/
void			(nml_f3_delete)(nml_f3bands* pT) {
  if (NULL != pT) {
    ++pT; nml_fv_delete(*pT);	/* T[+1]			*/
    --pT; nml_fv_delete(*pT);	/* T[ 0]			*/
    --pT; nml_fv_delete(*pT);	/* T[-1]			*/
    --pT; nml_free(pT);
    }
  }

/* Reconstructor							*/
nml_f3bands*	(nml_f3_resize)(nml_f3bands* *ppT,
    nml_extent n) {
  if (NULL != ppT) {
    (nml_f3_delete)(*ppT);
    *ppT = (nml_f3_new)(n);
    return *ppT;
    }
  else {
    nml_message(
	"In function nml_f3_resize(nml_f3bands**, nml_extent):\n"
	"invalid pointer to nml_f3bands*.");
    return NULL;
    }
  }

