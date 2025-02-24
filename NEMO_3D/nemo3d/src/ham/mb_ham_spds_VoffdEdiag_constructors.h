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
$Header:
*****************************************************************************/

/*
 * Inline code that returns the off-diagonal matrix elements
 * as well as the additional corrections to the diagonal elements
 * due to strain.
 */

/* STa_ */
#define get_STa_STc(Vmat, Ecorr) \
	    r1 = Vmat  \
		=     ststs; \
	    r2	= pVp_ststs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stst; \
	    Ecorr += r3;
#define get_STa_Sc(Vmat, Ecorr) \
	    r1 = Vmat  \
		=     stass; \
	    r2	= pVp_stass; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stas; \
	    Ecorr += r3;
#define get_STa_Pxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*     stapcs; \
	    r2	= l0*pVp_stapcs; \
 \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stapc; \
	    Ecorr += r3;
#define get_STa_Pyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*     stapcs; \
	    r2	= m0*pVp_stapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stapc; \
	    Ecorr += r3;
#define get_STa_Pzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*     stapcs; \
	    r2	= n0*pVp_stapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stapc; \
	    Ecorr += r3;
#define get_STa_Dxyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*      stadcs; \
	    r2	= sqr3*l0*m0*pVp_stadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stadc; \
	    Ecorr += r3;
#define get_STa_Dyzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*      stadcs; \
	    r2	= sqr3*m0*n0*pVp_stadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stadc; \
	    Ecorr += r3;
#define get_STa_Dzxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*n* l*     stadcs; \
	    r2	= sqr3*n0*l0*pVp_stadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stadc; \
	    Ecorr += r3;
#define get_STa_Dx2my2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2- m2)*     stadcs; \
	    r2	= 0.5*sqr3*(l02-m02)*pVp_stadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stadc; \
	    Ecorr += r3;
#define get_STa_Dz2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= (n2 -0.5*(l2 +m2)) *    stadcs; \
	    r2	= (n02-0.5*(l02+m02))*pVp_stadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stadc; \
	    Ecorr += r3;
/* Sa_ */
#define get_Sa_STc(Vmat, Ecorr) \
	    r1 = Vmat  \
		=     sasts; \
	    r2	= pVp_sasts; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sast; \
	    Ecorr += r3;
#define get_Sa_Sc(Vmat, Ecorr) \
	    r1 = Vmat  \
		=     sss; \
	    r2	= pVp_sss; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_ss; \
	    Ecorr += r3;
#define get_Sa_Pxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*     sapcs; \
	    r2	= l0*pVp_sapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sapc; \
	    Ecorr += r3;
#define get_Sa_Pyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*     sapcs; \
	    r2	= m0*pVp_sapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sapc; \
	    Ecorr += r3;
#define get_Sa_Pzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*     sapcs; \
	    r2	= n0*pVp_sapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sapc; \
	    Ecorr += r3;
#define get_Sa_Dxyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l* m*     sadcs; \
	    r2	= sqr3*l0*m0*pVp_sadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sadc; \
	    Ecorr += r3;
#define get_Sa_Dyzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m *n *    sadcs; \
	    r2	= sqr3*m0*n0*pVp_sadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sadc; \
	    Ecorr += r3;
#define get_Sa_Dzxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l *n *    sadcs; \
	    r2	= sqr3*l0*n0*pVp_sadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sadc; \
	    Ecorr += r3;
#define get_Sa_Dx2my2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2 -m2) *    sadcs; \
	    r2	= 0.5*sqr3*(l02-m02)*pVp_sadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sadc; \
	    Ecorr += r3;
#define get_Sa_Dz2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= (n2 -0.5*(l2 +m2)) *    sadcs; \
	    r2	= (n02-0.5*(l02+m02))*pVp_sadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sadc; \
	    Ecorr += r3;
/* Pxa_ */
#define get_Pxa_STc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*     stcpas; \
	    r2	= l0*pVp_stcpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcpa; \
	    Ecorr += r3;
#define get_Pxa_Sc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*     scpas; \
	    r2	= l0*pVp_scpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scpa; \
	    Ecorr += r3;
#define get_Pxa_Pxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l2*         pps  \
		+ (1-l2)*     ppp; \
	    r2	= l02*    pVp_pps  \
		+ (1-l02)*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pxa_Pyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*m*      pps  \
		- l*m*      ppp; \
	    r2	= l0*m0*pVp_pps  \
		- l0*m0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pxa_Pzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*n*      pps  \
		- l*n*      ppp; \
	    r2	= l0*n0*pVp_pps  \
		- l0*n0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pxa_Dxyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l2*m*       padcs  \
		+ m*(1-2*l2)*      padcp; \
	    r2	= sqr3*l02*m0* pVp_padcs  \
		+ m0*(1-2*l02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pxa_Dyzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*       padcs  \
		- 2*l*m*n*          padcp; \
	    r2	= sqr3*l0*m0*n0*pVp_padcs  \
		- 2*l0*m0*n0*   pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pxa_Dzxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l2*n*       padcs  \
		+ n*(1-2*l2)*      padcp; \
	    r2	= sqr3*l02*n0* pVp_padcs  \
		+ n0*(1-2*l02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pxa_Dx2my2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*l*(l2-m2)*       padcs  \
		+ l*(1-l2+m2)*              padcp; \
	    r2	= 0.5*sqr3*l0*(l02-m02)*pVp_padcs  \
		+ l0*(1-l02+m02)*       pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pxa_Dz2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*(n2-0.5*(l2+m2))*        padcs  \
		- sqr3*l*n2*                 padcp; \
	    r2	= l0*(n02-0.5*(l02+m02))*pVp_padcs  \
		- sqr3*l0*n02*           pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
/* Pya_ */
#define get_Pya_STc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*     stcpas; \
	    r2	= m0*pVp_stcpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcpa; \
	    Ecorr += r3;
#define get_Pya_Sc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*     scpas; \
	    r2	= m0*pVp_scpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scpa; \
	    Ecorr += r3;
#define get_Pya_Pxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*m*      pps  \
		- l*m*      ppp; \
	    r2	= l0*m0*pVp_pps  \
		- l0*m0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pya_Pyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m2*         pps  \
		+ (1-m2)*     ppp; \
	    r2	= m02*    pVp_pps  \
		+ (1-m02)*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pya_Pzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*n*      pps  \
		- m*n*      ppp; \
	    r2	= m0*n0*pVp_pps  \
		- m0*n0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pya_Dxyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m2*       padcs  \
		+ l*(1-2*m2)*      padcp; \
	    r2	= sqr3*l0*m02* pVp_padcs  \
		+ l0*(1-2*m02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pya_Dyzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*n*m2*       padcs  \
		+ n*(1-2*m2)*      padcp; \
	    r2	= sqr3*n0*m02* pVp_padcs  \
		+ n0*(1-2*m02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pya_Dzxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*       padcs  \
		- 2*l*m*n*          padcp; \
	    r2	= sqr3*l0*m0*n0*pVp_padcs  \
		- 2*l0*m0*n0*   pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pya_Dx2my2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*m*(l2-m2)*       padcs  \
		- m*(1+l2-m2)*              padcp; \
	    r2	= 0.5*sqr3*m0*(l02-m02)*pVp_padcs  \
		- m0*(1+l02-m02)*       pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pya_Dz2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*(n2-0.5*(l2+m2))*        padcs  \
		- sqr3*m*n2*                 padcp; \
	    r2	= m0*(n02-0.5*(l02+m02))*pVp_padcs  \
		- sqr3*m0*n02*           pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
/* Pza_ */
#define get_Pza_STc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*     stcpas; \
	    r2	= n0*pVp_stcpas;  \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcpa;  \
	    Ecorr += r3;
#define get_Pza_Sc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*     scpas; \
	    r2	= n0*pVp_scpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scpa; \
	    Ecorr += r3;
#define get_Pza_Pxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*n*      pps  \
		- l*n*      ppp; \
	    r2	= l0*n0*pVp_pps  \
		- l0*n0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pza_Pyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*n*      pps  \
		- m*n*      ppp; \
	    r2	= m0*n0*pVp_pps  \
		- m0*n0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pza_Pzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n2*pps  \
		+ (1-n2)*ppp; \
	    r2	= n02*pVp_pps  \
		+ (1-n02)*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pza_Dxyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*padcs  \
		- 2*l*m*n*padcp; \
	    r2	= sqr3*l0*m0*n0*pVp_padcs  \
		- 2*l0*m0*n0*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pza_Dyzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n2*padcs  \
		+ m*(1-2*n2)*padcp; \
	    r2	= sqr3*m0*n02*pVp_padcs  \
		+ m0*(1-2*n02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pza_Dzxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*n2*padcs  \
		+ l*(1-2*n2)*padcp; \
	    r2	= sqr3*l0*n02*pVp_padcs  \
		+ l0*(1-2*n02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pza_Dx2my2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*n*(l2-m2)*padcs  \
		- n*(l2-m2)*padcp; \
	    r2	= 0.5*sqr3*n0*(l02-m02)*pVp_padcs  \
		- n0*(l02-m02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Pza_Dz2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*(n2-0.5*(l2+m2))*padcs  \
		+ sqr3*n*(l2+m2)*padcp; \
	    r2	= n0*(n02-0.5*(l02+m02))*pVp_padcs  \
		+ sqr3*n0*(l02+m02)*pVp_padcp;  \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc;  \
	    Ecorr += r3;
/* Dxya_ */
#define get_Dxya_STc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*stcdas; \
	    r2	= sqr3*l0*m0*pVp_stcdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcda; \
	    Ecorr += r3;
#define get_Dxya_Sc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*scdas; \
	    r2	= sqr3*l0*m0*pVp_scdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scda; \
	    Ecorr += r3;
#define get_Dxya_Pxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l2*m*pcdas  \
		+ m*(1-2*l2)*pcdap; \
	    r2	= sqr3*l02*m0*pVp_pcdas  \
		+ m0*(1-2*l02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dxya_Pyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m2*pcdas  \
		+ l*(1-2*m2)*pcdap; \
	    r2	= sqr3*l0*m02*pVp_pcdas  \
		+ l0*(1-2*m02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dxya_Pzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*pcdas  \
		- 2*l*m*n*pcdap; \
	    r2	= sqr3*l0*m0*n0*pVp_pcdas  \
		- 2*l0*m0*n0*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dxya_Dxyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*l2*m2*dds  \
		+ (l2+m2-4*l2*m2)*ddp  \
		+ (n2+l2*m2)*ddd; \
	    r2	= 3*l02*m02*pVp_dds  \
		+ (l02+m02-4*l02*m02)*pVp_ddp  \
		+ (n02+l02*m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dxya_Dyzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*l*m2*n*dds  \
		+ l*n*(1-4*m2)*ddp  \
		+ l*n*(m2-1)*ddd; \
	    r2	= 3*l0*m02*n0*pVp_dds  \
		+ l0*n0*(1-4*m02)*pVp_ddp  \
		+ l0*n0*(m02-1)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dxya_Dzxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*l2*m*n*dds  \
		+ m*n*(1-4*l2)*ddp  \
		+ m*n*(l2-1)*ddd; \
	    r2	= 3*l02*m0*n0*pVp_dds  \
		+ m0*n0*(1-4*l02)*pVp_ddp  \
		+ m0*n0*(l02-1)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dxya_Dx2my2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*l*m*(l2-m2)*dds  \
		+ 2*l*m*(m2-l2)*ddp  \
		+ 0.5*l*m*(l2-m2)*ddd; \
	    r2	= 1.5*l0*m0*(l02-m02)*pVp_dds  \
		+ 2*l0*m0*(m02-l02)*pVp_ddp  \
		+ 0.5*l0*m0*(l02-m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dxya_Dz2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*(n2-0.5*(l2+m2))*dds  \
		- sqr3*2*l*m*n2*ddp  \
		+ 0.5*sqr3*l*m*(1+n2)*ddd; \
	    r2	= sqr3*l0*m0*(n02-0.5*(l02+m02))*pVp_dds  \
		- sqr3*2*l0*m0*n02*pVp_ddp  \
		+ 0.5*sqr3*l0*m0*(1+n02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
/* Dyza_ */
#define get_Dyza_STc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*stcdas; \
	    r2	= sqr3*m0*n0*pVp_stcdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcda; \
	    Ecorr += r3;
#define get_Dyza_Sc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*scdas; \
	    r2	= sqr3*m0*n0*pVp_scdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scda; \
	    Ecorr += r3;
#define get_Dyza_Pxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*pcdas  \
		- 2*l*m*n*pcdap; \
	    r2	= sqr3*l0*m0*n0*pVp_pcdas  \
		- 2*l0*m0*n0*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dyza_Pyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*n*m2*pcdas  \
		+ n*(1-2*m2)*pcdap; \
	    r2	= sqr3*n0*m02*pVp_pcdas  \
		+ n0*(1-2*m02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dyza_Pzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n2*pcdas  \
		+ m*(1-2*n2)*pcdap; \
	    r2	= sqr3*m0*n02*pVp_pcdas  \
		+ m0*(1-2*n02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dyza_Dxyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*l*m2*n*dds  \
		+ l*n*(1-4*m2)*ddp  \
		+ l*n*(m2-1)*ddd; \
	    r2	= 3*l0*m02*n0*    pVp_dds  \
		+ l0*n0*(1-4*m02)*pVp_ddp  \
		+ l0*n0*(m02-1)*  pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dyza_Dyzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*m2*n2*dds  \
		+ (m2+n2-4*m2*n2)*ddp  \
		+ (l2+m2*n2)*ddd; \
	    r2	= 3*m02*n02*          pVp_dds  \
		+ (m02+n02-4*m02*n02)*pVp_ddp  \
		+ (l02+m02*n02)*      pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dyza_Dzxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*m*n2*l*dds  \
		+ m*l*(1-4*n2)*ddp  \
		+ m*l*(n2-1)*ddd; \
	    r2	= 3*m0*n02*l0*    pVp_dds  \
		+ m0*l0*(1-4*n02)*pVp_ddp  \
		+ m0*l0*(n02-1)*  pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dyza_Dx2my2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*m*n*(l2-m2)*dds  \
		- m*n*(1+2*(l2-m2))*ddp  \
		+ m*n*(1+0.5*(l2-m2))*ddd; \
	    r2	= 1.5*m0*n0*(l02-m02)*pVp_dds  \
		- m0*n0*(1+2*(l02-m02))*pVp_ddp  \
		+ m0*n0*(1+0.5*(l02-m02))*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dyza_Dz2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*m*n*(l2+m2-n2)*ddp  \
		- 0.5*sqr3*m*n*(l2+m2)*ddd; \
	    r2	= sqr3*m0*n0*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*m0*n0*(l02+m02-n02)*pVp_ddp  \
		- 0.5*sqr3*m0*n0*(l02+m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
/* Dzxa_ */
#define get_Dzxa_STc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*n*l*      stcdas; \
	    r2	= sqr3*n0*l0*pVp_stcdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcda; \
	    Ecorr += r3;
#define get_Dzxa_Sc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*n*l*      scdas; \
	    r2	= sqr3*n0*l0*pVp_scdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scda; \
	    Ecorr += r3;
#define get_Dzxa_Pxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l2*n*pcdas  \
		+ n*(1-2*l2)*pcdap; \
	    r2	= sqr3*l02*n0*pVp_pcdas  \
		+ n0*(1-2*l02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dzxa_Pyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*pcdas  \
		- 2*l*m*n*pcdap; \
	    r2	= sqr3*l0*m0*n0*pVp_pcdas  \
		- 2*l0*m0*n0*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dzxa_Pzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*n2*pcdas  \
		+ l*(1-2*n2)*pcdap; \
	    r2	= sqr3*l0*n02*pVp_pcdas  \
		+ l0*(1-2*n02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dzxa_Dxyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*l2*m*n*dds  \
		+ m*n*(1-4*l2)*ddp  \
		+ m*n*(l2-1)*ddd; \
	    r2	= 3*l02*m0*n0*pVp_dds  \
		+ m0*n0*(1-4*l02)*pVp_ddp  \
		+ m0*n0*(l02-1)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dzxa_Dyzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*n2*l*m*dds  \
		+ l*m*(1-4*n2)*ddp  \
		+ l*m*(n2-1)*ddd; \
	    r2	= 3*n02*l0*m0*pVp_dds  \
		+ l0*m0*(1-4*n02)*pVp_ddp  \
		+ l0*m0*(n02-1)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dzxa_Dzxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*n2*l2*dds  \
		+ (n2+l2-4*n2*l2)*ddp  \
		+ (m2+n2*l2)*ddd; \
	    r2	= 3*n02*l02*pVp_dds  \
		+ (n02+l02-4*n02*l02)*pVp_ddp  \
		+ (m02+n02*l02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dzxa_Dx2my2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*n*l*(l2-m2)*dds  \
		+ n*l*(1-2*(l2-m2))*ddp  \
		- n*l*(1-0.5*(l2-m2))*ddd; \
	    r2	= 1.5*n0*l0*(l02-m02)*pVp_dds  \
		+ n0*l0*(1-2*(l02-m02))*pVp_ddp  \
		- n0*l0*(1-0.5*(l02-m02))*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dzxa_Dz2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*n*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*l*n*(l2+m2-n2)*ddp  \
		- 0.5*sqr3*l*n*(l2+m2)*ddd; \
	    r2	= sqr3*l0*n0*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*l0*n0*(l02+m02-n02)*pVp_ddp  \
		- 0.5*sqr3*l0*n0*(l02+m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
/* Dx2my2a_ */
#define get_Dx2my2a_STc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2-m2)*      stcdas; \
	    r2	= 0.5*sqr3*(l02-m02)*pVp_stcdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcda; \
	    Ecorr += r3;
#define get_Dx2my2a_Sc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2-m2)*      scdas; \
	    r2	= 0.5*sqr3*(l02-m02)*pVp_scdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scda; \
	    Ecorr += r3;
#define get_Dx2my2a_Pxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*l*(l2-m2)*       pcdas  \
		+ l*(1-l2+m2)*              pcdap; \
	    r2	= 0.5*sqr3*l0*(l02-m02)*pVp_pcdas  \
		+ l0*(1-l02+m02)*       pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dx2my2a_Pyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*m*(l2-m2)*       pcdas  \
		- m*(1+l2-m2)*              pcdap; \
	    r2	= 0.5*sqr3*m0*(l02-m02)*pVp_pcdas  \
		- m0*(1+l02-m02)*       pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dx2my2a_Pzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*n*(l2-m2)*       pcdas  \
		- n*(l2-m2)*                pcdap; \
	    r2	= 0.5*sqr3*n0*(l02-m02)*pVp_pcdas  \
		- n0*(l02-m02)*         pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dx2my2a_Dxyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*l*m*(l2-m2)*dds  \
		+ 2*l*m*(m2-l2)*ddp  \
		+ 0.5*l*m*(l2-m2)*ddd; \
	    r2	= 1.5*l0*m0*(l02-m02)*pVp_dds  \
		+ 2*l0*m0*(m02-l02)*pVp_ddp  \
		+ 0.5*l0*m0*(l02-m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dx2my2a_Dyzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*m*n*(l2-m2)*dds  \
		- m*n*(1+2*(l2-m2))*ddp  \
		+ m*n*(1+0.5*(l2-m2))*ddd; \
	    r2	= 1.5*m0*n0*(l02-m02)*pVp_dds  \
		- m0*n0*(1+2*(l02-m02))*pVp_ddp  \
		+ m0*n0*(1+0.5*(l02-m02))*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dx2my2a_Dzxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*n*l*(l2-m2)*dds  \
		+ n*l*(1-2*(l2-m2))*ddp  \
		- n*l*(1-0.5*(l2-m2))*ddd; \
	    r2	= 1.5*n0*l0*(l02-m02)*pVp_dds  \
		+ n0*l0*(1-2*(l02-m02))*pVp_ddp  \
		- n0*l0*(1-0.5*(l02-m02))*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dx2my2a_Dx2my2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.75*(l2-m2)*(l2-m2)*dds  \
		+ (l2+m2-(l2-m2)*(l2-m2))*ddp  \
		+ (n2+0.25*(l2-m2)*(l2-m2))*ddd; \
	    r2	= 0.75*(l02-m02)*(l02-m02)*pVp_dds  \
		+ (l02+m02-(l02-m02)*(l02-m02))*pVp_ddp  \
		+ (n02+0.25*(l02-m02)*(l02-m02))*pVp_ddd;  \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dx2my2a_Dz2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2-m2)*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*n2*(m2-l2)*ddp  \
		+ 0.25*sqr3*(1+n2)*(l2-m2)*ddd; \
	    r2	= 0.5*sqr3*(l02-m02)*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*n02*(m02-l02)*pVp_ddp  \
		+ 0.25*sqr3*(1+n02)*(l02-m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
/* Dz2a_ */
#define get_Dz2a_STc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= (n2-0.5*(l2+m2))*stcdas; \
	    r2	= (n02-0.5*(l02+m02))*pVp_stcdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcda; \
	    Ecorr += r3;
#define get_Dz2a_Sc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= (n2-0.5*(l2+m2))*scdas; \
	    r2	= (n02-0.5*(l02+m02))*pVp_scdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scda; \
	    Ecorr += r3;
#define get_Dz2a_Pxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*(n2-0.5*(l2+m2))*pcdas  \
		- sqr3*l*n2*pcdap; \
	    r2	= l0*(n02-0.5*(l02+m02))*pVp_pcdas  \
		- sqr3*l0*n02*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dz2a_Pyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*(n2-0.5*(l2+m2))*pcdas  \
		- sqr3*m*n2*pcdap; \
	    r2	= m0*(n02-0.5*(l02+m02))*pVp_pcdas  \
		- sqr3*m0*n02*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dz2a_Pzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*(n2-0.5*(l2+m2))*pcdas  \
		+ sqr3*n*(l2+m2)*pcdap; \
	    r2	= n0*(n02-0.5*(l02+m02))*pVp_pcdas  \
		+ sqr3*n0*(l02+m02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Dz2a_Dxyc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*(n2-0.5*(l2+m2))*dds  \
		- sqr3*2*l*m*n2*ddp  \
		+ 0.5*sqr3*l*m*(1+n2)*ddd; \
	    r2	= sqr3*l0*m0*(n02-0.5*(l02+m02))*pVp_dds  \
		- sqr3*2*l0*m0*n02*pVp_ddp  \
		+ 0.5*sqr3*l0*m0*(1+n02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dz2a_Dyzc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*m*n*(l2+m2-n2)*ddp  \
		- 0.5*sqr3*m*n*(l2+m2)*ddd; \
	    r2	= sqr3*m0*n0*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*m0*n0*(l02+m02-n02)*pVp_ddp  \
		- 0.5*sqr3*m0*n0*(l02+m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dz2a_Dzxc(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*n*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*l*n*(l2+m2-n2)*ddp  \
		- 0.5*sqr3*l*n*(l2+m2)*ddd; \
	    r2	= sqr3*l0*n0*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*l0*n0*(l02+m02-n02)*pVp_ddp  \
		- 0.5*sqr3*l0*n0*(l02+m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dz2a_Dx2my2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2-m2)*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*n2*(m2-l2)*ddp  \
		+ 0.25*sqr3*(1+n2)*(l2-m2)*ddd; \
	    r2	= 0.5*sqr3*(l02-m02)*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*n02*(m02-l02)*pVp_ddp  \
		+ 0.25*sqr3*(1+n02)*(l02-m02)*pVp_ddd;  \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dz2a_Dz2c(Vmat, Ecorr) \
	    r1 = Vmat  \
		= (n2-0.5*(l2+m2))*(n2-0.5*(l2+m2))*          dds  \
		+ 3*n2*(l2+m2)*                               ddp  \
		+ 0.75*(l2+m2)*(l2+m2)*                       ddd; \
	    r2	= (n02-0.5*(l02+m02))*(n02-0.5*(l02+m02))*pVp_dds  \
		+ 3*n02*(l02+m02)*                        pVp_ddp  \
		+ 0.75*(l02+m02)*(l02+m02)*               pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
/* STc_ */
#define get_STc_STa(Vmat, Ecorr) \
	    r1 = Vmat  \
		=     ststs; \
	    r2	= pVp_ststs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stst; \
	    Ecorr += r3;
#define get_STc_Sa(Vmat, Ecorr) \
	    r1 = Vmat  \
		=     sasts; \
	    r2	= pVp_sasts; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sast; \
	    Ecorr += r3;
#define get_STc_Pxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*     stcpas; \
	    r2	= l0*pVp_stcpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcpa; \
	    Ecorr += r3;
#define get_STc_Pya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*     stcpas; \
	    r2	= m0*pVp_stcpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcpa; \
	    Ecorr += r3;
#define get_STc_Pza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*     stcpas; \
	    r2	= n0*pVp_stcpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcpa; \
	    Ecorr += r3;
#define get_STc_Dxya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*      stcdas; \
	    r2	= sqr3*l0*m0*pVp_stcdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcda; \
	    Ecorr += r3;
#define get_STc_Dyza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*      stcdas; \
	    r2	= sqr3*m0*n0*pVp_stcdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcda; \
	    Ecorr += r3;
#define get_STc_Dzxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*n*l*      stcdas; \
	    r2	= sqr3*n0*l0*pVp_stcdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcda; \
	    Ecorr += r3;
#define get_STc_Dx2my2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2-m2)*      stcdas; \
	    r2	= 0.5*sqr3*(l02-m02)*pVp_stcdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcda; \
	    Ecorr += r3;
#define get_STc_Dz2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= (n2-0.5*(l2+m2))*       stcdas; \
	    r2	= (n02-0.5*(l02+m02))*pVp_stcdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stcda; \
	    Ecorr += r3;
/* Sc_ */
#define get_Sc_STa(Vmat, Ecorr) \
	    r1 = Vmat  \
		=     stass; \
	    r2	= pVp_stass; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stas; \
	    Ecorr += r3;
#define get_Sc_Sa(Vmat, Ecorr) \
	    r1 = Vmat  \
		=     sss; \
	    r2	= pVp_sss; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_ss; \
	    Ecorr += r3;
#define get_Sc_Pxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*     scpas; \
	    r2	= l0*pVp_scpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scpa; \
	    Ecorr += r3;
#define get_Sc_Pya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*     scpas; \
	    r2	= m0*pVp_scpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scpa; \
	    Ecorr += r3;
#define get_Sc_Pza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*     scpas; \
	    r2	= n0*pVp_scpas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scpa; \
	    Ecorr += r3;
#define get_Sc_Dxya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*      scdas; \
	    r2	= sqr3*l0*m0*pVp_scdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scda; \
	    Ecorr += r3;
#define get_Sc_Dyza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*      scdas; \
	    r2	= sqr3*m0*n0*pVp_scdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scda; \
	    Ecorr += r3;
#define get_Sc_Dzxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*n*      scdas; \
	    r2	= sqr3*l0*n0*pVp_scdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scda; \
	    Ecorr += r3;
#define get_Sc_Dx2my2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2-m2)*      scdas; \
	    r2	= 0.5*sqr3*(l02-m02)*pVp_scdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scda; \
	    Ecorr += r3;
#define get_Sc_Dz2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= (n2-0.5*(l2+m2))*       scdas; \
	    r2	= (n02-0.5*(l02+m02))*pVp_scdas; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_scda; \
	    Ecorr += r3;
/* Pxc_ */
#define get_Pxc_STa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*     stapcs; \
	    r2	= l0*pVp_stapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stapc; \
	    Ecorr += r3;
#define get_Pxc_Sa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*     sapcs; \
	    r2	= l0*pVp_sapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sapc; \
	    Ecorr += r3;
#define get_Pxc_Pxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l2*         pps  \
		+ (1-l2)*     ppp; \
	    r2	= l02*    pVp_pps  \
		+ (1-l02)*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pxc_Pya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*m*      pps  \
		- l*m*      ppp; \
	    r2	= l0*m0*pVp_pps  \
		- l0*m0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pxc_Pza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*n*      pps  \
		- l*n*      ppp; \
	    r2	= l0*n0*pVp_pps  \
		- l0*n0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pxc_Dxya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l2*m*pcdas  \
		+ m*(1-2*l2)*pcdap; \
	    r2	= sqr3*l02*m0*pVp_pcdas  \
		+ m0*(1-2*l02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pxc_Dyza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*pcdas  \
		- 2*l*m*n*pcdap; \
	    r2	= sqr3*l0*m0*n0*pVp_pcdas  \
		- 2*l0*m0*n0*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pxc_Dzxa(Vmat, Ecorr) \
	    r1 = Vmat \
		 = sqr3*l2*n*pcdas  \
		+ n*(1-2*l2)*pcdap; \
	    r2 = sqr3*l02*n0*pVp_pcdas  \
		+ n0*(1-2*l02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pxc_Dx2my2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*l*(l2-m2)*pcdas  \
		+ l*(1-l2+m2)*pcdap; \
	    r2	= 0.5*sqr3*l0*(l02-m02)*pVp_pcdas  \
		+ l0*(1-l02+m02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pxc_Dz2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*(n2-0.5*(l2+m2))*pcdas  \
		- sqr3*l*n2*pcdap; \
	    r2	= l0*(n02-0.5*(l02+m02))*pVp_pcdas  \
		- sqr3*l0*n02*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
/* Pyc_ */
#define get_Pyc_STa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*stapcs; \
	    r2	= m0*pVp_stapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stapc; \
	    Ecorr += r3;
#define get_Pyc_Sa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*sapcs; \
	    r2	= m0*pVp_sapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sapc; \
	    Ecorr += r3;
#define get_Pyc_Pxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*m*pps  \
		- l*m*ppp; \
	    r2	= l0*m0*pVp_pps  \
		- l0*m0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pyc_Pya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m2*pps  \
		+ (1-m2)*ppp; \
	    r2	= m02*pVp_pps  \
		+ (1-m02)*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pyc_Pza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*n*pps  \
		- m*n*ppp; \
	    r2	= m0*n0*pVp_pps  \
		- m0*n0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pyc_Dxya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m2*pcdas  \
		+ l*(1-2*m2)*pcdap; \
	    r2	= sqr3*l0*m02*pVp_pcdas  \
		+ l0*(1-2*m02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pyc_Dyza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*n*m2*pcdas  \
		+ n*(1-2*m2)*pcdap; \
	    r2	= sqr3*n0*m02*pVp_pcdas  \
		+ n0*(1-2*m02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pyc_Dzxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*pcdas  \
		- 2*l*m*n*pcdap; \
	    r2	= sqr3*l0*m0*n0*pVp_pcdas  \
		- 2*l0*m0*n0*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pyc_Dx2my2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*m*(l2-m2)*pcdas  \
		- m*(1+l2-m2)*pcdap; \
	    r2	= 0.5*sqr3*m0*(l02-m02)*pVp_pcdas  \
		- m0*(1+l02-m02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pyc_Dz2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*(n2-0.5*(l2+m2))*pcdas  \
		- sqr3*m*n2*pcdap; \
	    r2	= m0*(n02-0.5*(l02+m02))*pVp_pcdas  \
		- sqr3*m0*n02*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
/* Pzc_ */
#define get_Pzc_STa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*stapcs; \
	    r2	= n0*pVp_stapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stapc; \
	    Ecorr += r3;
#define get_Pzc_Sa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*sapcs; \
	    r2	= n0*pVp_sapcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sapc; \
	    Ecorr += r3;
#define get_Pzc_Pxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*n*pps  \
		- l*n*ppp; \
	    r2	= l0*n0*pVp_pps  \
		- l0*n0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pzc_Pya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*n*pps  \
		- m*n*ppp; \
	    r2	= m0*n0*pVp_pps  \
		- m0*n0*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pzc_Pza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n2*pps  \
		+ (1-n2)*ppp; \
	    r2	= n02*pVp_pps  \
		+ (1-n02)*pVp_ppp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pp; \
	    Ecorr += r3;
#define get_Pzc_Dxya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*pcdas  \
		- 2*l*m*n*pcdap; \
	    r2	= sqr3*l0*m0*n0*pVp_pcdas  \
		- 2*l0*m0*n0*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pzc_Dyza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n2*pcdas  \
		+ m*(1-2*n2)*pcdap; \
	    r2	= sqr3*m0*n02*pVp_pcdas  \
		+ m0*(1-2*n02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pzc_Dzxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*n2*pcdas  \
		+ l*(1-2*n2)*pcdap; \
	    r2	= sqr3*l0*n02*pVp_pcdas  \
		+ l0*(1-2*n02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pzc_Dx2my2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*n*(l2-m2)*pcdas  \
		- n*(l2-m2)*pcdap; \
	    r2	= 0.5*sqr3*n0*(l02-m02)*pVp_pcdas  \
		- n0*(l02-m02)*pVp_pcdap; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda; \
	    Ecorr += r3;
#define get_Pzc_Dz2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*(n2-0.5*(l2+m2))*pcdas  \
		+ sqr3*n*(l2+m2)*pcdap; \
	    r2	= n0*(n02-0.5*(l02+m02))*pVp_pcdas  \
		+ sqr3*n0*(l02+m02)*pVp_pcdap;  \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_pcda;  \
	    Ecorr += r3;
/* Dxyc_ */
#define get_Dxyc_STa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*stadcs; \
	    r2	= sqr3*l0*m0*pVp_stadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stadc; \
	    Ecorr += r3;
#define get_Dxyc_Sa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*sadcs; \
	    r2	= sqr3*l0*m0*pVp_sadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sadc; \
	    Ecorr += r3;
#define get_Dxyc_Pxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l2*m*padcs  \
		+ m*(1-2*l2)*padcp; \
	    r2	= sqr3*l02*m0*pVp_padcs  \
		+ m0*(1-2*l02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dxyc_Pya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m2*padcs  \
		+ l*(1-2*m2)*padcp; \
	    r2	= sqr3*l0*m02*pVp_padcs  \
		+ l0*(1-2*m02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dxyc_Pza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*padcs  \
		- 2*l*m*n*padcp; \
	    r2	= sqr3*l0*m0*n0*pVp_padcs  \
		- 2*l0*m0*n0*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dxyc_Dxya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*l2*m2*dds  \
		+ (l2+m2-4*l2*m2)*ddp  \
		+ (n2+l2*m2)*ddd; \
	    r2	= 3*l02*m02*pVp_dds  \
		+ (l02+m02-4*l02*m02)*pVp_ddp  \
		+ (n02+l02*m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dxyc_Dyza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*l*m2*n*dds  \
		+ l*n*(1-4*m2)*ddp  \
		+ l*n*(m2-1)*ddd; \
	    r2	= 3*l0*m02*n0*pVp_dds  \
		+ l0*n0*(1-4*m02)*pVp_ddp  \
		+ l0*n0*(m02-1)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dxyc_Dzxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*l2*m*n*dds  \
		+ m*n*(1-4*l2)*ddp  \
		+ m*n*(l2-1)*ddd; \
	    r2	= 3*l02*m0*n0*pVp_dds  \
		+ m0*n0*(1-4*l02)*pVp_ddp  \
		+ m0*n0*(l02-1)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dxyc_Dx2my2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*l*m*(l2-m2)*dds  \
		+ 2*l*m*(m2-l2)*ddp  \
		+ 0.5*l*m*(l2-m2)*ddd; \
	    r2	= 1.5*l0*m0*(l02-m02)*pVp_dds  \
		+ 2*l0*m0*(m02-l02)*pVp_ddp  \
		+ 0.5*l0*m0*(l02-m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dxyc_Dz2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*(n2-0.5*(l2+m2))*dds  \
		- sqr3*2*l*m*n2*ddp  \
		+ 0.5*sqr3*l*m*(1+n2)*ddd; \
	    r2	= sqr3*l0*m0*(n02-0.5*(l02+m02))*pVp_dds  \
		- sqr3*2*l0*m0*n02*pVp_ddp  \
		+ 0.5*sqr3*l0*m0*(1+n02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
/* Dyzc_ */
#define get_Dyzc_STa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*stadcs; \
	    r2	= sqr3*m0*n0*pVp_stadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stadc; \
	    Ecorr += r3;
#define get_Dyzc_Sa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*sadcs; \
	    r2	= sqr3*m0*n0*pVp_sadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sadc; \
	    Ecorr += r3;
#define get_Dyzc_Pxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*padcs  \
		- 2*l*m*n*padcp; \
	    r2	= sqr3*l0*m0*n0*pVp_padcs  \
		- 2*l0*m0*n0*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dyzc_Pya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*n*m2*padcs  \
		+ n*(1-2*m2)*padcp; \
	    r2	= sqr3*n0*m02*pVp_padcs  \
		+ n0*(1-2*m02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dyzc_Pza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n2*padcs  \
		+ m*(1-2*n2)*padcp; \
	    r2	= sqr3*m0*n02*pVp_padcs  \
		+ m0*(1-2*n02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dyzc_Dxya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*l*m2*n*dds  \
		+ l*n*(1-4*m2)*ddp  \
		+ l*n*(m2-1)*ddd; \
	    r2	= 3*l0*m02*n0*pVp_dds  \
		+ l0*n0*(1-4*m02)*pVp_ddp  \
		+ l0*n0*(m02-1)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dyzc_Dyza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*m2*n2*dds  \
		+ (m2+n2-4*m2*n2)*ddp  \
		+ (l2+m2*n2)*ddd; \
	    r2	= 3*m02*n02*pVp_dds  \
		+ (m02+n02-4*m02*n02)*pVp_ddp  \
		+ (l02+m02*n02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dyzc_Dzxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*m*n2*l*dds  \
		+ m*l*(1-4*n2)*ddp  \
		+ m*l*(n2-1)*ddd; \
	    r2	= 3*m0*n02*l0*pVp_dds  \
		+ m0*l0*(1-4*n02)*pVp_ddp  \
		+ m0*l0*(n02-1)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dyzc_Dx2my2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*m*n*(l2-m2)*dds  \
		- m*n*(1+2*(l2-m2))*ddp  \
		+ m*n*(1+0.5*(l2-m2))*ddd; \
	    r2	= 1.5*m0*n0*(l02-m02)*pVp_dds  \
		- m0*n0*(1+2*(l02-m02))*pVp_ddp  \
		+ m0*n0*(1+0.5*(l02-m02))*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dyzc_Dz2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*m*n*(l2+m2-n2)*ddp  \
		- 0.5*sqr3*m*n*(l2+m2)*ddd; \
	    r2	= sqr3*m0*n0*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*m0*n0*(l02+m02-n02)*pVp_ddp  \
		- 0.5*sqr3*m0*n0*(l02+m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
/* Dzxc_ */
#define get_Dzxc_STa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*n*l*stadcs; \
	    r2	= sqr3*n0*l0*pVp_stadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stadc; \
	    Ecorr += r3;
#define get_Dzxc_Sa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*n*l*sadcs; \
	    r2	= sqr3*n0*l0*pVp_sadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sadc; \
	    Ecorr += r3;
#define get_Dzxc_Pxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l2*n*padcs  \
		+ n*(1-2*l2)*padcp; \
	    r2	= sqr3*l02*n0*pVp_padcs  \
		+ n0*(1-2*l02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dzxc_Pya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*n*padcs  \
		- 2*l*m*n*padcp; \
	    r2	= sqr3*l0*m0*n0*pVp_padcs  \
		- 2*l0*m0*n0*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dzxc_Pza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*n2*padcs  \
		+ l*(1-2*n2)*padcp; \
	    r2	= sqr3*l0*n02*pVp_padcs  \
		+ l0*(1-2*n02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dzxc_Dxya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*l2*m*n*dds  \
		+ m*n*(1-4*l2)*ddp  \
		+ m*n*(l2-1)*ddd; \
	    r2	= 3*l02*m0*n0*pVp_dds  \
		+ m0*n0*(1-4*l02)*pVp_ddp  \
		+ m0*n0*(l02-1)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dzxc_Dyza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*n2*l*m*dds  \
		+ l*m*(1-4*n2)*ddp  \
		+ l*m*(n2-1)*ddd; \
	    r2	= 3*n02*l0*m0*pVp_dds  \
		+ l0*m0*(1-4*n02)*pVp_ddp  \
		+ l0*m0*(n02-1)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dzxc_Dzxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 3*n2*l2*dds  \
		+ (n2+l2-4*n2*l2)*ddp  \
		+ (m2+n2*l2)*ddd; \
	    r2	= 3*n02*l02*pVp_dds  \
		+ (n02+l02-4*n02*l02)*pVp_ddp  \
		+ (m02+n02*l02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dzxc_Dx2my2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*n*l*(l2-m2)*dds  \
		+ n*l*(1-2*(l2-m2))*ddp  \
		- n*l*(1-0.5*(l2-m2))*ddd; \
	    r2	= 1.5*n0*l0*(l02-m02)*pVp_dds  \
		+ n0*l0*(1-2*(l02-m02))*pVp_ddp  \
		- n0*l0*(1-0.5*(l02-m02))*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dzxc_Dz2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*n*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*l*n*(l2+m2-n2)*ddp  \
		- 0.5*sqr3*l*n*(l2+m2)*ddd; \
	    r2	= sqr3*l0*n0*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*l0*n0*(l02+m02-n02)*pVp_ddp  \
		- 0.5*sqr3*l0*n0*(l02+m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
/* Dx2my2c_ */
#define get_Dx2my2c_STa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2-m2)*stadcs; \
	    r2	= 0.5*sqr3*(l02-m02)*pVp_stadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stadc; \
	    Ecorr += r3;
#define get_Dx2my2c_Sa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2-m2)*sadcs; \
	    r2	= 0.5*sqr3*(l02-m02)*pVp_sadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sadc; \
	    Ecorr += r3;
#define get_Dx2my2c_Pxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*l*(l2-m2)*padcs  \
		+ l*(1-l2+m2)*padcp; \
	    r2	= 0.5*sqr3*l0*(l02-m02)*pVp_padcs  \
		+ l0*(1-l02+m02)*       pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dx2my2c_Pya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*m*(l2-m2)*padcs  \
		- m*(1+l2-m2)*padcp; \
	    r2	= 0.5*sqr3*m0*(l02-m02)*pVp_padcs  \
		- m0*(1+l02-m02)*       pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dx2my2c_Pza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*n*(l2-m2)*padcs  \
		- n*(l2-m2)*padcp; \
	    r2	= 0.5*sqr3*n0*(l02-m02)*pVp_padcs  \
		- n0*(l02-m02)*         pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dx2my2c_Dxya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*l*m*(l2-m2)*dds  \
		+ 2*l*m*(m2-l2)*ddp  \
		+ 0.5*l*m*(l2-m2)*ddd; \
	    r2	= 1.5*l0*m0*(l02-m02)*pVp_dds  \
		+ 2*l0*m0*(m02-l02)*pVp_ddp  \
		+ 0.5*l0*m0*(l02-m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dx2my2c_Dyza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*m*n*(l2-m2)*dds  \
		- m*n*(1+2*(l2-m2))*ddp  \
		+ m*n*(1+0.5*(l2-m2))*ddd; \
	    r2	= 1.5*m0*n0*(l02-m02)*pVp_dds  \
		- m0*n0*(1+2*(l02-m02))*pVp_ddp  \
		+ m0*n0*(1+0.5*(l02-m02))*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dx2my2c_Dzxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 1.5*n*l*(l2-m2)*dds  \
		+ n*l*(1-2*(l2-m2))*ddp  \
		- n*l*(1-0.5*(l2-m2))*ddd; \
	    r2	= 1.5*n0*l0*(l02-m02)*pVp_dds  \
		+ n0*l0*(1-2*(l02-m02))*pVp_ddp  \
		- n0*l0*(1-0.5*(l02-m02))*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dx2my2c_Dx2my2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.75*(l2-m2)*(l2-m2)*dds  \
		+ (l2+m2-(l2-m2)*(l2-m2))*ddp  \
		+ (n2+0.25*(l2-m2)*(l2-m2))*ddd; \
	    r2	= 0.75*(l02-m02)*(l02-m02)*pVp_dds  \
		+ (l02+m02-(l02-m02)*(l02-m02))*pVp_ddp  \
		+ (n02+0.25*(l02-m02)*(l02-m02))*pVp_ddd;  \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd;  \
	    Ecorr += r3;
#define get_Dx2my2c_Dz2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2-m2)*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*n2*(m2-l2)*ddp  \
		+ 0.25*sqr3*(1+n2)*(l2-m2)*ddd; \
	    r2	= 0.5*sqr3*(l02-m02)*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*n02*(m02-l02)*pVp_ddp  \
		+ 0.25*sqr3*(1+n02)*(l02-m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
/* Dz2c_ */
#define get_Dz2c_STa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= (n2-0.5*(l2+m2))*stadcs; \
	    r2	= (n02-0.5*(l02+m02))*pVp_stadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_stadc; \
	    Ecorr += r3;
#define get_Dz2c_Sa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= (n2-0.5*(l2+m2))*sadcs; \
	    r2	= (n02-0.5*(l02+m02))*pVp_sadcs; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_sadc; \
	    Ecorr += r3;
#define get_Dz2c_Pxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= l*(n2-0.5*(l2+m2))*padcs  \
		- sqr3*l*n2*padcp; \
	    r2	= l0*(n02-0.5*(l02+m02))*pVp_padcs  \
		- sqr3*l0*n02*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dz2c_Pya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= m*(n2-0.5*(l2+m2))*padcs  \
		- sqr3*m*n2*padcp; \
	    r2	= m0*(n02-0.5*(l02+m02))*pVp_padcs  \
		- sqr3*m0*n02*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dz2c_Pza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= n*(n2-0.5*(l2+m2))*padcs  \
		+ sqr3*n*(l2+m2)*padcp; \
	    r2	= n0*(n02-0.5*(l02+m02))*pVp_padcs  \
		+ sqr3*n0*(l02+m02)*pVp_padcp; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_padc; \
	    Ecorr += r3;
#define get_Dz2c_Dxya(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*m*(n2-0.5*(l2+m2))*dds  \
		- sqr3*2*l*m*n2*ddp  \
		+ 0.5*sqr3*l*m*(1+n2)*ddd; \
	    r2	= sqr3*l0*m0*(n02-0.5*(l02+m02))*pVp_dds  \
		- sqr3*2*l0*m0*n02*pVp_ddp  \
		+ 0.5*sqr3*l0*m0*(1+n02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dz2c_Dyza(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*m*n*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*m*n*(l2+m2-n2)*ddp  \
		- 0.5*sqr3*m*n*(l2+m2)*ddd; \
	    r2	= sqr3*m0*n0*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*m0*n0*(l02+m02-n02)*pVp_ddp  \
		- 0.5*sqr3*m0*n0*(l02+m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dz2c_Dzxa(Vmat, Ecorr) \
	    r1 = Vmat  \
		= sqr3*l*n*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*l*n*(l2+m2-n2)*ddp  \
		- 0.5*sqr3*l*n*(l2+m2)*ddd; \
	    r2	= sqr3*l0*n0*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*l0*n0*(l02+m02-n02)*pVp_ddp  \
		- 0.5*sqr3*l0*n0*(l02+m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;
#define get_Dz2c_Dx2my2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= 0.5*sqr3*(l2-m2)*(n2-0.5*(l2+m2))*dds  \
		+ sqr3*n2*(m2-l2)*ddp  \
		+ 0.25*sqr3*(1+n2)*(l2-m2)*ddd; \
	    r2	= 0.5*sqr3*(l02-m02)*(n02-0.5*(l02+m02))*pVp_dds  \
		+ sqr3*n02*(m02-l02)*pVp_ddp  \
		+ 0.25*sqr3*(1+n02)*(l02-m02)*pVp_ddd;  \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		pVpE_dd;  \
	    Ecorr += r3;
#define get_Dz2c_Dz2a(Vmat, Ecorr) \
	    r1 = Vmat  \
		= (n2-0.5*(l2+m2))*(n2-0.5*(l2+m2))*dds  \
		+ 3*n2*(l2+m2)*ddp  \
		+ 0.75*(l2+m2)*(l2+m2)*ddd; \
	    r2	= (n02-0.5*(l02+m02))*(n02-0.5*(l02+m02))*pVp_dds  \
		+ 3*n02*(l02+m02)*pVp_ddp  \
		+ 0.75*(l02+m02)*(l02+m02)*pVp_ddd; \
	    r3  = (r2*r2-r1*r1)*cstrain/ \
		  pVpE_dd; \
	    Ecorr += r3;

