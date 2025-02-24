c /*****************************************************************************
c The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
c Copyright (C) 2002 California Institute of Technology (Caltech)
c
c This application is free software, which you can redistribute and/or modify
c under the terms of the GNU Lesser General Public License as published by the
c Free Software Foundation; either version 2.1 of the License, or (at your
c option) any later version.
c
c This library is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c Lesser General Public License for more details.
c
c You should have received a copy of the GNU Lesser General Public License
c along with this library; see the file COPYING. If not, write to the
c Free Software Foundation, Inc.,
c 59 Temple Place, Suite 330,
c Boston, MA  02111-1307  USA
c
c For additional information, please contact
c   Gerhard Klimeck (gekco@jpl.nasa.gov)
c   Fabiano Oyafuso (fabiano@jpl.nasa.gov)
c
c Written by:  Chris Bowen
c              Gerhard Klimeck
c              Fabiano Oyafuso
c              Seungwon Lee
c              Olga Lazarenkova
c              Hook Hua
c
c This product includes software developed by the Apache Software Foundation
c (http://www.apache.org/).
c
c *****************************************************************************
c $Header: /repo/nemo3d/src/base/blas_run3d_f.f,v 1.2 2003/10/08 16:16:22 hook Exp $
c *****************************************************************************/
c
	subroutine zgemvt_f( y, a, alpha, beta, x, n)

	integer n
        complex*16 alpha, beta
        complex*16 y(n), x(n)
	complex*16 a(n,n)

        call zgemv('T', n, n, alpha, a, n, x, 1, beta, y, 1)
 
	return 
  
	end


	subroutine zdotc_f( a, x, y, n )

	integer n
	complex*16 a
	complex*16 x(n), y(n)

	a = zdotc(n, x, 1, y, 1)
	
	return

	end
