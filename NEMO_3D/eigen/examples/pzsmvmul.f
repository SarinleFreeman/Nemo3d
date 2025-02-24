c *****************************************************************************
c The JPL Parallel Eigensolvers package.
c Copyright (C) 2002 California Institute of Technology (Caltech)
c
c This file is part of
c The Jet Propulsion Laboratory (JPL) Parallel Eigensolvers package.
c
c This library is free software, which you can redistribute and/or modify
c under the terms of the GNU Lesser General Public License as published by the
c Free Software Foundation; either version 2.1 of the License, or (at your
c option) any later version.
c
c This library is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty
c of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c See the GNU General Public License for more details.
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
c   E. Robert Tisdale (E.Robert.Tisdale@jpl.nasa.gov)
c
c Written by: Chris Bowen
c             Gerhard Klimeck
c             Fabiano Oyafuso
c             E. Robert Tisdale
c
c *****************************************************************************
c
c
c     parallel matrix-vector multiplication subroutine
c
c     The matrix used is the 2 dimensional discrete Laplacian
c     on the unit square with zero Dirichlet boundary condition.
c
c     Computes w <-- OP*v,
c     where OP is the nx*nx by nx*nx block tridiagonal matrix
c
c		| T -I		|
c		|-I  T -I	|
c	   OP = |   -I  T	|
c		|	 ...  -I|
c		|	    -I T|
c
c     The subroutine TV is called to compute y <-- T*x.

c
      subroutine av(comm, nloc, nx, mv_buf, v, w)
c
c	.. MPI Declarations ...
	include	'mpif.h'
	integer	comm, nprocs, myid, ierr, status(MPI_STATUS_SIZE)
	integer		nloc, nx, np, j, lo, next, prev
	complex*16	v(nloc), w(nloc), mv_buf(nx), one
	parameter	(one = (1.0, 0.0))
	external	zaxpy

	call MPI_COMM_RANK(comm, myid, ierr)
	call MPI_COMM_SIZE(comm, nprocs, ierr)
c
	np = nloc/nx
	call tv(nx, v(1), w(1))
	call zaxpy(nx, -one, v(nx+1), 1, w(1), 1)
c
	if (np .gt. 2) then
	  do 10 j = 2, np-1
	    lo = (j-1)*nx
	    call tv(nx, v(lo+1), w(lo+1))
	    call zaxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
	    call zaxpy(nx, -one, v(lo+nx+1), 1, w(lo+1), 1)
  10	    continue
	  end if
c
	if (np .gt. 1) then
	  lo = (np-1)*nx
	  call tv(nx, v(lo+1), w(lo+1))
	  call zaxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
	  end if
c
	next = myid + 1
	prev = myid - 1
	if (myid .lt. nprocs-1) then
	  call mpi_send(v((np-1)*nx+1), nx, MPI_DOUBLE_COMPLEX,
     &	      next, myid+1, comm, ierr)
	  end if
	if (myid .gt. 0) then
	  call mpi_recv(mv_buf, nx, MPI_DOUBLE_COMPLEX,
     &	      prev, myid, comm, status, ierr)
	  call zaxpy(nx, -one, mv_buf, 1, w(1), 1)
	  end if
c
	if (myid .gt. 0) then
	  call mpi_send(v(1), nx, MPI_DOUBLE_COMPLEX,
     &	      prev, myid-1, comm, ierr)
	  end if
	if (myid .lt. nprocs-1) then
	  call mpi_recv(mv_buf, nx, MPI_DOUBLE_COMPLEX,
     &	      next, myid, comm, status, ierr)
	  call zaxpy(nx, -one, mv_buf, 1, w(lo+1), 1)
	  end if
c
	return
	end

c
      subroutine tv(nx, x, y)
c
	integer		nx, j
	complex*16	x(nx), y(nx), dd, dl, du
c
	complex*16	one
	parameter	(one = (1.0, 0.0))
c
c	Compute the matrix-vector multiplication y <-- T*x
c	where T is a nx by nx tridiagonal matrix with DD on the diagonal,
c	DL on the subdiagonal, and DU on the superdiagonal.
c
c
	dd = (4.0, 0.0)
	dl = -one
	du = -one
c
	y(1) =  dd*x(1) + du*x(2)
	do 10 j = 2, nx-1
	  y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1)
 10	  continue
	y(nx) =  dl*x(nx-1) + dd*x(nx)
	return
	end

