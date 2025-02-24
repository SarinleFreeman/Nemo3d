      program pzsdrv1
c
c     Message Passing Layer: MPI
c
c     Example program to illustrate the idea of reverse communication
c     in regular mode for a standard symmetric eigenvalue problem.
c
c     We implement example one of ex-sym.doc in SRC directory
c
c\Example-1
c     ... Suppose we want to solve A*x = lambda*x in regular mode,
c         where A is derived from the central difference discretization
c         of the 2-dimensional Laplacian on the unit square with
c         zero Dirichlet boundary condition.
c     ... OP = A  and  B = I.
c     ... Assume "call av(comm, nloc, nx, mv_buf, x, y)" computes y = A*x.
c     ... Use mode 1 of PZNAUPD.
c
c\BeginLib
c
c\Routines called:
c     pznaupd   Parallel ARPACK reverse communication interface routine.
c     pzneupd   Parallel ARPACK routine
c               that returns Ritz values and (optionally) Ritz vectors.
c     pdznorm2  Parallel version of Level 1 BLAS
c               that computes the norm of a complex vector.
c     zaxpy     Level 1 BLAS that computes y <- alpha*x + y.
c     av        Distributed matrix-vector multiplication routine
c               that computes A*x.
c     tv        Matrix vector multiplication routine that computes T*x,
c               where T is a tridiagonal matrix.  It is used in routine av.
c
c\Author
c     Richard Lehoucq
c     Danny Sorensen
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Parallel Modifications
c     Kristi Maschhoff
c
c\Revision history:
c     Starting Point: Complex Code FILE: ndrv1.F   SID: 2.1
c
c\SCCS Information:
c FILE: ndrv1.F   SID: 1.1   DATE OF SID: 8/13/96   RELEASE: 1
c
c\Remarks
c     1. None
c
c\EndLib
c
c---------------------------------------------------------------------------
c
        include 'mpif.h'
        include 'debug.h'
        include 'stat.h'
c
c       %---------------%
c       | MPI INTERFACE |
c       %---------------%
c
        integer         comm, myid, nprocs, rc, nloc
c
c       %-----------------------------%
c       | Define maximum dimensions   |
c       | for all arrays.             |
c       | MAXNLOC:Maximum dimension   |
c       |         of the A allowed.   |
c       | MAXNEV: Maximum NEV allowed |
c       | MAXNCV: Maximum NCV allowed |
c       %-----------------------------%
c
        integer         processors, blocks
        parameter       (processors = PROCESSORS, blocks = BLOCKS)
        integer         maxnloc, maxnev, maxncv, ldv
c       parameter       (maxnloc=256, maxnev=12, maxncv=30, ldv=maxnloc)
        parameter       (maxnloc=processors*blocks*blocks,
     &                  maxnev=VALUES, maxncv=VECTORS, ldv=maxnloc)
c
c       %--------------%
c       | Local Arrays |
c       %--------------%
c
        integer         iparam(11), ipntr(14)
        logical         select(maxncv)
        complex*16      ax(maxnloc), d(maxncv),
     &                  v(ldv,maxncv), workd(3*maxnloc),
     &                  workev(3*maxncv), resid(maxnloc),
     &                  workl(3*maxncv*maxncv+5*maxncv)
        double precision
     &                  rwork(maxncv), rd(maxncv,3)
c
c       %---------------%
c       | Local Scalars |
c       %---------------%
c
        character       bmat*1, which*2
        integer         ido, n, nx, nev, ncv, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode
        complex*16      sigma
        double precision
     &                  tol
        logical         rvec
c
c       %----------------------------------------------%
c       | Local Buffers needed for MPI communication |
c       %----------------------------------------------%
c
        complex*16      mv_buf(maxnloc)
c
        double precision
     &                  zero
        parameter       (zero = 0.0)
c
c       %-----------------------------%
c       | BLAS & LAPACK routines used |
c       %-----------------------------%
c
        double precision
     &                  pdznorm2
        external        pdznorm2, zaxpy
      integer   multiplication
c
c       %-----------------------%
c       | Executable Statements |
c       %-----------------------%
c
        call MPI_INIT(ierr)
        comm = MPI_COMM_WORLD
        call MPI_COMM_RANK(comm, myid, ierr)
        call MPI_COMM_SIZE(comm, nprocs, ierr)
c
        ndigit = -3
        logfil = +6
        mcaupd = +1
c
c       %----------------------------------------------------%
c       | The number NX is the number of interior points     |
c       | in the discretization of the 2-dimensional         |
c       | Laplacian on the unit square with zero Dirichlet   |
c       | boundary condition.  The number N(=NX*NX) is the   |
c       | dimension of the matrix.  A standard eigenvalue    |
c       | problem is solved (BMAT = 'I'). NEV is the number  |
c       | of eigenvalues to be approximated.  The user can   |
c       | modify NEV, NCV, WHICH to solve problems of        |
c       | different sizes, and to get different parts of the |
c       | spectrum.  However, The following conditions must  |
c       | be satisfied:                                      |
c       |                   N <= MAXNLOC,                    |
c       |                 NEV <= MAXNEV,                     |
c       |             NEV + 2 <= NCV <= MAXNCV               |
c       %----------------------------------------------------%
c
c       nx      = 10
        nx      = processors*blocks
        n       = nx*nx
c       nev     = 4
c       ncv     = 20
        nev     = maxnev
        ncv     = maxncv
c
c       %--------------------------------------%
c       | Set up distribution of data to nodes |
c       %--------------------------------------%
c
        nloc = (nx/nprocs)*nx
        if (mod(nx, nprocs) .gt. myid) nloc = nloc + nx
        if (nloc .lt. ncv) ncv = nloc
c
        if (nloc .gt. maxnloc) then
          print *, ' ERROR with _NDRV1: NLOC is greater than MAXNLOC '
          go to 9000
        else if (nev .gt. maxnev) then
          print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
          go to 9000
        else if (ncv .gt. maxncv) then
          print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
          go to 9000
          end if
        bmat  = 'I'
c       which = 'SM'
        which = 'LM'
c
c       %--------------------------------------------------%
c       | The work array WORKL is used in PZNAUPD as       |
c       | workspace.  Its dimension LWORKL is set as       |
c       | illustrated below.  The parameter TOL determines |
c       | the stopping criterion.  If TOL <= 0, machine    |
c       | precision is used.  The variable IDO is used for |
c       | reverse communication and is initially set to 0. |
c       | Setting INFO = 0 indicates that a random vector  |
c       | is generated in PZNAUPD to start the Arnoldi     |
c       | iteration.                                       |
c       %--------------------------------------------------%
c
        lworkl  = 3*ncv**2+5*ncv
c       tol     = zero
        tol     = 1.0D-8
        ido     = 0
        info    = 0
c
c       %---------------------------------------------------%
c       | This program uses exact shifts with respect to    |
c       | the current Hessenberg matrix (IPARAM(1) = 1).    |
c       | IPARAM(3) specifies                               |
c       | the maximum number of Arnoldi iterations allowed. |
c       | Mode 1 of PZNAUPD is used (IPARAM(7) = 1).        |
c       | All these options can be changed by the user.     |
c       | For details, see the documentation in ZNAUPD.     |
c       %---------------------------------------------------%
c
        ishfts = 1
c       maxitr = 300
        maxitr = MAXITR
        mode   = 1
c
        iparam(1) = ishfts
        iparam(3) = maxitr
        iparam(7) = mode
c
c       %-------------------------------------------%
c       | M A I N   L O O P (Reverse communication) |
c       %-------------------------------------------%
c
      multiplication = 0
 10     continue
c
c         %---------------------------------------------%
c         | Repeatedly call the routine PZNAUPD         |
c         | and take actions indicated by parameter IDO |
c         | until either convergence is indicated       |
c         | or maxitr has been exceeded.                |
c         %---------------------------------------------%
c
          call pznaupd(comm, ido, bmat, nloc, which,
     &          nev, tol, resid, ncv, v, ldv, iparam, ipntr,
     &          workd, workl, lworkl, rwork, info)
c
      multiplication = multiplication + 1
      if ((0 .eq. myid) .and. (0 .eq. mod(multiplication, 1024))) then
        print *, multiplication, ' = multiplication'
        end if
          if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %------------------------------------------------%
c           | Perform matrix-vector multiplication           |
c           |              y <-- OP*x                        |
c           | The user should supply his/her own             |
c           | matrix-vector multiplication routine here      |
c           | that takes workd(ipntr(1)) as the input vector |
c           | and returns the result to workd(ipntr(2)).     |
c           %------------------------------------------------%
c
            call av(comm, nloc, nx, mv_buf,
     &          workd(ipntr(1)), workd(ipntr(2)))
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call PZNAUPD again. |
c           %-----------------------------------------%
c
            go to 10
            end if
c
c       %----------------------------------------%
c       | Either we have convergence or there is |
c       | an error.                              |
c       %----------------------------------------%
c
        if (info .lt. 0) then
c
c         %--------------------------%
c         | Error message. Check the |
c         | documentation in PZNAUPD.|
c         %--------------------------%
c
c         if (myid .eq. 0) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation in _naupd'
            print *, ' '
c           end if
c
        else
c
c         %-------------------------------------------%
c         | No fatal errors occurred.                 |
c         | Post-Process using PZNEUPD.               |
c         |                                           |
c         | Computed eigenvalues may be extracted.    |
c         |                                           |
c         | Eigenvectors may also be computed now if  |
c         | desired.  (indicated by rvec = .true.)    |
c         %-------------------------------------------%
c
          rvec = .true.
c
          call pzneupd(comm, rvec, 'A', select, d, v, ldv, sigma,
     &          workev, bmat, nloc, which, nev, tol, resid, ncv, v, ldv,
     &          iparam, ipntr, workd, workl, lworkl, rwork, info)
c
c         %----------------------------------------------%
c         | Eigenvalues are returned in the one          |
c         | dimensional array D.  The corresponding      |
c         | eigenvectors are returned in the first NCONV |
c         | (=IPARAM(5)) columns of the two dimensional  |
c         | array V if requested.  Otherwise, an         |
c         | orthogonal basis for the invariant subspace  |
c         | corresponding to the eigenvalues in D is     |
c         | returned in V.                               |
c         %----------------------------------------------%
c
          if (info .ne. 0) then
c
c           %------------------------------------%
c           | Error condition:                   |
c           | Check the documentation of PZNEUPD.|
c           %------------------------------------%
c
            if (myid .eq. 0) then
              print *, ' '
              print *, ' Error with _neupd, info = ', info
              print *, ' Check the documentation of _neupd.'
              print *, ' '
              end if
c
          else
c
            nconv = iparam(5)
            do 20 j = 1, nconv
c
c             %---------------------------%
c             | Compute the residual norm |
c             |                           |
c             |   ||  A*x - lambda*x ||   |
c             |                           |
c             | for the NCONV accurately  |
c             | computed eigenvalues and  |
c             | eigenvectors.  (iparam(5) |
c             | indicates how many are    |
c             | accurate to the requested |
c             | tolerance)                |
c             %---------------------------%
c
              call av(comm, nloc, nx, mv_buf, v(1,j), ax)
              call zaxpy(nloc, -d(j), v(1,j), 1, ax, 1)
              rd(j, 1) = dble(d(j))
              rd(j, 2) = dimag(d(j))
              rd(j, 3) = pdznorm2(comm, nloc, ax, 1)
c
 20           continue
c
c           %-----------------------------%
c           | Display computed residuals. |
c           %-----------------------------%
c
            call pdmout(comm, 6, nconv, 3, rd, maxncv, -6,
     &          'Ritz values (Real, Imag) and direct residuals')
            end if
c
c         %-------------------------------------------%
c         | Print additional convergence information. |
c         %-------------------------------------------%
c
          if (myid .eq. 0)then
            if (info .eq. 1) then
              print *, ' '
              print *, ' Maximum number of iterations reached.'
              print *, ' '
            else if (info .eq. 3) then
              print *, ' '
              print *, ' No shifts could be applied during implicit
     &          Arnoldi update, try increasing NCV.'
              print *, ' '
              end if
c
            print *, ' '
            print *, '_NDRV1'
            print *, '====== '
            print *, ' '
            print *, ' Size of the matrix is ', n
            print *, ' The number of processors is ', nprocs
            print *, ' The number of Ritz values requested is ', nev
            print *, ' The number of Arnoldi vectors generated',
     &          ' (NCV) is ', ncv
            print *, ' What portion of the spectrum: ', which
            print *, ' The number of converged Ritz values is ', nconv
            print *, ' The number of Implicit Arnoldi update',
     &          ' iterations taken is ', iparam(3)
            print *, ' The number of OP*x is ', iparam(9)
            print *, ' The convergence criterion is ', tol
            print *, ' '
            end if
c
          end if
c
c       %----------------------------%
c       | Done with program pzndrv1. |
c       %----------------------------%
c
 9000   continue
c
c
c       %-------------------------%
c       | Release resources MPI   |
c       %-------------------------%
c
        call MPI_FINALIZE(rc)
c
        end
c
c     parallel matrix-vector multiplication subroutine
c
c     The matrix used is the 2 dimensional discrete Laplacian
c     on the unit square with zero Dirichlet boundary condition.
c
c     Computes w <-- OP*v,
c     where OP is the nx*nx by nx*nx block tridiagonal matrix.
c
c               | T -I          |
c               |-I  T -I       |
c          OP = |   -I  T       |
c               |        ...  -I|
c               |           -I T|
c
c     The subroutine TV is called to computed y <-- T*x.
c
      subroutine av(comm, nloc, nx, mv_buf, v, w)
c
c       .. MPI Declarations ...
        include 'mpif.h'
        integer comm, nprocs, myid, ierr, status(MPI_STATUS_SIZE)
        integer         nloc, nx, np, j, lo, next, prev
        complex*16      v(nloc), w(nloc), mv_buf(nx), one
        parameter       (one = (1.0D+0, 0.0D+0))
        external        zaxpy

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
  10        continue
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
     &          next, myid+1, comm, ierr)
          end if
        if (myid .gt. 0) then
          call mpi_recv(mv_buf, nx, MPI_DOUBLE_COMPLEX,
     &          prev, myid, comm, status, ierr)
          call zaxpy(nx, -one, mv_buf, 1, w(1), 1)
          end if
c
        if (myid .gt. 0) then
          call mpi_send(v(1), nx, MPI_DOUBLE_COMPLEX,
     &          prev, myid-1, comm, ierr)
          end if
        if (myid .lt. nprocs-1) then
          call mpi_recv(mv_buf, nx, MPI_DOUBLE_COMPLEX,
     &          next, myid, comm, status, ierr)
          call zaxpy(nx, -one, mv_buf, 1, w(lo+1), 1)
          end if
c
        return
        end
c
      subroutine tv(nx, x, y)
c
        integer         nx, j
        complex*16      x(nx), y(nx), dd, dl, du
c
        complex*16      one
        parameter       (one = (1.0D+0, 0.0D+0))
c
c       Compute the matrix-vector multiplication y <-- T*x
c       where T is a nx by nx tridiagonal matrix with DD on the diagonal,
c       DL on the subdiagonal, and DU on the superdiagonal.
c
c
        dd = (4.0D+0, 0.0D+0)
        dl = -one
        du = -one
c
        y(1) =  dd*x(1) + du*x(2)
        do 10 j = 2, nx-1
          y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1)
 10       continue
        y(nx) =  dl*x(nx-1) + dd*x(nx)
        return
        end

