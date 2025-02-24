module pctracemin_functions
    use mpi	
    use tracemin_types
    implicit none 
    contains


    subroutine test_zgtracemin(M,ai,Y,E,R,ne,rank,np,comm,info) 
      type(zmatrix)                                 :: M
      type(alginfo),                  intent(inout) :: ai    !structure containing different info for the algorithm 
      double complex,  dimension(:,:),intent(inout) :: Y     !neig eigenvectros of our matrix
      double precision,dimension(:),  intent(inout) :: E     !neig eigenvalues of our matrix
      double complex,  dimension(:,:),intent(inout) :: R     !residual of the solution 
      integer,         dimension(:),  intent(in)    :: ne    !number of elements on every processor 
      integer,                        intent(in)    :: rank  !rank of the processor (in sequential version not used)
      integer,                        intent(in)    :: np    !number of processors  (in sequential version not used)
      integer,                        intent(in)    :: comm  !MPI communicator
      integer,                        intent(inout) :: info  !variable that reports errors
      !auxiliary variables for the zgtracemin routine
      double complex,  dimension(:,:),allocatable   :: Q     !temporary variable, when matrix-vector multiply is requested it is Q=M*Q where M= A or B
      double complex,  dimension(:,:),allocatable   :: Qt    !temporary variable, with dimensions of transpose of Qt
      double complex,  dimension(:,:),allocatable   :: H     !temporary variable, contains Y'*A*Y
      double complex,  dimension(:,:),allocatable   :: Hloc  !temporary variable, contains Y'*A*Y
      double complex,  dimension(:,:),allocatable   :: G     !temporary variable, contains Y'*A*Y
      double complex,  dimension(:),  allocatable   :: Zwork !temporary array of size 2*neig-1 by 1
      double precision,dimension(:),  allocatable   :: RWork !temporary array of size 3*neig-2 by 1
      integer,         dimension(:),  allocatable   :: IWork !temporary array of size neig
      !auxiliary variables for the map_chebyshev_shift routine
      double complex,  dimension(:,:),allocatable   :: S     !temporary variable 
      double complex,  dimension(:,:),allocatable   :: Vn    !temporary variable
      double complex,  dimension(:,:),allocatable   :: Vp    !temporary variable
      double complex,  dimension(:,:),allocatable   :: T     !temporary variable
      !simple variables
      integer                                       :: neig,n,nloc,i,j,mcsinfo,einfo,mpierror
      double precision                              :: gamma, ytloc, yyloc, yt, yy      
      double complex                                :: zdotc 


      neig = ai%neig
      n    = M%sizem  
      nloc = ne(rank+1)/ai%neig 
      gamma=  min(TWO/((ai%c-ai%a)*(ai%c-ai%b)), TWO/((ai%d-ai%a)*(ai%d-ai%b)))
      !write(*,*) rank,'allocate temporary variables needed by zgtracemin'
      allocate(Q    (nloc,neig))
      allocate(Qt   (neig,nloc))
      allocate(H    (neig,neig)) 
      allocate(Hloc (neig,neig)) 
      allocate(G    (neig,neig)) 
      allocate(ZWork(max(2*neig-1,1)))
      allocate(RWork(max(3*neig-2,1)))
      allocate(IWork(neig))   
      !write(*,*) rank,'allocate temporary variables needed by map_chebyshev'
      allocate(S    (nloc,neig))
      allocate(Vn   (nloc,neig))
      allocate(Vp   (nloc,neig))
      allocate(T    (nloc,neig)) 

      !write(*,*) rank,'perform the map_chebyshev_shift multiply'
      call zcopy(neig*nloc,Y,1,Q,1) 
      call zcopy(neig*nloc,Q,1,S,1)
      call pcfmv(M,Q,Vn,neig,nloc,ne,rank,np)
      do j=1,4*ai%k
         call map_chebyshev_shift(n,nloc,neig,Vn,Q,R,Vp,S,ai%k,gamma,ai%a,ai%b,j,mcsinfo)
         if (mcsinfo == 0) then
            write(*,*) rank, '--- map chebyshev shift done ---'
            exit
         end if
         if ((mcsinfo == 2) .OR. (mcsinfo == 1)) then
            call pcfmv(M,Vn,T,neig,nloc,ne,rank,np)
            call zcopy(neig*nloc,T,1,Vn,1)
         end if
         if (mcsinfo < 0) then
            write(*,*) '--- Error in the map chebyshev shift ---'
            info = mcsinfo
            exit
         end if
      end do

      !write(*,*) rank,'----- start zgtracemin -----'
      do i=1,4*ai%maxi !every iteration requires two matrix-vector multiplies (the algorithm will stop when needed by itself)
         call zgtracemin(n,neig,Y,E,Q,Qt,H,Hloc,G,ZWork,RWork,IWork,R,ai%tol,ai%maxi,i,nloc,rank,np,comm,info)
         if (info == 0) then
            write(*,*) '--- Algorithm Converged --- (iteration =',i,')'
            exit
         end if         
         if ((info == 1) .OR. (info == 2)) then
            !write(*,*) rank,'perform map_chebyshev_shift multiply' 
            call zcopy(neig*nloc,Q,1,S,1)
            call pcfmv(M,Q,Vn,neig,nloc,ne,rank,np)
            do j=1,4*ai%k
               call map_chebyshev_shift(n,nloc,neig,Vn,Q,R,Vp,S,ai%k,gamma,ai%a,ai%b,j,mcsinfo)
               if (mcsinfo == 0) then   
                  write(*,*) rank,'--- map chebyshev shift done ---'
                  exit
               end if
               if ((mcsinfo == 1) .OR. (mcsinfo == 2)) then
                  call pcfmv(M,Vn,T,neig,nloc,ne,rank,np)
                  call zcopy(neig*nloc,T,1,Vn,1)
               end if
               if (mcsinfo < 0) then
                 write(*,*) '--- Error in the map chebyshev shift ---'
                 info = mcsinfo
                 exit 
               end if
            end do
         end if
         if (info == 3) then
            write(*,*) '--- Algorithm reached maximum # of iterations --- (iteration =',i,')'
            exit
         end if
         if (info < 0) then
            write(*,*) '--- Error in the zgtracemin routine --- (iteration =',i,'; info =',info,')'
            exit 
         end if
      end do
      ai%maxi=i

      !write(*,*) rank,'A orthonormalize Y'
      einfo=1
      call extract_eigenpairs(Y,E,neig,nloc,H,Hloc,ZWork,RWork,T,rank,comm,einfo)
      call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
      call extract_eigenpairs(Y,E,neig,nloc,H,Hloc,ZWork,RWork,T,rank,comm,einfo)


      !write(*,*) rank, 'find the eigenvalues lambda using lambda=xt*A*x/xt*x where x is the eigenvector and symbol "t" is conjugate transpose'
      call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
      do i=1,neig
         ytloc= dble(zdotc(nloc,Y(1,i),1,T(1,i),1))
         yyloc= dble(zdotc(nloc,Y(1,i),1,Y(1,i),1))
         call MPI_Allreduce (ytloc, yt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
         call MPI_Allreduce (yyloc, yy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
         E(i)=yt/yy
      end do

      !write(*,*) rank, 'deallocate temporary variables'
      deallocate(Q)
      deallocate(Qt)
      deallocate(H)
      deallocate(Hloc)
      deallocate(G)
      deallocate(ZWork)
      deallocate(RWork)
      deallocate(IWork)
      deallocate(S)
      deallocate(Vn)
      deallocate(Vp)
      deallocate(T)
    end subroutine test_zgtracemin
 
    !We want to find p smallest eigenpairs of generalized eigenvalue problem I*x=lambda*A*x with Hermitian positive definite A.
    !The routine below is an implementation of generalized tracemin algorithm using complex arithmetic for the above problem. 
    !It uses reverse communication and needs to be called in an outside loop, the algorithm will stop if it converged or if it did maxi iterations.
    !The algorithm assumes that the first time its is called Q=B*Y, info=1 and that nobody touches Y,E,Q,Qt,H,info between successive calls (unless requested through info)  
    subroutine zgtracemin(n,neig,Y,E,Q,Qt,H,Hloc,G,ZWork,RWork,IWork,R,tol,maxi,ncall,nloc,rank,np,comm,info)
      integer,                        intent(in)   :: n      !order of the matrix 
      integer,                        intent(in)   :: neig   !number of eigenvalues we want to know
      double complex,  dimension(:,:),intent(inout):: Y      !neig eigenvectros of our matrix      
      double precision,dimension(:),  intent(inout):: E      !neig eigenvalues of our matrix 
      double complex,  dimension(:,:),intent(inout):: Q      !temporary variable, when matrix-vector multiply is requested it is Q=M*Q where M= A or B
      double complex,  dimension(:,:),intent(inout):: Qt     !temporary variable, with dimensions of transpose of Q
      double complex,  dimension(:,:),intent(inout):: H      !temporary variable, contains Y'*A*Y
      double complex,  dimension(:,:),intent(inout):: Hloc   !temporary variable, contains Y'*A*Y
      double complex,  dimension(:,:),intent(inout):: G      !temporary variable, contains Y'*A*Y
      double complex,  dimension(:),  intent(in)   :: Zwork  !temporary array of size 2*neig-1 by 1
      double precision,dimension(:),  intent(in)   :: RWork  !temporary array of size 3*neig-2 by 1
      integer,         dimension(:),  intent(in)   :: IWork  !temporary array of size neig 
      double complex,  dimension(:,:),intent(out)  :: R      !vector containing the residual = Y-lambda*B*Y
      double precision,               intent(in)   :: tol    !required tolerance    
      integer,                        intent(in)   :: maxi   !maximum number of iterations
      integer,                        intent(in)   :: ncall  !number of calls to the routine so far (every iteration takes two calls, starts with 1) 
      integer,                        intent(in)   :: nloc   !number of elements of a vector on the current processor
      integer,                        intent(in)   :: rank   !rank of the processor 
      integer,                        intent(in)   :: np     !number of processors  
      integer,                        intent(in)   :: comm   !MPI communicator
      integer,                        intent(inout):: info   !on exit 0  - done, 
                                                             !        1,2- perform matrix vector multiply Q=B*Q
                                                             !        3  - maximum number of iterations reached
                                                             !        <0 - error    
      integer                                      :: locinfo!used to report LAPACK errors  
      integer                                      :: i, converged, mpierror
      double precision                             :: alfa, dznrm2, Rnrmloc, Rnrm 


      !write(*,*) rank,'info=',info 
      if ((info == 1) .OR. (ncall == 1))  then
         !write(*,*) rank,'compute residual of current set of eigenvectors Y and eigenvalues E'
         call zcopy(nloc*neig,Q,1,R,1)
         converged = TRUE
         do i=1,neig 
            call zscal(nloc,dcmplx(E(i),ZERO),R(1,i),1)   
            call zaxpy(nloc,ZMONE,Y(1,i),1,R(1,i),1)
            call zscal(nloc,ZMONE,R(1,i),1)
            Rnrmloc = dznrm2(nloc,R(1,i),1)
            Rnrm    = ZERO
            call MPI_Allreduce (Rnrmloc, Rnrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
            if (tol < Rnrm) then
               converged = FALSE
            end if 
         end do
!!         if ((converged == TRUE) .OR. (ncall > 2*maxi)) then
         if ((converged == TRUE) .OR. (ncall > maxi)) then !!new
            if (converged == TRUE) then
               info=0
            else
               info=3
            end if
            return
         else
            !write(*,*) rank,'Step 2 - B-orthonormalization (the assumption is that the first time routine is called Q=B*Y)'   
            call zgemm('C','N',neig,neig,nloc,ZONE,Y,nloc,Q,nloc,ZZERO,Hloc,neig)

!!          call MPI_Allreduce (Hloc(1:neig*neig,1), H(1:neig*neig,1), neig*neig, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpierror)
            call MPI_Allreduce (Hloc(1:neig*neig,1), G(1:neig*neig,1), neig*neig, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpierror) !!new

!!          call zheev('V','L',neig,H,neig,E,ZWork,max(2*neig-1,1),RWork,locinfo)
            call zheev('V','L',neig,G,neig,E,ZWork,max(2*neig-1,1),RWork,locinfo) !!new
            if (locinfo /= 0) then
               !write(*,*) 'Error in LAPACK function zheev (info = ',locinfo,') - step 2 of tracemin' 
               info = locinfo 
               return
            end if
            do i=1,neig
               alfa = ONE/dsqrt(E(i))
!!             call zscal(neig,dcmplx(alfa,ZERO),H(1,i),1) !notice that B is hermitian positive definite (and Y is complex), so eigenvalues of H=Y'*B*Y must be real 
               call zscal(neig,dcmplx(alfa,ZERO),G(1,i),1) !!new
            end do
!!          call zgemm('N','N',nloc,neig,neig,ZONE,Y,nloc,H,neig,ZZERO,Q,nloc)
!!          call zcopy(nloc*neig,Q,1,Y,1)
            call zgemm('N','N',nloc,neig,neig,ZONE,Y,nloc,G,neig,ZZERO,R,nloc) !!new
            call zcopy(nloc*neig,R,1,Y,1)                                      !!new
            !write(*,*) rank,'Step 3 - form H'
!!          call zgemm('C','N',neig,neig,nloc,ZONE,Y,nloc,Q,nloc,ZZERO,Hloc,neig)
            call zgemm('C','N',neig,neig,nloc,ZONE,Y,nloc,R,nloc,ZZERO,Hloc,neig) !!new
            call MPI_Allreduce (Hloc(1:neig*neig,1), H(1:neig*neig,1), neig*neig, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpierror)
            call zheev('V','L',neig,H,neig,E,ZWork,max(2*neig-1,1),RWork,locinfo)
            !write(*,*) rank,'Step 4 done'
            if (locinfo /= 0) then
               !write(*,*) 'Error in LAPACK function zheev (info = ',locinfo,') - step 3 of tracemin'
               info = locinfo
               return 
            end if
            !write(*,*) rank,'Step 5 - form Y that is section of A*x=lambda*B*x'       
!!          call zgemm('N','N',nloc,neig,neig,ZONE,Y,nloc,H,neig,ZZERO,Q,nloc)
            !write(*,*) 'extra steps to for BY=V Omega^-1/2 W'        !!new
            call zgemm('N','N',nloc,neig,neig,ZONE,Q,nloc,G,neig,ZZERO,R,nloc) !!new
            call zgemm('N','N',nloc,neig,neig,ZONE,R,nloc,H,neig,ZZERO,Q,nloc) !!new 

!!            info = 2
!!            return          
!!         end if
!!      end if

!!      if (info == 2) then
         !write(*,*) rank,'Step 6 - update Y'
         call zcopy(nloc*neig,Q,1,Y,1)
         call zgemm('C','N',neig,neig,nloc,ZONE,Y,nloc,Q,nloc,ZZERO,Hloc,neig)
         call MPI_Allreduce (Hloc(1:neig*neig,1), H(1:neig*neig,1), neig*neig, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpierror)
         !in principle you could use hermitian positive definite solve, but there is no option to transpose the matrix there 
         Qt=transpose(Q)                
         call zgetrf(neig,neig,H,neig,IWork,locinfo)
         if (locinfo /= 0) then
            !write(*,*) 'Error in LAPACK function zgetrf (info = ',locinfo,') - step 6 of tracemin'
            info = locinfo
            return
         end if
         call zgetrs('T',neig,nloc,H,neig,IWork,Qt,neig,locinfo)
         if (locinfo /= 0) then
            !write(*,*) 'Error in LAPACK function zgetrs (info = ',locinfo,') - step 6 of tracemin'
            info = locinfo
            return
         end if
         Q=transpose(Qt)                    
         call zcopy(nloc*neig,Q,1,Y,1)
         info = 1         
         return
      end if

      !extra  !!new
      end if  !!new

    end subroutine zgtracemin

    !The below routine performs a map Q(x)=1-gamma*(x-a)*(x-b) then applies chebyshev polynomial to it and finally shifts by 2.
    !It also assumes that nobody touched none of the parameters between successive calls (unless requested to do so through info).
    subroutine map_chebyshev_shift(n,nloc,nrhs,Vn,Vc,Vch,Vp,S,k,gamma,a,b,ncall,info)
      integer,                        intent(in)    :: n         !order of the matrix 
      integer,                        intent(in)    :: nloc      !number of elements on the current processor
      integer,                        intent(in)    :: nrhs      !number of vectors
      double complex, dimension(:,:), intent(inout) :: Vn        !B times set of p current vectors (input: vector B*U)
      double complex, dimension(:,:), intent(inout) :: Vc        !set of p current vectors         (input: vector U by which we are multiplying,output: result of multiply) 
      double complex, dimension(:,:), intent(inout) :: Vch       !set of p temporary vectors
      double complex, dimension(:,:), intent(inout) :: Vp        !temporary array for set of p current vectors
      double complex, dimension(:,:), intent(inout) :: S         !set of original vectors used to perform the shift (input: copy of vector U, never modified in the routine)
      integer,                        intent(in)    :: k         !degree of chebyshev polynomial 
      double precision,               intent(in)    :: gamma     !parameter used in the mapping Q(x)=1-gamma*(x-a)*(x-b)
      double precision,               intent(in)    :: a         !left bound of the interval containing p desired eigenvalues 
      double precision,               intent(in)    :: b         !right bound of the interval containing p desired eigenvalues
      integer,                        intent(in)    :: ncall     !number of times routine was called (every iteration takes two calls, starts with 1)
      integer,                        intent(inout) :: info      !on exit 0  - done
                                                                 !        1,2- compute Vn=B*Vn  
                                                                 !        <0 - error 
    

      if (ncall > 2*k) then
         info=0
         call zaxpy(nloc*nrhs,ZTWO,S,1,Vc,1)
         return 
      else
         if ((info ==  1) .OR. (ncall == 1)) then
            call zaxpy(nloc*nrhs,dcmplx(MONE*b,ZERO),Vc,1,Vn,1)
            call zcopy(nloc*nrhs,Vn,1,Vch,1)
            info=2 
            return  
         end if
         if (info == 2) then
            call zaxpy(nloc*nrhs,dcmplx(MONE*a,ZERO),Vch,1,Vn,1)
            if (ncall == 2) then 
               call zscal(nloc*nrhs,dcmplx(MONE*gamma,ZERO),Vn,1)
               call zaxpy(nloc*nrhs,ZONE,Vc,1,Vn,1)
            else
               call zscal(nloc*nrhs,dcmplx(MTWO*gamma,ZERO),Vn,1)
               call zaxpy(nloc*nrhs,ZTWO,Vc,1,Vn,1)
               call zaxpy(nloc*nrhs,ZMONE,Vp,1,Vn,1)
            end if
            
            call zcopy(nloc*nrhs,Vc,1,Vp,1)
            call zcopy(nloc*nrhs,Vn,1,Vc,1)
            info=1
            return 
         end if
      end if        
    end subroutine map_chebyshev_shift

    subroutine extract_eigenpairs(Y,E,neig,nloc,H,Hloc,ZWork,RWork,T,rank,comm,info)
      double complex,  dimension(:,:),intent(inout) :: Y     !eigenvectors   
      double precision,dimension(:),  intent(inout) :: E     !eigenvalues
      integer,                        intent(in)    :: neig  !number of eiegnpairs    
      integer,                        intent(in)    :: nloc  !number of vector elements on each processor
      double complex,  dimension(:,:),intent(inout) :: H     !temporary array of size neig by neig                                                                    
      double complex,  dimension(:,:),intent(inout) :: Hloc  !temporary array of size neig by neig
      double complex,  dimension(:),  intent(inout) :: Zwork !temporary array of size 2*neig-1 by 1
      double precision,dimension(:),  intent(inout) :: RWork !temporary array of size 3*neig-2 by 1
      double complex,  dimension(:,:),intent(inout) :: T     !on second call T=A*Y
      integer,                        intent(in)    :: rank  !rank of the current processor
      integer,                        intent(in)    :: comm  !MPI communicator
      integer,                        intent(inout) :: info  !on entry first is one then is two, on exit reports errors 
      double precision                              :: alfa
      integer                                       :: i,mpierror
      
      if (info == 1) then
         !compute Omega   
         call zgemm('C','N',neig,neig,nloc,ZONE,Y,nloc,Y,nloc,ZZERO,Hloc,neig)
         call MPI_Allreduce (Hloc(1:neig*neig,1), H(1:neig*neig,1), neig*neig, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpierror)
         !compute H=Omega^-1/2 * X'*A*X * Omega^-1/2 
         do i=1,neig
            alfa = ONE/dsqrt(dble(H(i,i)))
            call zscal(nloc,dcmplx(alfa,ZERO),Y(1,i),1)
         end do

         !return to do the outside multiply
         info=2
         return
      end if
      !call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
      if (info == 2) then
         call zgemm('C','N',neig,neig,nloc,ZONE,Y,nloc,T,nloc,ZZERO,Hloc,neig)
         call MPI_Allreduce (Hloc(1:neig*neig,1), H(1:neig*neig,1), neig*neig, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpierror)
         !find eigenvectors of H
         call zheev('V','L',neig,H,neig,E,ZWork,max(2*neig-1,1),RWork,info)
         !form new Y by multiplying scaled Y by eigenvectors of H   
         call zgemm('N','N',nloc,neig,neig,ZONE,Y,nloc,H,neig,ZZERO,T,nloc)
         call zcopy(nloc*neig,T,1,Y,1)
         
         !done
         info=0
         return
      end if

      if ((info /= 1) .AND. (info /=2)) then
         write(*,*) rank,'ERROR: info is not 1, neither 2'
         return
      end if
    end subroutine extract_eigenpairs

    !parallel complex sparse hermitian matrix vector multiply (full matrix is stored)
    subroutine pcfmv(M,p,q,nrhs,nloc,ne,rank,np)
      type(zmatrix),                  intent(in)    :: M     !matrix
      double complex, dimension(:,:), intent(in)    :: p     !set of vectors p
      double complex, dimension(:,:), intent(out)   :: q     !set of vectors q
      integer,                        intent(in)    :: nrhs  !number of vectors
      integer,                        intent(in)    :: nloc  !number of elements on current processor (#of rows in p and q)
      integer,         dimension(:),  intent(in)    :: ne    !number of elements on every processor
      integer,                        intent(in)    :: rank  !rank of the processor (in sequential version not used)
      integer,                        intent(in)    :: np    !number of processors  (in sequential version not used)
      integer                                       :: i,j,k,r,c, i_low,i_high,mpierror
      double complex                                :: v
      double complex, dimension(:,:), allocatable   :: qloc  !set of temporary vectors q that will be allocated and deallocated here


      i_low =1          
      i_high=0          
      do k=0,np-1
         allocate(qloc(ne(k+1)/nrhs,nrhs))         
         qloc(1:ne(k+1)/nrhs,1:nrhs)=ZZERO
         i_high= i_high+ ne(k+1)/nrhs

         do i=1,M%rsize
            r = M%row(i)
            c = M%column(i)
            v = M%coefficient(i)
            if ((i_low <= r) .AND. (r <= i_high)) then
               !write(*,*) rank,k, 'v',v,'r',r,'c',c
               do j=1,nrhs
                  qloc(r-i_low+1,j) = qloc(r-i_low+1,j) + v*p(c,j)               
               end do
               !alternative way of doing the above loops
               !call zaxpy(nrhs,v,p(c,1),ne(rank+1)/nrhs,qloc(r-i_low+1,1),ne(k+1)/nrhs)         
            end if
         end do

         call MPI_Reduce(qloc(1:ne(k+1),1),q(1:nloc*nrhs,1),ne(k+1),MPI_DOUBLE_COMPLEX,MPI_SUM,k,MPI_COMM_WORLD,mpierror)
         i_low = i_low + ne(k+1)/nrhs
         deallocate(qloc)
      end do
    end subroutine pcfmv

    subroutine read_matrix_from_file(M, ai)
      type(zmatrix),   intent(inout) :: M                 !matrix
      type(alginfo),   intent(inout) :: ai                !structure containing different info needed by the algorithm
      integer                        :: fdescriptorA      !file descriptor
      character(len=256)             :: fpartitionA       !file name
      integer                        :: k
      double precision               :: rp, ip, shift 
      
      fpartitionA='matrixA.in'
      fdescriptorA=11        
      
      !write(*,*) 'read matrix A'
      open(fdescriptorA,file=fpartitionA,status='old')
      read(fdescriptorA,*) M%rsize, M%sizem, M%sizen, M%bsize, shift, ai%neig, ai%a, ai%b, ai%c, ai%d, ai%k, ai%maxi, ai%tol
      !write(*,*) M%rsize, M%sizem,  M%sizen, M%bsize, shift, ai%neig, ai%a, ai%b, ai%c, ai%d, ai%k, ai%maxi, ai%tol
      ai%a = ai%a + shift
      ai%b = ai%b + shift
      ai%c = ai%c + shift
      ai%d = ai%d + shift
      
      !this allocation has to be done here since before we don't know M%rsize
      allocate(M%row(1:M%rsize))
      allocate(M%column(1:M%rsize))
      allocate(M%coefficient(1:M%rsize))
      
      !assuming the input file values are in the following format: "i j (a,b)"
      !do k=1,M%rsize
      !   read(fdescriptorA,*) M%row(k), M%column(k), M%coefficient(k)
      !end do
      !assuming the input file values are in the following format: "i j a b"
      do k=1,M%rsize
         read(fdescriptorA,*) M%row(k), M%column(k), rp, ip 
         if (M%row(k) == M%column(k)) then
            M%coefficient(k) = dcmplx(rp+shift,ip)
         else
            M%coefficient(k) = dcmplx(rp,ip)
         end if
      end do
      
      !check the matrix
      !write(*,*) 'A.row',        M%row(1:10)
      !write(*,*) 'A.column',     M%column(1:10)
      !write(*,*) 'A.coefficient',M%coefficient(1:10)
      
      close(fdescriptorA)
    end subroutine read_matrix_from_file

    subroutine print_matrix(A,m,n)
      double complex, dimension(:,:), intent(in) :: A
      integer,                        intent(in) :: m
      integer,                        intent(in) :: n      
      integer                                    :: i,j   

      write(*,*) '--- print matrix ---'
      do i=1,m 
         !I'm not sure how to print complex numbers (the below statement doesn't work) 
         !write(*,"(F10.6)"), (A(i,j), j=1,n)  
      end do
    end subroutine print_matrix

    subroutine compute_matrix_splitting(M, MS, np, op)
      type(zmatrix),                 intent(in) :: M             !matrix
      type(matrix_splitting),        intent(out):: MS            !matrix splitting
      integer,                       intent(in) :: np            !number of processors
      integer,                       intent(in) :: op            !op - 0 row splitting, op - 1 column splitting (any op different from 0 is treated as 1) 
      integer                                   :: i,j,s,sz 
       
       
      allocate(MS%nvector (np))
      allocate(MS%svector (np))
      allocate(MS%evector (np))
      allocate(MS%isvector(np))
      allocate(MS%ievector(np))
      allocate(MS%ne      (np))

      !compute number, starting and ending vectors for each processor
      if(op == 0) then
         sz=M%sizem
      else
         sz=M%sizen
      end if
      do i=1,np
         MS%nvector(i) = sz/np
         MS%svector(i) = (i-1)*MS%nvector(i)+1
         MS%evector(i) = MS%svector(i)+MS%nvector(i)-1
      end do
      MS%nvector(np) = sz - MS%svector(np) + 1
      MS%evector(np) = sz
      
      !compute index of starting vector for each processor
      j=1
      do i=1,M%rsize
         if (op == 0) then
            s=M%row(i)
         else
            s=M%column(i)
         end if
         if((j <= np) .AND. (s == MS%svector(j))) then
            MS%isvector(j) = i
            j=j+1
         end if
      end do
       
      do i=1,np-1
         MS%ievector(i) = MS%isvector(i+1)-1
         MS%ne(i)       = MS%ievector(i)  -MS%isvector(i)+1
      end do
      MS%ievector(np) = M%rsize
      MS%ne(np)       = MS%ievector(np)-MS%isvector(np)+1
    end subroutine compute_matrix_splitting

    subroutine column2row_major(M)
      type(zmatrix),                 intent(inout):: M             !matrix
      integer,          dimension(:),allocatable  :: Index, ITemp
      double complex,   dimension(:),allocatable  :: ZTemp
      integer                                     :: i,j
      
      allocate(Index(M%rsize))
      allocate(ITemp(M%rsize))
      allocate(ZTemp(M%rsize))
      
      !sort row indecis
      call sort2(M%row,Index,ITemp,M%rsize)
      do i=1,M%rsize
         ITemp(i) = M%row(Index(i))
      end do      
      do i=1,M%rsize
         M%row(i)=ITemp(i)
      end do

      !sort accordingly column indecis and coefficient values
      do i=1,M%rsize
         j=Index(i)
         ITemp(i) = M%column(j)
         ZTemp(i) = M%coefficient(j)
      end do
      !update columns and coefficients
      do i=1,M%rsize
         M%column(i)     = ITemp(i)
         M%coefficient(i)= ZTemp(i)
      end do

      deallocate(Index)
      deallocate(ITemp)
      deallocate(ZTemp)
    end subroutine column2row_major

    subroutine row2column_major(M)
      type(zmatrix),                 intent(inout):: M             !matrix
      integer,          dimension(:),allocatable  :: Index, ITemp
      double complex,   dimension(:),allocatable  :: ZTemp
      integer                                     :: i,j
       
       allocate(Index(M%rsize))
       allocate(ITemp(M%rsize))
       allocate(ZTemp(M%rsize))
       
       !sort column indecis
       call sort2(M%column,Index,ITemp,M%rsize)
       do i=1,M%rsize
          ITemp(i) = M%column(Index(i))
       end do
       do i=1,M%rsize
          M%column(i)=ITemp(i)
       end do

       !sort accordingly row indecis and coefficient values
       do i=1,M%rsize
          j=Index(i)
          ITemp(i) = M%row(j)
          ZTemp(i) = M%coefficient(j)
       end do
       !update row and coefficients
       do i=1,M%rsize
          M%row(i)        = ITemp(i)
          M%coefficient(i)= ZTemp(i)
       end do

       deallocate(Index)
       deallocate(ITemp)
       deallocate(ZTemp)
     end subroutine row2column_major

     subroutine allocate_matrix_splitting(MS,np)
       type(matrix_splitting),        intent(out):: MS            !matrix splitting
       integer,                       intent(in) :: np            !number of processors

       allocate(MS%nvector (np))
       allocate(MS%svector (np))
       allocate(MS%evector (np))
       allocate(MS%isvector(np))
       allocate(MS%ievector(np))
       allocate(MS%ne      (np))
     end subroutine allocate_matrix_splitting

     subroutine MPI_Bcast_alginfo(ai,root,mpierror)
       type(alginfo),  intent(inout) :: ai
       integer,        intent(in)    :: root
       integer,        intent(out)   :: mpierror

       call MPI_Bcast(ai%neig,    1, MPI_INTEGER,         root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(ai%a,       1, MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(ai%b,       1, MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(ai%c,       1, MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(ai%d,       1, MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(ai%k,       1, MPI_INTEGER,         root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(ai%maxi,    1, MPI_INTEGER,         root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(ai%tol,     1, MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,mpierror)
     end subroutine MPI_Bcast_alginfo

     subroutine MPI_Bcast_matrix_splitting(MS,root,np,mpierror)
       type(matrix_splitting), intent(inout) :: MS 
       integer,                intent(in)    :: root
       integer,                intent(in)    :: np
       integer,                intent(out)   :: mpierror

       call MPI_Bcast(MS%nvector, np,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)       
       call MPI_Bcast(MS%svector, np,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(MS%evector, np,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(MS%isvector,np,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(MS%ievector,np,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)
       call MPI_Bcast(MS%ne,      np,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)
    end subroutine MPI_Bcast_matrix_splitting
       
end module pctracemin_functions
