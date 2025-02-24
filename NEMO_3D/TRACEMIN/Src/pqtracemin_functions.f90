module pqtracemin_functions
    use mpi
    use tracemin_types
    implicit none 
    contains

    subroutine test_ztracemin(M,ai,Y,E,R,ne,rank,np,comm,info) 
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
      double complex,  dimension(:,:),allocatable   :: H     !temporary variable, contains Y'*A*Y
      double complex,  dimension(:,:),allocatable   :: Hloc  !temporary variable, contains Y'*A*Y
      double complex,  dimension(:),  allocatable   :: Zwork !temporary array of size 2*neig-1 by 1
      double precision,dimension(:),  allocatable   :: RWork !temporary array of size 3*neig-2 by 1
      integer,         dimension(:),  allocatable   :: IWork !temporary array of size neig
      !auxiliary variables for the modified CG routine
      double complex,  dimension(:,:),allocatable   :: X     !temporary variable 
      double complex,  dimension(:,:),allocatable   :: Qcg   !temporary variable
      double complex,  dimension(:,:),allocatable   :: G     !temporary variable
      double complex,  dimension(:,:),allocatable   :: Gloc  !temporary variable
      double complex,  dimension(:,:),allocatable   :: P     !temporary variable
      double complex,  dimension(:,:),allocatable   :: Rcg   !temporary variable
      double complex,  dimension(:),  allocatable   :: gamma_old  !temporary variable
      double precision,dimension(:),  allocatable   :: rcgtol!relative cg tolerance
      double complex,  dimension(:,:),allocatable   :: T     !temporary variable
      double precision                              :: cgtol !cg tolerance
      !simple variables
      integer                                       :: neig,n,nloc,i,j,cginfo,cgmaxi,mpierror      
      double complex                                :: zdotc 
      double precision                              :: ytloc, yyloc, yt, yy


      !write(*,*) rank,'setup some variables'
      neig = ai%neig
      n    = M%sizem  
      nloc = ne(rank+1)/ai%neig 
      cgtol= ONE
      !write(*,*) rank,'allocate temporary variables needed by zgtracemin'
      allocate(Q    (nloc,neig))
      allocate(H    (neig,neig))
      allocate(Hloc (neig,neig))
      allocate(ZWork(max(2*neig-1,1)))
      allocate(RWork(max(3*neig-2,1)))
      allocate(IWork(neig))   
      !write(*,*) rank,'allocate temporary variables needed by modified CG'
      allocate(X        (nloc,neig))
      allocate(Qcg      (nloc,neig))
      allocate(G        (neig,neig))
      allocate(Gloc     (neig,neig))
      allocate(P        (nloc,neig))
      allocate(Rcg      (nloc,neig))
      allocate(gamma_old(neig))
      allocate(rcgtol(neig))
      write(*,*) rank,'allocate temporary variables'
      allocate(T        (nloc,neig))

      write(*,*) rank,'start ztracemin (every iteration requires one matrix-vector multiply and one solve done through CG)' 
      !call cshmv(M,Y,T,neig)
      !call cshmv(M,T,Q,neig)
      call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
      call pcfmv(M,T,Q,neig,nloc,ne,rank,np)

      do i=1,2*ai%maxi+1 
         call ztracemin(n,neig,Y,E,Q,H,Hloc,ZWork,RWork,IWork,R,ai%tol,ai%maxi,i,nloc,rank,np,comm,info)
         if(info == 1) then
            write(*,*) rank,'--- ztracemin has converged --- (iteration =',i/2,'; info =',info,')'
            exit
         end if
         if (info == 2) then
            write(*,*) rank,'Perform a solve'
            !setup CG parameters
            !call cshmv(M,Y,T,neig)
            !call cshmv(M,T,Qcg,neig)
            call pcfmv(M,Y,T,  neig,nloc,ne,rank,np)
            call pcfmv(M,T,Qcg,neig,nloc,ne,rank,np)

            cgmaxi=min(n,1000)
            if (mod(i,10) == 0) then
               cgtol = cgtol/10
               if (cgtol < ai%tol) then
                  cgtol=ai%tol
               end if
            end if
            !start CG
            do j=1,cgmaxi+2 !notice that one iteration is for setup
               call modified_conjugate_gradient(n,neig,X,Y,Qcg,G,Gloc,P,Rcg,gamma_old,rcgtol,j,cgmaxi,cgtol,ai%tol,nloc,rank,np,comm,cginfo)
               !write(*,*) 'rcgtol', rcgtol
               if (cginfo == 1) then
                  write(*,*) rank,'modified CG converged at iteration',j
                  exit
               end if
               if (cginfo == 2) then
                  !write(*,*) 'Perform two matrix-vector multiplies (CG)'
                  !call cshmv(M,P,T,neig)
                  !call cshmv(M,T,Qcg,neig)
                  call pcfmv(M,P,T,  neig,nloc,ne,rank,np)
                  call pcfmv(M,T,Qcg,neig,nloc,ne,rank,np)
               end if
               if (cginfo == 3) then
                  write(*,*) rank,'modified CG has reached maximum number of iterations at iteration',j
                  exit
               end if
               if (cginfo < 0) then
                  write(*,*) rank,'error in modified CG, at iteration',j
                  exit
               end if
            end do
            !copy the solution X=Y-D back into Y and compute Q=A*Y
            call zcopy(nloc*neig,X,1,Y,1)
            !call cshmv(M,Y,T,neig)
            !call cshmv(M,T,Q,neig)
            call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
            call pcfmv(M,T,Q,neig,nloc,ne,rank,np)
         end if
         if (info == 3) then
            write(*,*) rank,'Perform two matrix-vector multiplies'
            !call cshmv(M,Y,T,neig)
            !call cshmv(M,T,Q,neig)
            call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
            call pcfmv(M,T,Q,neig,nloc,ne,rank,np)
         end if
         if (info == 4) then
            write(*,*) rank,'--- Maximum number of iterations reached ---(iteration =',i/2,'; info =',info,')'
            exit
         end if
         if (info < 0) then
            write(*,*) rank,'--- Error in the ztracemin routine --- (iteration =',i/2,'; info =',info,')'
            exit 
         end if
         write(*,*) rank,'done with do loop iteration -->', i
      end do
      ai%maxi=i/2 

      !find the eigenvalues lambda using lambda=x'*A*x/x'*x where x is the eigenvector and symbol "'" is conjugate transpose
      !call cshmv(M,Y,T,neig)
      !do i=1,neig
      !   E(i)=dble(zdotc(n,Y(1,i),1,T(1,i),1))/dble(zdotc(n,Y(1,i),1,Y(1,i),1))
      !end do
      call pcfmv(M,Y,T,neig,nloc,ne,rank,np)
      do i=1,neig
         ytloc= dble(zdotc(nloc,Y(1,i),1,T(1,i),1))
         yyloc= dble(zdotc(nloc,Y(1,i),1,Y(1,i),1))
         call MPI_Allreduce (ytloc, yt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
         call MPI_Allreduce (yyloc, yy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
         E(i)=yt/yy
      end do
      !write(*,*) 'E', E
      !write(*,*) 'Y',Y 

      !refinement


      !deallocate temporary variables
      deallocate(Q)
      deallocate(H)
      deallocate(Hloc)
      deallocate(ZWork)
      deallocate(RWork)
      deallocate(IWork)
      deallocate(X)
      deallocate(Qcg)
      deallocate(G)
      deallocate(Gloc)
      deallocate(P) 
      deallocate(Rcg)
      deallocate(gamma_old)
    end subroutine test_ztracemin
 
    !We want to find p smallest eigenpairs of standard eigenvalue problem A*x=lambda*x with Hermitian positive definite A.
    !The routine below is an implementation of tracemin algorithm using complex arithmetic for the above problem. 
    !It uses reverse communication and needs to be called in an outside loop, the algorithm will stop if it converged or if it did maxi iterations.
    !The algorithm assumes that the first time its is called Q=B*Y, info=1 and that nobody touches Y,E,Q,Qt,H,info between successive calls (unless requested through info)  
    subroutine ztracemin(n,neig,Y,E,Q,H,Hloc,ZWork,RWork,IWork,R,tol,maxi,ncall,nloc,rank,np,comm,info)
      integer,                        intent(in)   :: n      !order of the matrix 
      integer,                        intent(in)   :: neig   !number of eigenvalues we want to know
      double complex,  dimension(:,:),intent(inout):: Y      !neig eigenvectros of our matrix      
      double precision,dimension(:),  intent(inout):: E      !neig eigenvalues of our matrix 
      double complex,  dimension(:,:),intent(inout):: Q      !temporary variable, on entrance Q=A*Y 
      double complex,  dimension(:,:),intent(inout):: H      !temporary variable, contains Y'*A*Y
      double complex,  dimension(:,:),intent(inout):: Hloc   !temporary variable, contains Y'*A*Y
      double complex,  dimension(:),  intent(in)   :: Zwork  !temporary array of size 2*neig-1 by 1
      double precision,dimension(:),  intent(in)   :: RWork  !temporary array of size 3*neig-2 by 1
      integer,         dimension(:),  intent(in)   :: IWork  !temporary array of size neig
      double complex,  dimension(:,:),intent(out)  :: R      !vector containing the residual = Y-lambda*B*Y
      double precision,               intent(in)   :: tol    !required tolerance    
      integer,                        intent(in)   :: maxi   !maximum number of iterations
      integer,                        intent(in)   :: ncall  !number of calls to the routine so far (every iteration takes two calls, starts with 1)  
      integer,                        intent(in)   :: nloc   !number of elements of a vector on the current processor (in sequential version not used)
      integer,                        intent(in)   :: rank   !rank of the processor                                   (in sequential version not used)
      integer,                        intent(in)   :: np     !number of processors                                    (in sequential version not used)
      integer,                        intent(in)   :: comm   !MPI communicator                                        (in sequential version not used)
      integer,                        intent(inout):: info   !on exit  1 - done
                                                             !         2 - perform solve of (I-P)A(I-P)D=(I-P)AY to find updated Y=Y-(I-P)D where P=YY'
                                                             !         3 - perform matrix vector multiply Q=A*Y
                                                             !         4 - maximum number of iterations reached
                                                             !        <0 - error    
      integer                                      :: locinfo!used to report LAPACK errors  
      integer                                      :: i, converged, mpierror
      double precision                             :: alfa, Rnrmloc, Rnrm
      double complex                               :: zdotc


      !write(*,*) rank,'inside ztracemin'
      if ((info == 2) .OR. (ncall == 1)) then
         !compute residual of current set of eigenvectors Y and eigenvalues E
         call zcopy(nloc*neig,Q,1,R,1)
         converged = TRUE
         do i=1,neig 
            call zaxpy(nloc,dcmplx(-E(i),ZERO),Y(1,i),1,R(1,i),1)
            Rnrmloc = dble(zdotc(nloc,R(1,i),1,R(1,i),1))
            call MPI_Allreduce(Rnrmloc,Rnrm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierror)
            if (tol < sqrt(Rnrm)) then  
               converged = FALSE
            end if 
         end do
         if ((converged == TRUE) .OR. (ncall > 2*maxi)) then
            if (converged == TRUE) then
               info=1
            else
               info=4
            end if
            return
         else
            !write(*,*) rank,'Y', Y
            !Step 2 - orthonormalize Y 
            call zgemm('C','N',neig,neig,nloc,ZONE,Y,nloc,Y,nloc,ZZERO,Hloc,neig)
            call MPI_Allreduce(Hloc(1:neig*neig,1),H(1:neig*neig,1),neig*neig,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpierror)
            !write(*,*) 'G',H
            call zheev('V','L',neig,H,neig,E,ZWork,max(2*neig-1,1),RWork,locinfo)
            !write(*,*) 'Hvec',H,'E',1/E
            if (locinfo /= 0) then
               write(*,*) 'Error in LAPACK function zheev (info = ',locinfo,') - step 2 of tracemin' 
               info = locinfo 
               return
            end if
            do i=1,neig
               alfa = ONE/sqrt(E(i))
               call zscal(neig,dcmplx(alfa,ZERO),H(1,i),1) 
            end do
            call zgemm('N','N',nloc,neig,neig,ZONE,Y,nloc,H,neig,ZZERO,Q,nloc)
            call zcopy(nloc*neig,Q,1,Y,1)
            info = 3
            return    !------------->do matvec outside Q=A*Y
         end if
      end if

      if (info == 3) then
         !Step 3 - form H
         call zgemm('C','N',neig,neig,nloc,ZONE,Y,nloc,Q,nloc,ZZERO,Hloc,neig)            
         call MPI_Allreduce(Hloc(1:neig*neig,1),H(1:neig*neig,1),neig*neig,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpierror)
         !write(*,*) 'H',H
         !Step 4 - find eigenvalues of H
         call zheev('V','L',neig,H,neig,E,ZWork,max(2*neig-1,1),RWork,locinfo)
         if (locinfo /= 0) then
            write(*,*) 'Error in LAPACK function zheev (info = ',locinfo,') - step 3 of tracemin'
            info = locinfo
            return 
         end if
         !Step 5 - form Y that is section of A*x=lambda*x       
         call zgemm('N','N',nloc,neig,neig,ZONE,Y,nloc,H,neig,ZZERO,Q,nloc)
         call zcopy(nloc*neig,Q,1,Y,1)
         !write(*,*) rank,'Y2',Y
         info = 2  
         return    !------------>do solve outside to find Y=Y-D      
      end if

    end subroutine ztracemin

    subroutine modified_conjugate_gradient(n,neig,X,Y,Q,G,Gloc,P,R,gamma_old,rtol,ncall,maxi,tol,gtol,nloc,rank,np,comm,info)
      integer,                          intent(in)    :: n         !order of the matrix
      integer,                          intent(in)    :: neig      !number of right hand sides
      double complex, dimension(:,:),   intent(inout) :: X         !array of size n    by neig, approximation to the solution 
      double complex, dimension(:,:),   intent(in)    :: Y         !array of size n    by neig, set of vetors used to perform (I-Y*Y')*Z multiplications             
      double complex, dimension(:,:),   intent(inout) :: Q         !array of size n    by neig, Q=A*Y on entrance and Q=A*P during iterations of CG                         
      double complex, dimension(:,:),   intent(inout) :: G         !array of size neig by neig, temporary matrix                         
      double complex, dimension(:,:),   intent(inout) :: Gloc      !array of size neig by neig, temporary matrix                         
      double complex, dimension(:,:),   intent(inout) :: P         !array of size n    by neig, search directions of CG                         
      double complex, dimension(:,:),   intent(inout) :: R         !array of size n    by neig, residual of CG             
      double complex, dimension(:),     intent(inout) :: gamma_old !array of size neig, contains gamma values from previous iteration of CG
      double precision, dimension(:),   intent(inout) :: rtol      !array of size neig,         relative tolerance computed inside this routine
      integer,                          intent(inout) :: ncall     !number of times this routine was called
      integer,                          intent(inout) :: maxi      !maximum number of CG iterations specified by the user   
      double precision,                 intent(in)    :: tol       !tolerance specified by the user
      double precision,                 intent(in)    :: gtol      !global tolerance for the tracemin specified by the user
      integer,                          intent(in)    :: nloc   !number of elements of a vector on the current processor (in sequential version not used)
      integer,                          intent(in)    :: rank   !rank of the processor                                   (in sequential version not used)
      integer,                          intent(in)    :: np     !number of processors                                    (in sequential version not used)
      integer,                          intent(in)    :: comm   !MPI communicator                                        (in sequential version not used)
      integer,                          intent(inout) :: info   !on exit  1 - done
                                                                !         2 - perform matrix vector multiply Q=A*P
                                                                !         3 - maximum number of iterations achieved 
                                                                !        <0 - error
      double complex                                  :: alfa(neig),beta(neig),tau(neig),gamma_new(neig), gamma_old_loc(neig),tau_loc(neig),gamma_new_loc(neig)        
      integer                                         :: i,j,flag, mpierror
      double complex                                  :: zdotc


      if (ncall > maxi+1) then !notice that one iteration is for setup
         info = 3
         return
      else
         if (ncall == 1) then
            !compute right hand side 
            call zcopy(nloc*neig,Q,1,R,1)
            call zgemm('C','N',neig,neig,nloc,ZONE,Y,nloc,R,nloc,ZZERO,Gloc,neig)
            call MPI_Allreduce(Gloc(1:neig*neig,1),G(1:neig*neig,1),neig*neig,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpierror)
            call zgemm('N','N',nloc,neig,neig,ZMONE,Y,nloc,G,neig,ZONE,R,nloc)                  

            !setup preliminary variables         
            call zcopy(nloc*neig,Y,1,X,1)
            call zcopy(nloc*neig,R,1,P,1)
            !write(*,*) 'p',P

            do j=1,neig
               gamma_old_loc(j) = zdotc(nloc,R(1,j),1,R(1,j),1) 
            end do
            call MPI_Allreduce(gamma_old_loc,gamma_old,neig,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpierror)
            do j=1,neig
               rtol(j)      = sqrt(dble(gamma_old(j)))*tol   !may be use just tol     
               !do not require tolerance smaller then the tracemin tolerance
               if (rtol(j) < gtol) then
                  rtol(j)   = gtol
               end if
            end do
         
            !check if you already have a solution
            flag = TRUE
            do j=1,neig
               if (sqrt(dble(gamma_old(j))) > rtol(j)) then
                  flag = FALSE
               end if
            end do
            if (flag == TRUE) then
               info = 1 !----------------> done
               return
            else
               info = 2 !----------------> do matvec
               return
            end if
         else
            !CG iteration
            do j=1,neig
               tau_loc(j)  = zdotc(nloc,Q(1,j),1,P(1,j),1)
            end do
            call MPI_Allreduce(tau_loc,tau,neig,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpierror)
            do j=1,neig
               if (tau(j) /= ZZERO) then 
                  alfa(j) = gamma_old(j)/tau(j)
               else
                  write(*,*) 'ERROR: pt*A*p is zero in the modified CG algorithm, the matrix is not s.p.d., try to recover by setting alfa=0, continuing...'  
                  alfa(j) = ZZERO
               end if
               call zaxpy(nloc,-alfa(j),P(1,j),1,X(1,j),1)   
            end do
            !write(*,*) 'alfa',alfa
            !write(*,*) 'x',X

            call zgemm('C','N',neig,neig,nloc,ZONE,Y,nloc,Q,nloc,ZZERO,Gloc,neig)
            call MPI_Allreduce(Gloc(1:neig*neig,1),G(1:neig*neig,1),neig*neig,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpierror)
            call zgemm('N','N',nloc,neig,neig,ZMONE,Y,nloc,G,neig,ZONE,Q,nloc)             
            do j=1,neig
               call zaxpy(nloc,-alfa(j),Q(1,j),1,R(1,j),1)
            end do
            !write(*,*) 'r',R

            do j=1,neig
               gamma_new_loc(j) = zdotc(nloc,R(1,j),1,R(1,j),1)
            end do
            call MPI_Allreduce(gamma_new_loc,gamma_new,neig,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpierror)
            do j=1,neig
               if (gamma_old(j) /= ZZERO) then
                  beta(j) = gamma_new(j)/gamma_old(j)
               else
                  write(*,*) 'ERROR: rt*r from the previous iteration is zero, should have stopped earlier, try to recover by setting beta=0, continuing...'
                  beta(j) = ZZERO
               end if
               gamma_old(j)= gamma_new(j) 
               call zscal(nloc,beta(j),P(1,j),1)               
            end do
            call zaxpy(nloc*neig,ZONE,R,1,P,1)  
            !write(*,*) 'beta',beta
            !write(*,*) 'p',P

            !check the convergence
            flag = TRUE 
            do j=1,neig
               if (sqrt(dble(gamma_new(j))) > rtol(j)) then 
                  flag = FALSE
               end if
            end do   
            if (flag == TRUE) then
               info = 1 !----------------> done
               return
            else
               info = 2 !----------------> do matvec
               return
            end if
         end if
      end if

      info = -1
      return
    end subroutine modified_conjugate_gradient 

    !parallel complex sparse hermitian matrix vector multiply (full matrix is stored)
    subroutine pcfmv(M,p,q,nrhs,nloc,ne,rank,np)
      type(zmatrix),                  intent(in)    :: M     !matrix
      double complex, dimension(:,:), intent(in)    :: p     !set of vectors p
      double complex, dimension(:,:), intent(out)   :: q     !set of vectors q
      integer,                        intent(in)    :: nrhs  !number of vectors
      integer,                        intent(in)    :: nloc  !number of elements on current processor (#of rows in p and q)
      integer,        dimension(:),   intent(in)    :: ne    !number of elements on every processor
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
      read(fdescriptorA,*) M%rsize, M%sizem, M%sizen, M%bsize, shift, ai%neig,  ai%maxi, ai%tol
      !write(*,*) M%rsize, M%sizem,  M%sizen, M%bsize, shift, ai%neig,  ai%maxi, ai%tol
      
      !this allocation has to be done here since before we don't know M%rsize
      allocate(M%row        (M%rsize))
      allocate(M%column     (M%rsize))
      allocate(M%coefficient(M%rsize))
      
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

end module pqtracemin_functions
