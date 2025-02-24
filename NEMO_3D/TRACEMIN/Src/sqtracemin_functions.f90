module sqtracemin_functions
    use tracemin_types
    !use zbicgstab
    implicit none 
    contains
 

    subroutine test_ztracemin(M,ai,Y,E,R,rank,np,comm,info) 
      type(zmatrix)                                 :: M
      type(alginfo),                  intent(inout) :: ai    !structure containing different info for the algorithm 
      double complex,  dimension(:,:),intent(inout) :: Y     !neig eigenvectros of our matrix
      double precision,dimension(:),  intent(inout) :: E     !neig eigenvalues of our matrix
      double complex,  dimension(:,:),intent(inout) :: R     !residual of the solution 
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
      integer                                       :: neig,n,i,j,cginfo,cgmaxi      
      double complex                                :: zdotc 
      !BiCGstab variables
      double complex,  dimension(:,:),allocatable   :: Work_bicgstab,RWork_bicgstab
      integer,         dimension(:),  allocatable   :: IWork_bicgstab
      integer                                       :: l,ldw,ldrw,bicgstabinfo,bicgstabmaxi,row_index,column_index
      integer                                       :: flag_refine


      write(*,*) 'setup some variables'
      neig = ai%neig
      n    = M%sizem  
      cgtol= ONE
      write(*,*) 'allocate temporary variables needed by zgtracemin'
      allocate(Q    (n,neig))
      allocate(H    (neig,neig))
      allocate(Hloc (neig,neig))
      allocate(ZWork(max(2*neig-1,1)))
      allocate(RWork(max(3*neig-2,1)))
      allocate(IWork(neig))   
      write(*,*) 'allocate temporary variables needed by modified CG'
      allocate(X        (n,neig))
      allocate(Qcg      (n,neig))
      allocate(G        (neig,neig))
      allocate(Gloc     (neig,neig))
      allocate(P        (n,neig))
      allocate(Rcg      (n,neig))
      allocate(gamma_old(neig))
      allocate(rcgtol(neig))
      write(*,*) 'allocate temporary variables'
      allocate(T        (n,neig))

      write(*,*) 'start ztracemin (every iteration requires one matrix-vector multiply and one solve done through CG)' 
      call cshmv(M,Y,T,neig)
      call cshmv(M,T,Q,neig)

      do i=1,2*ai%maxi+1 
         call ztracemin(n,neig,Y,E,Q,H,Hloc,ZWork,RWork,IWork,R,ai%tol,ai%maxi,i,0,rank,np,comm,info)
         if(info == 1) then
            write(*,*) '--- ztracemin has converged --- (iteration =',i/2,'; info =',info,')'
            exit
         end if
         if (info == 2) then
            write(*,*) 'Perform a solve'
            !setup CG parameters
            call cshmv(M,Y,T,neig)
            call cshmv(M,T,Qcg,neig)
            cgmaxi=min(n,1000)
            if (mod(i,10) == 0) then
               cgtol = cgtol/10
               if (cgtol < ai%tol) then
                  cgtol=ai%tol
               end if
            end if
            !start CG
            do j=1,cgmaxi+2 !notice that one iteration is for setup
               call modified_conjugate_gradient(n,neig,X,Y,Qcg,G,Gloc,P,Rcg,gamma_old,rcgtol,j,cgmaxi,cgtol,ai%tol,0,rank,np,comm,cginfo)
               !write(*,*) 'rcgtol', rcgtol
               if (cginfo == 1) then
                  write(*,*) 'modified CG converged at iteration',j
                  exit
               end if
               if (cginfo == 2) then
                  !write(*,*) 'Perform two matrix-vector multiplies (CG)'
                  call cshmv(M,P,T,neig)
                  call cshmv(M,T,Qcg,neig)
               end if
               if (cginfo == 3) then
                  write(*,*) 'modified CG has reached maximum number of iterations at iteration',j
                  exit
               end if
               if (cginfo < 0) then
                  write(*,*) 'error in modified CG, at iteration',j
                  exit
               end if
            end do
            !copy the solution X=Y-D back into Y and compute Q=A*Y
            call zcopy(n*neig,X,1,Y,1)
            call cshmv(M,Y,T,neig)
            call cshmv(M,T,Q,neig)
         end if
         if (info == 3) then
            write(*,*) 'Perform two matrix-vector multiplies'
            call cshmv(M,Y,T,neig)
            call cshmv(M,T,Q,neig)
         end if
         if (info == 4) then
            write(*,*) '--- Maximum number of iterations reached ---(iteration =',i/2,'; info =',info,')'
            exit
         end if
         if (info < 0) then
            write(*,*) '--- Error in the ztracemin routine --- (iteration =',i/2,'; info =',info,')'
            exit 
         end if
         write(*,*) 'done with do loop iteration -->', i
      end do
      ai%maxi=i/2 

      !find the eigenvalues lambda using lambda=x'*A*x/x'*x where x is the eigenvector and symbol "'" is conjugate transpose
      call cshmv(M,Y,T,neig)
      do i=1,neig
         E(i)=dble(zdotc(n,Y(1,i),1,T(1,i),1))/dble(zdotc(n,Y(1,i),1,Y(1,i),1))
      end do
      !write(*,*) 'E', E
      !write(*,*) 'Y',Y 

      !refinement
      !flag_refine = 1
      !if (flag_refine == 1) then
      !   bicgstabmaxi = min(n,1000)
      !   l   = 2
      !   ldw = n*(3+2*(l+1))
      !   ldrw= (l+1)*(3+2*(l+1))
      !   allocate(Work_bicgstab (n,3+2*(l+1)))
      !   allocate(RWork_bicgstab(l+1,3+2*(l+1)))
      !   allocate(IWork_bicgstab(l+1))
      !   do i=1,neig 
      !      !write(*,*) 'compute matrix A-lambdaI' 
      !      do j=1,M%rsize
      !         row_index    = M%row(j)
      !         column_index = M%column(j)
      !         if (row_index == column_index) then
      !            if (i /= 1) then 
      !               M%coefficient(j) =  M%coefficient(j) - E(i) + E(i-1)
      !            else
      !               M%coefficient(j) =  M%coefficient(j) - E(i) 
      !            end if
      !         end if
      !      end do
      !    
      !      !write(*,*) 'form rhs T = -(I-YYc)*(A-a*I)*Yi'  
      !      call cshmv(M,Y(:,i:i),T(:,i:i),1)
      !      call zgemm('C','N',neig,1,n,ZONE, Y,n,T(1,i),n,ZZERO,H(1,i),neig)                              
      !      call zgemm('N','N',n,1,neig,ZONE,Y,n,H(1,i),neig,ZMONE,T(1,i),n)
      !
      !      !write(*,*) 'solve system (I-YYc)*(A-a*I)*(I-YYc)Qi = Ti for Qi using BiCGstab'
      !      call zbistbl(M,Y,neig, l,n,Q(:,i),T(:,i),matvec_bicgstab,precondition_bicgstab,ai%tol,bicgstabmaxi, &
      !                   Work_bicgstab,ldw,Rwork_bicgstab,ldrw,IWork_bicgstab, Qcg(:,1:1),H(:,1), bicgstabinfo)
      !      if (bicgstabinfo == 0) then
      !         write(*,*) 'BiCGstab has converged, in',bicgstabmaxi,'iterations, info=',bicgstabinfo
      !         !write(*,*) 'solution:',Q(:,i)
      !      else
      !         if(bicgstabinfo < 0) then
      !            write(*,*) 'ERROR in BiCGstab, no update done, info=',bicgstabinfo
      !         else
      !            if (bicgstabinfo == 1) then
      !               write(*,*) 'BiCGstab did not converge, info=',bicgstabinfo
      !            else
      !               write(*,*) 'BiCGstab broke down, no update done, info=',bicgstabinfo
      !            end if
      !         end if
      !         Q(:,i)=ZZERO
      !      end if
      !   end do
      !
      !   !write(*,*) 'return matrix to its original form and deallocate temporary vectors'
      !   do j=1,M%rsize
      !      row_index    = M%row(j)
      !      column_index = M%column(j)
      !      if (row_index == column_index) then
      !         M%coefficient(j) =  M%coefficient(j) + E(neig)
      !      end if
      !   end do
      !   deallocate(Work_bicgstab)
      !   deallocate(RWork_bicgstab)
      !   deallocate(IWork_bicgstab)
      !	
      !   !write(*,*) 'compute D=(I-YYc)*Q and update eigenvectors Y=Y+D and recompute eigenvalues E'
      !   call zgemm('C','N',neig,neig,n,ZONE,Y,n,Q,n,ZZERO,H,neig)
      !   call zgemm('N','N',n,neig,neig,ZMONE,Y,n,H,neig,ZONE,Q,n)
      !   call zaxpy(n*neig,ZONE,Q,1,Y,1)
      !   call cshmv(M,Y,T,neig)
      !   do i=1,neig
      !      E(i)=dble(zdotc(n,Y(1,i),1,T(1,i),1))/dble(zdotc(n,Y(1,i),1,Y(1,i),1))
      !   end do
      !end if
      !end refinement


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
      integer                                      :: i, converged
      double precision                             :: alfa, dznrm2 



      if ((info == 2) .OR. (ncall == 1)) then
         !compute residual of current set of eigenvectors Y and eigenvalues E
         call zcopy(n*neig,Q,1,R,1)
         converged = TRUE
         do i=1,neig 
            call zaxpy(n,dcmplx(-E(i),ZERO),Y(1,i),1,R(1,i),1)
            !write(*,*) 'resid norm', dznrm2(n,R(1,i),1)
            if (tol < dznrm2(n,R(1,i),1)) then
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
            !write(*,*) 'Y', Y
            !Step 2 - orthonormalize Y 
            call zgemm('C','N',neig,neig,n,ZONE,Y,n,Y,n,ZZERO,H,neig)
            call zheev('V','L',neig,H,neig,E,ZWork,max(2*neig-1,1),RWork,locinfo)
            if (locinfo /= 0) then
               write(*,*) 'Error in LAPACK function zheev (info = ',locinfo,') - step 2 of tracemin' 
               info = locinfo 
               return
            end if
            do i=1,neig
               alfa = ONE/sqrt(E(i))
               call zscal(neig,dcmplx(alfa,ZERO),H(1,i),1) 
            end do
            call zgemm('N','N',n,neig,neig,ZONE,Y,n,H,neig,ZZERO,Q,n)
            call zcopy(n*neig,Q,1,Y,1)
            info = 3
            return    !------------->do matvec outside Q=A*Y
         end if
      end if

      if (info == 3) then
         !Step 3 - form H
         call zgemm('C','N',neig,neig,n,ZONE,Y,n,Q,n,ZZERO,H,neig)            
         !write(*,*) 'H',H
         !Step 4 - find eigenvalues of H
         call zheev('V','L',neig,H,neig,E,ZWork,max(2*neig-1,1),RWork,locinfo)
         if (locinfo /= 0) then
            write(*,*) 'Error in LAPACK function zheev (info = ',locinfo,') - step 3 of tracemin'
            info = locinfo
            return 
         end if
         !Step 5 - form Y that is section of A*x=lambda*x       
         call zgemm('N','N',n,neig,neig,ZONE,Y,n,H,neig,ZZERO,Q,n)
         call zcopy(n*neig,Q,1,Y,1)
         !write(*,*) 'Y2',Y
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
      double complex                                  :: gamma_new(neig), alfa(neig), beta(neig), tau(neig)        
      integer                                         :: i,j,flag
      double complex                                  :: zdotc


      if (ncall > maxi+1) then !notice that one iteration is for setup
         info = 3
         return
      else
         if (ncall == 1) then
            !compute right hand side 
            call zcopy(n*neig,Q,1,R,1)
            call zgemm('C','N',neig,neig,n,   ZONE, Y,n,R,n,   ZZERO,G,neig)       
            call zgemm('N','N',n,   neig,neig,ZMONE,Y,n,G,neig,ZONE, R,n)                  

            !setup preliminary variables         
            call zcopy(n*neig,Y,1,X,1)
            call zcopy(n*neig,R,1,P,1)
            !write(*,*) 'p',P

            do j=1,neig
               gamma_old(j) = zdotc(n,R(1,j),1,R(1,j),1) 
               !rtol(j)     = tol
               rtol(j)      = sqrt(dble(gamma_old(j)))*tol     
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
               tau(j)  = zdotc(n,Q(1,j),1,P(1,j),1)
               if (tau(j) /= ZZERO) then 
                  alfa(j) = gamma_old(j)/tau(j)
               else
                  write(*,*) 'ERROR: pt*A*p is zero in the modified CG algorithm, the matrix is not s.p.d., try to recover by setting alfa=0, continuing...'  
                  alfa(j) = ZZERO
               end if
               call zaxpy(n,-alfa(j),P(1,j),1,X(1,j),1)   
            end do
            !write(*,*) 'alfa',alfa
            !write(*,*) 'x',X

            call zgemm('C','N',neig,neig,n,   ZONE, Y,n,Q,n,   ZZERO,G,neig)
            call zgemm('N','N',n,   neig,neig,ZMONE,Y,n,G,neig,ZONE, Q,n)             
            do j=1,neig
               call zaxpy(n,-alfa(j),Q(1,j),1,R(1,j),1)
            end do
            !write(*,*) 'r',R

            do j=1,neig
               gamma_new(j) = zdotc(n,R(1,j),1,R(1,j),1)
               if (gamma_old(j) /= ZZERO) then
                  beta(j) = gamma_new(j)/gamma_old(j)
               else
                  write(*,*) 'ERROR: rt*r from the previous iteration is zero, should have stopped earlier, try to recover by setting beta=0, continuing...'
                  beta(j) = ZZERO
               end if
               gamma_old(j)= gamma_new(j) 
               call zscal(n,beta(j),P(1,j),1)               
            end do
            call zaxpy(n*neig,ZONE,R,1,P,1)  
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
    
    !complex sparse hermitian matrix vector multiply
    subroutine cshmv(M,p,q,nrhs)
      type(zmatrix),                  intent(in)    :: M         !matrix
      double complex, dimension(:,:), intent(in)    :: p         !set of vectors p
      double complex, dimension(:,:), intent(out)   :: q         !set of vectors q
      integer,                        intent(in)    :: nrhs      !number of vectors
      integer                                       :: i,k,r,c
      double complex                                :: v
      
      !write(*,*) '---------p-----------',p
      q(1:M%sizem,1:nrhs)=ZZERO
      do k=1,M%rsize
         r = M%row(k)
         c = M%column(k)
         v = M%coefficient(k)
         !write(*,*),'r',r,'c',c,'v',v
         do i=1,nrhs
            q(r,i) = q(r,i) + v*p(c,i)
            if (r /= c) then
               q(c,i) = q(c,i) + dconjg(v)*p(r,i)
            end if
         end do
         !alternative way of doing the above loops
         !call zaxpy(nrhs,v,p(c,1),M%sizem,q(r,1),M%sizem)
         !if (r /= c) then
         !  call zaxpy(nrhs,dconjg(v),p(r,1),M%sizem,q(c,1),M%sizem)
         !end if
      end do
      !write(*,*) '--------------q------------',q

    end subroutine cshmv

    !complex sparse hermitian matrix vector multiply by transpose of the matrix M
    subroutine cshmvt(M,p,q,nrhs)
      type(zmatrix),                  intent(in)    :: M         !matrix
      double complex, dimension(:,:), intent(in)    :: p         !set of vectors p
      double complex, dimension(:,:), intent(out)   :: q         !set of vectors q
      integer,                        intent(in)    :: nrhs      !number of vectors
      integer                                       :: i,k,r,c
      double complex                                :: v


      q(1:M%sizem,1:nrhs)=ZZERO
      do k=1,M%rsize
         c = M%row(k)
         r = M%column(k)
         v = M%coefficient(k)
         do i=1,nrhs
            q(r,i) = q(r,i) + v*p(c,i)
            if (r /= c) then
               q(c,i) = q(c,i) + dconjg(v)*p(r,i)
            end if
         end do
         !alternative way of doing the above loops
         !call zaxpy(nrhs,v,p(c,1),M%sizem,q(r,1),M%sizem)
         !if (r /= c) then
         !  call zaxpy(nrhs,dconjg(v),p(r,1),M%sizem,q(c,1),M%sizem)
         !end if
      end do
    end subroutine cshmvt

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
      read(fdescriptorA,*) M%rsize, M%sizem, M%bsize, shift, ai%neig, ai%maxi, ai%tol
      !write(*,*) M%rsize, M%sizem, M%bsize, shift, ai%neig, ai%maxi, ai%tol
      
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

      !write(*,*) 'Row',        M%row
      !write(*,*) 'Column',     M%column
      !write(*,*) 'Coefficient',M%coefficient

      close(fdescriptorA)
    end subroutine read_matrix_from_file


    subroutine matvec_bicgstab(M,Y,neig, n,Pv,Qv, Tv,Hv)
      type(zmatrix)                               :: M
      !!double complex, dimension(:,:),intent(in) :: Y
      double complex                              :: Y(n,neig)
      integer,                       intent(in)   :: neig
      integer,                       intent(in)   :: n
      !!double complex, dimension(:),  intent(in) :: Pv 
      !!double complex, dimension(:,:),intent(out):: Qv !the second dimension is one
      !!double complex, dimension(:,:)            :: Tv !the second dimension is one
      !!double complex, dimension(:)              :: Hv
      double complex                              :: Pv(n)
      double complex                              :: Qv(n,1)
      double complex                              :: Tv(n,1)
      double complex                              :: Hv(neig)

      
      !write(*,*) 'compute Qv=(I-YYc)*M*(I-YYc)*Pv'
      !write(*,*) 'Msz',M%sizem
      !write(*,*) 'neig',neig
      !write(*,*) 'n',n
      !write(*,*) 'Y',Y
      !write(*,*) 'Pv',Pv
      !write(*,*) 'Qv',Qv
      !write(*,*) 'Tv',Tv(1,1)
      !write(*,*) 'Hv',Hv

      call zcopy(n,Pv,1,Tv,1)                                      
      call zgemm('C','N',neig,1,n,ZONE,Y,n,Tv,n,ZZERO,Hv,neig)     
      call zgemm('N','N',n,1,neig,ZMONE,Y,n,Hv,neig,ZONE,Tv,n)
      call cshmv(M,Tv,Qv,1)
      call zgemm('C','N',neig,1,n,ZONE, Y,n,Qv,n,ZZERO,Hv,neig)
      call zgemm('N','N',n,1,neig,ZMONE,Y,n,Hv,neig,ZONE,Qv,n)
    end subroutine matvec_bicgstab
    
    subroutine precondition_bicgstab(n,x)
      integer,                      intent(in)    :: n
      !!double complex, dimension(:), intent(inout) :: x
      double complex                              :: x(n)

      !there is no preconditioning currently so this routine returns unmodified x
    end subroutine precondition_bicgstab

end module sqtracemin_functions
