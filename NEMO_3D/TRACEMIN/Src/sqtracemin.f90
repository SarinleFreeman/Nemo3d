program sqtracemin
      use tracemin_types
      use sqtracemin_functions
      implicit none

      !----- Varibales Declaration -----
      !general variables
      type(zmatrix)                              :: M                     !matrix B in generalized eigenvalue problem I*x=lambda*B*x 
      type(alginfo)                              :: ai                    !structure containing different information for the algorithm 
      integer                                    :: iseed(4)      
      double complex,dimension(:,:), allocatable :: result_v, R
      double precision,dimension(:), allocatable :: result_e, Rnrm, iRnrm 
      double precision                           :: t, t2, dznrm2, timer
      integer                                    :: i, info    
      

      write(*,*) 'read matrix'
      call read_matrix_from_file(M, ai)
  
      write(*,*) 'generate random eigenvector guess and initialize eigenvalues to zero'
      allocate(result_v(M%sizem,ai%neig))
      allocate(result_e(ai%neig))          
      allocate(R       (M%sizem,ai%neig))
      iseed(1)=2 
      iseed(2)=17 
      iseed(3)=51 
      iseed(4)=101
      call zlarnv(1, iseed, M%sizem*ai%neig, result_v)        
      result_e(1:ai%neig)=ZERO

      write(*,*) 'call tracemin'       
      !t= timer()      
      call test_ztracemin(M,ai,result_v,result_e,R, 0,0,0,info)  
      !t2= timer()

      write(*,*) 'print results'
      !compute norm of residual of tracemin algorithm
      allocate(Rnrm(ai%neig)) 
      do i=1,ai%neig 
        Rnrm(i) = dznrm2(M%sizem,R(1,i),1)                   
      end do
      !compute norm of residual after refinement of real A*x=lambda*x eigenvalue problem
      allocate(iRnrm(ai%neig))
      call cshmv(M,result_v,R,ai%neig)
      do i=1,ai%neig
         call zaxpy(M%sizem,dcmplx(MONE*result_e(i),ZERO),result_v(1,i),1,R(1,i),1)
         iRnrm(i) = dznrm2(M%sizem,R(1,i),1)
      end do
      write(*,*) 'Smallest eigenpairs are:'
      do i=1,ai%neig
        write(*,*) 'eigenvalue', result_e(i)
        !write(*,*) 'eigenvector', result_v(:,i) 
      end do
      write(*,*) 'Norm of the residual: ', Rnrm
      write(*,*) 'Norm of the improved, real residual: ', iRnrm
      write(*,*) '# of iterations taken:', ai%maxi
      write(*,*) 'Time taken:', t2-t
  
      !deallocate unnecessary arrays
      deallocate(M%row)
      deallocate(M%column)
      deallocate(M%coefficient)
      deallocate(Rnrm)
      deallocate(R)
      deallocate(result_v)
      deallocate(result_e) 
end program sqtracemin
