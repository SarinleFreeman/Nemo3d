program pctracemin
      use tracemin_types
      use pctracemin_functions
      implicit none
 
      !----- Varibales Declaration -----
      !general variables
      type(zmatrix)                              :: M,Mp                  !matrix M in generalized eigenvalue problem I*x=lambda*M*x 
      type(alginfo)                              :: ai                    !structure containing different information for the algorithm 
      integer                                    :: iseed(4)              !variable used as a seed for random vector generation
      double complex,dimension(:,:), allocatable :: result_v, R, ig, ig_t, result_v_t       !resulting eigenvectors, residual and initial guess 
      double precision,dimension(:), allocatable :: result_e, Rnrm, iRnrm !resulting eigenvalues and norm of the residual  
      double precision                           :: t, t2, timer, dznrm2  !time variables, timer function, complex vector norm function
      integer                                    :: i, info               
      !MPI variables
      integer                                    :: mpierror, status(MPI_STATUS_SIZE)
      integer, parameter                         :: root = 0              !usually the processor with rank 0    
      integer                                    :: rank, np              !rank of the processor, number of processors  
      integer,         dimension(:), allocatable :: ne,disp               !number of elements and displacement on every processor (array used in MPI_scatterv) 
      type(matrix_splitting)                     :: MS                    !matrix splitting
      double precision                           :: alpha

      
      !write(*,*) rank, 'MPI initializations (rank is not set yet)'       
      call MPI_Init(mpierror)
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, mpierror)
      call MPI_Comm_size(MPI_COMM_WORLD, np, mpierror) 
            
      !write(*,*) rank, 'read matrix'
      if (rank == root) then
         call read_matrix_from_file(M, ai)                 !read the matrix from a file
         call row2column_major(M)                          !transform matrix into column major order 
         call compute_matrix_splitting(M,MS,np,1)          !compute matrix column splitting across processors (last parameter=1 indicates column splitting)           
      else                                          
         call allocate_matrix_splitting(MS,np)               
      end if
       
      !write(*,*) rank, 'Broadcast the control information'
      call MPI_Bcast_alginfo(ai,root,mpierror)
      call MPI_Bcast_matrix_splitting(MS,root,np,mpierror) 
      
      !write(*,*) rank, 'Broadcast the whole matrix'
      !!call MPI_Bcast(M%size, 1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)
      !!call MPI_Bcast(M%rsize,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)
      !!call MPI_Bcast(M%bsize,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)      
      !!if (rank /= root) then
      !!   allocate(M%row        (M%rsize))
      !!   allocate(M%column     (M%rsize))
      !!   allocate(M%coefficient(M%rsize))
      !!end if
      !!call MPI_Bcast(M%row,        M%rsize,MPI_INTEGER,       root,MPI_COMM_WORLD,mpierror)
      !!call MPI_Bcast(M%column,     M%rsize,MPI_INTEGER,       root,MPI_COMM_WORLD,mpierror)
      !!call MPI_Bcast(M%coefficient,M%rsize,MPI_DOUBLE_COMPLEX,root,MPI_COMM_WORLD,mpierror)

      !write(*,*) rank, 'Spread the matrix across processors'
      Mp%bsize=0                                
      Mp%sizem=M%sizem 
      call MPI_Bcast  (Mp%sizem,  1,MPI_INTEGER,                       root,MPI_COMM_WORLD,mpierror)
      call MPI_scatter(MS%nvector,1,MPI_INTEGER,Mp%sizen,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)
      call MPI_scatter(MS%ne,     1,MPI_INTEGER,Mp%rsize,1,MPI_INTEGER,root,MPI_COMM_WORLD,mpierror)      
      allocate(Mp%row        (Mp%rsize))
      allocate(Mp%column     (Mp%rsize))
      allocate(Mp%coefficient(Mp%rsize))
      allocate(disp(np))
      do i=1,np
         disp(i)=MS%isvector(i)-1
      end do      
      call MPI_scatterv(M%row,        MS%ne,disp,MPI_INTEGER       ,Mp%row        ,MS%ne(rank+1),MPI_INTEGER       ,root,MPI_COMM_WORLD,mpierror)
      call MPI_scatterv(M%column,     MS%ne,disp,MPI_INTEGER       ,Mp%column     ,MS%ne(rank+1),MPI_INTEGER       ,root,MPI_COMM_WORLD,mpierror)
      call MPI_scatterv(M%coefficient,MS%ne,disp,MPI_DOUBLE_COMPLEX,Mp%coefficient,MS%ne(rank+1),MPI_DOUBLE_COMPLEX,root,MPI_COMM_WORLD,mpierror)
      Mp%column(1:MS%ne(rank+1)) = Mp%column(1:MS%ne(rank+1)) - MS%svector(rank+1) + 1          !adjust the column index
      call column2row_major(Mp)		
      !write(*,*) rank, '-->Mp%sizem',       Mp%sizem
      !write(*,*) rank, '-->Mp%sizen',       Mp%sizen
      !write(*,*) rank, '-->Mp%rsize',       Mp%rsize
      !write(*,*) rank, '-->Mp%bsize',       Mp%bsize
      !write(*,*) rank, '-->Mp%row',         Mp%row
      !write(*,*) rank, '-->Mp%column',      Mp%column
      !write(*,*) rank, '-->Mp%coefficient', Mp%coefficient
      
      !write(*,*) rank,'create and scatter random initial guess'       
      allocate(ig  (Mp%sizem,ai%neig))   
      allocate(ig_t(ai%neig,Mp%sizem))   
      if (rank == root) then
         iseed(1)=2 
         iseed(2)=17 
         iseed(3)=51 
         iseed(4)=101
         call zlarnv(1, iseed, Mp%sizem*ai%neig, ig)        
      end if
      allocate(ne  (np))
      do i=1,np-1
         ne(i)=(Mp%sizem/np)*ai%neig
      end do
      ne(np)=(Mp%sizem-(Mp%sizem/np)*(np-1))*ai%neig
      disp(1)=0
      do i=2,np
         disp(i)=ne(i-1)
      end do
      allocate(result_v  (ne(rank+1)/ai%neig,ai%neig))
      allocate(result_v_t(ai%neig,ne(rank+1)/ai%neig))
      ig_t = transpose(ig)
      call MPI_scatterv(ig_t(1:Mp%sizem*ai%neig,1),ne,disp,MPI_DOUBLE_COMPLEX,result_v_t(1:ne(rank+1),1),ne(rank+1),MPI_DOUBLE_COMPLEX,root,MPI_COMM_WORLD,mpierror) 
      result_v=transpose(result_v_t)

      !write(*,*) rank, 'allocate space for eigenvalues and residual and set them to zero'      
      allocate(result_e(ai%neig))          
      allocate(R(ne(rank+1)/ai%neig,ai%neig))
      result_e(1:ai%neig)                     = ZERO
      R       (1:ne(rank+1)/ai%neig,1:ai%neig)=ZZERO
      
      !write(*,*) rank,'call generalized complex tracemin' 
      t= MPI_wtime()      
      call test_zgtracemin(Mp,ai,result_v,result_e,R, ne,rank,np,MPI_COMM_WORLD,info)  
      t2=MPI_wtime()

      !write(*,*) rank,'Print results'
      allocate(Rnrm(ai%neig)) 
      do i=1,ai%neig 
         alpha  =dznrm2(ne(rank+1)/ai%neig,R(1,i),1)
         call MPI_Allreduce(alpha, Rnrm(i), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)                  
      end do


      !compute norm of residual after refinement of real A*x=lambda*x eigenvalue problem
      allocate(iRnrm(ai%neig))
      !call cshmv(M,result_v,R,ai%neig)
      call pcfmv(Mp,result_v,R,ai%neig,ne(rank+1)/ai%neig,ne,rank,np)
      do i=1,ai%neig
         call zaxpy(ne(rank+1)/ai%neig,dcmplx(MONE*result_e(i),ZERO),result_v(1,i),1,R(1,i),1)
         alpha = dznrm2(ne(rank+1)/ai%neig,R(1,i),1)
         call MPI_Allreduce(alpha, iRnrm(i), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)                  
      end do




      !write(*,*) 'Smallest eigenpairs are:'
      if (rank == root) then
         do i=1,ai%neig
            write(*,*) rank,'eigenvalue', result_e(i)
            !write(*,*) rank,'eigenvector', result_v(:,i) 
         end do
         write(*,*) 'Norm of the residual: ', Rnrm
         write(*,*) 'Norm of the real residual: ', iRnrm
         write(*,*) '# of iterations taken:', ai%maxi
         write(*,*) 'Time taken:', t2-t
      end if

      !write(*,*) rank,'deallocate unnecessary arrays'
      if (rank == root) then
         deallocate(M%row)
         deallocate(M%column)
         deallocate(M%coefficient)
      end if
      deallocate(disp)
      deallocate(ne)
      deallocate(Rnrm)
      deallocate(R)
      deallocate(result_v)
      deallocate(result_e) 
      deallocate(result_v_t)
      
      !write(*,*) rank,'finalize MPI'  
      call MPI_Finalize(mpierror)  

end program pctracemin
