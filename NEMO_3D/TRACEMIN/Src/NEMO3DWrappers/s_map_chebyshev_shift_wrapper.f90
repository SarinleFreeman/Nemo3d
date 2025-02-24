subroutine map_chebyshev_shift_wrapper(n,nloc,nrhs,Vn,Vc,Vch,Vp,S,k,gamma,a,b,ncall,info)
  use tracemin_types
  use sctracemin_functions
  integer,                              intent(in)    :: n         !order of the matrix 
  integer,                              intent(in)    :: nloc      !number of elements on the current processor
  integer,                              intent(in)    :: nrhs      !number of vectors
  double complex, dimension(nloc,nrhs), intent(inout) :: Vn        !B times set of p current vectors (input: vector B*U)
  double complex, dimension(nloc,nrhs), intent(inout) :: Vc        !set of p current vectors         (input: vector U by which we are multiplying,output: result of multiply) 
  double complex, dimension(nloc,nrhs), intent(inout) :: Vch       !set of p temporary vectors
  double complex, dimension(nloc,nrhs), intent(inout) :: Vp        !temporary array for set of p current vectors
  double complex, dimension(nloc,nrhs), intent(inout) :: S         !set of original vectors used to perform the shift (input: copy of vector U, never modified in the routine)
  integer,                              intent(in)    :: k         !degree of chebyshev polynomial 
  double precision,                     intent(in)    :: gamma     !parameter used in the mapping Q(x)=1-gamma*(x-a)*(x-b)
  double precision,                     intent(in)    :: a         !left bound of the interval containing p desired eigenvalues 
  double precision,                     intent(in)    :: b         !right bound of the interval containing p desired eigenvalues
  integer,                              intent(in)    :: ncall     !number of times routine was called (every iteration takes two calls, starts with 1)
  integer,                              intent(inout) :: info      !on exit 0  - done
                                                             !        1,2- compute Vn=B*Vn  
                                                             !        <0 - error 
  !write(*,*) 'IN THE WRAPPER3'
  !write(*,*) 'n,nloc,nrhs',n,nloc,nrhs
  !write(*,*) 'k,gamma,a,b,ncall,info',k,gamma,a,b,ncall,info
  !write(*,*) 'Vn(1,1)', Vn(1,1)
  !write(*,*) 'Vc(1,1)', Vc(1,1)
  !write(*,*) 'IN THE WRAPPER AFTER THE PRINT'

  call map_chebyshev_shift(n,nloc,nrhs,Vn,Vc,Vch,Vp,S,k,gamma,a,b,ncall,info)
end subroutine map_chebyshev_shift_wrapper
