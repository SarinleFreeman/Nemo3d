subroutine modified_conjugate_gradient_wrapper(n,neig,X,Y,Q,G,Gloc,P,R,gamma_old,rtol,ncall,maxi,tol,gtol,nloc,rank,np,comm,info)
  use tracemin_types
  use sqtracemin_functions
  integer,                                 intent(in)    :: n         !order of the matrix
  integer,                                 intent(in)    :: neig      !number of right hand sides
  double complex,  dimension(nloc,neig),   intent(inout) :: X         !array of size n    by neig, approximation to the solution 
  double complex,  dimension(nloc,neig),   intent(in)    :: Y         !array of size n    by neig, set of vetors used to perform (I-Y*Y')*Z multiplications             
  double complex,  dimension(nloc,neig),   intent(inout) :: Q         !array of size n    by neig, Q=A*Y on entrance and Q=A*P during iterations of CG                        
  double complex,  dimension(neig,neig),   intent(inout) :: G         !array of size neig by neig, temporary matrix                         
  double complex,  dimension(neig,neig),   intent(inout) :: Gloc      !array of size neig by neig, temporary matrix                         
  double complex,  dimension(nloc,neig),   intent(inout) :: P         !array of size n    by neig, search directions of CG                         
  double complex,  dimension(nloc,neig),   intent(inout) :: R         !array of size n    by neig, residual of CG             
  double complex,  dimension(neig),        intent(inout) :: gamma_old !array of size neig, contains gamma values from previous iteration of CG
  double precision,dimension(neig),        intent(inout) :: rtol      !array of size neig,         relative tolerance computed inside this routine
  integer,                                 intent(inout) :: ncall     !number of times this routine was called
  integer,                                 intent(inout) :: maxi      !maximum number of CG iterations specified by the user   
  double precision,                        intent(in)    :: tol       !tolerance specified by the user
  double precision,                        intent(in)    :: gtol      !global tolerance for the tracemin specified by the user
  integer,                                 intent(in)    :: nloc   !number of elements of a vector on the current processor (in sequential version not used)
  integer,                                 intent(in)    :: rank   !rank of the processor                                   (in sequential version not used)
  integer,                                 intent(in)    :: np     !number of processors                                    (in sequential version not used)
  integer,                                 intent(in)    :: comm   !MPI communicator                                        (in sequential version not used)
  integer,                                 intent(inout) :: info   !on exit  1 - done
                                                                   !         2 - perform matrix vector multiply Q=A*P
                                                                   !         3 - maximum number of iterations achieved 
                                                                   !        <0 - error


  call modified_conjugate_gradient(n,neig,X,Y,Q,G,Gloc,P,R,gamma_old,rtol,ncall,maxi,tol,gtol,nloc,rank,np,comm,info)

end subroutine modified_conjugate_gradient_wrapper
