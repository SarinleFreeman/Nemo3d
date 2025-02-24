subroutine pcextract_eigenpairs_wrapper(Y,E,neig,nloc,H,Hloc,ZWork,RWork,T,rank,comm,info)
  use tracemin_types
  use sctracemin_functions
  double complex,  dimension(nloc,neig),      intent(inout):: Y     !eigenvectors   
  double precision,dimension(neig),           intent(inout):: E     !eigenvalues
  integer,                                    intent(in)   :: neig  !number of eiegnpairs    
  integer,                                    intent(in)   :: nloc  !number of vector elements on each processor
  double complex,  dimension(neig,neig),      intent(inout):: H     !temporary array of size neig by neig
  double complex,  dimension(neig,neig),      intent(inout):: Hloc  !temporary array of size neig by neig
  double complex,  dimension(max(1,2*neig-1)),intent(inout):: Zwork !temporary array of size 2*neig-1 by 1
  double precision,dimension(max(1,3*neig-2)),intent(inout):: RWork !temporary array of size 3*neig-2 by 1
  double complex,  dimension(nloc,neig),      intent(inout):: T     !on second call T=A*Y
  integer,                                    intent(in)   :: rank  !rank of the current processor
  integer,                                    intent(in)   :: comm  !MPI communicator
  integer,                                    intent(inout):: info  !on entry first is one then is two, on exit reports errors 
      
  call extract_eigenpairs(Y,E,neig,nloc,H,Hloc,ZWork,RWork,T,rank,comm,info)
end subroutine pcextract_eigenpairs_wrapper
