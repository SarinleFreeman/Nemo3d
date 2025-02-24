subroutine ztracemin_wrapper(n,neig,Y,E,Q,H,Hloc,ZWork,RWork,IWork,R,tol,maxi,ncall,nloc,rank,np,comm,info)
    use tracemin_types
    use sqtracemin_functions
    integer,                                    intent(in)   :: n      !order of the matrix 
    integer,                                    intent(in)   :: neig   !number of eigenvalues we want to know
    double complex,  dimension(nloc,neig),      intent(inout):: Y      !neig eigenvectros of our matrix      
    double precision,dimension(neig),           intent(inout):: E      !neig eigenvalues of our matrix 
    double complex,  dimension(nloc,neig),      intent(inout):: Q      !temporary variable, when matrix-vector multiply is requested it is Q=M*Q where M= A or B
    double complex,  dimension(neig,neig),      intent(inout):: H      !temporary variable, contains Y'*A*Y
    double complex,  dimension(neig,neig),      intent(inout):: Hloc   !temporary variable, contains Y'*A*Y
    double complex,  dimension(max(1,2*neig-1)),intent(in)   :: Zwork  !temporary array of size 2*neig-1 by 1
    double precision,dimension(max(1,3*neig-2)),intent(in)   :: RWork  !temporary array of size 3*neig-2 by 1
    integer,         dimension(neig),           intent(in)   :: IWork  !temporary array of size neig 
    double complex,  dimension(nloc,neig),      intent(out)  :: R      !vector containing the residual = Y-lambda*B*Y
    double precision,                           intent(in)   :: tol    !required tolerance    
    integer,                                    intent(in)   :: maxi   !maximum number of iterations
    integer,                                    intent(in)   :: ncall  !number of calls to the routine so far (every iteration takes two calls, starts with 1) 
    integer,                                    intent(in)   :: nloc   !number of elements of a vector on the current processor
    integer,                                    intent(in)   :: rank   !rank of the processor 
    integer,                                    intent(in)   :: np     !number of processors  
    integer,                                    intent(in)   :: comm   !MPI communicator
    integer,                                    intent(inout):: info   !on exit  1 - done
                                                                       !         2 - perform solve of (I-P)A(I-P)D=(I-P)AY to find updated Y=Y-(I-P)D where P=YY'
                                                                       !         3 - perform matrix vector multiply Q=A*Y
                                                                       !         4 - maximum number of iterations reached
                                                                       !        <0 - error    


    call ztracemin(n,neig,Y,E,Q,H,Hloc,ZWork,RWork,IWork,R,tol,maxi,ncall,nloc,rank,np,comm,info)

end subroutine ztracemin_wrapper
