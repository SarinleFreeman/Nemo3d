module tracemin_types
      !use mpi
      implicit none

      double precision, parameter :: ZERO    = 0.0D+0 
      double precision, parameter :: ONE     = 1.0D+0 
      double precision, parameter :: TWO     = 2.0D+0 
      double precision, parameter :: MONE    =-1.0D+0  
      double precision, parameter :: MTWO    =-2.0D+0  
      double precision, parameter :: EPSILON = 1.0D-16 

      double complex,   parameter :: ZZERO   = (0.0D+0, 0.0D+0)
      double complex,   parameter :: ZONE    = (1.0D+0, 0.0D+0)
      double complex,   parameter :: ZTWO    = (2.0D+0, 0.0D+0)
      double complex,   parameter :: ZMONE   =(-1.0D+0, 0.0D+0)
      double complex,   parameter :: ZMTWO   =(-2.0D+0, 0.0D+0) 

      integer :: TRUE = 1
      integer :: FALSE= 0 

      !this structure is used to hold a sparse matrix
      type dmatrix
        double precision,dimension(:),allocatable :: coefficient
        integer,         dimension(:),allocatable :: row, column 
        integer                                   :: size,rsize, bsize !size of the matrix, representation size, size of the block (not necessarily used)
      end type dmatrix

      type zmatrix
        double complex,  dimension(:),allocatable :: coefficient
        integer,         dimension(:),allocatable :: row, column
        integer                                   :: sizem,sizen, rsize, bsize            
      end type zmatrix

      !splitting of the matrix 
      type matrix_splitting
        integer, dimension(:), allocatable :: nvector  !the number of vectors=[rows|columns] belonging to each processor
        integer, dimension(:), allocatable :: svector  !the starting vector on each processor
        integer, dimension(:), allocatable :: evector  !the ending   vector on each processor 
        integer, dimension(:), allocatable :: isvector !the index into array [M%row|M%column] where the first element of set of vectors belonging to each processor is located
        integer, dimension(:), allocatable :: ievector !the index into array [M%row|M%column] where the last  element of set of vectors belonging to each processor is located
        integer, dimension(:), allocatable :: ne       !the # of elements on every processor 
      end type matrix_splitting

      type alginfo
         integer          :: neig !number of eigenvalues you are looking for
         double precision :: a    !the interval to which neig eigenvalues you are looking for belong to [a,b] 
         double precision :: b 
         double precision :: c    !the interval to which all eigenvalues belong to [c,d]
         double precision :: d         
         integer          :: k    !degree of chebyshev polynomial to be used      
         integer          :: maxi !maximum number of iterations
         double precision :: tol  !tolerance
      end type alginfo

end module tracemin_types
