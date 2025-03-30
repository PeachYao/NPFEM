! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                       .
! .                             N P F E M                                 .
! .                                                                       .
! .     Nodal Position Finite Element Method PROGRAM IN FORTRAN 90        .
! .     Adapted from STAP (KJ Bath) and STAP90 (Xiong Zhang)              .
! .                                                                       .
! .     Fuzhen Yao, (2025)                                                .
! .     Department of Mechanical Engineering, York University             .
! .                                                                       .
! . . . . . . . . . . . . . .  . . .  . . . . . . . . . . . . . . . . . . .

! .  Define global variables

module GLOBALS

   integer, parameter :: IELMNT=1	! Unit storing element data 
   integer, parameter :: ILOAD=2	! Unit storing load vectors
   integer, parameter :: INODE=3	! Unit storing node vectors
   integer, parameter :: IIN=5		! Unit used for input
   integer, parameter :: IOUT=6		! Unit used for output

   integer :: NUMNP		! Total number of nodal points
						! = 0 : Program stop
   integer :: NEQ		! Number of equations
   integer :: NWK		! Number of matrix elements
   integer :: MK		! Maximum half bandwidth

   integer :: IND		! Solution phase indicator
						!   1 - Read and generate element information
						!   2 - Assemble structure stiffness matrix
						!   3 - Stress calculations
   integer :: NPAR(10)	! Element group control data
						!   NPAR(1) - Element type
						!             1 : Brick element
						!             2 : Beam element   
						!   NPAR(2) - Number of elements
						!   NPAR(3) - Number of different sets of material and 
						!             cross-sectional  constants
   						!   NPAR(4) - Consider Gravity Or Not
						!   NPAR(5) - Numder of Nodes
   integer :: NUMEG		! Total number of element groups, > 0

   integer :: MODEX		! Solution mode: 0 - data check only;  1 -  execution                                   

   real :: TIM(5)		! Timing information
   character*80 :: HED	! Master heading information for use in labeling the output

   integer :: NFIRST
   integer :: NLAST
   integer :: NUMEST
   integer :: MIDEST
   integer :: MAXEST

   integer :: NG
   
   logical :: PARDISODOOR =  .FALSE.
   logical :: DYNANALYSIS =  .TRUE.

end module GLOBALS
