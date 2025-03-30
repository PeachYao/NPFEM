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
subroutine prepare_MassMatrix
    use globals
    use memallocate
    implicit none

    CALL MEMALLOC(10,"M    ",NWK,ITWO)
    if (PARDISODOOR) then
        call memalloc(9,"MrInd",NEQ+1,1)
        call memalloc(8,"Mcolm",NWK,1)
        IA(NP(8):NP(8)+NWK) = IA(NP(5):NP(5)+NWK)
        IA(NP(9):NP(9)+NEQ+1) = IA(NP(2):NP(2)+NEQ+1)
    end if
end subroutine prepare_MassMatrix 
    
    
SUBROUTINE transmatrix(EulerAngles, R24)
  !--------------------------------------------------------------
  ! Build a 24x24 block-diagonal rotation matrix for 8-node brick element
  ! Input:
  !   EulerAngles(3) - Rotation angles [alpha, beta, gamma] in radians
  !                    (rotation about Z, then Y, then X axes)
  ! Output:
  !   R24(24,24)     - Block-diagonal rotation matrix (8 blocks of 3x3)
  !--------------------------------------------------------------

  IMPLICIT NONE
  REAL(8), INTENT(IN)  :: EulerAngles(3)
  REAL(8), INTENT(OUT) :: R24(24,24)

  REAL(8) :: alpha, beta, gamma
  REAL(8) :: R3(3,3)
  INTEGER :: i, j, blkIdx

  ! Unpack Euler angles
  alpha = EulerAngles(1)  ! rotation about Z
  beta  = EulerAngles(2)  ! rotation about Y
  gamma = EulerAngles(3)  ! rotation about X

  ! Clear the full rotation matrix
  R24 = 0.0D0

  !--------------------------------------------------------------
  ! Compute the basic 3x3 rotation matrix using ZYX convention
  !--------------------------------------------------------------

  R3(1,1) = COS(alpha)*COS(beta)
  R3(1,2) = COS(alpha)*SIN(beta)*SIN(gamma) - SIN(alpha)*COS(gamma)
  R3(1,3) = COS(alpha)*SIN(beta)*COS(gamma) + SIN(alpha)*SIN(gamma)

  R3(2,1) = SIN(alpha)*COS(beta)
  R3(2,2) = SIN(alpha)*SIN(beta)*SIN(gamma) + COS(alpha)*COS(gamma)
  R3(2,3) = SIN(alpha)*SIN(beta)*COS(gamma) - COS(alpha)*SIN(gamma)

  R3(3,1) = -SIN(beta)
  R3(3,2) = COS(beta)*SIN(gamma)
  R3(3,3) = COS(beta)*COS(gamma)

  !--------------------------------------------------------------
  ! Copy R3 into each of the 8 block-diagonal positions in R24
  !--------------------------------------------------------------

  DO blkIdx = 0, 7
     DO i = 1, 3
        DO j = 1, 3
           R24(blkIdx*3 + i, blkIdx*3 + j) = R3(i,j)
        END DO
     END DO
  END DO

END SUBROUTINE transmatrix
