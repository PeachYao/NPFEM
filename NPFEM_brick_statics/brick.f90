﻿! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
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
    
    
subroutine brick
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the brick element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    use globals
    use memallocate
    
    implicit none
    integer :: NumberOfElements, NumberOfMaterials, ElementGroupSize
    integer :: N(8) !Pointers
    
    NPAR(5) = 8
    NPAR(4) = 0
    
    NumberOfElements = NPAR(2)
    NumberOfMaterials = NPAR(3)
    
    !Preallocate Memory
    
!pointer lists
! Calculate the pointer to the arrays in the element group data
! N101: E(NumberOfMaterials)
! N102: v(NumberOfMaterials)
! N103: Density(NumberOfMaterials)
! N104: LM(3*NPAR(5),NumberOfElements)
! N105: PositionData(3*NPAR(5),NumberOfElements)
! N106: MaterialData(NumberOfElements)
! N107: Node Connection Matrix
    

    N(1) = 0
    N(2) = N(1)+NumberOfMaterials*ITWO
    N(3) = N(2)+NumberOfMaterials*ITWO
    
    if (DYNANALYSIS) then
        N(4) = N(3)+NumberOfMaterials*ITWO
    else
        N(4) = N(3)
    end if
    
    N(5) = N(4) + 3*NPAR(5)*NumberOfElements
    N(6) = N(5) + 3*NPAR(5)*NumberOfElements*ITWO
    N(7) = N(6) + NumberOfElements
    N(8) = N(7) + NPAR(5)*NPAR(2)
    
    MIDEST = N(8)
    
    if (IND .EQ. 1) then
        call MemAlloc(11,"ELEGP",MIDEST,1)
    end if
    NFIRST = NP(11)
    N(:) = N(:) + NFIRST
    NLAST  = N(8)
    
    call C3D8 (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
                  A(N(1)),A(N(2)),A(N(3)),A(N(4)),A(N(5)),A(N(6)),A(N(7)))
    
    !Reuse DA(NP(4)) at Solution Phase 3 as displacement U
    return
end subroutine brick

subroutine C3D8 (ID,X,Y,Z,U,MHT,E, PoissonRatio, Density, LM, PositionData, MaterialData, Node)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   C3D8 element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    use globals
    use MemAllocate
    use MathKernel
    USE MKL_PARDISO
    
    implicit none
    integer ::  ElementShapeNodes, NumberOfMaterials, NumberOfElements
    integer ::  QuadratureOrder     = 2                     ! NPAR(6) -- Element Load Nodes
    INTEGER ::  ID(3,NUMNP), MHT(NEQ), MaterialData(NPAR(2))
    integer ::  MaterialType, MaterialComp, ND, L, N, i, j, LM(24,NPAR(2)), ElementType, ind0, iprint, k, &
                ind1, ind2, Node(NPAR(2),NPAR(5))
    real(8) ::  X(NUMNP), Y(NUMNP), Z(NUMNP), U(NEQ), &
                DetJ(2,2,2), E(NPAR(3)), PoissonRatio(NPAR(3)), ElementDisp(24)
    real(8) ::  BMatrix(6, 3*NPAR(5)), PositionData(3*NPAR(5), NPAR(2)), DMatrix(6,6), &
                Transformed(3), W(2), Weight(2,2,2), GaussianPts(2), GaussianCollection(3, NPAR(2)*2**3), &
                StressCollection(6,NPAR(2)*2**3), M(3*NPAR(5),3*NPAR(5)), Rho
    real(8) ::  Young, v, S(3*NPAR(5),3*NPAR(5)), GaussianPtsPosit(3,2**3), Strain(6,2**3), Stress(6,2**3), &
                Density(NPAR(3)), NMatrix(3,3*NPAR(5)), NormalVec(3), Point(3*NPAR(5),3*NPAR(5))
    REAL(8) ::  I24(24,24),Q24(24,24),EulerAngles(3),R24(24,24),FKe(24)
                
    ElementType         = NPAR(1)
    NumberOfElements    = NPAR(2)
    NumberOfMaterials   = NPAR(3)
    ElementShapeNodes   = NPAR(5)                           !NPAR(5)=8
    ND = 24
    
    SELECT CASE (IND)
    CASE(1)
        WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I10,/,   &
                   '     EQ.1, 3D HEXAHEDRAL',/,       &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I10,/)") ElementType, NumberOfElements
        IF (NumberOfMaterials.EQ.0) NumberOfMaterials=1
        WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I10,/)") NumberOfMaterials
        
        WRITE(IOUT,'(A)') '  SET       YOUNG''S MODULUS    POISSON RATIO     DENSITY (RHO)'
        WRITE(IOUT,'(A)') ' NUMBER         E (Pa)               ν              ρ (kg/m^3)'
        WRITE(IOUT,'(A)') '---------------------------------------------------------------'

        if (DYNANALYSIS) then
            DO I=1,NumberOfMaterials
            READ (IIN,'(I10,3F10.0)') N,E(N), PoissonRatio(N), Density(N)      ! Read material information for 3D Homogeneous
            WRITE(IOUT,'(I8,4X,E13.5,4X,E13.6,4X,E13.6)') N, E(N), PoissonRatio(N), Density(N)
            END DO
        else
            DO I=1,NumberOfMaterials
                READ (IIN,'(I10,2F10.0)') N,E(N), PoissonRatio(N)      ! Read material information for 3D Homogeneous
                WRITE (IOUT,"(I8,4X,E13.5,4X,E13.6)") N,E(N), PoissonRatio(N)
            END DO
        end if
        WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT        |------------------------- NODES -------------------------|       MATERIAL',/,   &
                      ' NUMBER-N        1       2       3       4       5       6       7       8       SET NUMBER')")
        
        N=0

            
        DO WHILE (N .NE. NumberOfElements)
            READ (IIN,'(11I5)') N,Node(N,1:ElementShapeNodes),MaterialType          ! Read in element information
    !       Save element information
            PositionData(1:ElementShapeNodes*3:3,N)=X(Node(N,:))        ! Coordinates of the element's nodes
            PositionData(2:ElementShapeNodes*3:3,N)=Y(Node(N,:))
            PositionData(3:ElementShapeNodes*3:3,N)=Z(Node(N,:))

            
            MaterialData(N) = MaterialType                                ! Material type

            DO L=1,ND
                LM(L,N)=0
            END DO
            DO L=1,3
                LM(L:ElementShapeNodes*3+L-2:3,N) = ID(L,Node(N,:))       ! Connectivity matrix
            END DO
            
    !       Update column heights and bandwidth        
            if (.NOT. PARDISODOOR) CALL COLHT (MHT,ND,LM(:,N))
            WRITE (IOUT,"(I7,5X,7(I7,1X),I7,4X,I10)") N,Node(N,1:ElementShapeNodes),MaterialType       
        END DO
        return
    
    CASE (2)
        ! NPFEM 的中间量定义
        ! 初始化为0
        Q24 = 0.0D0 
        I24 = 0.0D0

        ! 构造I和Q
        DO j = 0, 7
            DO i = 1, 3
            Q24(j*3 + i, i) = 1.0D0
            END DO
        END DO
        DO i = 1, 24
           I24(i,i) = 1.0D0
        END DO
        
        ! 设置欧拉角（单位：弧度）
        EulerAngles = (/ 0.0D0, 0.0D0, 0.0D0 /)

        ! 调用子程序生成旋转矩阵
        CALL transmatrix(EulerAngles, R24)
  
        
        CALL GaussianMask(GaussianPts, W, QuadratureOrder)
        do i = 1, QuadratureOrder
            do j = 1, QuadratureOrder
                do k = 1, QuadratureOrder
                    Weight(i,j,k)  = W(i)*W(j)*W(k)
                end do
            end do
        end do
        
        MaterialComp = -1
        do N = 1, NumberOfElements
            S(:,:) = 0
            MaterialType = MaterialData(N)
            if (MaterialType .NE. MaterialComp) then
                Young        = E(MaterialType)
                v            = PoissonRatio(MaterialType)
                call GetDMat(Young, v, DMatrix)
                IF (DYNANALYSIS .EQV. .TRUE.) Rho = Density(MaterialType)
                MaterialComp = MaterialType
            end if
            
            !DetJ = reshape(Jacobian((n-1)*QuadratureOrder**3+1 : n*QuadratureOrder**3), &
            !               (/QuadratureOrder, QuadratureOrder, QuadratureOrder/))
            DetJ(:, :, :)=0
            do i = 1, QuadratureOrder
                do j = 1, QuadratureOrder
                    do k = 1, QuadratureOrder
                        Transformed   = (/GaussianPts(i), GaussianPts(j), GaussianPts(k)/)
                        CALL HexB(BMatrix, DetJ(i,j,k), ElementShapeNodes, Transformed, &
                                  (/PositionData(1:ElementShapeNodes*3:3,N), &
                                    PositionData(2:ElementShapeNodes*3:3,N), &
                                    PositionData(3:ElementShapeNodes*3:3,N)/))
                        Point   = matmul(matmul(transpose(BMatrix),DMatrix),BMatrix)
                        S = S + (Weight(i,j,k)*DetJ(i,j,k))*Point
                        IF (DYNANALYSIS .EQV. .TRUE.) then
                            call HexN (NMatrix, ElementShapeNodes, Transformed)             !Initialize N Matrix for mass assembly
                            Point = Rho*matmul(transpose(NMatrix),NMatrix)
                            M = M + (Weight(i,j,k)*DetJ(i,j,k))*Point
                        END IF
                    end do
                end do
            end do
            
            ! Important part of NPFEM
            S = matmul(transpose(R24), matmul(transpose(I24 - Q24), matmul(S, matmul(I24 - Q24, R24))))
            FKe = matmul(transpose(R24), matmul(transpose(I24 - Q24), matmul(S, matmul(R24, PositionData(:,N)))))
            !! check the stiffness matrix
            !if (N .EQ. 1) then
            !    DO i = 1, 24
            !        WRITE(IOUT, '(24(1X,E12.5))') FKe(i)
            !    END DO
            !end if
            
            !write(*,*) "S",S
            
            !if(pardisodoor) then
            !    call pardiso_addban(DA(NP(3)),IA(NP(2)),IA(NP(5)),S,LM(:,N),ND)
            !    if (DYNANALYSIS) CALL pardiso_addban(DA(NP(10)),IA(NP(9)), IA(NP(8)),M,LM(:,N),ND)
            !else
                CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(:,N),ND)                              !Assemble Stiffness Matrix
                CALL ADDRHS (DA(NP(6)),FKe,LM(:,N),ND)                                      !Assemble equivalent nodal force vector
                IF (DYNANALYSIS) CALL ADDBAN (DA(NP(10)),IA(NP(2)),M,LM(:,N),ND)            !Assemble Mass Matrix
            !end if
          
        end do
        
    CASE (3)
        
        IPRINT=0
        call GaussianMask(GaussianPts, W, QuadratureOrder)
        DO N=1,NumberOfElements
            IPRINT=IPRINT + 1
            IF (IPRINT.GT.50) IPRINT=1
            IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                           'E L E M E N T  G R O U P',I4,//,   &
            '  ELEMENT',13X, 'COORDINSTES',19X, 'Sigma_xx',9X, 'Sigma_yy',9X,'Sigma_zz',9X, &
            'Sigma_xy',9X,'Sigma_yz',9X,'Sigma_xz', /, &
            '  NUMBER', 8X,'X',10X,'Y',10X,'Z')") NG
            MaterialType = MaterialData(N)
            Young        = E(MaterialType)
            v            = PoissonRatio(MaterialType)
            call GetDMat(Young, v, DMatrix)
            
            ElementDisp(:) = 0
            
            do i = 1,ND
                if (LM(i,N) .NE. 0) ElementDisp(i)  =  U(LM(i,N))
            end do
            
            ind0 = 1
            do i = 1,QuadratureOrder
                do j = 1,QuadratureOrder
                    do k = 1,QuadratureOrder
                        Transformed   = (/GaussianPts(i), GaussianPts(j), GaussianPts(k)/)
                        call HexN (NMatrix, ElementShapeNodes, Transformed)
                        GaussianPtsPosit(:,ind0) = matmul(reshape(PositionData(:,N), (/3,ElementShapeNodes/)), &
                                                         NMatrix(1, 1:3*ElementShapeNodes:3))
                        CALL HexB (BMatrix, DetJ(i,j,k), ElementShapeNodes, Transformed, &
                                  (/PositionData(1:ElementShapeNodes*3-1:3,N), &
                                    PositionData(2:ElementShapeNodes*3  :3,N), &
                                    PositionData(3:ElementShapeNodes*3+1:3,N)/))
                        Strain(:,ind0) = matmul(BMatrix, ElementDisp)             
                        Stress(:,ind0) = matmul(DMatrix, Strain(:,ind0))
                        
                        ind0 = ind0 + 1
                    end do
                end do
            end do
            
            write (IOUT,"(I6,3(3X, F10.4),6(4X, E13.6),/,7(6X, 3(3X, F10.4),6(4X, E13.6),/))") &
                               N, (GaussianPtsPosit(:,I), Stress(:,I), I=1,QuadratureOrder**3)
            ind1 = (N-1)*QuadratureOrder**3+1
            ind2 = N*QuadratureOrder**3
            GaussianCollection (:,ind1:ind2) = GaussianPtsPosit
            StressCollection (:,ind1:ind2) = Stress
        END DO

        !call PostProcessor(ElementType, 3, PositionData, &
        !                   Node, QuadratureOrder**3, GaussianCollection, StressCollection, U)
                           
                
    END SELECT

end subroutine C3D8


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              Calculate Shape Function Matrix             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HexN (NMatrix, ElementShapeNodes, Transformed)
    implicit none
    integer ::  ElementShapeNodes, i
    real(8) ::  NMatrix(3, 3*ElementShapeNodes), Transformed(3), N(ElementShapeNodes)
    
    select case (ElementShapeNodes)
    case (8)
        N(1) = (1-Transformed(1))*(1-Transformed(2))*(1-Transformed(3))/8
        N(2) = (1+Transformed(1))*(1-Transformed(2))*(1-Transformed(3))/8
        N(3) = (1+Transformed(1))*(1+Transformed(2))*(1-Transformed(3))/8
        N(4) = (1-Transformed(1))*(1+Transformed(2))*(1-Transformed(3))/8
        N(5) = (1-Transformed(1))*(1-Transformed(2))*(1+Transformed(3))/8
        N(6) = (1+Transformed(1))*(1-Transformed(2))*(1+Transformed(3))/8
        N(7) = (1+Transformed(1))*(1+Transformed(2))*(1+Transformed(3))/8
        N(8) = (1-Transformed(1))*(1+Transformed(2))*(1+Transformed(3))/8
        NMatrix = reshape((/((/N(i), 0D0, 0D0, 0D0, N(i), 0D0, 0D0, 0D0, N(i)/),i=1,ElementShapeNodes)/), &
                          shape(NMatrix))
    end select
end subroutine HexN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       Calculate Gradiient of Shape Function Matrix       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HexB (BMatrix, DetJ, ElementShapeNodes, Transformed, Original)
    USE MathKernel
    implicit none
    integer, intent(in) ::  ElementShapeNodes
    real(8) ::  BMatrix(6, 3*ElementShapeNodes), DetJ, Transformed(3), Original(ElementShapeNodes, 3)
    real(8) ::  GradN(3,ElementShapeNodes), J(3,3), InvMatJ(3,3), DerivN(3,ElementShapeNodes)
    real(8) ::  Bx(ElementShapeNodes), By(ElementShapeNodes), Bz(ElementShapeNodes), xi, eta, zta
    integer ::  i
    logical ::  OK_Flag
    
    select case (ElementShapeNodes)
    case (8)
        xi      = Transformed(1)
        eta     = Transformed(2)
        zta     = Transformed(3)
        GradN   = 0.125*reshape((/-(1-eta)*(1-zta), -(1-xi)*(1-zta), -(1-xi)*(1-eta), &
                                   (1-eta)*(1-zta), -(1+xi)*(1-zta), -(1+xi)*(1-eta), &
                                   (1+eta)*(1-zta),  (1+xi)*(1-zta), -(1+xi)*(1+eta), &
                                  -(1+eta)*(1-zta),  (1-xi)*(1-zta), -(1-xi)*(1+eta), &
                                  -(1-eta)*(1+zta), -(1-xi)*(1+zta),  (1-xi)*(1-eta), &
                                   (1-eta)*(1+zta), -(1+xi)*(1+zta),  (1+xi)*(1-eta), &
                                   (1+eta)*(1+zta),  (1+xi)*(1+zta),  (1+xi)*(1+eta), &
                                  -(1+eta)*(1+zta),  (1-xi)*(1+zta),  (1-xi)*(1+eta)/), &
                                shape(GradN))
        J       = matmul(GradN, Original)
        DetJ    = Det(J, 3)
        call InvMat3(J, InvMatJ, OK_Flag)
        if (OK_Flag .EQV. .FALSE.) STOP "***ERROR*** Derivative Of Shape Function Is SINGULAR"
        DerivN  = matmul(InvMatJ,GradN)
        Bx(:)   = DerivN(1,:)
        By(:)   = DerivN(2,:)
        Bz(:)   = DerivN(3,:)
        BMatrix = reshape((/((/Bx(i), 0D0, 0D0, By(i), 0D0, Bz(i), &
                             0D0, By(i), 0D0, Bx(i), Bz(i), 0D0, &
                             0D0, 0D0, Bz(i), 0D0, By(i), Bx(i)/), &
                             i = 1, ElementShapeNodes)/), Shape(BMatrix))
    end select
end subroutine HexB