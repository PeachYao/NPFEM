subroutine pardiso_crop(A, rowIndex, columns)
    use GLOBALS, only : NWK, NEQ
    real(8)    :: A(NWK)                  ! 稀疏矩阵的值数组
    integer    :: rowIndex(NEQ+1)         ! 行指针（CSR 格式）
    integer    :: columns(NWK)           ! 每个非零元素对应的列索引
    integer    :: i, j, k, temp           ! 循环变量与计数器

    ! 初始化
    j = 0           ! 新数组的当前位置指针
    k = 2           ! 当前正在处理的行号，从第2行开始
    temp = 0        ! 当前行的有效非零元素计数器

    do i = 1, NWK
        if (i == rowIndex(k)) then
            ! 到达新的一行，更新 rowIndex(k)
            rowIndex(k) = rowIndex(k-1) + temp
            k = k + 1
            temp = 0
        end if

        if (abs(A(i)) .gt. 1.0d-5) then
            ! 如果是“有效”的非零元素，就保留下来
            j = j + 1
            A(j) = A(i)
            columns(j) = columns(i)
            temp = temp + 1
        end if
    end do

    ! 更新最后一行的 rowIndex
    rowIndex(NEQ + 1) = rowIndex(NEQ) + temp
    NWK = j           ! 更新当前矩阵的实际存储长度
end subroutine

subroutine pardiso_addban(A, rowIndex, columns, S, LM, ND)
    use GLOBALS, only : NWK, NEQ
    implicit none

    real(8)    :: A(NWK), S(ND,ND)              ! A: 全局稀疏矩阵值，S: 单元刚度矩阵
    integer    :: rowIndex(NEQ+1), columns(NWK) ! CSR 结构
    integer    :: LM(ND)                        ! 单元自由度映射到全局自由度
    integer    :: ND                            ! 单元自由度个数
    integer    :: I, J, II, JJ, K, KK, KKK      ! 循环变量、索引变量
    integer    :: tempJ                         ! 可选临时变量

    do J = 1, ND
        JJ = LM(J)
        if (JJ > 0) then
            do I = J, ND
                II = LM(I)
                if (II > 0) then
                    ! 将 (I,J) 项添加到 (II,JJ) 或 (JJ,II)，取对称位置
                    KK  = min(II, JJ)
                    KKK = max(II, JJ)
                    ! 在稀疏结构中找到匹配位置并累加 S(I,J)
                    loop: do K = rowIndex(KK), rowIndex(KK+1) - 1
                        if (columns(K) == KKK) then
                            A(K) = A(K) + S(I,J)
                            exit loop
                        end if
                    end do loop
                end if
            end do
        end if
    end do

    return
end subroutine

subroutine pardiso_solver(K, V, rowIndex, columns)
    use mkl_pardiso                     ! Intel MKL PARDISO 接口
    use GLOBALS, only : neq, nwk, tim, IOUT
    implicit none

    real(8)    :: K(NWK), V(NEQ)                 ! 系数矩阵和右端向量
    integer    :: rowIndex(NEQ+1), columns(NWK)  ! CSR 格式

    ! 内部变量
    type(MKL_PARDISO_HANDLE) :: pt(64)           ! PARDISO 内部内存指针
    integer :: maxfct, mnum, mtype, phase, nrhs, error, msglvl, error1
    integer :: iparm(64), idum(1), i
    real(8) :: x(NEQ), ddum(1)

    ! 初始化控制参数
    iparm = 0
    iparm(1)  = 1   ! 自定义控制参数
    iparm(2)  = 3   ! 使用 METIS 重排序
    iparm(4)  = 0   ! 使用直接解法
    iparm(8)  = 2   ! 最多2步迭代修正
    iparm(10) = 13  ! 设定主元扰动值为1e-13
    iparm(11) = 1   ! 非对称比例缩放
    iparm(13) = 0   ! 关闭最大加权匹配
    iparm(60) = 1   ! 使用磁盘存储（如内存不足）

    ! 初始化指针
    do i = 1, 64
        pt(i)%DUMMY = 0
    end do

    ! 解方程组 phase=13（仅求解）
    mtype   = 2      ! 矩阵类型：对称不定
    nrhs    = 1
    maxfct  = 1
    mnum    = 1
    msglvl  = 1
    error   = 0
    phase   = 13

    call pardiso(pt, maxfct, mnum, mtype, phase, neq, K, rowIndex, columns, &
                 idum, nrhs, iparm, msglvl, V, x, error)

    if (error /= 0) then
        write(IOUT, *) 'PARDISO 求解错误，错误码: ', error
        stop
    end if

    ! 输出解向量
    V = x

    ! 清理内存
    phase = -1
    call pardiso(pt, maxfct, mnum, mtype, phase, neq, ddum, idum, idum, &
                 idum, nrhs, iparm, msglvl, ddum, ddum, error1)

    if (error1 /= 0) then
        write(IOUT, *) 'PARDISO 释放内存时出错，错误码: ', error1
        stop
    end if

end subroutine
