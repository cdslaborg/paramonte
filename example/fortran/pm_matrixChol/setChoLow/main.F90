program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK, RKS, RKD, RKH
    use pm_io, only: display_type
    use pm_matrixChol, only: setChoLow

    implicit none

    real(RKH), allocatable :: Matrix_RKH(:,:), ChoDia_RKH(:)
    real(RKD), allocatable :: Matrix_RKD(:,:), ChoDia_RKD(:)
    real(RKS), allocatable :: Matrix_RKS(:,:), ChoDia_RKS(:)

    logical(LK) :: failed
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Construct the Cholesky lower-triangle: RKH
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Matrix_RKH = reshape(   [ 1._RKH, 0._RKH, 0._RKH &
                            , 0._RKH, 4._RKH, 0._RKH &
                            , 2._RKH, 0._RKH, 8._RKH ], shape = [3,3])

    call disp%show("Matrix Upper Triangle")
    call disp%show( Matrix_RKH )

    allocate(ChoDia_RKH(size(Matrix_RKH, dim = 1)))
    call setChoLow(Matrix_RKH, ChoDia_RKH, size(ChoDia_RKH, 1, IK))
    if (ChoDia_RKH(1) <= 0) error stop

    call disp%show("Matrix Upper Triangle / Cholesky Lower Triangle")
    call disp%show( Matrix_RKH )

    call disp%show("Cholesky Diagonal")
    call disp%show( ChoDia_RKH )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Construct the Cholesky lower-triangle: RKD
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Matrix_RKD = real(Matrix_RKH, kind = RKD)

    call disp%skip()
    call disp%show("Matrix Upper Triangle")
    call disp%show( Matrix_RKD )

    allocate(ChoDia_RKD(size(Matrix_RKD, dim = 1)))
    call setChoLow(Matrix_RKD, ChoDia_RKD, size(ChoDia_RKD, 1, IK))
    if (ChoDia_RKD(1) <= 0) error stop

    call disp%show("Matrix Upper Triangle / Cholesky Lower Triangle")
    call disp%show( Matrix_RKD )

    call disp%show("Cholesky Diagonal")
    call disp%show( ChoDia_RKD )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Construct the Cholesky lower-triangle: RKS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Matrix_RKS = real(Matrix_RKH, kind = RKS)

    call disp%skip()
    call disp%show("Matrix Upper Triangle")
    call disp%show( Matrix_RKS )

    allocate(ChoDia_RKS(size(Matrix_RKS, dim = 1)))
    call setChoLow(Matrix_RKS, ChoDia_RKS, size(ChoDia_RKS, 1, IK))
    if (ChoDia_RKS(1) <= 0) error stop

    call disp%show("Matrix Upper Triangle / Cholesky Lower Triangle")
    call disp%show( Matrix_RKS )

    call disp%show("Cholesky Diagonal")
    call disp%show( ChoDia_RKS )

end program example