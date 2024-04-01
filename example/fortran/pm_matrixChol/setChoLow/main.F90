program example

    use pm_kind, only: SK
    use pm_kind, only: IK, LK, RKS, RKD, RKH
    use pm_io, only: display_type
    use pm_matrixChol, only: setChoLow

    implicit none

    real(RKH), allocatable :: matrix_RKH(:,:), chodia_RKH(:)
    real(RKD), allocatable :: matrix_RKD(:,:), chodia_RKD(:)
    real(RKS), allocatable :: matrix_RKS(:,:), chodia_RKS(:)

    logical(LK) :: failed
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Construct the Cholesky lower-triangle: RKH
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matrix_RKH = reshape(   [ 1._RKH, 0._RKH, 0._RKH &
                            , 0._RKH, 4._RKH, 0._RKH &
                            , 2._RKH, 0._RKH, 8._RKH ], shape = [3,3])

    call disp%show("matrix Upper Triangle")
    call disp%show( matrix_RKH )

    allocate(chodia_RKH(size(matrix_RKH, dim = 1)))
    call setChoLow(matrix_RKH, chodia_RKH, size(chodia_RKH, 1, IK))
    if (chodia_RKH(1) <= 0) error stop

    call disp%show("matrix Upper Triangle / Cholesky Lower Triangle")
    call disp%show( matrix_RKH )

    call disp%show("Cholesky Diagonal")
    call disp%show( chodia_RKH )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Construct the Cholesky lower-triangle: RKD
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matrix_RKD = real(matrix_RKH, kind = RKD)

    call disp%skip()
    call disp%show("matrix Upper Triangle")
    call disp%show( matrix_RKD )

    allocate(chodia_RKD(size(matrix_RKD, dim = 1)))
    call setChoLow(matrix_RKD, chodia_RKD, size(chodia_RKD, 1, IK))
    if (chodia_RKD(1) <= 0) error stop

    call disp%show("matrix Upper Triangle / Cholesky Lower Triangle")
    call disp%show( matrix_RKD )

    call disp%show("Cholesky Diagonal")
    call disp%show( chodia_RKD )

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Construct the Cholesky lower-triangle: RKS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matrix_RKS = real(matrix_RKH, kind = RKS)

    call disp%skip()
    call disp%show("matrix Upper Triangle")
    call disp%show( matrix_RKS )

    allocate(chodia_RKS(size(matrix_RKS, dim = 1)))
    call setChoLow(matrix_RKS, chodia_RKS, size(chodia_RKS, 1, IK))
    if (chodia_RKS(1) <= 0) error stop

    call disp%show("matrix Upper Triangle / Cholesky Lower Triangle")
    call disp%show( matrix_RKS )

    call disp%show("Cholesky Diagonal")
    call disp%show( chodia_RKS )

end program example