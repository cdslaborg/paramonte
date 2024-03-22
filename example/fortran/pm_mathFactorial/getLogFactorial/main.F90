program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKS, RKD, RKH
    use pm_io, only: display_type
    use pm_mathFactorial, only: getLogFactorial

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RKH)               :: gamIncLow_RKH, x_RKH, shape_RKH
    real(RKD)               :: gamIncLow_RKD, x_RKD, shape_RKD
    real(RKS)               :: gamIncLow_RKS, x_RKS, shape_RKS
    logical(LK)             :: failed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("exp(getLogFactorial(x = 3._RKS))")
    call disp%show( exp(getLogFactorial(x = 3._RKS)) )
    call disp%skip()

    call disp%skip()
    call disp%show("exp(getLogFactorial(x = 4._RKD))")
    call disp%show( exp(getLogFactorial(x = 4._RKD)) )
    call disp%skip()

    call disp%skip()
    call disp%show("exp(getLogFactorial(x = 5._RKH))")
    call disp%show( exp(getLogFactorial(x = 5._RKH)) )
    call disp%skip()

    call disp%skip()
    call disp%show("exp(getLogFactorial(x = real([3., 5., 7., 10.],RKS)))")
    call disp%show( exp(getLogFactorial(x = real([3., 5., 7., 10.],RKS))) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Lower Incomplete Gamma function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        real(RKS), allocatable :: WholeNumber(:)
        integer :: fileUnit, i

        WholeNumber = getRange(start = 0_IK, stop = 30_IK, step = 5_IK)
        open(newunit = fileUnit, file = "getLogFactorial.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (WholeNumber(i), exp(getLogFactorial(WholeNumber(i))), i = 1, size(WholeNumber))
        close(fileUnit)

    end block

end program example