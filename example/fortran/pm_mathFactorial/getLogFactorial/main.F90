program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_mathFactorial, only: getLogFactorial

    implicit none

    integer(IK) , parameter :: NP = 1000_IK
    real(RK128)             :: gamIncLow_RK128, x_RK128, shape_RK128
    real(RK64 )             :: gamIncLow_RK64 , x_RK64 , shape_RK64
    real(RK32 )             :: gamIncLow_RK32 , x_RK32 , shape_RK32
    logical(LK)             :: failed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("exp(getLogFactorial(x = 3._RK32))")
    call disp%show( exp(getLogFactorial(x = 3._RK32)) )
    call disp%skip()

    call disp%skip()
    call disp%show("exp(getLogFactorial(x = 4._RK64))")
    call disp%show( exp(getLogFactorial(x = 4._RK64)) )
    call disp%skip()

    call disp%skip()
    call disp%show("exp(getLogFactorial(x = 5._RK128))")
    call disp%show( exp(getLogFactorial(x = 5._RK128)) )
    call disp%skip()

    call disp%skip()
    call disp%show("exp(getLogFactorial(x = real([3., 5., 7., 10.],RK32)))")
    call disp%show( exp(getLogFactorial(x = real([3., 5., 7., 10.],RK32))) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Lower Incomplete Gamma function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arrayRange, only: getRange
        real(RK32), allocatable :: WholeNumber(:)
        integer :: fileUnit, i

        WholeNumber = getRange(start = 0_IK, stop = 30_IK, step = 5_IK)
        open(newunit = fileUnit, file = "getLogFactorial.RK.txt")
        write(fileUnit,"(2(g0,:,' '))") (WholeNumber(i), exp(getLogFactorial(WholeNumber(i))), i = 1, size(WholeNumber))
        close(fileUnit)

    end block

end program example