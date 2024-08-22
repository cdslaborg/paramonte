program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKH ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distBeta, only: setBetaCDF
    use pm_mathBeta, only: getLogBeta

    implicit none

    real(RKG) :: cdf(4)
    logical(LK) :: signed
    integer(IK) :: info(4)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    signed = .false.

    call disp%skip()
    call disp%show("call setBetaCDF(cdf(1), x = 0._RKG, alpha = 2._RKG, beta = 3._RKG, logFuncBeta = getLogBeta(2._RKG, 3._RKG), signed = signed, info = info(1))")
                    call setBetaCDF(cdf(1), x = 0._RKG, alpha = 2._RKG, beta = 3._RKG, logFuncBeta = getLogBeta(2._RKG, 3._RKG), signed = signed, info = info(1))
    call disp%show("if (info(1) /= 0) error stop")
                    if (info(1) /= 0) error stop
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setBetaCDF(cdf(1), x = .5_RKG, alpha = 2._RKG, beta = 3._RKG, logFuncBeta = getLogBeta(2._RKG, 3._RKG), signed = signed, info = info(1))")
                    call setBetaCDF(cdf(1), x = .5_RKG, alpha = 2._RKG, beta = 3._RKG, logFuncBeta = getLogBeta(2._RKG, 3._RKG), signed = signed, info = info(1))
    call disp%show("if (info(1) /= 0) error stop")
                    if (info(1) /= 0) error stop
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setBetaCDF(cdf(1), x = 1._RKG, alpha = 2._RKG, beta = 3._RKG, logFuncBeta = getLogBeta(2._RKG, 3._RKG), signed = signed, info = info(1))")
                    call setBetaCDF(cdf(1), x = 1._RKG, alpha = 2._RKG, beta = 3._RKG, logFuncBeta = getLogBeta(2._RKG, 3._RKG), signed = signed, info = info(1))
    call disp%show("if (info(1) /= 0) error stop")
                    if (info(1) /= 0) error stop
    call disp%show("cdf(1)")
    call disp%show( cdf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setBetaCDF(cdf, x = [0._RKG, 0.1_RKG, 0.5_RKG, 1._RKG], alpha = 2._RKG, beta = 3._RKG, logFuncBeta = getLogBeta(2._RKG, 3._RKG), signed = signed, info = info)")
                    call setBetaCDF(cdf, x = [0._RKG, 0.1_RKG, 0.5_RKG, 1._RKG], alpha = 2._RKG, beta = 3._RKG, logFuncBeta = getLogBeta(2._RKG, 3._RKG), signed = signed, info = info)
    call disp%show("if (any(info /= 0)) error stop")
                    if (any(info /= 0)) error stop
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()

    call disp%skip()
    call disp%show("call setBetaCDF(cdf, x = [0._RKG, .1_RKG, 0.5_RKG, 1._RKG], alpha = [0.1_RKG, .5_RKG, 1._RKG, 10._RKG], beta = 3._RKG, logFuncBeta = getLogBeta([0.1_RKG, .5_RKG, 1._RKG, 10._RKG], 3._RKG), signed = signed, info = info)")
                    call setBetaCDF(cdf, x = [0._RKG, .1_RKG, 0.5_RKG, 1._RKG], alpha = [0.1_RKG, .5_RKG, 1._RKG, 10._RKG], beta = 3._RKG, logFuncBeta = getLogBeta([0.1_RKG, .5_RKG, 1._RKG, 10._RKG], 3._RKG), signed = signed, info = info)
    call disp%show("if (any(info /= 0)) error stop")
                    if (any(info /= 0)) error stop
    call disp%show("cdf")
    call disp%show( cdf )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block

        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 1000_IK
        real(RKG) :: cdf(4), X(NP)
        integer :: fileUnit, i

        call setLinSpace(X, 0._RKG, 1._RKG)
        open(newunit = fileUnit, file = "setBetaCDF.RK.txt")
        do i = 1, NP
            call setBetaCDF(cdf, X(i), [.5_RKG, 5._RKG, .5_RKG, 5._RKG], [.5_RKG, 1.0_RKG, 5.0_RKG, 10._RKG], getLogBeta([.5_RKG, 5._RKG, .5_RKG, 5._RKG], [.5_RKG, 1.0_RKG, 5.0_RKG, 10._RKG]), signed = signed, info = info)
            if (any(info /= 0)) error stop
            write(fileUnit, "(*(g0,:,' '))") X(i), cdf
        end do
        close(fileUnit)

    end block

end program example