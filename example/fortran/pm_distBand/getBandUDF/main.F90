program example

    use pm_kind, only: SK, IK
    use pm_distBand, only: getBandUDF
    use pm_distBand, only: getBandZeta
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: getLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 999
    real :: udf(NP), point(NP), alpha, beta, ebreak

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(point, x1 = 0.01, x2 = 10.)

    call disp%skip()
    call disp%show("point(1)")
    call disp%show( point(1) )
    call disp%show("alpha = -.5; beta = -1.5; ebreak = 100.")
                    alpha = -.5; beta = -1.5; ebreak = 100.
    call disp%show("udf(1) = getBandUDF(point(1), alpha, beta, ebreak)")
                    udf(1) = getBandUDF(point(1), alpha, beta, ebreak)
    call disp%show("udf(1)")
    call disp%show( udf(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("point(1)")
    call disp%show( point(1) )
    call disp%show("alpha = -.5; beta = -1.5; ebreak = 100.")
                    alpha = -.5; beta = -1.5; ebreak = 100.
    call disp%show("udf(1) = getBandUDF(point(1), alpha, beta, ebreak, zeta = getBandZeta(alpha, beta, ebreak), invEfold = (alpha - beta) / ebreak)")
                    udf(1) = getBandUDF(point(1), alpha, beta, ebreak, zeta = getBandZeta(alpha, beta, ebreak), invEfold = (alpha - beta) / ebreak)
    call disp%show("udf(1)")
    call disp%show( udf(1) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example udf array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "getBandUDF.RK.txt")
        do i = 1, NP
            write(fileUnit,"(*(g0,:,' '))") point(i), getBandUDF(point(i), [.5, 1.5, -.5, +1.5], -[.5, 1.0, 2., 3.], [.5, 1.0, 2.0, 5.])
        end do
        close(fileUnit)
    end block

end program example