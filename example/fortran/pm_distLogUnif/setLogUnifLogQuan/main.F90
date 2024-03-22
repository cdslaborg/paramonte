program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distLogUnif, only: getLogUnifPDFNF
    use pm_distLogUnif, only: setLogUnifLogQuan

    implicit none

    real :: logx(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Inverse CDF (Quantile) Function of the LogUniform distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogUnifLogQuan(logx(1), cdf = 0.5, logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))")
                    call setLogUnifLogQuan(logx(1), cdf = 0.5, logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogUnifLogQuan(logx(1:3), cdf = [.3, .4, .25], logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))")
                    call setLogUnifLogQuan(logx(1:3), cdf = [.3, .4, .25], logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: logMinX(2), logMaxX(2), logx(2), CDF(2000)
        integer(IK) :: fileUnit, i
        call setLinSpace(CDF, x1 = 0., x2 = 1.)
        logMinX = log([3., 2.0])
        logMaxX = log([7., 10.])
        open(newunit = fileUnit, file = "setLogUnifLogQuan.RK.txt")
        do i = 1, size(CDF)
            call setLogUnifLogQuan(logx, CDF(i), logMinX, getLogUnifPDFNF(logMinX, logMaxX))
            write(fileUnit, "(*(g0,:,', '))") CDF(i), exp(logx)
        end do
        close(fileUnit)
    end block

end program example