program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_distLogUnif, only: getLogUnifPDFNF
    use pm_distLogUnif, only: setLogUnifLogRand

    implicit none

    real :: logx(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random value(s) from the LogUniform distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogUnifLogRand(logx(1), urand = getUnifRand(0., 1.), logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))")
                    call setLogUnifLogRand(logx(1), urand = getUnifRand(0., 1.), logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))
    call disp%show("logx(1)")
    call disp%show( logx(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogUnifLogRand(logx(1:3), urand = getUnifRand(0., 1.), logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))")
                    call setLogUnifLogRand(logx(1:3), urand = getUnifRand(0., 1.), logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))
    call disp%show("logx(1:3)")
    call disp%show( logx(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: logMinX(2), logMaxX(2), LogRand(2), UnifRand(2000)
        integer(IK) :: fileUnit, i
        call setUnifRand(UnifRand)
        logMinX = log([3., 2.0])
        logMaxX = log([7., 10.])
        open(newunit = fileUnit, file = "setLogUnifLogRand.RK.txt")
        do i = 1, size(UnifRand)
            call setLogUnifLogRand(LogRand, UnifRand(i), logMinX, getLogUnifPDFNF(logMinX, logMaxX))
            write(fileUnit, "(*(g0,:,', '))") exp(LogRand)
        end do
        close(fileUnit)
    end block

end program example