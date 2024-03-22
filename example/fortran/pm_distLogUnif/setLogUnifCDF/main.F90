program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distLogUnif, only: setLogUnifCDF
    use pm_distLogUnif, only: getLogUnifPDFNF

    implicit none

    real :: CDF(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Cumulative Distribution Function (CDF) of the LogUniform distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogUnifCDF(CDF(1), logx = log(3.), logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))")
                    call setLogUnifCDF(CDF(1), logx = log(3.), logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))
    call disp%show("CDF(1)")
    call disp%show( CDF(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("call setLogUnifCDF(CDF(1:3), logx = log([3., 4., 2.5]), logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))")
                    call setLogUnifCDF(CDF(1:3), logx = log([3., 4., 2.5]), logMinX = log(2.), pdfnf = getLogUnifPDFNF(logMinX = log(2.), logMaxX = log(5.)))
    call disp%show("CDF(1:3)")
    call disp%show( CDF(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: logMinX(2), logMaxX(2), CDF(2), logx(2000)
        integer(IK) :: fileUnit, i, j
        call setLinSpace(logx, x1 = log(1.), x2 = log(13.))
        logMinX = log([3., 2.0])
        logMaxX = log([7., 10.])
        open(newunit = fileUnit, file = "setLogUnifCDF.RK.txt")
        do i = 1, size(logx)
            do j = 1, size(CDF)
                if (logMaxX(j) <= logx(i)) then
                    CDF(j) = 1.
                elseif (logMinX(j) <= logx(i)) then
                    call setLogUnifCDF(CDF(j), logx(i), logMinX(j), getLogUnifPDFNF(logMinX(j), logMaxX(j)))
                else
                    CDF(j) = 0.
                end if
            end do
            write(fileUnit, "(*(g0,:,', '))") exp(logx(i)), CDF
        end do
        close(fileUnit)
    end block

end program example