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
        real :: LogMinX(2), LogMaxX(2), CDF(2), LogX(2000)
        integer(IK) :: fileUnit, i, j
        call setLinSpace(LogX, x1 = log(1.), x2 = log(13.))
        LogMinX = log([3., 2.0])
        LogMaxX = log([7., 10.])
        open(newunit = fileUnit, file = "setLogUnifCDF.RK.txt")
        do i = 1, size(LogX)
            do j = 1, size(CDF)
                if (LogMaxX(j) <= LogX(i)) then
                    CDF(j) = 1.
                elseif (LogMinX(j) <= LogX(i)) then
                    call setLogUnifCDF(CDF(j), LogX(i), LogMinX(j), getLogUnifPDFNF(LogMinX(j), LogMaxX(j)))
                else
                    CDF(j) = 0.
                end if
            end do
            write(fileUnit, "(*(g0,:,', '))") exp(LogX(i)), CDF
        end do
        close(fileUnit)
    end block

end program example