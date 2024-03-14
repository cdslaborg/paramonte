program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_arraySpace, only: getLinSpace
    use pm_distMultiNorm, only: getMultiNormLogPDFNF

    implicit none

    integer(IK) :: info
    real, allocatable :: logPDFNF(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the normalization factor of the MultiNormorm distribution PDF.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF = getMultiNormLogPDFNF(ndim = 3_IK, logSqrtDetInvCov = [-3., .1, 2.])")
                    logPDFNF = getMultiNormLogPDFNF(ndim = 3_IK, logSqrtDetInvCov = [-3., .1, 2.])
    call disp%show("logPDFNF")
    call disp%show( logPDFNF )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF = getMultiNormLogPDFNF(invCov = reshape([1., .5, .5, 1.], [2, 2]))")
                    logPDFNF = getMultiNormLogPDFNF(invCov = reshape([1., .5, .5, 1.], [2, 2]))
    call disp%show("logPDFNF")
    call disp%show( logPDFNF )
    call disp%skip()

    call disp%skip()
    call disp%show("logPDFNF = getMultiNormLogPDFNF(reshape([1., .5, .5, 1.], [2, 2]), info)")
                    logPDFNF = getMultiNormLogPDFNF(reshape([1., .5, .5, 1.], [2, 2]), info)
    call disp%show("logPDFNF")
    call disp%show( logPDFNF )
    call disp%show("info")
    call disp%show( info )
    call disp%show("if (info /= 0) error stop")
                    if (info /= 0) error stop
    call disp%skip()

end program example