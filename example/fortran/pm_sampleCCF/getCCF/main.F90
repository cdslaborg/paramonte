program example

    use pm_kind, only: SK, IK, LK, RK
    use pm_io, only: getErrTableWrite
    use pm_io, only: display_type
    use pm_io, only: getFormat
    use pm_sampleCCF, only: getCCF
    use pm_sampleCCF, only: stdscale
    use pm_arrayRange, only: getRange
    use pm_distUnif, only: getUnifRand
    use pm_distNorm, only: getNormLogPDF
    use pm_arrayPad, only: getPaddedr
    use pm_arrayFill, only: getFilled
    use pm_sampleShift, only: getShifted
    use pm_arrayResize, only: setResized
    use pm_arraySpace, only: getLinSpace
    use pm_sampleNorm, only: getNormed
    use pm_sampleMean, only: getMean
    use pm_sampleVar, only: getVar

    implicit none

    type(display_type) :: disp
    character(:), allocatable :: format
    integer(IK) :: nsam, shift, itry, ntry = 1
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the cross-correlation of two samples.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported.
        integer(IK), allocatable :: lag(:)
        real(TKC), allocatable :: f(:), g(:)
        real(TKC), allocatable :: range(:), ccf(:)
        call disp%skip()
        call disp%show("nsam = 41; shift = nsam")
                        nsam = 41; shift = nsam
        call disp%show("range = getLinSpace(-10., 10., nsam)")
                        range = getLinSpace(-10., 10., nsam)
        call disp%show("f = sin(range)")
                        f = sin(range)
        call disp%show("f")
        call disp%show( f )
        call disp%show("g = cos(range)")
                        g = cos(range)
        call disp%show("g")
        call disp%show( g )
        call disp%show("ccf = getCCF(f, g)")
                        ccf = getCCF(f, g)
        call disp%show("size(ccf)")
        call disp%show( size(ccf) )
        call disp%show("ccf")
        call disp%show( ccf )
        call disp%show("lag = getRange(-nsam + 1, nsam - 1)")
                        lag = getRange(-nsam + 1, nsam - 1)
        call disp%show("ccf = getCCF(f, g, lag)")
                        ccf = getCCF(f, g, lag)
        call disp%show("size(ccf)")
        call disp%show( size(ccf) )
        call disp%show("ccf")
        call disp%show( ccf )
        call disp%show("if (0 /= getErrTableWrite(SK_'getCCF.crd.sin.cos.RK.txt', reshape([range, f, g], [nsam, 3_IK]), header = SK_'crd,f,g')) error stop 'ccf outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'getCCF.crd.sin.cos.RK.txt', reshape([range, f, g], [nsam, 3_IK]), header = SK_'crd,f,g')) error stop 'ccf outputting failed.'
        call disp%show("if (0 /= getErrTableWrite(SK_'getCCF.ccf.sin.cos.RK.txt', reshape([real(lag, TKC), ccf], [size(lag), 2]), header = SK_'lag,ccf')) error stop 'ccf outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'getCCF.ccf.sin.cos.RK.txt', reshape([real(lag, TKC), ccf], [size(lag), 2]), header = SK_'lag,ccf')) error stop 'ccf outputting failed.'
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported.
        complex(TKC), allocatable :: f(:), g(:)
        complex(TKC), allocatable :: range(:), ccf(:)
        format = getFormat(deliml = '', subsep = SK_'', delimr = 'i', subcount = 2_IK, signed = .true._LK)
        call disp%skip()
        call disp%show("nsam = 41; shift = nsam")
                        nsam = 41; shift = nsam
        call disp%show("range = cmplx(getLinSpace(-10., 10., nsam), getLinSpace(10., -10., nsam), TKC)")
                        range = cmplx(getLinSpace(-10., 10., nsam), getLinSpace(10., -10., nsam), TKC)
        call disp%show("f = sin(range)")
                        f = sin(range)
        call disp%show("f")
        call disp%show( f , format = format )
        call disp%show("g = cos(range)")
                        g = cos(range)
        call disp%show("g")
        call disp%show( g , format = format )
        call disp%show("ccf = getCCF(f, g, norm = stdscale)")
                        ccf = getCCF(f, g, norm = stdscale)
        call disp%show("ccf")
        call disp%show( ccf , format = format )
        call disp%skip()
    end block

end program example