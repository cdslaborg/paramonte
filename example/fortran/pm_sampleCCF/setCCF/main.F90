program example

    use pm_kind, only: SK, IK, LK, RK
    use pm_io, only: getErrTableWrite
    use pm_io, only: display_type
    use pm_io, only: getFormat
    use pm_fftpack, only: allocatable
    use pm_fftpack, only: getFactorFFT
    use pm_sampleCCF, only: setCCF
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

    logical(LK) :: inf
    type(display_type) :: disp
    character(:), allocatable :: format
    integer(IK) :: nsam, itry, ntry = 1
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the cross-correlation of two samples.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported.
        integer(IK), allocatable :: factor(:), lag(:)
        real(TKC), allocatable :: f(:), g(:), fcopy(:), gcopy(:), ccf(:), coef(:), range(:)
        real(TKC), parameter :: ZERO = 0._TKC
        call disp%skip()
        call disp%show("nsam = 41")
                        nsam = 41
        call disp%show("range = getLinSpace(-10., 10., nsam)")
                        range = getLinSpace(-10., 10., nsam)
        call disp%show("f = sin(range)")
                        f = sin(range)
        call disp%show("f = getPaddedr(getNormed(f, -getMean(f), 1._TKC / getVar(f)), nsam - 1, ZERO)")
                        f = getPaddedr(getNormed(f, -getMean(f), 1._TKC / getVar(f)), nsam - 1, ZERO)
        call disp%show("f")
        call disp%show( f )
        call disp%show("g = cos(range)")
                        g = cos(range)
        call disp%show("g = getPaddedr(getNormed(g, -getMean(g), 1._TKC / getVar(g)), nsam - 1, ZERO)")
                        g = getPaddedr(getNormed(g, -getMean(g), 1._TKC / getVar(g)), nsam - 1, ZERO)
        call disp%show("g")
        call disp%show( g )
        call disp%show("fcopy = f; gcopy = g")
                        fcopy = f; gcopy = g
        call disp%show("call setResized(ccf, size(f, 1, IK))")
                        call setResized(ccf, size(f, 1, IK))
        call disp%show("factor = getFactorFFT(f, coef, allocatable)")
                        factor = getFactorFFT(f, coef, allocatable)
        call disp%show("call setCCF(factor, coef, fcopy, gcopy, ccf, inf)")
                        call setCCF(factor, coef, fcopy, gcopy, ccf, inf)
        call disp%show("ccf = merge(fcopy, gcopy, inf) / size(ccf)")
                        ccf = merge(fcopy, gcopy, inf) / size(ccf)
        call disp%show("size(ccf)")
        call disp%show( size(ccf) )
        call disp%show("ccf")
        call disp%show( ccf )
        call disp%show("lag = getRange(-nsam + 1_IK, nsam - 1_IK)")
                        lag = getRange(-nsam + 1_IK, nsam - 1_IK)
        call disp%show("if (0 /= getErrTableWrite(SK_'setCCF.crd.sin.cos.RK.txt', reshape([range, f(1:nsam), g(1:nsam)], [nsam, 3_IK]), header = SK_'crd,f,g')) error stop 'ccf outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setCCF.crd.sin.cos.RK.txt', reshape([range, f(1:nsam), g(1:nsam)], [nsam, 3_IK]), header = SK_'crd,f,g')) error stop 'ccf outputting failed.'
        call disp%show("if (0 /= getErrTableWrite(SK_'setCCF.ccf.sin.cos.RK.txt', reshape([real(lag, TKC), ccf], [size(lag), 2]), header = SK_'lag,ccf')) error stop 'ccf outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setCCF.ccf.sin.cos.RK.txt', reshape([real(lag, TKC), ccf], [size(lag), 2]), header = SK_'lag,ccf')) error stop 'ccf outputting failed.'
        call disp%skip()
    end block

end program example