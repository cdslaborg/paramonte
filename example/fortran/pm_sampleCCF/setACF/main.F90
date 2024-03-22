program example

    use pm_kind, only: SK, IK, LK, RK
    use pm_io, only: getErrTableWrite
    use pm_io, only: display_type
    use pm_io, only: getFormat
    use pm_fftpack, only: allocatable
    use pm_fftpack, only: getFactorFFT
    use pm_sampleCCF, only: setACF
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
        real(TKC), allocatable :: seq(:), f(:), g(:), acf(:), coef(:), range(:)
        real(TKC), parameter :: ZERO = 0._TKC
        call disp%skip()
        call disp%show("nsam = 41")
                        nsam = 41
        call disp%show("range = getLinSpace(-10., 10., nsam)")
                        range = getLinSpace(-10., 10., nsam)
        call disp%show("seq = sin(range)")
                        seq = sin(range)
        call disp%show("seq = getPaddedr(getShifted(seq, -getMean(seq)), nsam - 1, ZERO)")
                        seq = getPaddedr(getShifted(seq, -getMean(seq)), nsam - 1, ZERO)
        call disp%show("seq")
        call disp%show( seq )
        call disp%show("call setResized(acf, size(seq, 1, IK))")
                        call setResized(acf, size(seq, 1, IK))
        call disp%show("factor = getFactorFFT(seq, coef, allocatable)")
                        factor = getFactorFFT(seq, coef, allocatable)
        call disp%show("f = seq")
                        f = seq
        call disp%show("call setACF(factor, coef, f, acf, inf)")
                        call setACF(factor, coef, f, acf, inf)
        call disp%show("acf = merge(f, acf, inf) / size(acf)")
                        acf = merge(f, acf, inf) / size(acf)
        call disp%show("size(acf)")
        call disp%show( size(acf) )
        call disp%show("acf")
        call disp%show( acf )
        call disp%show("f = seq; g = seq")
                        f = seq; g = seq
        call disp%show("call setCCF(factor, coef, f, g, acf, inf) ! for comparison with the above.")
                        call setCCF(factor, coef, f, g, acf, inf) ! for comparison with the above.
        call disp%show("acf = merge(f, g, inf) / size(acf)")
                        acf = merge(f, g, inf) / size(acf)
        call disp%show("size(acf)")
        call disp%show( size(acf) )
        call disp%show("acf")
        call disp%show( acf )
        call disp%show("lag = getRange(-nsam + 1_IK, nsam - 1_IK)")
                        lag = getRange(-nsam + 1_IK, nsam - 1_IK)
        call disp%show("if (0 /= getErrTableWrite(SK_'setACF.crd.sin.RK.txt', reshape([range, seq], [nsam, 2_IK]), header = SK_'crd,f')) error stop 'acf outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setACF.crd.sin.RK.txt', reshape([range, seq], [nsam, 2_IK]), header = SK_'crd,f')) error stop 'acf outputting failed.'
        call disp%show("if (0 /= getErrTableWrite(SK_'setACF.acf.sin.RK.txt', reshape([real(lag, TKC), acf], [size(lag), 2]), header = SK_'lag,acf')) error stop 'acf outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'setACF.acf.sin.RK.txt', reshape([real(lag, TKC), acf], [size(lag), 2]), header = SK_'lag,acf')) error stop 'acf outputting failed.'
        call disp%skip()
    end block

end program example