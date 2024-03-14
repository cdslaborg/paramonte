program example

    use pm_kind, only: SK, IK, LK, RK
    use pm_io, only: getErrTableWrite
    use pm_io, only: display_type
    use pm_io, only: getFormat
    use pm_sampleCCF, only: getACF
    use pm_sampleCCF, only: getCCF
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
    integer(IK) :: nsam, shift, itry, ntry = 1
    character(:), allocatable :: format
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the cross-correlation of two samples.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported, e.g., RK32, RK64, RK128.
        integer(IK), allocatable :: lag(:)
        real(TKC), allocatable :: f(:), g(:)
        real(TKC), allocatable :: range(:), acf(:)
        call disp%skip()
        call disp%show("nsam = 41; shift = nsam")
                        nsam = 41; shift = nsam
        call disp%show("range = getLinSpace(-10., 10., nsam)")
                        range = getLinSpace(-10., 10., nsam)
        call disp%show("f = sin(range)")
                        f = sin(range)
        call disp%show("f")
        call disp%show( f )
        call disp%show("acf = getACF(f)")
                        acf = getACF(f)
        call disp%show("size(acf)")
        call disp%show( size(acf) )
        call disp%show("acf")
        call disp%show( acf )
        call disp%show("g = f")
                        g = f
        call disp%show("acf = getCCF(f, g, lag = getRange(1 - nsam, nsam - 1_IK)) ! same as above.")
                        acf = getCCF(f, g, lag = getRange(1 - nsam, nsam - 1_IK)) ! same as above.
        call disp%show("size(acf)")
        call disp%show( size(acf) )
        call disp%show("acf")
        call disp%show( acf )
        call disp%show("acf = getCCF(f, g, lag = getRange(0_IK, nsam - 1_IK)) ! same as above.")
                        acf = getCCF(f, g, lag = getRange(0_IK, nsam - 1_IK)) ! same as above.
        call disp%show("size(acf)")
        call disp%show( size(acf) )
        call disp%show("acf")
        call disp%show( acf )
        call disp%show("lag = getRange(-nsam + 1_IK, nsam - 1_IK)")
                        lag = getRange(-nsam + 1_IK, nsam - 1_IK)
        call disp%show("acf = getACF(f, lag)")
                        acf = getACF(f, lag)
        call disp%show("size(acf)")
        call disp%show( size(acf) )
        call disp%show("acf")
        call disp%show( acf )
        call disp%show("if (0 /= getErrTableWrite(SK_'getACF.crd.sin.RK.txt', reshape([range, f], [nsam, 2_IK]), header = SK_'crd,f')) error stop 'acf outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'getACF.crd.sin.RK.txt', reshape([range, f], [nsam, 2_IK]), header = SK_'crd,f')) error stop 'acf outputting failed.'
        call disp%show("if (0 /= getErrTableWrite(SK_'getACF.acf.sin.RK.txt', reshape([real(lag, TKC), acf], [size(lag), 2]), header = SK_'lag,acf')) error stop 'acf outputting failed.'")
                        if (0 /= getErrTableWrite(SK_'getACF.acf.sin.RK.txt', reshape([real(lag, TKC), acf], [size(lag), 2]), header = SK_'lag,acf')) error stop 'acf outputting failed.'
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKS ! All other real types are also supported, e.g., RK32, RK64, RK128.
        integer(IK), allocatable :: lag(:)
        complex(TKC), allocatable :: f(:), g(:), acf(:), range(:)
        call disp%skip()
        call disp%show("nsam = 41; shift = nsam")
                        nsam = 41; shift = nsam
        call disp%show("range = getLinSpace((-10._TKC, +10._TKC), (+10._TKC, -10._TKC), nsam)")
                        range = getLinSpace((-10._TKC, +10._TKC), (+10._TKC, -10._TKC), nsam)
        call disp%show("f = sin(range)")
                        f = sin(range)
        call disp%show("f")
        call disp%show( f )
        call disp%show("acf = getACF(f)")
                        acf = getACF(f)
        call disp%show("size(acf)")
        call disp%show( size(acf) )
        call disp%show("acf")
        call disp%show( acf )
        call disp%show("acf%re")
        call disp%show( acf%re )
        call disp%show("g = f")
                        g = f
        call disp%show("acf = getCCF(f, g, lag = getRange(0_IK, nsam - 1_IK)) ! same as above.")
                        acf = getCCF(f, g, lag = getRange(0_IK, nsam - 1_IK)) ! same as above.
        call disp%show("size(acf)")
        call disp%show( size(acf) )
        call disp%show("acf")
        call disp%show( acf )
        call disp%show("acf%re")
        call disp%show( acf%re )
        call disp%show("lag = getRange(-nsam + 1_IK, nsam - 1_IK)")
                        lag = getRange(-nsam + 1_IK, nsam - 1_IK)
        call disp%show("acf = getACF(f, lag)")
                        acf = getACF(f, lag)
        call disp%show("size(acf)")
        call disp%show( size(acf) )
        call disp%show("acf")
        call disp%show( acf )
        call disp%show("acf%re")
        call disp%show( acf%re )
        call disp%skip()
    end block

end program example