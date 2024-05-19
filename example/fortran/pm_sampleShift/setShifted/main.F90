program example

    use pm_kind, only: SK, IK, LK
    use pm_sampleMean, only: getMean
    use pm_sampleShift, only: transHerm
    use pm_sampleShift, only: setShifted
    use pm_arraySpace, only: getLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK) :: dim
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Shift a 1D sample to have a zero mean.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        real(TKG), allocatable :: sample(:)
        call disp%show("sample = getLinSpace(x1 = 0., x2 = 10., count = 11_IK)")
                        sample = getLinSpace(x1 = 0., x2 = 10., count = 11_IK)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getMean(sample)")
        call disp%show( getMean(sample) )
        call disp%show("call setShifted(sample, -getMean(sample))")
                        call setShifted(sample, -getMean(sample))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getMean(sample)")
        call disp%show( getMean(sample) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        complex(TKG), allocatable :: sample(:)
        call disp%show("sample = cmplx(getLinSpace(x1 = 0., x2 = 10., count = 11_IK), -getLinSpace(x1 = 0., x2 = 10., count = 11_IK), TKG)")
                        sample = cmplx(getLinSpace(x1 = 0., x2 = 10., count = 11_IK), -getLinSpace(x1 = 0., x2 = 10., count = 11_IK), TKG)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getMean(sample)")
        call disp%show( getMean(sample) )
        call disp%show("call setShifted(sample, -getMean(sample))")
                        call setShifted(sample, -getMean(sample))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getMean(sample)")
        call disp%show( getMean(sample) )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Shift a 2D sample to have a zero mean along the second dimension.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        real(TKG), allocatable :: sample(:,:)
        call disp%show("sample = reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5])")
                        sample = reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5])
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("dim = 2")
                        dim = 2
        call disp%show("getMean(sample, dim)")
        call disp%show( getMean(sample, dim) )
        call disp%show("call setShifted(sample, dim, -getMean(sample, dim))")
                        call setShifted(sample, dim, -getMean(sample, dim))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getMean(sample, dim)")
        call disp%show( getMean(sample, dim) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        complex(TKG), allocatable :: sample(:,:)
        call disp%show("sample = cmplx(reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]), -reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]))")
                        sample = cmplx(reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]), -reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("dim = 2")
                        dim = 2
        call disp%show("getMean(sample, dim)")
        call disp%show( getMean(sample, dim) )
        call disp%show("call setShifted(sample, dim, -getMean(sample, dim))")
                        call setShifted(sample, dim, -getMean(sample, dim))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getMean(sample, dim)")
        call disp%show( getMean(sample, dim) )
        call disp%skip()
    end block

end program example