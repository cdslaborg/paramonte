program example

    use pm_kind, only: SK, IK, LK
    use pm_sampleVar, only: getVar
    use pm_sampleMean, only: getMean
    use pm_sampleScale, only: transHerm
    use pm_sampleScale, only: getScaled
    use pm_sampleShift, only: setShifted
    use pm_arraySpace, only: getLinSpace
    use pm_io, only: display_type

    implicit none

    integer(IK) :: dim
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Scale a 1D zero-mean sample to have a unit variance.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        real(TKG), allocatable :: sample(:)
        call disp%show("sample = getLinSpace(x1 = 0., x2 = 10., count = 11_IK)")
                        sample = getLinSpace(x1 = 0., x2 = 10., count = 11_IK)
        call disp%show("call setShifted(sample, -getMean(sample))")
                        call setShifted(sample, -getMean(sample))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample)")
        call disp%show( getVar(sample) )
        call disp%show("sample = getScaled(sample, 1._TKG / sqrt(getVar(sample)))")
                        sample = getScaled(sample, 1._TKG / sqrt(getVar(sample)))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample)")
        call disp%show( getVar(sample) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        complex(TKG), allocatable :: sample(:)
        call disp%show("sample = cmplx(getLinSpace(x1 = 0., x2 = 10., count = 11_IK), -getLinSpace(x1 = 0., x2 = 10., count = 11_IK), TKG)")
                        sample = cmplx(getLinSpace(x1 = 0., x2 = 10., count = 11_IK), -getLinSpace(x1 = 0., x2 = 10., count = 11_IK), TKG)
        call disp%show("call setShifted(sample, -getMean(sample))")
                        call setShifted(sample, -getMean(sample))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample)")
        call disp%show( getVar(sample) )
        call disp%show("sample = getScaled(sample, 1._TKG / sqrt(getVar(sample)))")
                        sample = getScaled(sample, 1._TKG / sqrt(getVar(sample)))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample)")
        call disp%show( getVar(sample) )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Scale a 2D zero-mean sample to have a unit variance.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        real(TKG), allocatable :: sample(:,:)
        call disp%show("sample = reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5])")
                        sample = reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5])
        call disp%show("dim = 2")
                        dim = 2
        call disp%show("call setShifted(sample, dim, -getMean(sample, dim))")
                        call setShifted(sample, dim, -getMean(sample, dim))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample, dim)")
        call disp%show( getVar(sample, dim) )
        call disp%show("sample = getScaled(sample, dim, 1._TKG / sqrt(getVar(sample, dim)))")
                        sample = getScaled(sample, dim, 1._TKG / sqrt(getVar(sample, dim)))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample, dim)")
        call disp%show( getVar(sample, dim) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        complex(TKG), allocatable :: sample(:,:)
        call disp%show("sample = cmplx(reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]), -reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]))")
                        sample = cmplx(reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]), -reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]))
        call disp%show("dim = 2")
                        dim = 2
        call disp%show("call setShifted(sample, dim, -getMean(sample, dim))")
                        call setShifted(sample, dim, -getMean(sample, dim))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample, dim)")
        call disp%show( getVar(sample, dim) )
        call disp%show("sample = getScaled(sample, dim, 1._TKG / sqrt(getVar(sample, dim)))")
                        sample = getScaled(sample, dim, 1._TKG / sqrt(getVar(sample, dim)))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample, dim)")
        call disp%show( getVar(sample, dim) )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Scale a 2D zero-mean sample to have a unit variance and transpose the result upon return.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        real(TKG), allocatable :: sample(:,:)
        call disp%show("sample = reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5])")
                        sample = reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5])
        call disp%show("dim = 2")
                        dim = 2
        call disp%show("call setShifted(sample, dim, -getMean(sample, dim))")
                        call setShifted(sample, dim, -getMean(sample, dim))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample, dim)")
        call disp%show( getVar(sample, dim) )
        call disp%show("sample = getScaled(sample, dim, 1._TKG / sqrt(getVar(sample, dim)), transHerm)")
                        sample = getScaled(sample, dim, 1._TKG / sqrt(getVar(sample, dim)), transHerm)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample, 3_IK - dim)")
        call disp%show( getVar(sample, 3_IK - dim) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKS ! All kinds are supported.
        complex(TKG), allocatable :: sample(:,:)
        call disp%show("sample = cmplx(reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]), -reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]))")
                        sample = cmplx(reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]), -reshape(getLinSpace(x1 = 1._TKG, x2 = 10._TKG, count = 10_IK), [2,5]))
        call disp%show("dim = 2")
                        dim = 2
        call disp%show("call setShifted(sample, dim, -getMean(sample, dim))")
                        call setShifted(sample, dim, -getMean(sample, dim))
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample, dim)")
        call disp%show( getVar(sample, dim) )
        call disp%show("sample = getScaled(sample, dim, 1._TKG / sqrt(getVar(sample, dim)), transHerm)")
                        sample = getScaled(sample, dim, 1._TKG / sqrt(getVar(sample, dim)), transHerm)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getVar(sample, 3_IK - dim)")
        call disp%show( getVar(sample, 3_IK - dim) )
        call disp%skip()
    end block

end program example