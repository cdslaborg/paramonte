program example

    use pm_kind, only: SK, IK, LK, RK
    use pm_kind, only: TKG => RKS ! All other real types are also supported.
    use pm_io, only: display_type
    use pm_sampleVar, only: setVar
    use pm_arraySpace, only: getLinSpace
    use pm_sampleCov, only: fweight, rweight
    use pm_sampleShift, only: getShifted
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_sampleMean, only: getMean
    use pm_sampleMean, only: setMean

    implicit none

    integer(IK) :: dim
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the variance of a 1-D sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG) :: var, mean
        real(TKG), allocatable :: sample(:)
        call disp%show("sample = getLinSpace(1._TKG, 9._TKG, 5_IK)")
                        sample = getLinSpace(1._TKG, 9._TKG, 5_IK)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("mean = getMean(sample)")
                        mean = getMean(sample)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("call setVar(var, mean, sample)")
                        call setVar(var, mean, sample)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, sample - mean)")
                        call setVar(var, sample - mean)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, mean, sample, dim = 1_IK)")
                        call setVar(var, mean, sample, dim = 1_IK)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, sample - mean, dim = 1_IK)")
                        call setVar(var, sample - mean, dim = 1_IK)
        call disp%show("var")
        call disp%show( var )
        call disp%skip()
    end block

    block
        real(TKG) :: var
        complex(TKG) :: mean
        complex(TKG), allocatable :: sample(:)
        call disp%show("sample = cmplx(getLinSpace(1., 9., 5_IK), -getLinSpace(1., 9., 5_IK), TKG)")
                        sample = cmplx(getLinSpace(1., 9., 5_IK), -getLinSpace(1., 9., 5_IK), TKG)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("mean = getMean(sample)")
                        mean = getMean(sample)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("call setVar(var, mean, sample)")
                        call setVar(var, mean, sample)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, sample - mean)")
                        call setVar(var, sample - mean)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, mean, sample, dim = 1_IK)")
                        call setVar(var, mean, sample, dim = 1_IK)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, sample - mean, dim = 1_IK)")
                        call setVar(var, sample - mean, dim = 1_IK)
        call disp%show("var")
        call disp%show( var )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the variance of a 1-D weighted sample.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG) :: var, weisum
        real(TKG), allocatable :: weight(:)
        real(TKG), allocatable :: sample(:)
        real(TKG) :: mean
        call disp%show("sample = getLinSpace(1._TKG, 9._TKG, 5_IK)")
                        sample = getLinSpace(1._TKG, 9._TKG, 5_IK)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("weight = getLinSpace(1._TKG, 9._TKG, size(sample, kind = IK))")
                        weight = getLinSpace(1._TKG, 9._TKG, size(sample, kind = IK))
        call disp%show("weight")
        call disp%show( weight )
        call disp%show("call setMean(mean, sample, weight, weisum)")
                        call setMean(mean, sample, weight, weisum)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("call setVar(var, mean, sample, weight, weisum)")
                        call setVar(var, mean, sample, weight, weisum)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, sample - mean, weight, weisum)")
                        call setVar(var, sample - mean, weight, weisum)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, mean, sample, 1_IK, weight, weisum)")
                        call setVar(var, mean, sample, 1_IK, weight, weisum)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, sample - mean, 1_IK, weight, weisum)")
                        call setVar(var, sample - mean, 1_IK, weight, weisum)
        call disp%show("var")
        call disp%show( var )
        call disp%skip()
    end block

    block
        real(TKG) :: var, weisum
        real(TKG), allocatable :: weight(:)
        complex(TKG), allocatable :: sample(:)
        complex(TKG) :: mean
        call disp%show("sample = cmplx(getLinSpace(1., 9., 5_IK), -getLinSpace(1., 9., 5_IK), TKG)")
                        sample = cmplx(getLinSpace(1., 9., 5_IK), -getLinSpace(1., 9., 5_IK), TKG)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("weight = getLinSpace(1._TKG, 9._TKG, size(sample, kind = IK))")
                        weight = getLinSpace(1._TKG, 9._TKG, size(sample, kind = IK))
        call disp%show("weight")
        call disp%show( weight )
        call disp%show("call setMean(mean, sample, weight, weisum)")
                        call setMean(mean, sample, weight, weisum)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("call setVar(var, mean, sample, weight, weisum)")
                        call setVar(var, mean, sample, weight, weisum)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, sample - mean, weight, weisum)")
                        call setVar(var, sample - mean, weight, weisum)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, mean, sample, 1_IK, weight, weisum)")
                        call setVar(var, mean, sample, 1_IK, weight, weisum)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var, sample - mean, 1_IK, weight, weisum)")
                        call setVar(var, sample - mean, 1_IK, weight, weisum)
        call disp%show("var")
        call disp%show( var )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the variance of a 2-D array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG), allocatable :: var(:), mean(:), sample(:,:)
        call disp%skip()
        call disp%show("sample = getUnifRand(1._TKG, 9._TKG, 4_IK, 5_IK)")
                        sample = getUnifRand(1._TKG, 9._TKG, 4_IK, 5_IK)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("mean = [getMean(sample)]")
                        mean = [getMean(sample)]
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("call setResized(var, size(mean, 1, IK))")
                        call setResized(var, size(mean, 1, IK))
        call disp%show("call setVar(var(1), mean(1), sample)")
                        call setVar(var(1), mean(1), sample)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var(1), sample - mean(1))")
                        call setVar(var(1), sample - mean(1))
        call disp%show("var")
        call disp%show( var )
        call disp%skip()
        do dim = 1, 2
            call disp%show("dim ! The observations axis.")
            call disp%show( dim )
            call disp%show("sample = getUnifRand(1._TKG, 9._TKG, 4_IK, 5_IK)")
                            sample = getUnifRand(1._TKG, 9._TKG, 4_IK, 5_IK)
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("mean = getMean(sample, dim)")
                            mean = getMean(sample, dim)
            call disp%show("mean")
            call disp%show( mean )
            call disp%show("call setResized(var, size(mean, 1, IK))")
                            call setResized(var, size(mean, 1, IK))
            call disp%show("call setVar(var, mean, sample, dim)")
                            call setVar(var, mean, sample, dim)
            call disp%show("var")
            call disp%show( var )
            call disp%show("call setVar(var, getShifted(sample, dim, -mean), dim)")
                            call setVar(var, getShifted(sample, dim, -mean), dim)
            call disp%show("var")
            call disp%show( var )
            call disp%skip()
        end do
    end block

    block
        real(TKG), allocatable :: var(:)
        complex(TKG), allocatable :: mean(:), sample(:,:)
        call disp%skip()
        call disp%show("sample = cmplx(getUnifRand(1., 9., 4_IK, 5_IK), -getUnifRand(1., 9., 4_IK, 5_IK), TKG)")
                        sample = cmplx(getUnifRand(1., 9., 4_IK, 5_IK), -getUnifRand(1., 9., 4_IK, 5_IK), TKG)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("mean = [getMean(sample)]")
                        mean = [getMean(sample)]
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("call setResized(var, size(mean, 1, IK))")
                        call setResized(var, size(mean, 1, IK))
        call disp%show("call setVar(var(1), mean(1), sample)")
                        call setVar(var(1), mean(1), sample)
        call disp%show("var")
        call disp%show( var )
        call disp%show("call setVar(var(1), sample - mean(1))")
                        call setVar(var(1), sample - mean(1))
        call disp%show("var")
        call disp%show( var )
        call disp%skip()
        do dim = 1, 2
            call disp%show("dim ! The observations axis.")
            call disp%show( dim )
            call disp%show("sample = getUnifRand(1._TKG, 9._TKG, 4_IK, 5_IK)")
                            sample = getUnifRand(1._TKG, 9._TKG, 4_IK, 5_IK)
            call disp%show("sample")
            call disp%show( sample )
            call disp%show("mean = getMean(sample, dim)")
                            mean = getMean(sample, dim)
            call disp%show("mean")
            call disp%show( mean )
            call disp%show("call setResized(var, size(mean, 1, IK))")
                            call setResized(var, size(mean, 1, IK))
            call disp%show("call setVar(var, mean, sample, dim)")
                            call setVar(var, mean, sample, dim)
            call disp%show("var")
            call disp%show( var )
            call disp%show("call setVar(var, getShifted(sample, dim, -mean), dim)")
                            call setVar(var, getShifted(sample, dim, -mean), dim)
            call disp%show("var")
            call disp%show( var )
            call disp%skip()
        end do
    end block

end program example