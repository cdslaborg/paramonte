program example

    use pm_kind, only: SK, IK, LK, RK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_arrayReverse, only: getReversed
    use pm_distCov, only: getCovRand
    use pm_arrayRange, only: getRange
    use pm_sampleCor, only: setRho
    use pm_sampleCor, only: uppDia
    use pm_sampleCor, only: lowDia
    use pm_sampleCor, only: upp
    use pm_sampleCor, only: low
    use pm_sampleMean, only: getMean
    use pm_sampleMean, only: setMean
    use pm_sampleShift, only: getShifted
    use pm_arraySpace, only: getLinSpace
    use pm_arrayResize, only: setResized
    use pm_arrayFill, only: getFilled
    use pm_io, only: getFormat

    implicit none

    type(display_type) :: disp
    integer(IK) :: itry, ntry = 10
    character(:), allocatable :: format
    real(RK), allocatable :: rweight(:), rho(:,:), frank(:,:), frankt(:,:)
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Spearman correlation matrix for a pair of time series.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKG => RKS ! All other real types are also supported.
        integer(IK) :: ndim, nsam
        real(RKG), allocatable :: sample(:,:)
        format = getFormat(mold = [0._RKG], ed = SK_"es", signed = .true._LK)
        call disp%show("ndim = 2; nsam = 10")
                        ndim = 2; nsam = 10
        call disp%show("sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample , format = format )
        call disp%show("call setResized(frank, shape(sample, IK))")
                        call setResized(frank, shape(sample, IK))
        call disp%show("call setResized(frankt, getReversed(shape(sample, IK)))")
                        call setResized(frankt, getReversed(shape(sample, IK)))
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, uppDia, frank, sample, dim = 2_IK)")
                        call setRho(rho, uppDia, frank, sample, dim = 2_IK)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, lowDia, frank, sample, dim = 2_IK)")
                        call setRho(rho, lowDia, frank, sample, dim = 2_IK)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%show("frank")
        call disp%show( frank , format = format )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, uppDia, frankt, transpose(sample), dim = 1_IK)")
                        call setRho(rho, uppDia, frankt, transpose(sample), dim = 1_IK)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%show("frankt")
        call disp%show( frankt , format = format )
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, lowDia, frankt, transpose(sample), dim = 1_IK)")
                        call setRho(rho, lowDia, frankt, transpose(sample), dim = 1_IK)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%show("frankt")
        call disp%show( frankt , format = format )
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("call setRho(rho(1,1), frank(1,:), frank(2,:), sample(1,:), sample(2,:))")
                        call setRho(rho(1,1), frank(1,:), frank(2,:), sample(1,:), sample(2,:))
        call disp%show("rho(1,1)")
        call disp%show( rho(1,1) , format = format )
        call disp%show("frank")
        call disp%show( frank , format = format )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Spearman correlation matrix for a weighted pair of time series.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_arrayVerbose, only: getVerbose
        use pm_kind, only: RKG => RKS ! All other real types are also supported.
        integer(IK) :: ndim, nsam
        integer(IK), allocatable :: iweight(:)
        real(RKG), allocatable :: sample(:,:)
        format = getFormat(mold = [0._RKG], ed = SK_"es", signed = .true._LK)
        call disp%show("ndim = 2; nsam = 10")
                        ndim = 2; nsam = 10
        call disp%show("sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample , format = format )
        call disp%show("call setResized(frank, shape(sample, IK))")
                        call setResized(frank, shape(sample, IK))
        call disp%show("call setResized(frankt, getReversed(shape(sample, IK)))")
                        call setResized(frankt, getReversed(shape(sample, IK)))
        call disp%show("iweight = getUnifRand(1, 10, nsam) ! integer-valued weights.")
                        iweight = getUnifRand(1, 10, nsam) ! integer-valued weights.
        call disp%show("iweight")
        call disp%show( iweight )
        call disp%show("rweight = iweight ! or real-valued weights.")
                        rweight = iweight ! or real-valued weights.
        call disp%show("iweight")
        call disp%show( iweight )

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Compute the correlation matrix with integer weights.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, uppDia, frank, sample, 2_IK, iweight)")
                        call setRho(rho, uppDia, frank, sample, 2_IK, iweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%show("frank")
        call disp%show( frank , format = format )
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, lowDia, frank, sample, 2_IK, iweight)")
                        call setRho(rho, lowDia, frank, sample, 2_IK, iweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%show("frank")
        call disp%show( frank , format = format )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, uppDia, frankt, transpose(sample), 1_IK, iweight)")
                        call setRho(rho, uppDia, frankt, transpose(sample), 1_IK, iweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, lowDia, frankt, transpose(sample), 1_IK, iweight)")
                        call setRho(rho, lowDia, frankt, transpose(sample), 1_IK, iweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%show("frankt")
        call disp%show( frankt , format = format )
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("call setRho(rho(1,1), frank(1,:), frank(2,:), sample(1,:), sample(2,:), iweight)")
                        call setRho(rho(1,1), frank(1,:), frank(2,:), sample(1,:), sample(2,:), iweight)
        call disp%show("rho(1,1)")
        call disp%show( rho(1,1) , format = format )
        call disp%show("frank")
        call disp%show( frank , format = format )
        call disp%skip()

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Compute the correlation matrix with    real weights.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, uppDia, frank, sample, 2_IK, rweight)")
                        call setRho(rho, uppDia, frank, sample, 2_IK, rweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%show("frank")
        call disp%show( frank , format = format )
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, lowDia, frank, sample, 2_IK, rweight)")
                        call setRho(rho, lowDia, frank, sample, 2_IK, rweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%show("frank")
        call disp%show( frank , format = format )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, uppDia, frankt, transpose(sample), 1_IK, rweight)")
                        call setRho(rho, uppDia, frankt, transpose(sample), 1_IK, rweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%show("frankt")
        call disp%show( frankt , format = format )
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("call setRho(rho, lowDia, frankt, transpose(sample), 1_IK, rweight)")
                        call setRho(rho, lowDia, frankt, transpose(sample), 1_IK, rweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%show("frankt")
        call disp%show( frankt , format = format )
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("call setRho(rho(1,1), frank(1,:), frank(2,:), sample(1,:), sample(2,:), rweight)")
                        call setRho(rho(1,1), frank(1,:), frank(2,:), sample(1,:), sample(2,:), rweight)
        call disp%show("rho(1,1)")
        call disp%show( rho(1,1) , format = format )
        call disp%show("frank")
        call disp%show( frank , format = format )
        call disp%skip()
    end block

end program example