program example

    use pm_kind, only: SK, IK, LK, RK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand
    use pm_arrayReverse, only: getReversed
    use pm_distCov, only: getCovRand
    use pm_arrayRange, only: getRange
    use pm_sampleCor, only: getRho
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
    real(RK), allocatable :: rho(:,:), rweight(:)
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Spearman correlation matrix for a pair of character-sequence time series.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: SKC => SK ! All other kinds are also supported.
        real(RK) :: rho
        integer(IK) :: nsam
        character(:,SKC), allocatable :: x, y
        call disp%show("nsam = 10")
                        nsam = 10
        call disp%show("x = getUnifRand(repeat('A', nsam), repeat('Z', nsam))")
                        x = getUnifRand(repeat('A', nsam), repeat('Z', nsam))
        call disp%show("x")
        call disp%show( x , deliml = SK_'''' )
        call disp%show("y = getUnifRand(repeat('A', nsam), repeat('Z', nsam))")
                        y = getUnifRand(repeat('A', nsam), repeat('Z', nsam))
        call disp%show("y")
        call disp%show( y , deliml = SK_'''' )
        call disp%show("rho = getRho(x, y)")
                        rho = getRho(x, y)
        call disp%show("rho")
        call disp%show( rho )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Spearman correlation matrix for a pair of character-valued time series.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: SKC => SK ! All other kinds are also supported.
        integer(IK) :: ndim, nsam
        character(2,SKC), allocatable :: sample(:,:)
        call disp%show("ndim = 2; nsam = 10")
                        ndim = 2; nsam = 10
        call disp%show("call setResized(sample, [ndim, nsam])")
                        call setResized(sample, [ndim, nsam])
        call disp%show("call setUnifRand(sample, SKC_'AA', SKC_'ZZ')")
                        call setUnifRand(sample, SKC_'AA', SKC_'ZZ')
        call disp%show("sample")
        call disp%show( sample , deliml = SK_'''' )
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("rho = getRho(sample, dim = 2_IK)")
                        rho = getRho(sample, dim = 2_IK)
        call disp%show("rho")
        call disp%show( rho )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("rho = getRho(transpose(sample), dim = 1_IK)")
                        rho = getRho(transpose(sample), dim = 1_IK)
        call disp%show("rho")
        call disp%show( rho )
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho(1,1) = getRho(sample(1,:), sample(2,:))")
                        rho(1,1) = getRho(sample(1,:), sample(2,:))
        call disp%show("rho(1,1)")
        call disp%show( rho(1,1) )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Spearman correlation matrix for a pair of integer-valued time series.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: IKC => IK ! All other kinds are also supported.
        integer(IK) :: ndim, nsam
        integer(IKC), allocatable :: sample(:,:)
        call disp%show("ndim = 2; nsam = 10")
                        ndim = 2; nsam = 10
        call disp%show("sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("rho = getRho(sample, dim = 2_IK)")
                        rho = getRho(sample, dim = 2_IK)
        call disp%show("rho")
        call disp%show( rho )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("rho = getRho(transpose(sample), dim = 1_IK)")
                        rho = getRho(transpose(sample), dim = 1_IK)
        call disp%show("rho")
        call disp%show( rho )
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho(1,1) = getRho(sample(1,:), sample(2,:))")
                        rho(1,1) = getRho(sample(1,:), sample(2,:))
        call disp%show("rho(1,1)")
        call disp%show( rho(1,1) )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the Spearman correlation matrix for a pair of real-valued time series.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKC => RKS ! All other real types are also supported.
        integer(IK) :: ndim, nsam
        real(RKC), allocatable :: sample(:,:)
        format = getFormat(mold = [0._RKC], ed = SK_"es", signed = .true._LK)
        call disp%show("ndim = 2; nsam = 10")
                        ndim = 2; nsam = 10
        call disp%show("sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample , format = format )
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("rho = getRho(sample, dim = 2_IK)")
                        rho = getRho(sample, dim = 2_IK)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("rho = getRho(transpose(sample), dim = 1_IK)")
                        rho = getRho(transpose(sample), dim = 1_IK)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho(1,1) = getRho(sample(1,:), sample(2,:))")
                        rho(1,1) = getRho(sample(1,:), sample(2,:))
        call disp%show("rho(1,1)")
        call disp%show( rho(1,1) , format = format )
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
        use pm_kind, only: RKC => RKS ! All other real types are also supported.
        integer(IK) :: ndim, nsam
        integer(IK), allocatable :: iweight(:)
        real(RKC), allocatable :: sample(:,:)
        format = getFormat(mold = [0._RKC], ed = SK_"es", signed = .true._LK)
        call disp%show("ndim = 2; nsam = 10")
                        ndim = 2; nsam = 10
        call disp%show("sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])")
                        sample = reshape(getUnifRand(1, 20, ndim * nsam), shape = [ndim, nsam], order = [2, 1])
        call disp%show("sample")
        call disp%show( sample , format = format )
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
        call disp%show("rho = getRho(sample, 2_IK, iweight)")
                        rho = getRho(sample, 2_IK, iweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("rho = getRho(transpose(sample), 1_IK, iweight)")
                        rho = getRho(transpose(sample), 1_IK, iweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho(1,1) = getRho(sample(1,:), sample(2,:), iweight)")
                        rho(1,1) = getRho(sample(1,:), sample(2,:), iweight)
        call disp%show("rho(1,1)")
        call disp%show( rho(1,1) , format = format )
        call disp%skip()

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Compute the correlation matrix with    real weights.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("rho = getRho(sample, 2_IK, rweight)")
                        rho = getRho(sample, 2_IK, rweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%skip()
        call disp%show("Compute the sample correlation along the first dimension.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%show("rho = getRho(transpose(sample), 1_IK, rweight)")
                        rho = getRho(transpose(sample), 1_IK, rweight)
        call disp%show("rho")
        call disp%show( rho , format = format )
        call disp%skip()
        call disp%show("rho = getFilled(0., ndim, ndim)")
                        rho = getFilled(0., ndim, ndim)
        call disp%skip()
        call disp%show("Compute the full sample correlation for a pair of time series.", deliml = SK_'''')
        call disp%skip()
        call disp%show("rho(1,1) = getRho(sample(1,:), sample(2,:), rweight)")
                        rho(1,1) = getRho(sample(1,:), sample(2,:), rweight)
        call disp%show("rho(1,1)")
        call disp%show( rho(1,1) , format = format )
        call disp%skip()
    end block

end program example