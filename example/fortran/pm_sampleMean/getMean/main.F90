program example

    use pm_kind, only: SK, IK
    use pm_kind, only: TKG => RKS ! All other real types are also supported.
    use pm_sampleMean, only: getMean
    use pm_arrayRange, only: getRange
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    integer(IK) :: idim, ndim, nsam
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the mean of a 1-D array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG) :: mean
        real(TKG), allocatable :: sample(:)
        call disp%skip()
        call disp%show("nsam = getUnifRand(1, 5)")
                        nsam = getUnifRand(1, 5)
        call disp%show("sample = getUnifRand(0., 1., nsam)")
                        sample = getUnifRand(0., 1., nsam)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getMean(sample)")
        call disp%show( getMean(sample) )
        call disp%show("mean = getMean(sample, dim = 1_IK)")
                        mean = getMean(sample, dim = 1_IK)
        call disp%show("mean")
        call disp%show( mean )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the mean of a 2-D array along a specific dimension.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG), allocatable :: mean(:)
        real(TKG), allocatable :: sample(:,:)
        call disp%skip()
        call disp%show("ndim = getUnifRand(1, 3); nsam = getUnifRand(1, 5)")
                        ndim = getUnifRand(1, 3); nsam = getUnifRand(1, 5)
        call disp%show("sample = getUnifRand(0., 1., ndim, nsam)")
                        sample = getUnifRand(0., 1., ndim, nsam)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("getMean(sample)")
        call disp%show( getMean(sample) )
        call disp%show("mean = getMean(sample, dim = 2_IK)")
                        mean = getMean(sample, dim = 2_IK)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("mean = getMean(transpose(sample), dim = 1_IK)")
                        mean = getMean(transpose(sample), dim = 1_IK)
        call disp%show("mean")
        call disp%show( mean )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the mean of a 1-D weighted array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG) :: mean
        real(TKG), allocatable :: sample(:)
        integer(IK), allocatable :: weight(:)
        call disp%skip()
        call disp%show("nsam = getUnifRand(1, 5)")
                        nsam = getUnifRand(1, 5)
        call disp%show("sample = getUnifRand(0., 1., nsam)")
                        sample = getUnifRand(0., 1., nsam)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("weight = getUnifRand(1, 9, nsam)")
                        weight = getUnifRand(1, 9, nsam)
        call disp%show("weight")
        call disp%show( weight )
        call disp%show("mean = getMean(sample)")
                        mean = getMean(sample)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("mean = getMean(sample, weight)")
                        mean = getMean(sample, weight)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("mean = getMean(sample, 1_IK, weight)")
                        mean = getMean(sample, 1_IK, weight)
        call disp%show("mean")
        call disp%show( mean )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the mean of a 1-D weighted array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG) :: mean
        real(TKG), allocatable :: sample(:)
        real(TKG), allocatable :: weight(:)
        call disp%skip()
        call disp%show("nsam = getUnifRand(1, 5)")
                        nsam = getUnifRand(1, 5)
        call disp%show("sample = getUnifRand(0., 1., nsam)")
                        sample = getUnifRand(0., 1., nsam)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("weight = getUnifRand(0., 1., nsam)")
                        weight = getUnifRand(0., 1., nsam)
        call disp%show("weight")
        call disp%show( weight )
        call disp%show("mean = getMean(sample)")
                        mean = getMean(sample)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("mean = getMean(sample, weight)")
                        mean = getMean(sample, weight)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("mean = getMean(sample, 1_IK, weight)")
                        mean = getMean(sample, 1_IK, weight)
        call disp%show("mean")
        call disp%show( mean )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the mean of a 2-D weighted array along a specific dimension.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG), allocatable :: mean(:)
        real(TKG), allocatable :: sample(:,:)
        real(TKG), allocatable :: weight(:)
        call disp%skip()
        call disp%show("ndim = getUnifRand(1, 3); nsam = getUnifRand(1, 5)")
                        ndim = getUnifRand(1, 3); nsam = getUnifRand(1, 5)
        call disp%show("sample = getUnifRand(0., 1., ndim, nsam)")
                        sample = getUnifRand(0., 1., ndim, nsam)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("weight = getUnifRand(0., 1., nsam)")
                        weight = getUnifRand(0., 1., nsam)
        call disp%show("weight")
        call disp%show( weight )
        call disp%show("getMean(sample, [(weight, idim = 1, ndim)])")
        call disp%show( getMean(sample, [(weight, idim = 1, ndim)]) )
        call disp%show("mean = getMean(sample, 2_IK, weight)")
                        mean = getMean(sample, 2_IK, weight)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("mean = getMean(transpose(sample), 1_IK, weight)")
                        mean = getMean(transpose(sample), 1_IK, weight)
        call disp%show("mean")
        call disp%show( mean )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Compute the mean of a multidimensional array by associating it with a 1D pointer.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer(IK) :: nslice
        real(TKG), allocatable :: mean
        real(TKG), allocatable, target :: sample(:,:,:)
        real(TKG), pointer :: samptr(:)
        call disp%skip()
        call disp%show("ndim = getUnifRand(1, 2); nsam = getUnifRand(1, 3); nslice = getUnifRand(1, 4)")
                        ndim = getUnifRand(1, 2); nsam = getUnifRand(1, 3); nslice = getUnifRand(1, 4)
        call disp%show("sample = getUnifRand(0., 1., ndim, nsam, nslice)")
                        sample = getUnifRand(0., 1., ndim, nsam, nslice)
        call disp%show("sample")
        call disp%show( sample )
        call disp%show("shape(sample)")
        call disp%show( shape(sample) )
        call disp%show("samptr(1:product(shape(sample))) => sample")
                        samptr(1:product(shape(sample))) => sample
        call disp%show("mean = getMean(samptr)")
                        mean = getMean(samptr)
        call disp%show("nullify(samptr)")
                        nullify(samptr)
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("mean = getMean(reshape(sample, [product(shape(sample))]))")
                        mean = getMean(reshape(sample, [product(shape(sample))]))
        call disp%show("mean")
        call disp%show( mean )
        call disp%skip()
    end block

end program example