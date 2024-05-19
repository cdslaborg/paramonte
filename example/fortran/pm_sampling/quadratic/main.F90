module logfunc
    use pm_kind, only: SK, IK, RKG => RKD ! all processor kinds are supported.
    integer(IK), parameter :: NDIM = 2
contains
    recursive function getLogFunc(state) result(logFunc)
        real(RKG) :: logFunc
        real(RKG), intent(in), contiguous :: state(:)
        if (size(state) /= 2) error stop "The input state vector size must be 2."
        logFunc = -sum(state**2)
    end function
end module logfunc

program example
    use logfunc, only: SK, IK, NDIM, getLogFunc
    use pm_sampling, only: getErrSampling, paradram_type, err_type
    type(paradram_type) :: sampler
    type(err_type) :: err
    sampler = paradram_type()
    sampler%randomSeed = 284651_IK
    sampler%outputFileName = SK_"./out/"
    sampler%parallelismNumThread = 4_IK
    sampler%proposalScale = SK_"4 * gelman"
    sampler%outputChainSize = max(50_IK, ndim + 1_IK)
    err = getErrSampling(sampler, getLogFunc, NDIM)
    if (err%occurred) error stop "sampler failed: "//err%msg
end program example