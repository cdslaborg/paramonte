module logfunc
    use pm_kind, only: IK, RKG => RKD ! all processor kinds are supported.
    integer(IK), parameter :: NDIM = 2
contains
    recursive function getLogFunc(state) result(logFunc)
        ! Return the negative natural logarithm of Himmelblau function evaluated at the input vector state.
        real(RKG), intent(in), contiguous :: state(:)
        real(RKG) :: logFunc
        if (size(state) /= 2) error stop "The input state vector size must be 2."
        logFunc = -log((state(1)**2 + state(2) - 11)**2 + (state(1) + state(2)**2 - 7)**2 + 0.1_RKG)
    end function
end module logfunc

program example
    use logfunc, only: NDIM, getLogFunc
    use pm_sampling, only: getErrSampling, paradram_type, err_type
    type(paradram_type) :: sampler
    type(err_type) :: err
    ! \bug
    ! gfortran bug requires setting sampler components explicitly.
    sampler = paradram_type()
    sampler%inputFile = "input.nml"
    sampler%outputFileName = "./out/himmelblau"
    err = getErrSampling(sampler, getLogFunc, NDIM)
    if (err%occurred) error stop "sampler failed: "//err%msg
end program example