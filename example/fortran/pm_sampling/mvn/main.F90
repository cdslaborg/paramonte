module logfunc
    use pm_kind, only: IK, RKG => RKS ! All real kinds are supported.
    implicit none
    integer(IK) , parameter :: NDIM = 4_IK
    real(RKG)   , parameter :: MEAN(NDIM) = [-6, -2, 2, 6] ! The mean vector of the MVN
    real(RKG)   , parameter :: LOG_INVERSE_SQRT_TWO_PI = log(1./sqrt(2.*acos(-1.))) ! MVN distribution coefficient: log(1/sqrt(2*Pi)^ndim).
    real(RKG)   , parameter :: COVMAT(NDIM, NDIM) = reshape([ +1.0, +0.5, +0.5, +0.5 &
                                                            , +0.5, +1.0, +0.5, +0.5 &
                                                            , +0.5, +0.5, +1.0, +0.5 &
                                                            , +0.5, +0.5, +0.5, +1.0 ], shape(COVMAT)) ! The covariance matrix of the MVN
    real(RKG)   , parameter :: INVCOV(NDIM, NDIM) = reshape([ +1.6, -0.4, -0.4, -0.4 &
                                                            , -0.4, +1.6, -0.4, -0.4 &
                                                            , -0.4, -0.4, +1.6, -0.4 &
                                                            , -0.4, -0.4, -0.4, +1.6 ], shape(INVCOV)) ! The inverse covariance matrix of the MVN
    real(RKG)   , parameter :: MVN_COEF = NDIM * LOG_INVERSE_SQRT_TWO_PI + 0.581575404902840 ! log(sqrt(det(INVCOV)))
contains
    function getLogFunc(state) result(logFunc)
        ! Return the negative natural logarithm of MVN distribution 
        ! evaluated at the input vector `state` of length `ndim`.
        real(RKG), intent(in), contiguous :: state(:)
        real(RKG) :: stateNormed(size(state))
        real(RKG) :: logFunc
        stateNormed = state - MEAN
        logFunc = MVN_COEF - 0.5_RKG * (dot_product(stateNormed, matmul(INVCOV, stateNormed)))
    end function
end module logfunc

program example
    use logfunc, only: RKG, NDIM, getLogFunc
    use pm_sampling, only: getErrSampling, paradram_type, err_type
    type(err_type) :: err
    err = getErrSampling(paradram_type(inputFile = 'input.nml'), getLogFunc, NDIM)
    if (err%occurred) error stop "sampler failed: "//err%msg
end