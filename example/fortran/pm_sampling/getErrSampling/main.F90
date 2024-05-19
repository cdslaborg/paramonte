program example

    use pm_kind, only: RKG => RK
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_sampling, only: getErrSampling, paradram_type
    use pm_matrixInit, only: getMatInit, uppLowDia
    use pm_arrayFill, only: getFilled
    use pm_sysPath, only: glob
    use pm_err, only: err_type

    implicit none

    type(err_type) :: err
    type(display_type) :: disp
    integer(IK), parameter :: NDIM = 4
    real(RKG), parameter :: MEAN(*) = 4 * [real(RKG) :: -3, -1, 1, 3]

    disp = display_type(file = SK_"main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Simple minimally-informed simulation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        call disp%skip()
        call disp%show("NDIM")
        call disp%show( NDIM )
        call disp%show("err = getErrSampling(paradram_type(parallelismMpiFinalizeEnabled = .false._LK), getLogFunc, NDIM)")
                        err = getErrSampling(paradram_type(parallelismMpiFinalizeEnabled = .false._LK), getLogFunc, NDIM)
        call disp%show("if (err%occurred) error stop err%msg")
                        if (err%occurred) error stop err%msg
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Optionally-informed simulation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        type(paradram_type) :: sampler
        call disp%skip()
        call disp%show("NDIM")
        call disp%show( NDIM )
       ! GNU gfortran yields a segfault with the following default constructor.
       !call disp%show("sampler = paradram_type(parallelismMpiFinalizeEnabled = .false._LK, outputFileName = 'mvn', outputChainSize = 100000, outputStatus = 'retry')")
       !                sampler = paradram_type(parallelismMpiFinalizeEnabled = .false._LK, outputFileName = 'mvn', outputChainSize = 100000, outputStatus = 'retry')
        call disp%show("sampler = paradram_type()")
                        sampler = paradram_type()
        call disp%show("sampler%parallelismMpiFinalizeEnabled = .false._LK; sampler%outputFileName = 'mvn'; sampler%outputChainSize = 100000; sampler%outputStatus = 'retry'")
                        sampler%parallelismMpiFinalizeEnabled = .false._LK; sampler%outputFileName = 'mvn'; sampler%outputChainSize = 100000; sampler%outputStatus = 'retry'
        call disp%show("err = getErrSampling(sampler, getLogFunc, NDIM)")
                        err = getErrSampling(sampler, getLogFunc, NDIM)
        call disp%show("if (err%occurred) error stop err%msg")
                        if (err%occurred) error stop err%msg
        call disp%show("glob(pattern = sampler%outputFileName)")
        call disp%show( glob(pattern = sampler%outputFileName) )
        call disp%skip()
    end block

contains

    recursive function getLogFunc(state) result(logFunc)
        use pm_distMultiNorm, only: getMultiNormLogPDF
        real(RKG), intent(in), contiguous :: state(:)
        real(RKG) :: logFunc
        integer :: i
        logFunc = getMultiNormLogPDF(state, mean = real(MEAN, RKG))
#if     OMP_ENABLED
        ! We are done. But kill some time for an illustration of parallel sampling.
        do i = 1, 500
            logFunc = logFunc + getMultiNormLogPDF(state, mean = 4 * [real(RKG) :: -3, -1, 1, 3]) * merge(1, -1, mod(i, 2) == 0)
        end do
#endif
    end function

    !subroutine setLogFunc(logFuncState)
    !    use pm_distMultiNorm, only: getMultiNormLogPDF
    !    real(RKG), intent(inout), contiguous :: logFuncState(0:,:)
    !    integer(IK) :: ithread
    !    do ithread = 1, size(logFuncState, 2)
    !        logFuncState(0, ithread) = getMultiNormLogPDF(logFuncState(1:, ithread), mean = real(MEAN, RKG))
    !    end do
    !end subroutine

end program example