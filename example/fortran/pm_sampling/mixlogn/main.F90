!>  \brief
!>  This code fits the BATSE T90 duration data with the mathematical model of Moharana and Piran (2018).
module logfunc

    use pm_kind, only: SK, IK, LK, RKC => RKD ! All real kinds are supported.
    use pm_mathLogAddExp, only: getLogAddExp
    use pm_mathLogSubExp, only: getLogSubExp
    use pm_distNorm, only: getNormLogPDF
    use pm_mathMinMax, only: getMinMax
    use pm_io, only: getErrTableWrite
    use pm_io, only: getErrTableRead
    use pm_io, only: display_type
    use pm_err, only: err_type
    use pm_io, only: disp

    implicit none

    real(RKC), parameter :: LARGE = sqrt(huge(0._RKC))

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! model

    type merger_type
        real(RKC) :: avg, std
    end type

    type collapsar_type
        real(RKC) :: avg, std
    end type

    type logmixfac_type
        real(RKC) :: merger, collapsar
    end type

    type param_type
        type(merger_type) :: merger
        type(collapsar_type) :: collapsar
        type(logmixfac_type) :: logmixfac
    end type

    ! data

    type vmm_type
        real(RKC) :: min, max
        real(RKC), allocatable :: val(:)
    end type

    type, extends(vmm_type) :: dur_type
        type(vmm_type) :: log
    end type

    type(dur_type) :: mc_dur ! log(T90) to fit.

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Read the data to fit.
    subroutine setInitModel()
        real(RKC), allocatable :: table(:,:)
        character(:, SK), allocatable :: file
        character(255, SK) :: iomsg
        integer(IK) :: stat
        file = "data.csv"
        stat = getErrTableRead(file, table, sep = SK_",", roff = 1_IK, iomsg = iomsg)
        if (stat /= 0) error stop "Failed to read data: "//trim(iomsg)
        mc_dur%val = table(:,4)
        mc_dur%log%val = log(mc_dur%val)
        mc_dur%min = minval(mc_dur%val, dim = 1)
        mc_dur%max = maxval(mc_dur%val, dim = 1)
        mc_dur%log%min = minval(mc_dur%log%val, dim = 1)
        mc_dur%log%max = maxval(mc_dur%log%val, dim = 1)
        call disp%skip()
        call disp%show("mc_dur%min")
        call disp%show( mc_dur%min )
        call disp%show("mc_dur%max")
        call disp%show( mc_dur%max )
        call disp%show("mc_dur%log%min")
        call disp%show( mc_dur%log%min )
        call disp%show("mc_dur%log%max")
        call disp%show( mc_dur%log%max )
        call disp%skip()
        stat = getErrTableWrite(file//SK_".out", mc_dur%log%val, header = SK_"logdur", sep = SK_",")
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental function getLogProbMix(logdur, param) result(logProb)
        real(RKC), intent(in) :: logdur
        type(param_type), intent(in) :: param
        real(RKC) :: logProb, logProbMerger, logProbCollapsar
        logProbMerger = param%logmixfac%merger + getNormLogPDF(logdur, mu = param%merger%avg, sigma = param%merger%std)
        logProbCollapsar = param%logmixfac%collapsar + getNormLogPDF(logdur, mu = param%collapsar%avg, sigma = param%collapsar%std)
        logProb = getLogAddExp(getMinMax(logProbMerger, logProbCollapsar))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogLike(state) result(logLike)

        real(RKC), intent(in), contiguous :: state(:)
        real(RKC) :: logLike, mixfrac
        type(param_type) :: param

        param%merger%avg = state(1)
        param%merger%std = exp(state(2))

        param%collapsar%avg = state(3)
        param%collapsar%std = exp(state(4))

        mixfrac = .5_RKC * (tanh(state(5)) + 1)
        param%logmixfac%merger = log(mixfrac)
        param%logmixfac%collapsar = log(1 - mixfrac)

        if (param%merger%avg < param%collapsar%avg) then
            logLike = sum(getLogProbMix(mc_dur%log%val, param))
        else
            logLike = -LARGE
        end if

    end function

end module logfunc

program example

    use logfunc
    use pm_sampling, only: getErrSampling, paradram_type
    type(paradram_type) :: sampler
    type(err_type) :: err

    call setInitModel()

    !call disp%show(getLogLike([real(RKC) :: log(.5), -2, log(10.), 2, 1, 0.5]))

    sampler = paradram_type()
    sampler%outputStatus ="retry"
    sampler%outputChainSize = 30000
    !sampler%proposalScaleFactor = "0.1 * gelman"
    sampler%outputFileName = "./out/mixlogn"
    sampler%domainAxisName= [ "log_merger_avg   " &
                            , "log_merger_std   " &
                            , "log_collapsar_avg" &
                            , "log_collapsar_std" &
                            , "atanh_2mixfracm1 " &
                            ]
    sampler%proposalStart = [ -.5 & ! log(merger%avg)
                            , +0. & ! log(merger%std)
                            , 3.5 & ! log(collapsar%avg)
                            , +0. & ! log(collapsar%std)
                            , +.3 & ! atanh(mixfrac * 2 - 1)
                            ]
    sampler%proposalCovMat = reshape(   [  1.5E-2,  5.4E-3,  3.8E-3, -2.2E-3,  4.0E-3 &
                                        ,  5.4E-3,  3.2E-3,  1.5E-3, -8.9E-4,  1.6E-3 &
                                        ,  3.8E-3,  1.5E-3,  1.9E-3, -8.1E-4,  1.2E-3 &
                                        , -2.2E-3, -8.9E-4, -8.1E-4,  9.8E-4, -7.7E-4 &
                                        ,  4.0E-3,  1.6E-3,  1.2E-3, -7.7E-4,  1.8E-3 ], shape = [size(sampler%proposalStart), size(sampler%proposalStart)])

    err = getErrSampling(sampler, getLogLike, size(sampler%proposalStart, 1, IK))
    if (err%occurred) error stop "sampler failed: "//err%msg

    block
        use pm_sampleMean, only: getMean
        integer(IK) :: stat, offset = 1
        real(RKC), allocatable :: table(:,:), param(:)
        stat = getErrTableRead(sampler%outputFileName//"_run1_pid1_sample.txt", table, roff = 1_IK)
        if (stat /= 0) error stop "Failed to read output sample."
        param = getMean(table, dim = 1_IK)
        call disp%show("[merger%avg, merger%std, collapsar%avg, collapsar%std, merger_fraction]")
        call disp%show( [exp(param(offset + 1)), param(offset + 2), exp(param(offset + 3)), param(offset + 4 : offset + 5), .5_RKC * (tanh(param(offset + 6)) + 1)] )
    end block

    block
        integer(IK) :: stat, offset = 7
        real(RKC), allocatable :: table(:,:)
        real(RKC) :: maxloglike, bic
        stat = getErrTableRead(sampler%outputFileName//"_run1_pid1_chain.txt", table, roff = 1_IK)
        if (stat /= 0) error stop "Failed to read output sample."
        maxloglike = maxval(table(:, offset), dim = 1)
        bic = size(sampler%proposalStart) * log(real(size(mc_dur%log%val), RKC)) - 2 * maxloglike
        call disp%show("[maxloglike, bic]")
        call disp%show( [maxloglike, bic] )
        ! result:
        ! [merger%avg, merger%std, collapsar%avg, collapsar%std, merger_fraction]
        ! +0.81737715750468454, +0.37408303769528523, +35.039670487904210, -0.12260201209700368E-1, -0.33565399808200586, +1.0000000000000000
        ! [maxloglike, bic]
        ! -4150.6243864147127, +8339.3547486434745
    end block

end