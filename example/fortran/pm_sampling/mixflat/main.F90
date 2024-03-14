!>  \brief
!>  This code fits the BATSE T90 duration data with the mathematical model of Moharana and Piran (2018).
module logfunc

    use pm_kind, only: SK, IK, LK, RKC => RKD ! All real kinds are supported.
    use pm_mathLogAddExp, only: getLogAddExp
    use pm_mathLogSubExp, only: getLogSubExp
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
        real(RKC) :: tb, logtb, alpha, alphap1
    end type

    type collapsar_type
        real(RKC) :: tb, logtb, alpha, alphap1, beta
    end type

    type logmixfac_type
        real(RKC) :: merger, collapsar
    end type

    type param_type
        type(merger_type) :: merger
        type(collapsar_type) :: collapsar
        type(logmixfac_type) :: logmixfac
    end type

    ! log of model integrals.

    type logmint_type
        real(RKC) :: merger, collapsar
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

    !>  \brief
    !>  Compute the merger model density.
    pure elemental function getLogDenMerger(logdur, merger) result(logDen)
        type(merger_type), intent(in) :: merger
        real(RKC), intent(in) :: logdur
        real(RKC) :: logDen
        if (logdur < merger%logtb) then
            logDen = logdur
        else
            logDen = merger%alphap1 * logdur - merger%alpha * merger%logtb
        end if
    end function

    !>  \brief
    !>  Compute the merger model integral.
    function getLogModelIntMerger(merger) result(logModelInt)
        real(RKC) :: logModelInt, logModelInt1, logModelInt2, small, large !, logDenom
        type(merger_type), intent(in) :: merger
        logModelInt1 = log(merger%tb - mc_dur%min)
        if (merger%alphap1 < 0._RKC) then
            logModelInt2 = -log(-merger%alphap1)
            small = merger%alphap1 * mc_dur%log%max
            large = merger%alphap1 * merger%logtb
        elseif (0._RKC < merger%alphap1) then
            logModelInt2 = -log(merger%alphap1)
            small = merger%alphap1 * merger%logtb
            large = merger%alphap1 * mc_dur%log%max
        else
            logModelInt2 = 0._RKC
            small = merger%logtb
            large = mc_dur%log%max
        end if
        logModelInt2 = logModelInt2 - merger%alpha * merger%logtb + getLogSubExp(smaller = small, larger = large)
        logModelInt = getLogAddExp(getMinMax(logModelInt1, logModelInt2))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Compute the collapsar model integral.
    pure elemental function getLogDenCollapsar(dur, logdur, collapsar) result(logDen)
        type(collapsar_type), intent(in) :: collapsar
        real(RKC), intent(in) :: dur, logdur
        real(RKC) :: logDen
        if (logdur < collapsar%logtb) then
            logDen = logdur
        else
            logDen = collapsar%alphap1 * logdur - collapsar%alpha * collapsar%logtb - collapsar%beta * (dur - collapsar%tb)
        end if
    end function

    !>  \brief
    !>  Compute the collapsar model integral.
    function getLogModelIntCollapsar(collapsar) result(logModelInt)
        use pm_quadPack, only: isFailedQuad, getQuadErr, weps
        use pm_mathLogAddExp, only: getLogAddExp
        use pm_mathMinMax, only: getMinMax
        logical(LK) :: failed
        character(255, SK) :: msg
        type(collapsar_type), intent(in) :: collapsar
        real(RKC) :: logModelInt, logModelInt1, logModelInt2, logNormFac
        logModelInt1 = log(collapsar%tb - mc_dur%min)
        logNormFac = max(collapsar%logtb, collapsar%alphap1 * mc_dur%log%max - collapsar%alpha * collapsar%logtb - collapsar%beta * (mc_dur%max - collapsar%tb))
        failed = isFailedQuad(getDenCollapsarNormed, lb = collapsar%logtb, ub = mc_dur%log%max, help = weps, integral = logModelInt2, msg = msg)
        if (failed) error stop "@getLogModelInt(): "//trim(msg)
        logModelInt2 = log(logModelInt2) + logNormFac
        logModelInt = getLogAddExp(getMinMax(logModelInt1, logModelInt2))
    contains
        function getDenCollapsarNormed(logdur) result(density)
            real(RKC), intent(in) :: logdur
            real(RKC) :: density
            density = exp(collapsar%alphap1 * logdur - collapsar%alpha * collapsar%logtb - collapsar%beta * (exp(logdur) - collapsar%tb) - logNormFac)
        end function
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure elemental function getLogProbMix(dur, logdur, param, logmint) result(logProb)
        use pm_mathMinMax, only: getMinMax
        real(RKC), intent(in) :: dur, logdur
        type(param_type), intent(in) :: param
        type(logmint_type), intent(in) :: logmint
        real(RKC) :: logProb, logProbMerger, logProbCollapsar
        logProbMerger = param%logmixfac%merger + getLogDenMerger(logdur, param%merger) - logmint%merger
        logProbCollapsar = param%logmixfac%collapsar + getLogDenCollapsar(dur, logdur, param%collapsar) - logmint%collapsar
        logProb = getLogAddExp(getMinMax(logProbMerger, logProbCollapsar))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogLike(state) result(logLike)

        real(RKC), intent(in), contiguous :: state(:)
        real(RKC) :: logLike, mixfrac
        type(logmint_type) :: logmint
        type(param_type) :: param
        integer(IK) :: idata

        param%merger%logtb = state(1)
        param%merger%tb = exp(param%merger%logtb)

        param%merger%alpha = state(2)
        param%merger%alphap1 = param%merger%alpha + 1

        param%collapsar%logtb = state(3)
        param%collapsar%tb = exp(param%collapsar%logtb)

        param%collapsar%alpha = state(4)
        param%collapsar%alphap1 = param%collapsar%alpha + 1

        param%collapsar%beta = state(5) ! 0.0156_RKC

        mixfrac = .5_RKC * (tanh(state(6)) + 1) ! state(6) !
        param%logmixfac%merger = log(mixfrac)
        param%logmixfac%collapsar = log(1 - mixfrac)

        if (param%merger%tb < param%collapsar%tb) then
            logmint%merger = getLogModelIntMerger(param%merger)
            logmint%collapsar = getLogModelIntCollapsar(param%collapsar)
            logLike = 0._RKC
            do idata = 1, size(mc_dur%val)
                logLike = logLike + getLogProbMix(mc_dur%val(idata), mc_dur%log%val(idata), param, logmint)
            end do
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
    sampler%outputChainSize = 30000
    sampler%outputFileName = "./out/mixflat"
    sampler%domainAxisName= [ "merger%logtb     " &
                            , "merger%alpha     " &
                            , "collapsar%logtb  " &
                            , "collapsar%alpha  " &
                            , "collapsar%beta   " &
                            , "atanh_2mixfracm1 " &
                            ]
    sampler%proposalStart = [ log(.37_RKC)  & ! merger%logtb
                            , -1.4_RKC      & ! merger%alpha
                            , log(21.3_RKC) & ! collapsar%logtb
                            , -.5_RKC       & ! collapsar%alpha
                            , 0.0156_RKC    & ! collapsar%beta
                            , .3_RKC        & ! atanh(mixfrac * 2 - 1)
                            ]
    sampler%proposalCovMat = reshape(   [  1.E-2, -3.E-3,  1.E-2, -7.E-3, -1.E-4, -6.E-4 &
                                        , -3.E-3,  3.E-3, -4.E-3,  6.E-3,  1.E-4,  1.E-3 &
                                        ,  1.E-2, -4.E-3,  1.E-1, -7.E-2, -8.E-4,  3.E-4 &
                                        , -7.E-3,  6.E-3, -7.E-2,  7.E-2,  8.E-4,  3.E-3 &
                                        , -1.E-4,  1.E-4, -8.E-4,  8.E-4,  1.E-5,  6.E-5 &
                                        , -6.E-4,  1.E-3,  3.E-4,  3.E-3,  6.E-5,  2.E-3 ], shape = [size(sampler%proposalStart), size(sampler%proposalStart)])
    sampler%outputStatus ="retry"
    !sampler%proposalScaleFactor = "0.1 * gelman"
    err = getErrSampling(sampler, getLogLike, size(sampler%proposalStart, 1, IK))
    if (err%occurred) error stop "sampler failed: "//err%msg

    block
        use pm_sampleMean, only: getMean
        integer(IK) :: stat, offset = 1
        real(RKC), allocatable :: table(:,:), param(:)
        stat = getErrTableRead(sampler%outputFileName//"_run1_pid1_sample.txt", table, roff = 1_IK)
        if (stat /= 0) error stop "Failed to read output sample."
        param = getMean(table, dim = 1_IK)
        call disp%show("[tbm, alpham, tbc, alphac, betac, merger_fraction]")
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
        ! [tbm, alpham, tbc, alphac, betac, merger_fraction]
        ! +0.36451513509216688, -1.4068429995163054, +22.325504212059144, -0.54223278907592198, +0.15472485166018178E-1, +0.34074631224225671
        ! [maxloglike, bic]
        ! -4154.3911188154070, +8354.5094086076733
    end block

end