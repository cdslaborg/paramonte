 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This include file contains the implementations of procedures of [pm_batse](@ref pm_batse).
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getCorrectionLogEffPPF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        real(RKC), parameter :: THRESH_ERFC_AMP_RKC = real(THRESH_ERFC_AMP, RKC)
        real(RKC), parameter :: THRESH_ERFC_AVG_RKC = real(THRESH_ERFC_AVG, RKC)
        real(RKC), parameter :: THRESH_ERFC_STD_INV_RKC = real(THRESH_ERFC_STD_INV, RKC)
        correction = THRESH_ERFC_AMP_RKC * erfc((logT90 - THRESH_ERFC_AVG_RKC) / THRESH_ERFC_STD_INV_RKC)
        ! + THRESH_ERFC_BASE ! adding this term will make the effective peak flux equivalent to PF1024ms

        !%%%%%%%%%%%%%%%%%%%
#elif   getLogEffPPF_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        logEffPPF  = logPPF64 - getCorrectionLogEffPPF(logT90)

        !%%%%%%%%%%%%%%%%%%%
#elif   getLogEffPPF_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        logEffPPF  = logPPF64 - getCorrectionLogEffPPF(logT90)

        !%%%%%%%%%%%%%%%%%
#elif   getLogPbol_ENABLED
        !%%%%%%%%%%%%%%%%%

        logPbol = logPF53 - getLogPF53(logEpk, 0._RKC)

        !%%%%%%%%%%%%%%%%%
#elif   getLogPF53_ENABLED
        !%%%%%%%%%%%%%%%%%

        real(RKC), parameter :: LOGPF53_MINUS_LOGPBOL_RKC = real(LOGPF53_MINUS_LOGPBOL, RKC)
        if (logEpk < -6.712165960423344_RKC) then
            logPF53 = logPbol   + LOGPF53_MINUS_LOGPBOL_RKC ! 11.328718657530706
        elseif (logEpk < 3.453877639491069_RKC) then
            logPF53 = logPbol   + 13.20790440362500600_RKC + logEpk * &
                                ( 0.309360000000000000_RKC + logEpk * &
                                ( 0.001980382837478830_RKC + logEpk * &
                                ( 0.000299892598248466_RKC + logEpk * &
                                ( 1.25602147173493e-05_RKC - logEpk   &
                                * 1.27171265917873e-05_RKC ))))
        elseif (logEpk < 5.756462732485115_RKC) then
            logPF53 = logPbol   + 4.400884836537660000_RKC + logEpk * &
                                ( 39.71039000000000000_RKC - logEpk * &
                                ( 41.95557432120050000_RKC - logEpk * &
                                ( 20.60525451895990000_RKC - logEpk * &
                                ( 5.510436247342930000_RKC - logEpk * &
                                ( 0.832525333390336000_RKC - logEpk * &
                                ( 0.067135977132248900_RKC - logEpk   &
                                * 0.002254876138523550_RKC ))))))
        elseif (logEpk < 9.210340371976184_RKC) then
            logPF53 = logPbol   + 6.451981585674900000_RKC + logEpk * &
                                ( 4.569070000000000000_RKC - logEpk * &
                                ( 0.837198158654537000_RKC - logEpk * &
                                ( 0.055416002698982300_RKC - logEpk   &
                                * 0.001219684856402480_RKC )))
        elseif (logEpk < 12.455573549219071_RKC) then
            logPF53 = logPbol   - 24.09731285126340000_RKC + logEpk * &
                                ( 26.70637000000000000_RKC - logEpk * &
                                ( 6.286981551320860000_RKC - logEpk * &
                                ( 0.667762738216888000_RKC - logEpk * &
                                ( 0.033549115287895400_RKC - logEpk   &
                                * 0.000651366755890191_RKC ))))
        else
            logPF53 = logPbol   + LOGPF53_MINUS_LOGPBOL_RKC
        end if
        !write(*,"(*(g0.13,:,', '))") "logEpk, logPbol, logPF53<0.0: ", logEpk, logPbol, logPF53

        !%%%%%%%%%%%%%%%%%%%
#elif   getLog10PF53_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        real(RKC), parameter :: LOG10PF53_MINUS_LOG10PBOL_RKC = real(LOG10PF53_MINUS_LOG10PBOL, RKC)
        if (log10epk < -2.915056638230699_RKC) then
            log10PF53 = log10pbol    + LOG10PF53_MINUS_LOG10PBOL_RKC ! 4.9200000000000000_RKC
        elseif (log10epk < 1.5_RKC) then
            log10PF53 = log10pbol   + 5.7361200000000000_RKC + log10epk * &
                                    ( 0.3093600000000000_RKC + log10epk * &
                                    ( 0.0045600000000000_RKC + log10epk * &
                                    ( 0.0015900000000000_RKC + log10epk * &
                                    ( 0.0001533360000000_RKC - log10epk   &
                                    * 0.0003574800000000_RKC ))))
        elseif (log10epk < 2.5_RKC) then
            log10PF53 = log10pbol   + 1.9112800000000000_RKC + log10epk * &
                                    ( 39.710390000000000_RKC - log10epk * &
                                    ( 96.606280000000000_RKC - log10epk * &
                                    ( 109.24696000000000_RKC - log10epk * &
                                    ( 67.271800000000000_RKC - log10epk * &
                                    ( 23.402390000000000_RKC - log10epk * &
                                    ( 4.3454400000000000_RKC - log10epk   &
                                    * 0.3360600000000000_RKC ))))))
        elseif (log10epk < 4._RKC) then
            log10PF53 = log10pbol   + 2.8020600000000000_RKC + log10epk * &
                                    ( 4.5690700000000000_RKC - log10epk * &
                                    ( 1.9277200000000000_RKC - log10epk * &
                                    ( 0.2938100000000000_RKC - log10epk   &
                                    * 0.0148900000000000_RKC )))
        elseif (log10epk < 5.4093868613659435_RKC) then
            log10PF53 = log10pbol   - 10.465330000000000_RKC + log10epk * &
                                    ( 26.706370000000000_RKC - log10epk * &
                                    ( 14.476310000000000_RKC - log10epk * &
                                    ( 3.5404100000000000_RKC - log10epk * &
                                    ( 0.4095700000000000_RKC - log10epk   &
                                    * 0.0183100000000000_RKC ))))
        else ! if (log10epk>=5.4093868613659435_RKC) then
            log10PF53 = log10pbol   + LOG10PF53_MINUS_LOG10PBOL_RKC
        end if

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif