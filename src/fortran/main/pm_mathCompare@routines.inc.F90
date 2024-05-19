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
!>  This include file contains procedure implementation of [pm_mathCompare](@ref pm_mathCompare).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 3:33 AM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        use pm_except, only: isNAN
        use pm_except, only: isInf, isInfNeg, isInfPos
        integer, parameter :: TKG = kind(x)
        real(TKG) :: reltol_def, abstol_def
        real(TKG) :: absDiff

        if (isNAN(x) .or. isNAN(y)) then
            close = .false._LK
            return
        end if

        if (.not. (isInf(x) .or. isInf(y))) then

            close = logical(x == y, LK)
            if (close) return

            if (present(reltol)) then
                CHECK_ASSERTION(__LINE__, 0._TKG <= reltol, \
                SK_"@isClose(): The condition `0. < reltol` must hold. reltol = "//getStr(reltol)) ! fpp
                reltol_def = reltol
            else
                reltol_def = epsilon(0._TKG)
            end if
            if (present(abstol)) then
                CHECK_ASSERTION(__LINE__, 0._TKG <= abstol, \
                SK_"@isClose(): The condition `0. <= abstol` must hold. abstol = "//getStr(abstol)) ! fpp
                abstol_def = abstol
            else
                abstol_def = tiny(0._TKG)
            end if
            absDiff = abs(y - x)
#if         isCloseReference_ENABLED
            close = logical(absDiff <= abs(reltol_def * x) .or. absDiff <= abstol_def, LK)
#elif       isCloseStrong_ENABLED
            close = logical((absDiff <= abs(reltol_def * x) .and. absDiff <= abs(reltol_def * y)) .or. absDiff <= abstol_def, LK)
#elif       isCloseWeak_ENABLED || isCloseDefault_ENABLED
            close = logical(absDiff <= abs(reltol_def * x) .or. absDiff <= abs(reltol_def * y) .or. absDiff <= abstol_def, LK)
#elif       isCloseMean_ENABLED
            close = logical(absDiff <= abs(reltol_def * 0.5_TKG * (x + y)) .or. absDiff <= abstol_def, LK)
#else
#error      "Unrecognized interface."
#endif
            return
        end if
        close = (isInfNeg(x) .and. isInfNeg(y)) .or. (isInfPos(x) .and. isInfPos(y))