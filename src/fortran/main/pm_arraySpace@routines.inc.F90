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
!>  This include file contains procedure implementation of [pm_arraySpace](@ref pm_arraySpace).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Sunday 3:33 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%
#if     getLinSpace_ENABLED
        !%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0 <= count, SK_"@getLinSpace(): The condition `0 <= count` must hold. count = "//getStr(count))
        call setLinSpace(linSpace, x1, x2, fopen, lopen)

        !%%%%%%%%%%%%%%%%%%
#elif   setLinSpace_ENABLED
        !%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define TYPE_OF_SPACE complex(TKG)
#elif   RK_ENABLED
#define TYPE_OF_SPACE real(TKG)
#else
#error  "Unrecognized interface."
#endif
        integer(IK) :: i, count
        logical(LK) :: fopen_def, lopen_def
        real(TKG), parameter :: HALF = 0.5_TKG
        TYPE_OF_SPACE :: stepSize, start
        count = size(linSpace, 1, IK)
        if (0_IK < count) then
            if (present(fopen)) then
                fopen_def = fopen
            else
                fopen_def = .false._LK
            end if
            if (present(lopen)) then
                lopen_def = lopen
            else
                lopen_def = .false._LK
            end if
            if (lopen_def) then
                stepSize = (x2 - x1) / count
                if (fopen_def) then
                    start = x1 + HALF * stepSize
                else
                    start = x1
                end if
            else
                if (fopen_def) then
                    stepSize = (x2 - x1) / count
                    start = x1 + stepSize
                elseif (count > 1_IK) then
                    stepSize = (x2 - x1) / (count - 1_IK)
                    start = x1
                elseif (count == 1_IK) then
                    linSpace = x1
                    return
                end if
                count = count - 1
            end if
            linSpace(0) = start
            do concurrent(i = 1 : count - 1)
                linSpace(i) = start + i * stepSize
            end do
            if (.not. lopen_def) linSpace(count) = x2
        end if
#undef  TYPE_OF_SPACE

        !%%%%%%%%%%%%%%%%%%
#elif   getLogSpace_ENABLED
        !%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, 0 <= count, SK_"@getLogSpace(): The condition `0 <= count` must hold. count = "//getStr(count))
        !check_assertion(__LINE__, 0. < base .and. base /= 1., SK_"@getLogSpace(): The condition `0. < base .and. base /= 1.` must hold. base = "//getStr(base))
        if (present(base)) then
            logSpace = base**getLinSpace(logx1, logx2, count, fopen, lopen)
        else
            logSpace = exp(getLinSpace(logx1, logx2, count, fopen, lopen))
        end if

        !%%%%%%%%%%%%%%%%%%
#elif   setLogSpace_ENABLED
        !%%%%%%%%%%%%%%%%%%

        !check_assertion(__LINE__, 0. < base .and. base /= 1., SK_"@getLogSpace(): The condition `0. < base .and. base /= 1.` must hold. base = "//getStr(base))
        call setLinSpace(logSpace, logx1, logx2, fopen, lopen)
        if (present(base)) then
            logSpace = base**logSpace
        else
            logSpace = exp(logSpace)
        end if
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif