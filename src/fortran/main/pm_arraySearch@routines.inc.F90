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
!>  This file contains the implementation details of the routines under the generic interfaces in [pm_arraySearch](@ref pm_arraySearch).
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the array indexing rules.
#if     D0_D0_ENABLED && SK_ENABLED
#define GET_INDEX(i) i:i+lenValueMinusOne
#define GET_SIZE len
#elif   D1_D0_ENABLED
#define GET_INDEX(i) i
#define GET_SIZE size
#else
#error  "Unrecognized interface."
#endif
        ! Define the comparison rules.
#if     CusCom_ENABLED
#define IS_LESS(a,b) isLess(a,b)
#elif   DefCom_ENABLED && CK_ENABLED
#define IS_LESS(a,b) a%re < b%re
#else
#define IS_LESS(a,b) a < b
#endif
        integer(IK) :: minbin, maxbin
#if     D0_D0_ENABLED && SK_ENABLED
#define MAXBIN GET_SIZE(array, kind = IK) - lenValueMinusOne
        integer(IK) :: lenValueMinusOne
        lenValueMinusOne = len(value, kind = IK) - 1_IK
#elif   D1_D0_ENABLED
#define MAXBIN GET_SIZE(array, kind = IK)
#else
#error  "Unrecognized interface."
#endif
        CHECK_ASSERTION(__LINE__, 0_IK < MAXBIN, SK_"@getBin(): The length of the input array must be non-zero: lenArray = "//getStr(MAXBIN))
        if (IS_LESS(value, array(GET_INDEX(1)))) then
            bin = 0_IK
        else
#if         D0_D0_ENABLED && SK_ENABLED
            maxbin = GET_SIZE(array, kind = IK) - lenValueMinusOne
#elif       D1_D0_ENABLED
            maxbin = GET_SIZE(array, kind = IK)
#else
#error      "Unrecognized interface."
#endif
            if (IS_LESS(value, array(GET_INDEX(maxbin)))) then
                minbin = 1_IK
                do
                    bin = minbin + (maxbin - minbin) / 2_IK
                    if (IS_LESS(value, array(GET_INDEX(bin)))) then
                        if (bin - minbin > 1_IK) then
                            maxbin = bin
                            cycle
                        end if
                        bin = bin - 1_IK
                        return
                    else
                        minbin = bin
                        if (maxbin - bin > 1_IK) then
                            minbin = bin
                            cycle
                        end if
                        return
                    end if
                end do
            else
                bin = GET_SIZE(array, kind = IK)
            end if
        end if

#undef  COMPONENT
#undef  GET_INDEX
#undef  GET_SIZE
#undef  IS_LESS
#undef  MAXBIN