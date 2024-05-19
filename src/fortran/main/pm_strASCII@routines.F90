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
!>  This file contains procedure implementations of [pm_strASCII](@ref pm_strASCII).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday 2:21 AM, August 30, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

submodule (pm_strASCII) routines ! LCOV_EXCL_LINE

#if CHECK_ENABLED
    use pm_err, only: getFine
    use pm_val2str, only: getStr
    use pm_err, only: setAsserted
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) \
call setAsserted(ASSERTION,getFine(__FILE__,LINE)//MODULE_NAME//MSG);
#else
#define CHECK_ASSERTION(LINE,ASSERTION,MSG) continue;
#endif

    implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLocSpace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getLocSpace_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getLocSpace_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getLocSpace_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getLocSpace_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getLocSpace_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLocSpace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getLocNonSpace_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getLocNonSpace_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getLocNonSpace_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getLocNonSpace_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getLocNonSpace_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getLocNonSpace_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getLocNonSpace_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCharDigit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isCharDigit_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isCharDigit_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isCharDigit_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isCharDigit_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isCharDigit_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isCharDigit_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrDigitAll_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrDigitAll_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrDigitAll_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrDigitAll_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrDigitAll_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrDigitAll_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrDigitAll_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrDigitAny_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrDigitAny_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrDigitAny_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrDigitAny_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrDigitAny_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrDigitAny_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrDigitAny_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrDigit_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrDigit_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrDigit_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrDigit_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrDigit_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrDigit_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrDigit_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrInteger_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrInteger_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrInteger_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrInteger_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrInteger_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrInteger_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrInteger_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrComplex_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrComplex_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrComplex_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrComplex_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrComplex_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrComplex_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrComplex_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrReal_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrReal_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrReal_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrReal_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrReal_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrReal_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrReal_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrNumber_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrNumber_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrNumber_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrNumber_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrNumber_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrNumber_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrNumber_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCharUpper_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isCharUpper_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isCharUpper_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isCharUpper_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isCharUpper_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isCharUpper_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isCharUpper_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCharLower_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isCharLower_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isCharLower_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isCharLower_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isCharLower_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isCharLower_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isCharLower_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrUpperAll_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrUpperAll_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrUpperAll_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrUpperAll_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrUpperAll_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrUpperAll_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrUpperAll_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrLowerAll_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrLowerAll_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrLowerAll_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrLowerAll_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrLowerAll_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrLowerAll_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrLowerAll_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrUpperAny_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrUpperAny_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrUpperAny_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrUpperAny_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrUpperAny_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrUpperAny_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrUpperAny_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrLowerAny_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrLowerAny_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrLowerAny_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrLowerAny_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrLowerAny_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrLowerAny_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrLowerAny_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrUpper_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrUpper_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrUpper_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrUpper_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrUpper_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrUpper_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrUpper_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrLower_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrLower_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrLower_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrLower_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrLower_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrLower_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrLower_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCharAlphaNum_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isCharAlphaNum_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isCharAlphaNum_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isCharAlphaNum_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isCharAlphaNum_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isCharAlphaNum_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isCharAlphaNum_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrAlphaNumAll_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrAlphaNumAll_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrAlphaNumAll_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrAlphaNumAll_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrAlphaNumAll_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrAlphaNumAll_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrAlphaNumAll_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrAlphaNumAny_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrAlphaNumAny_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrAlphaNumAny_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrAlphaNumAny_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrAlphaNumAny_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrAlphaNumAny_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrAlphaNumAny_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrAlphaNum_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrAlphaNum_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrAlphaNum_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrAlphaNum_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrAlphaNum_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrAlphaNum_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrAlphaNum_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isCharAlpha_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isCharAlpha_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isCharAlpha_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isCharAlpha_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isCharAlpha_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isCharAlpha_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isCharAlpha_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrAlphaAll_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrAlphaAll_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrAlphaAll_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrAlphaAll_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrAlphaAll_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrAlphaAll_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrAlphaAll_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrAlphaAny_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrAlphaAny_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrAlphaAny_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrAlphaAny_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrAlphaAny_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrAlphaAny_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrAlphaAny_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define isStrAlpha_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure isStrAlpha_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure isStrAlpha_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure isStrAlpha_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure isStrAlpha_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure isStrAlpha_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef isStrAlpha_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCharUpper_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getCharUpper_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getCharUpper_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getCharUpper_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getCharUpper_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getCharUpper_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCharUpper_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCharUpper_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setCharUpper_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setCharUpper_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setCharUpper_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setCharUpper_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setCharUpper_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCharUpper_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getCharLower_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getCharLower_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getCharLower_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getCharLower_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getCharLower_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getCharLower_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getCharLower_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setCharLower_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setCharLower_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setCharLower_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setCharLower_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setCharLower_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setCharLower_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setCharLower_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getStrUpper_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getStrUpper_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getStrUpper_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getStrUpper_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getStrUpper_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getStrUpper_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getStrUpper_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setStrUpper_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setStrUpper_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setStrUpper_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setStrUpper_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setStrUpper_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setStrUpper_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setStrUpper_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getStrLower_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getStrLower_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getStrLower_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getStrLower_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getStrLower_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getStrLower_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getStrLower_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setStrLower_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setStrLower_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setStrLower_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setStrLower_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setStrLower_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setStrLower_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setStrLower_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getStrQuoted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getStrQuoted_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getStrQuoted_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getStrQuoted_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getStrQuoted_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getStrQuoted_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getStrQuoted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setStrQuoted_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setStrQuoted_SK5
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setStrQuoted_SK4
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setStrQuoted_SK3
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setStrQuoted_SK2
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setStrQuoted_SK1
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setStrQuoted_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define getAsciiFromEscaped_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure getAsciiFromEscapedNew_SK5
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure getAsciiFromEscapedNew_SK4
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure getAsciiFromEscapedNew_SK3
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure getAsciiFromEscapedNew_SK2
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure getAsciiFromEscapedNew_SK1
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef getAsciiFromEscaped_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define setAsciiFromEscaped_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define Rep_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setAsciiFromEscapedRep_SK5
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setAsciiFromEscapedRep_SK4
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setAsciiFromEscapedRep_SK3
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setAsciiFromEscapedRep_SK2
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setAsciiFromEscapedRep_SK1
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef Rep_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define New_ENABLED 1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define SK_ENABLED 1

#if SK5_ENABLED
    module procedure setAsciiFromEscapedNew_SK5
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK5
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK4_ENABLED
    module procedure setAsciiFromEscapedNew_SK4
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK4
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK3_ENABLED
    module procedure setAsciiFromEscapedNew_SK3
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK3
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK2_ENABLED
    module procedure setAsciiFromEscapedNew_SK2
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK2
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#if SK1_ENABLED
    module procedure setAsciiFromEscapedNew_SK1
        use pm_mathNumSys, only: getDecimal
        use pm_kind, only: SKG => SK1
#include "pm_strASCII@routines.inc.F90"
    end procedure
#endif

#undef SK_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef New_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef setAsciiFromEscaped_ENABLED

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef CHECK_ASSERTION

end submodule routines