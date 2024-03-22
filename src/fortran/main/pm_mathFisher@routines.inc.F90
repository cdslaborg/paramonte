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
!>  This include file contains implementations of the procedures in module [pm_mathFisher](@ref pm_mathFisher).<br>
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%
#if     getFisher_ENABLED
        !%%%%%%%%%%%%%%%%

#if     FDD_ENABLED
        CHECK_ASSERTION(__LINE__, -1 < val .and. val < +1, SK_"@getFisher(): The condition `-1 < val .and. val < +1` must hold. val = "//getStr(val))
        fisherz = atanh(val)
#elif   FLU_ENABLED
        CHECK_ASSERTION(__LINE__, lb < val .and. val < ub, SK_"@getFisher(): The condition `lb < val .and. val < ub` must hold. val, lb, ub = "//getStr(val, lb, ub))
        fisherz = atanh(2 * (val - lb) / (ub - lb) - 1._RKC)
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%
#elif   getFisherInv_ENABLED
        !%%%%%%%%%%%%%%%%%%%

#if     FDD_ENABLED
        val = tanh(fisherz)
#elif   FLU_ENABLED
        CHECK_ASSERTION(__LINE__, lb < ub, SK_"@getFisherInv(): The condition `lb < ub` must hold. lb, ub = "//getStr(lb, ub))
        val = (.5_RKC * tanh(fisherz) + .5_RKC) * (ub - lb) + lb
#else
#error  "Unrecognized interface."
#endif

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif