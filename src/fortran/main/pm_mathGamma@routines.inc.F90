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
!>  This include file contains implementations of the procedures in module [pm_mathGamma](@ref pm_mathGamma).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Sunday 11:23 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getGammaIncLow_ENABLED || getGammaIncUpp_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        integer(IK) :: info
#if     getGammaIncLow_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@getGammaIncLow()"
        real(RKG) :: gammaIncUpp
#elif   getGammaIncUpp_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = SK_"@getGammaIncUpp()"
        real(RKG) :: gammaIncLow
#else
#error  "Unrecognized interface."
#endif
        call setGammaInc(gammaIncLow, gammaIncUpp, x, kappa, info)
        if (info < 0_IK) error stop SK_"@file::"//__FILE__//SK_"@line::"//getStr(__LINE__)//&! LCOV_EXCL_LINE
        MODULE_NAME//PROCEDURE_NAME//SK_": The Incomplete Gamma function failed to converge." ! LCOV_EXCL_LINE

        !%%%%%%%%%%%%%%%%%%
#elif   setGammaInc_ENABLED
        !%%%%%%%%%%%%%%%%%%

        call setGammaIncGil(gammaIncLow, gammaIncUpp, x, kappa, info)

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif