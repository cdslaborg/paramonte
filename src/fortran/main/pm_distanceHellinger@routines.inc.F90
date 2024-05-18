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
!>  This include file contains the implementation of procedures in [pm_distanceMahal](@ref pm_distanceMahal).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getDisHellSq_ENABLED && PMF_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CHECK_ASSERTION(__LINE__, size(p, 1, IK) > 0, SK_": The condition `size(p) > 0` must hold. size(p) = "//getStr(size(p, 1, IK)))
        CHECK_ASSERTION(__LINE__, all(0 <= p) .and. all(p <= 1), SK_": The condition `all(0 <= p) .and. all(p <= 1)` must hold. p = "//getStr(p))
        CHECK_ASSERTION(__LINE__, all(0 <= q) .and. all(q <= 1), SK_": The condition `all(0 <= q) .and. all(q <= 1)` must hold. q = "//getStr(q))
        CHECK_ASSERTION(__LINE__, size(p, 1, IK) == size(q, 1, IK), SK_": The condition `size(p) == size(q)` must hold. size(p), size(q) = "//getStr([size(p, 1, IK), size(q, 1, IK)]))
        CHECK_ASSERTION(__LINE__, abs(1 - sum(p)) < sqrt(epsilon(0._RKC)), SK_": The condition `abs(1 - sum(p)) < sqrt(epsilon(0._RKC))` must hold. sum(p), sqrt(epsilon(0._RKC)) = "//getStr([sum(p), sqrt(epsilon(0._RKC))]))
        CHECK_ASSERTION(__LINE__, abs(1 - sum(q)) < sqrt(epsilon(0._RKC)), SK_": The condition `abs(1 - sum(q)) < sqrt(epsilon(0._RKC))` must hold. sum(q), sqrt(epsilon(0._RKC)) = "//getStr([sum(q), sqrt(epsilon(0._RKC))]))
        hellsq = 1._RKC - sum(sqrt(p * q))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   getDisHellSq_ENABLED && PRC_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(*, SK), parameter :: PROCEDURE_NAME = MODULE_NAME//SK_"@getDisHellSq(): "
        logical(LK) :: failed_def
        real(RKC) :: lb_def, ub_def
        if (present(lb)) then
            lb_def = lb
        else
            lb_def = getInfNeg(lb_def)
        end if
        if (present(ub)) then
            ub_def = ub
        else
            ub_def = getInfPos(ub_def)
        end if
        failed_def = isFailedQuad(getBhatCoef, lb_def, ub_def, weps, hellsq, msg = msg)
        if (present(failed)) then
            failed = failed_def
        elseif (failed_def) then
            error stop PROCEDURE_NAME//trim(msg)
        end if
        hellsq = 1._RKC - hellsq

    contains

        function getBhatCoef(x) result(coef)
            real(RKC), intent(in) :: x
            real(RKC) :: coef
            coef = sqrt(p(x) * q(x))
        end function
#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif