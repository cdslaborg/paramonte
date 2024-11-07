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
!>  This include file contains procedure implementation of the generic interfaces of [pm_mathRootTest](@ref pm_mathRootTest).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     func1_typer_ENABLED
        self%lb = getOption(-3._RKG, lb)
        self%ub = getOption(+3._RKG, ub)
        self%root = [-2._RKG, -1._RKG, 0._RKG, 1._RKG, 2._RKG]
        CHECK_ASSERTION(__LINE__, all(self%lb <= self%root), SK_"@getFunc1(): The condition `all(self%lb <= self%root)` must hold. self%lb, self%ub = "//getStr([self%lb, self%ub]))
        CHECK_ASSERTION(__LINE__, all(self%root <= self%ub), SK_"@getFunc1(): The condition `all(self%root <= self%ub)` must hold. self%lb, self%ub = "//getStr([self%lb, self%ub]))
        self%desc = "Func1_type: an algebraic integrand of the form f(x) = x * (x**2 + 1) * (x**2 + 4) for x in (lb, ub) with real roots: "//getStr(self%root)
#elif   getFunc1_ENABLED
        CHECK_ASSERTION(__LINE__, self%lb <= x .and. x <= self%ub, SK_"@getFunc1(): The condition `self%lb <= x .and. x <= self%ub` must hold. self%lb, x, self%ub = "//getStr([self%lb, x, self%ub]))
        func = x * (x**2 + 1._RKG) * (x**2 + 4._RKG)
#else
#error  "Unrecognized interface."
#endif
