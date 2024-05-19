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
!>  This include file contains implementations of the procedures in module [pm_except](@ref pm_except).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Sunday 11:23 PM, September 19, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define zero.
#if     IK_ENABLED
        integer(IKG), parameter :: ZERO = 0_IKG, MINMAX = huge(0_IKG)
#elif   CK_ENABLED
        use pm_complexCompareAny, only: operator(<)
        complex(CKG), parameter :: ZERO = (0._CKG, 0._CKG), MINMAX = huge(0._CKG)
#elif   RK_ENABLED
        real(RKG)   , parameter :: ZERO = 0._RKG, MINMAX = huge(0._RKG)
#else
#error  "Unrecognized interface."
#endif
        !%%%%%%%%%%%%%%%%%%%%%%
#if     isAddOutflowPos_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        outflow = logical(ZERO < a .and. ZERO < b, LK)
        if (outflow) outflow = logical(+MINMAX - a < b, LK)

        !%%%%%%%%%%%%%%%%%%%%%%
#elif   isAddOutflowNeg_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%

        outflow = logical(a < ZERO .and. b < ZERO, LK)
        if (outflow) outflow = logical(b < -MINMAX - a, LK)

        !%%%%%%%%%%%%%%%%%%%
#elif   isAddOutflow_ENABLED
        !%%%%%%%%%%%%%%%%%%%

        outflow = isAddOutflowPos(a, b) .or. isAddOutflowNeg(a, b)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   RK_ENABLED && (getInfPos_ENABLED || setInfPos_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        infPos = ieee_value(infPos, ieee_positive_inf)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   RK_ENABLED && (getInfNeg_ENABLED || setInfNeg_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        infNeg = ieee_value(infNeg, ieee_negative_inf)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   CK_ENABLED && (getInfPos_ENABLED || setInfPos_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        infPos%re = ieee_value(infPos%re, ieee_positive_inf)
        infPos%im = ieee_value(infPos%im, ieee_positive_inf)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   CK_ENABLED && (getInfNeg_ENABLED || setInfNeg_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        infNeg%re = ieee_value(infNeg%re, ieee_negative_inf)
        infNeg%im = ieee_value(infNeg%im, ieee_negative_inf)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   RK_ENABLED && isInfPos_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        infPos = .not. (ieee_is_finite(x) .or. ieee_is_negative(x))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   RK_ENABLED && isInfNeg_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        infNeg = ieee_is_negative(x) .and. .not. ieee_is_finite(x)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   CK_ENABLED && isInfPos_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        infPos = isInfPos(x%re) .or. isInfPos(x%im)
        !infPos = .not. ( (ieee_is_finite(x%re) .or. ieee_is_negative(x%re)) .and. (ieee_is_finite(x%im) .or. ieee_is_negative(x%im)) )

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   CK_ENABLED && isInfNeg_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        infNeg = isInfNeg(x%re) .or. isInfNeg(x%im)
        !infNeg = (ieee_is_negative(x%re) .and. .not. ieee_is_finite(x%re)) .or. (ieee_is_negative(x%im) .and. .not. ieee_is_finite(x%im))

        !%%%%%%%%%%%%
#elif   isInf_ENABLED
        !%%%%%%%%%%%%

        inf = isInfPos(x) .or. isInfNeg(x)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   CK_ENABLED && (getNAN_ENABLED || setNAN_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nan%re = ieee_value(nan%re, ieee_quiet_nan)
        nan%im = ieee_value(nan%im, ieee_quiet_nan)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   RK_ENABLED && (getNAN_ENABLED || setNAN_ENABLED)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !real, parameter :: NAN = transfer(2143289344_IK, 1.) This works only on i686 and x86_64 arch.
        nan = ieee_value(nan, ieee_quiet_nan)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isNAN_ENABLED && IEEE_ENABLED && CK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        isNotANumber = ieee_is_nan(x%re) .or. ieee_is_nan(x%im)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isNAN_ENABLED && IEEE_ENABLED && RK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        isNotANumber = ieee_is_nan(x)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif   isNAN_ENABLED && XNEQ_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        isNotANumber = logical(x /= xcopy, LK)

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif