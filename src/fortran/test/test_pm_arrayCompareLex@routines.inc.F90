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
!>  This include file contains the implementations of the tests of procedures of [test_pm_arrayCompareLex](@ref test_pm_arrayCompareLex).
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_llt_ENABLED
#define COMPARES_WITH .llt.
        character(*, SK), parameter :: PROCEDURE_NAME = "isllt()"
#elif   test_lle_ENABLED
#define COMPARES_WITH .lle.
        character(*, SK), parameter :: PROCEDURE_NAME = "islle()"
#elif   test_lge_ENABLED
#define COMPARES_WITH .lge.
        character(*, SK), parameter :: PROCEDURE_NAME = "islge()"
#elif   test_lgt_ENABLED
#define COMPARES_WITH .lgt.
        character(*, SK), parameter :: PROCEDURE_NAME = "islgt()"
#else
#error  "Unrecognized interface."
#endif

        logical(LK) :: result, result_ref

#if     test_llt_D0_D0_SK_ENABLED || test_lle_D0_D0_SK_ENABLED || test_lge_D0_D0_SK_ENABLED || test_lgt_D0_D0_SK_ENABLED
#define test_D0_D0_SK_ENABLED 1
        character(:,SKC), allocatable   :: array1, array2
#elif   test_llt_D1_D1_SK_ENABLED || test_lle_D1_D1_SK_ENABLED || test_lge_D1_D1_SK_ENABLED || test_lgt_D1_D1_SK_ENABLED
#define test_D1_D1_SK_ENABLED 1
        character(2,SKC), allocatable   :: array1(:), array2(:)
#elif   test_llt_D1_D1_LK_ENABLED || test_lle_D1_D1_LK_ENABLED || test_lge_D1_D1_LK_ENABLED || test_lgt_D1_D1_LK_ENABLED
#define test_D1_D1_LK_ENABLED 1
        logical(LKC)    , allocatable   :: array1(:), array2(:)
#elif   test_llt_D1_D1_IK_ENABLED || test_lle_D1_D1_IK_ENABLED || test_lge_D1_D1_IK_ENABLED || test_lgt_D1_D1_IK_ENABLED
#define test_D1_D1_IK_ENABLED 1
        integer(IKC)    , allocatable   :: array1(:), array2(:)
#elif   test_llt_D1_D1_CK_ENABLED || test_lle_D1_D1_CK_ENABLED || test_lge_D1_D1_CK_ENABLED || test_lgt_D1_D1_CK_ENABLED
#define test_D1_D1_CK_ENABLED 1
        complex(CKC)    , allocatable   :: array1(:), array2(:)
#elif   test_llt_D1_D1_RK_ENABLED || test_lle_D1_D1_RK_ENABLED || test_lge_D1_D1_RK_ENABLED || test_lgt_D1_D1_RK_ENABLED
#define test_D1_D1_RK_ENABLED 1
        real(RKC)       , allocatable   :: array1(:), array2(:)
#else
#error  "Unrecognized interface."
#endif

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_D0_D0_SK_ENABLED
        array1 = SKC_"A"
        array2 = SKC_"A"
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKC_"AA"]
        array2 = [SKC_"AA"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKC]
        array2 = [.false._LKC]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKC]
        array2 = [0_IKC]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKC]
        array2 = [0._RKC]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKC, 0._CKC)]
        array2 = [(0._CKC, 0._CKC)]
#endif

#if     test_llt_ENABLED
        result_ref = .false._LK
#elif   test_lle_ENABLED
        result_ref = .true._LK
#elif   test_lge_ENABLED
        result_ref = .true._LK
#elif   test_lgt_ENABLED
        result_ref = .false._LK
#endif

        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_" must correctly lexically compare two equivalent arrays of length 1.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_D0_D0_SK_ENABLED
        array1 = SKC_"A"
        array2 = SKC_"B"
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKC_"AA"]
        array2 = [SKC_"BB"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKC]
        array2 = [.true._LKC]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKC]
        array2 = [1_IKC]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKC]
        array2 = [1._RKC]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKC, 0._CKC)]
        array2 = [(1._CKC, 1._CKC)]
#endif

#if     test_llt_ENABLED
        result_ref = .true._LK
#elif   test_lle_ENABLED
        result_ref = .true._LK
#elif   test_lge_ENABLED
        result_ref = .false._LK
#elif   test_lgt_ENABLED
        result_ref = .false._LK
#endif

        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_" must correctly perform lexical comparison when array1 is less than array2 both of length 1.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_D0_D0_SK_ENABLED
        array1 = SKC_"A"
        array2 = SKC_"A "
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKC_"AA"]
        array2 = [SKC_"AA", SKC_"AA"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKC]
        array2 = [.false._LKC, .false._LKC]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKC]
        array2 = [0_IKC, 0_IKC]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKC]
        array2 = [0._RKC, 0._RKC]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKC, 0._CKC)]
        array2 = [(0._CKC, 0._CKC), (0._CKC, 0._CKC)]
#endif

#if     test_llt_ENABLED
        result_ref = .true._LK
#elif   test_lle_ENABLED
        result_ref = .true._LK
#elif   test_lge_ENABLED
        result_ref = .false._LK
#elif   test_lgt_ENABLED
        result_ref = .false._LK
#endif

        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_" must correctly perform lexical comparison when array1 of length 1 is less than array2 of length 2.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_D0_D0_SK_ENABLED
        array1 = SKC_"A"
        array2 = SKC_"A "
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKC_"AA"]
        array2 = [SKC_"AA", SKC_"AA"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKC]
        array2 = [.false._LKC, .false._LKC]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKC]
        array2 = [0_IKC, 0_IKC]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKC]
        array2 = [0._RKC, 0._RKC]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKC, 0._CKC)]
        array2 = [(0._CKC, 0._CKC), (0._CKC, 0._CKC)]
#endif

#if     test_llt_ENABLED
        result_ref = .true._LK
#elif   test_lle_ENABLED
        result_ref = .true._LK
#elif   test_lge_ENABLED
        result_ref = .false._LK
#elif   test_lgt_ENABLED
        result_ref = .false._LK
#endif

        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_" must correctly perform lexical comparison when array1 of length 1 is less than array2 of length 2.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_D0_D0_SK_ENABLED
        array1 = SKC_"A "
        array2 = SKC_"A"
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKC_"AA", SKC_"AA"]
        array2 = [SKC_"AA"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKC, .false._LKC]
        array2 = [.false._LKC]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKC, 0_IKC]
        array2 = [0_IKC]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKC, 0._RKC]
        array2 = [0._RKC]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKC, 0._CKC), (0._CKC, 0._CKC)]
        array2 = [(0._CKC, 0._CKC)]
#endif

#if     test_llt_ENABLED
        result_ref = .false._LK
#elif   test_lle_ENABLED
        result_ref = .false._LK
#elif   test_lge_ENABLED
        result_ref = .true._LK
#elif   test_lgt_ENABLED
        result_ref = .true._LK
#endif

        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_" must correctly perform lexical comparison when array1 of length 2 is more than array2 of length 1.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_D0_D0_SK_ENABLED
        array1 = SKC_"AA"
        array2 = SKC_"AB"
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKC_"AA"]
        array2 = [SKC_"AA", SKC_"AB"]
#elif   test_D1_D1_LK_ENABLED 
        array1 = [.false._LKC, .false._LKC]
        array2 = [.false._LKC, .true._LKC]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKC, 0_IKC]
        array2 = [0_IKC, 1_IKC]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKC, 0._RKC]
        array2 = [0._RKC, 1._RKC]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKC, 0._CKC), (0._CKC, 0._CKC)]
        array2 = [(0._CKC, 0._CKC), (0._CKC, 1._CKC)]
#endif

#if     test_llt_ENABLED
        result_ref = .true._LK
#elif   test_lle_ENABLED
        result_ref = .true._LK
#elif   test_lge_ENABLED
        result_ref = .false._LK
#elif   test_lgt_ENABLED
        result_ref = .false._LK
#endif

        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_" must correctly perform lexical comparison when array1 of length 2 is less than array2 of length 2.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_D0_D0_SK_ENABLED
        array1 = SKC_"AB"
        array2 = SKC_"AA"
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKC_"AA", SKC_"AB"]
        array2 = [SKC_"AA"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKC, .true._LKC]
        array2 = [.false._LKC, .false._LKC]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKC, 1_IKC]
        array2 = [0_IKC, 0_IKC]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKC, 1._RKC]
        array2 = [0._RKC, 0._RKC]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKC, 0._CKC), (0._CKC, 1._CKC)]
        array2 = [(0._CKC, 0._CKC), (0._CKC, 0._CKC)]
#endif

#if     test_llt_ENABLED
        result_ref = .false._LK
#elif   test_lle_ENABLED
        result_ref = .false._LK
#elif   test_lge_ENABLED
        result_ref = .true._LK
#elif   test_lgt_ENABLED
        result_ref = .true._LK
#endif

        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_" must correctly perform lexical comparison when array1 of length 2 is less than array2 of length 2.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            result = array1 COMPARES_WITH array2
            assertion = assertion .and. result .eqv. result_ref
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "array1     ", array1
                write(test%disp%unit,"(*(g0,:,', '))") "array2     ", array2
                write(test%disp%unit,"(*(g0,:,', '))") "result     ", result
                write(test%disp%unit,"(*(g0,:,', '))") "result_ref ", result_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef test_D0_D0_SK_ENABLED
#undef test_D1_D1_SK_ENABLED
#undef test_D1_D1_LK_ENABLED
#undef test_D1_D1_IK_ENABLED
#undef test_D1_D1_CK_ENABLED
#undef test_D1_D1_RK_ENABLED
#undef COMPARES_WITH
