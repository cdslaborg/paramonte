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
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

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
        character(:,SKG), allocatable   :: array1, array2
#elif   test_llt_D1_D1_SK_ENABLED || test_lle_D1_D1_SK_ENABLED || test_lge_D1_D1_SK_ENABLED || test_lgt_D1_D1_SK_ENABLED
#define test_D1_D1_SK_ENABLED 1
        character(2,SKG), allocatable   :: array1(:), array2(:)
#elif   test_llt_D1_D1_LK_ENABLED || test_lle_D1_D1_LK_ENABLED || test_lge_D1_D1_LK_ENABLED || test_lgt_D1_D1_LK_ENABLED
#define test_D1_D1_LK_ENABLED 1
        logical(LKG)    , allocatable   :: array1(:), array2(:)
#elif   test_llt_D1_D1_IK_ENABLED || test_lle_D1_D1_IK_ENABLED || test_lge_D1_D1_IK_ENABLED || test_lgt_D1_D1_IK_ENABLED
#define test_D1_D1_IK_ENABLED 1
        integer(IKG)    , allocatable   :: array1(:), array2(:)
#elif   test_llt_D1_D1_CK_ENABLED || test_lle_D1_D1_CK_ENABLED || test_lge_D1_D1_CK_ENABLED || test_lgt_D1_D1_CK_ENABLED
#define test_D1_D1_CK_ENABLED 1
        complex(CKG)    , allocatable   :: array1(:), array2(:)
#elif   test_llt_D1_D1_RK_ENABLED || test_lle_D1_D1_RK_ENABLED || test_lge_D1_D1_RK_ENABLED || test_lgt_D1_D1_RK_ENABLED
#define test_D1_D1_RK_ENABLED 1
        real(RKG)       , allocatable   :: array1(:), array2(:)
#else
#error  "Unrecognized interface."
#endif

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     test_D0_D0_SK_ENABLED
        array1 = SKG_"A"
        array2 = SKG_"A"
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKG_"AA"]
        array2 = [SKG_"AA"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKG]
        array2 = [.false._LKG]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKG]
        array2 = [0_IKG]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKG]
        array2 = [0._RKG]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKG, 0._CKG)]
        array2 = [(0._CKG, 0._CKG)]
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
        array1 = SKG_"A"
        array2 = SKG_"B"
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKG_"AA"]
        array2 = [SKG_"BB"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKG]
        array2 = [.true._LKG]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKG]
        array2 = [1_IKG]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKG]
        array2 = [1._RKG]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKG, 0._CKG)]
        array2 = [(1._CKG, 1._CKG)]
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
        array1 = SKG_"A"
        array2 = SKG_"A "
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKG_"AA"]
        array2 = [SKG_"AA", SKG_"AA"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKG]
        array2 = [.false._LKG, .false._LKG]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKG]
        array2 = [0_IKG, 0_IKG]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKG]
        array2 = [0._RKG, 0._RKG]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKG, 0._CKG)]
        array2 = [(0._CKG, 0._CKG), (0._CKG, 0._CKG)]
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
        array1 = SKG_"A"
        array2 = SKG_"A "
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKG_"AA"]
        array2 = [SKG_"AA", SKG_"AA"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKG]
        array2 = [.false._LKG, .false._LKG]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKG]
        array2 = [0_IKG, 0_IKG]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKG]
        array2 = [0._RKG, 0._RKG]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKG, 0._CKG)]
        array2 = [(0._CKG, 0._CKG), (0._CKG, 0._CKG)]
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
        array1 = SKG_"A "
        array2 = SKG_"A"
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKG_"AA", SKG_"AA"]
        array2 = [SKG_"AA"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKG, .false._LKG]
        array2 = [.false._LKG]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKG, 0_IKG]
        array2 = [0_IKG]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKG, 0._RKG]
        array2 = [0._RKG]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKG, 0._CKG), (0._CKG, 0._CKG)]
        array2 = [(0._CKG, 0._CKG)]
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
        array1 = SKG_"AA"
        array2 = SKG_"AB"
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKG_"AA"]
        array2 = [SKG_"AA", SKG_"AB"]
#elif   test_D1_D1_LK_ENABLED 
        array1 = [.false._LKG, .false._LKG]
        array2 = [.false._LKG, .true._LKG]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKG, 0_IKG]
        array2 = [0_IKG, 1_IKG]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKG, 0._RKG]
        array2 = [0._RKG, 1._RKG]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKG, 0._CKG), (0._CKG, 0._CKG)]
        array2 = [(0._CKG, 0._CKG), (0._CKG, 1._CKG)]
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
        array1 = SKG_"AB"
        array2 = SKG_"AA"
#elif   test_D1_D1_SK_ENABLED
        array1 = [SKG_"AA", SKG_"AB"]
        array2 = [SKG_"AA"]
#elif   test_D1_D1_LK_ENABLED
        array1 = [.false._LKG, .true._LKG]
        array2 = [.false._LKG, .false._LKG]
#elif   test_D1_D1_IK_ENABLED
        array1 = [0_IKG, 1_IKG]
        array2 = [0_IKG, 0_IKG]
#elif   test_D1_D1_RK_ENABLED
        array1 = [0._RKG, 1._RKG]
        array2 = [0._RKG, 0._RKG]
#elif   test_D1_D1_CK_ENABLED
        array1 = [(0._CKG, 0._CKG), (0._CKG, 1._CKG)]
        array2 = [(0._CKG, 0._CKG), (0._CKG, 0._CKG)]
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
