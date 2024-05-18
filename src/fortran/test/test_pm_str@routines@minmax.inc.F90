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
!>  This module contains implementations of the tests of the procedures under the generic interfaces in [pm_str](@ref pm_str).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMaxLoc_SK_ENABLED
#define GEN_VAL getMaxLoc
#elif   getMaxVal_SK_ENABLED
#define GEN_VAL getMaxVal
#elif   getMinLoc_SK_ENABLED
#define GEN_VAL getMinLoc
#elif   getMinVal_SK_ENABLED
#define GEN_VAL getMinVal
#else
#error  "Unrecognized Interface."
#endif

        character(:, SK), allocatable   :: Str(:)
#if     getMaxLoc_SK_ENABLED || getMinLoc_SK_ENABLED
        integer(IK)     , allocatable   :: Out(:), Out_ref(:)
#elif   getMaxVal_SK_ENABLED || getMinVal_SK_ENABLED
        character(1, SK), allocatable   :: Out(:), Out_ref(:)
#else
#error  "Unrecognized Interface."
#endif
        logical(LK), allocatable :: Mask(:)
        integer(IK) :: i

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

        Str = [SK_"paramontethe", SK_"theparamonte"]
#if     getMaxLoc_SK_ENABLED
        Out_ref = [8_IK, 1_IK]
#elif   getMaxVal_SK_ENABLED
        Out_ref = [SK_"t", SK_"t"]
#elif   getMinLoc_SK_ENABLED
        Out_ref = [2_IK, 5_IK]
#elif   getMinVal_SK_ENABLED
        Out_ref = [SK_"a", SK_"a"]
#endif
        allocate(Out, mold = Out_ref)

        Out(1) = GEN_VAL(Str(1))
        assertion = assertion .and. Out(1) == Out_ref(1)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string.")

        Out(2) = GEN_VAL(Str(2))
        assertion = assertion .and. Out(2) == Out_ref(2)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string.")

        Out = GEN_VAL(Str)
        assertion = assertion .and. all(Out == Out_ref)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper vector value for an input string vector.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()

        Str = [SK_"paramontethe", SK_"theparamonte"]
#if     getMaxLoc_SK_ENABLED
        Out_ref = [10_IK, 1_IK]
        Mask = [( .true._LK, i = 1, len(str,IK) )]
        Mask(8) = .false._LK
#elif   getMaxVal_SK_ENABLED
        Out_ref = [SK_"t", SK_"t"]
        Mask = [( .true._LK, i = 1, len(str,IK) )]
        Mask(8) = .false._LK
#elif   getMinLoc_SK_ENABLED
        Out_ref = [4_IK, 5_IK]
        Mask = [( .true._LK, i = 1, len(str,IK) )]
        Mask(2) = .false._LK
#elif   getMinVal_SK_ENABLED
        Out_ref = [SK_"a", SK_"a"]
        Mask = [( .true._LK, i = 1, len(str,IK) )]
        Mask(2) = .false._LK
#endif
        allocate(Out, mold = Out_ref)

        Out(1) = GEN_VAL(Str(1), Mask)
        assertion = assertion .and. Out(1) == Out_ref(1)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string with `Mask` present.")

        Out(2) = GEN_VAL(Str(2), Mask)
        assertion = assertion .and. Out(2) == Out_ref(2)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string with `Mask` present.")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMaxLoc_SK_ENABLED || getMinLoc_SK_ENABLED

        call reset()

        Str = [SK_"paramontethe", SK_"theparamonte"]
#if     getMaxLoc_SK_ENABLED
        Out_ref = [8_IK, 1_IK]
#elif   getMinLoc_SK_ENABLED
        Out_ref = [2_IK, 5_IK]
#endif
        allocate(Out, mold = Out_ref)

        Out(1) = GEN_VAL(Str(1), back = .false._LK)
        assertion = assertion .and. Out(1) == Out_ref(1)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string with `back = .false._LK`.")

        Out(2) = GEN_VAL(Str(2), back = .false._LK)
        assertion = assertion .and. Out(2) == Out_ref(2)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string with `back = .false._LK`.")

        Out = GEN_VAL(Str, back = .false._LK)
        assertion = assertion .and. all(Out == Out_ref)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper vector value for an input string vector with `back = .false._LK`.")

#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMaxLoc_SK_ENABLED || getMinLoc_SK_ENABLED

        call reset()

        Str = [SK_"paramontethe", SK_"theparamonte"]
#if     getMaxLoc_SK_ENABLED
        Out_ref = [10_IK, 1_IK]
        Mask = [( .true._LK, i = 1, len(str,IK) )]
        Mask(8) = .false._LK
#elif   getMinLoc_SK_ENABLED
        Out_ref = [4_IK, 5_IK]
        Mask = [( .true._LK, i = 1, len(str,IK) )]
        Mask(2) = .false._LK
#endif
        allocate(Out, mold = Out_ref)

        Out(1) = GEN_VAL(Str(1), Mask, back = .false._LK)
        assertion = assertion .and. Out(1) == Out_ref(1)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string with `Mask` and `back = .false._LK`.")

        Out(2) = GEN_VAL(Str(2), back = .false._LK)
        assertion = assertion .and. Out(2) == Out_ref(2)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string with `Mask` and `back = .false._LK`.")

#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMaxLoc_SK_ENABLED || getMinLoc_SK_ENABLED

        call reset()

        Str = [SK_"paramontethe", SK_"theparamonte"]
#if     getMaxLoc_SK_ENABLED
        Out_ref = [10_IK, 11_IK]
#elif   getMinLoc_SK_ENABLED
        Out_ref = [4_IK, 7_IK]
#endif
        allocate(Out, mold = Out_ref)

        Out(1) = GEN_VAL(Str(1), back = .true._LK)
        assertion = assertion .and. Out(1) == Out_ref(1)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string with `back = .true._LK`.")

        Out(2) = GEN_VAL(Str(2), back = .true._LK)
        assertion = assertion .and. Out(2) == Out_ref(2)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string with `back = .true._LK`.")

        Out = GEN_VAL(Str, back = .true._LK)
        assertion = assertion .and. all(Out == Out_ref)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper vector value for an input string vector with `back = .true._LK`.")

#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getMaxLoc_SK_ENABLED || getMinLoc_SK_ENABLED

        call reset()

        Str = [SK_"paramontethe", SK_"theparamonte"]
#if     getMaxLoc_SK_ENABLED
        Out_ref = [8_IK, 11_IK]
        Mask = [( .true._LK, i = 1, len(str,IK) )]
        Mask(10) = .false._LK
#elif   getMinLoc_SK_ENABLED
        Out_ref = [2_IK, 7_IK]
        Mask = [( .true._LK, i = 1, len(str,IK) )]
        Mask(4) = .false._LK
#endif
        allocate(Out, mold = Out_ref)

        Out(1) = GEN_VAL(Str(1), Mask, back = .true._LK)
        assertion = assertion .and. Out(1) == Out_ref(1)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string with `Mask` and `back = .true._LK`.")

        Out(2) = GEN_VAL(Str(2), back = .true._LK)
        assertion = assertion .and. Out(2) == Out_ref(2)
        call report()
        call test%assert(assertion, PROCEDURE_NAME//SK_"must output the proper scalar value for an input scalar string with `Mask` and `back = .true._LK`.")

#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(str)) deallocate(str)
            if (allocated(Out)) deallocate(Out)
            if (allocated(Out_ref)) deallocate(Out_ref)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))")   "Str      ", Str
                write(test%disp%unit,"(*(g0,:,', '))")   "Out      ", Out
                write(test%disp%unit,"(*(g0,:,', '))")   "Out_ref  ", Out_ref
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  GEN_VAL