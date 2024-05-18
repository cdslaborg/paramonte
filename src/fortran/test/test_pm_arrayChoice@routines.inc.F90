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
!>  This file contains procedure implementations of [test_pm_arrayChoice](@ref test_pm_arrayChoice).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     SK_ENABLED && D0_ENABLED
#define GET_SIZE(Array) len(Array, kind = IK)
#define GET_SLICE(i) i:i
#elif   D1_ENABLED
#define GET_SIZE(Array) size(Array, 1, IK)
#define GET_SLICE(i) i
#elif   !D2_ENABLED
#error  "Unrecognized interface."
#endif
#if     SK_ENABLED && D0_ENABLED
        character(:,SKC)    , allocatable   :: choices, set
        character(1,SKC)    , parameter     :: lb = "A", ub = "Z"
       !character(1,SKC)                    :: choice
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKC)    , allocatable   :: choices(:), set(:)
        character(2,SKC)    , parameter     :: lb = "AA", ub = "AZ"
       !character(2,SKC)                    :: choice
#elif   IK_ENABLED && D1_ENABLED
        integer(IKC)        , allocatable   :: choices(:), set(:)
        integer(IKC)        , parameter     :: lb = 0, ub = 9
       !integer(IKC)                        :: choice
#elif   LK_ENABLED && D1_ENABLED
        logical(LKC)        , allocatable   :: choices(:), set(:)
        logical(LKC)        , parameter     :: lb = .false., ub = .true.
       !logical(LKC)                        :: choice
#elif   CK_ENABLED && D1_ENABLED
        complex(CKC)        , allocatable   :: choices(:), set(:)
        complex(CKC)        , parameter     :: lb = (-9._CKC, 0._CKC), ub = (0._CKC, +9._CKC)
       !complex(CKC)                        :: choice
#elif   RK_ENABLED && D1_ENABLED
        real(RKC)           , allocatable   :: choices(:), set(:)
        real(RKC)           , parameter     :: lb = 0._RKC, ub = 9._RKC
       !real(RKC)                           :: choice
#else
#error  "Unrecognized interface."
#endif
        logical(LK) :: unique
        logical(LK) :: rngfUsed
        integer(IK) :: itry, csize
        type(display_type) :: disp
        type(xoshiro256ssw_type) :: rngx
        assertion = .true._LK
        
        rngx = xoshiro256ssw_type()

        do itry = 1, 100

            call setResized(set, getUnifRand(1_IK, 9_IK))
            call setUnifRand(set, lb, ub)
            set = getUnique(set)
            csize = getUnifRand(0_IK, 2 * GET_SIZE(set))
            unique = csize <= GET_SIZE(set)
            rngfUsed = getUnifRand()

            call setResized(choices, csize)
#if         getChoice_ENABLED
            choices = getChoice(set, csize)
#elif       setChoice_ENABLED
            if (rngfUsed) then
                call setChoice(rngf, choices, set)
            else
                call setChoice(rngx, choices, set)
            end if
#else
#error      "Unrecognized interface."
#endif
            assertion = assertion .and. GET_SIZE(choices) == csize
            call report()
            call test%assert(assertion, SK_"The size of the output `choices` must be set correctly.", int(__LINE__, IK))

#if         getChoice_ENABLED
            choices = getChoice(set, csize)
#elif       setChoice_ENABLED
            if (rngfUsed) then
                call setChoice(rngf, choices, set)
            else
                call setChoice(rngx, choices, set)
            end if
#endif
            assertion = assertion .and. (choices .allin. set)
            call report()
            call test%assert(assertion, SK_"All output choices must be members of the input set.", int(__LINE__, IK))

#if         getChoice_ENABLED
            choices = getChoice(set, csize, unique)
#elif       setChoice_ENABLED
            if (rngfUsed) then
                call setChoice(rngf, choices, set, unique)
            else
                call setChoice(rngx, choices, set, unique)
            end if
#endif
            assertion = assertion .and. (choices .allin. set)
            call report(unique)
            call test%assert(assertion, SK_"All output choices must be members of the input set.", int(__LINE__, IK))
            if (unique) then
                assertion = assertion .and. isUniqueAll(choices)
                call report()
                call test%assert(assertion, SK_"All output choices must be members of the input set.", int(__LINE__, IK))
            end if

            csize = 1
            call setResized(choices, csize)
#if         getChoice_ENABLED
            choices(GET_SLICE(1)) = getChoice(set)
#elif       setChoice_ENABLED
            if (rngfUsed) then
                call setChoice(rngf, choices(GET_SLICE(1)), set)
            else
                call setChoice(rngx, choices(GET_SLICE(1)), set)
            end if
#endif
            assertion = assertion .and. (choices .allin. set)
            call report()
            call test%assert(assertion, SK_"The output scalar choice must be member of the input set.", int(__LINE__, IK))

        end do

    contains

        subroutine report(unique)
            logical(LK), intent(in), optional :: unique
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%skip
                call disp%show("set")
                call disp%show( set )
                call disp%show("choices")
                call disp%show( choices )
                call disp%show("csize")
                call disp%show( csize )
                call disp%show("present(unique)")
                call disp%show( present(unique) )
                if (present(unique)) then
                call disp%show("unique")
                call disp%show( unique )
                end if
#if             setChoice_ENABLED
                call disp%show("rngfUsed")
                call disp%show( rngfUsed )
#endif
                call disp%skip
                ! LCOV_EXCL_STOP
            end if
        end subroutine

#undef GET_SLICE
#undef IS_EQUAL
#undef GET_SIZE
#undef ALL