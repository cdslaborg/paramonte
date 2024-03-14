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
!>  This file contains procedure implementations of [test_pm_arrayChange](@ref test_pm_arrayChange).
!>
!>  \fintest
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
#define GET_DIFF(X,Y) ichar(X, IK) - ichar(Y, IK)
        integer(IK)         , parameter     :: START_STEP = 1_IK
        character(:,SKC)    , allocatable   :: choices, set
        character(1,SKC)    , parameter     :: lb = "A", ub = "Z"
        character(1,SKC)                    :: start, finit
        integer(IK)                         :: step
#elif   IK_ENABLED && D1_ENABLED
#define GET_DIFF(X,Y) X - Y
        integer(IKC)        , parameter     :: START_STEP = 1_IKC
        integer(IKC)        , allocatable   :: choices(:), set(:)
        integer(IKC)        , parameter     :: lb = 1   , ub = 9
        integer(IKC)                        :: start, finit
        integer(IKC)                        :: step
#elif   RK_ENABLED && D1_ENABLED
#define GET_DIFF(X,Y) X - Y
        real(RKC)           , parameter     :: START_STEP = 1._RKC
        real(RKC)           , allocatable   :: choices(:), set(:)
        real(RKC)           , parameter     :: lb = 1._RKC, ub = 9._RKC
        real(RKC)                           :: start, finit
        real(RKC)                           :: step
#else
#error  "Unrecognized interface."
#endif
        logical(LK) :: unique
        logical(LK) :: rngfUsed
        integer(IK) :: itry, csize
        type(display_type) :: disp
        type(xoshiro256ssw_type) :: rngx
        rngx = xoshiro256ssw_type()
        assertion = .true._LK

        do itry = 1, 100

            call setUnifRand(start, lb, ub)
            call setUnifRand(finit, lb, ub)
            if (finit < start) then
                step = -getUnifRand(START_STEP, START_STEP + GET_DIFF(start,finit))
            else
                step = +getUnifRand(START_STEP, START_STEP + GET_DIFF(finit,start))
            end if
            set = getRange(start, finit, step)
            csize = getUnifRand(0_IK, 2 * GET_SIZE(set))
            unique = csize <= GET_SIZE(set)
            rngfUsed = getUnifRand()

            call setResized(choices, csize)
#if         getChange_ENABLED
            choices = getChange(csize, start, finit, step)
#elif       setChange_ENABLED
            if (rngfUsed) then
                call setChange(rngf, choices, start, finit, step)
            else
                call setChange(rngx, choices, start, finit, step)
            end if
#else
#error      "Unrecognized interface."
#endif
            assertion = assertion .and. GET_SIZE(choices) == csize
            call report()
            call test%assert(assertion, SK_"The size of the output `choices` must be set correctly.", int(__LINE__, IK))

#if         getChange_ENABLED
            choices = getChange(csize, start, finit, step)
#elif       setChange_ENABLED
            if (rngfUsed) then
                call setChange(rngf, choices, start, finit, step)
            else
                call setChange(rngx, choices, start, finit, step)
            end if
#endif
            assertion = assertion .and. (choices .allin. set)
            call report()
            call test%assert(assertion, SK_"All output choices must be members of the input set.", int(__LINE__, IK))

#if         getChange_ENABLED
            choices = getChange(csize, start, finit, step, unique)
#elif       setChange_ENABLED
            if (rngfUsed) then
                call setChange(rngf, choices, start, finit, step, unique)
            else
                call setChange(rngx, choices, start, finit, step, unique)
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
#if             setChange_ENABLED
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
#undef GET_DIFF
#undef ALL