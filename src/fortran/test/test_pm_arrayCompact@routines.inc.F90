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
!>  This module contains implementations of the tests of the procedures under the generic interfaces
!>  [getCompact](@ref pm_arrayCompact::getCompact),
!>  [setCompact](@ref pm_arrayCompact::setCompact).
!>
!>  \todo
!>  \phigh The tests in this file still benefit from expansion and improvement.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     SK_ENABLED && D0_ENABLED
#define GET_SIZE(Array, DIM) len(Array, kind = IK)
#define GET_SLICE(i) i:i
#elif   D1_ENABLED
#define GET_SIZE(Array, DIM) size(Array, DIM, IK)
#define GET_SLICE(i) i
#elif   !D2_ENABLED
#error  "Unrecognized interface."
#endif
        integer(IK)                         :: i, dim, csize, csize_ref
#if     SK_ENABLED && D0_ENABLED
        character(:,SKG)    , allocatable   :: ArrayVerbose, ArrayCompact, ArrayCompact_ref
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG)    , allocatable   :: ArrayVerbose(:), ArrayCompact(:), ArrayCompact_ref(:)
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG)        , allocatable   :: ArrayVerbose(:), ArrayCompact(:), ArrayCompact_ref(:)
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG)        , allocatable   :: ArrayVerbose(:), ArrayCompact(:), ArrayCompact_ref(:)
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG)        , allocatable   :: ArrayVerbose(:), ArrayCompact(:), ArrayCompact_ref(:)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG)           , allocatable   :: ArrayVerbose(:), ArrayCompact(:), ArrayCompact_ref(:)
#elif   SK_ENABLED && D2_ENABLED
        character(2,SKG)    , allocatable   :: ArrayVerbose(:,:), ArrayCompact(:,:), ArrayCompact_ref(:,:)
#elif   IK_ENABLED && D2_ENABLED
        integer(IKG)        , allocatable   :: ArrayVerbose(:,:), ArrayCompact(:,:), ArrayCompact_ref(:,:)
#elif   LK_ENABLED && D2_ENABLED
        logical(LKG)        , allocatable   :: ArrayVerbose(:,:), ArrayCompact(:,:), ArrayCompact_ref(:,:)
#elif   CK_ENABLED && D2_ENABLED
        complex(CKG)        , allocatable   :: ArrayVerbose(:,:), ArrayCompact(:,:), ArrayCompact_ref(:,:)
#elif   RK_ENABLED && D2_ENABLED
        real(RKG)           , allocatable   :: ArrayVerbose(:,:), ArrayCompact(:,:), ArrayCompact_ref(:,:)
#else
#error  "Unrecognized interface."
#endif
        integer(IK)         , allocatable   :: Weight_ref(:)
#if     setCompact_ENABLED
        integer(IK)         , allocatable   :: Weight(:)
#elif   !getCompact_ENABLED
#error  "Unrecognized interface."
#endif

        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%
#if     D0_ENABLED || D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        dim = 1

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED && D0_ENABLED
        allocate(character(0,SKG) :: ArrayCompact_ref)
#else
        allocate(ArrayCompact_ref(0))
#endif
        call setUnifRand(ArrayCompact_ref)
        ArrayCompact_ref = getUnique(ArrayCompact_ref)
        csize_ref = GET_SIZE(ArrayCompact_ref, dim)
        allocate(Weight_ref(csize_ref))
        call setUnifRand(Weight_ref, 1_IK, 10_IK)
        ArrayVerbose = getVerbose(ArrayCompact_ref, Weight_ref, sum(Weight_ref, mask = Weight_ref > 0_IK))

#if     getCompact_ENABLED
        ArrayCompact = getCompact(ArrayVerbose)
        csize = GET_SIZE(ArrayCompact, dim)
#elif   setCompact_ENABLED
        ArrayCompact = ArrayVerbose
        allocate(Weight(GET_SIZE(ArrayCompact, dim)))
        call setCompact(ArrayCompact, Weight, csize)
#endif

        assertion = assertion .and. csize == csize_ref
        call report()
        call test%assert(assertion, SK_"The compact size of an empty verbose array must be properly set.", int(__LINE__, IK))

        assertion = .true._LK
        do i = 1, csize
            assertion = assertion .and. ArrayCompact(GET_SLICE(i)) IS_EQUAL ArrayCompact_ref(GET_SLICE(i))
            call report()
            call test%assert(assertion, SK_"An empty verbose array must be properly condensed.", int(__LINE__, IK))
#if         setCompact_ENABLED
            assertion = assertion .and. Weight(i) == Weight_ref(i)
            call report()
            call test%assert(assertion, SK_"The Weight of an empty verbose array must be properly set.", int(__LINE__, IK))
#endif
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call reset()
#if     SK_ENABLED && D0_ENABLED
        allocate(character(10,SKG) :: ArrayCompact_ref)
#else
        allocate(ArrayCompact_ref(10))
#endif
        call setUnifRand(ArrayCompact_ref)
        ArrayCompact_ref = getUnique(ArrayCompact_ref)
        csize_ref = GET_SIZE(ArrayCompact_ref, dim)
        allocate(Weight_ref(csize_ref))
        call setUnifRand(Weight_ref, 1_IK, 10_IK)
        ArrayVerbose = getVerbose(ArrayCompact_ref, Weight_ref, sum(Weight_ref, mask = Weight_ref > 0_IK))

#if     getCompact_ENABLED
        ArrayCompact = getCompact(ArrayVerbose)
        csize = GET_SIZE(ArrayCompact, dim)
#elif   setCompact_ENABLED
        ArrayCompact = ArrayVerbose
        allocate(Weight(GET_SIZE(ArrayCompact, dim)))
        call setCompact(ArrayCompact, Weight, csize)
#endif

        assertion = assertion .and. csize == csize_ref
        call report()
        call test%assert(assertion, SK_"The compact size of a non-empty verbose array must be properly set.", int(__LINE__, IK))

        assertion = .true._LK
        do i = 1, csize
            assertion = assertion .and. ArrayCompact(GET_SLICE(i)) IS_EQUAL ArrayCompact_ref(GET_SLICE(i))
            call report()
            call test%assert(assertion, SK_"A non-empty verbose array must be properly condensed.", int(__LINE__, IK))
#if         setCompact_ENABLED
            assertion = assertion .and. Weight(i) == Weight_ref(i)
            call report()
            call test%assert(assertion, SK_"The Weight of a non-empty verbose array must be properly set.", int(__LINE__, IK))
#endif
        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%
#elif   D2_ENABLED
        !%%%%%%%%%

        do dim = 1, 2

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
            allocate(ArrayCompact_ref(0,0))
            call setUnifRand(ArrayCompact_ref)
            csize_ref = size(ArrayCompact_ref, dim)
            allocate(Weight_ref(csize_ref))
            call setUnifRand(Weight_ref, 1_IK, 10_IK)
            ArrayVerbose = getVerbose(ArrayCompact_ref, Weight_ref, sum(Weight_ref, mask = Weight_ref > 0_IK), dim)

#if         getCompact_ENABLED
            ArrayCompact = getCompact(ArrayVerbose, dim)
            csize = size(ArrayCompact, dim)
#elif       setCompact_ENABLED
            ArrayCompact = ArrayVerbose
            allocate(Weight(size(ArrayCompact, dim)))
            call setCompact(ArrayCompact, Weight, csize, dim)
#endif

            assertion = assertion .and. csize == csize_ref
            call report()
            call test%assert(assertion, SK_"The compact size of an empty verbose array must be properly set.", int(__LINE__, IK))

            assertion = .true._LK
            do i = 1, csize
                assertion = assertion .and. all(ArrayCompact IS_EQUAL ArrayCompact_ref)
                call report()
                call test%assert(assertion, SK_"An empty verbose array must be properly condensed.", int(__LINE__, IK))
#if             setCompact_ENABLED
                assertion = assertion .and. Weight(i) == Weight_ref(i)
                call report()
                call test%assert(assertion, SK_"The Weight of an empty verbose array must be properly set.", int(__LINE__, IK))
#endif
            end do

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()
#if         SK_ENABLED
            ArrayVerbose = reshape( [ "AB", "CD", "EF" &
                                    , "AB", "CD", "EF" &
                                    , "GH", "IJ", "KL" &
                                    , "GH", "IJ", "KL" &
                                    , "GH", "IJ", "KL" &
                                    , "MN", "OP", "QR" &
                                    ], shape = [6, 3], order = [2, 1])
#elif       IK_ENABLED
            ArrayVerbose = int(reshape( [ +1, +2, +3 &
                                        , +1, +2, +3 &
                                        , +4, +5, +6 &
                                        , +4, +5, +6 &
                                        , +4, +5, +6 &
                                        , +7, +8, +9 &
                                        ], shape = [6, 3], order = [2, 1]), IKG)
#elif       LK_ENABLED
            ArrayVerbose = logical(reshape( [ .false., .false., .true. &
                                            , .false., .false., .true. &
                                            , .true., .false., .false. &
                                            , .true., .false., .false. &
                                            , .true., .false., .false. &
                                            , .false., .true., .false. &
                                            ], shape = [6, 3], order = [2, 1]), LKG)
#elif       CK_ENABLED
            ArrayVerbose = cmplx( reshape(  [ (+1, -1), (+2, -2), (+3, -3) &
                                            , (+1, -1), (+2, -2), (+3, -3) &
                                            , (+4, -4), (+5, -5), (+6, -6) &
                                            , (+4, -4), (+5, -5), (+6, -6) &
                                            , (+4, -4), (+5, -5), (+6, -6) &
                                            , (+7, -7), (+8, -8), (+9, -9) &
                                            ], shape = [6, 3], order = [2, 1]), kind = CKG)
#elif       RK_ENABLED
            ArrayVerbose = real( reshape(   [ +1, +2, +3 &
                                            , +1, +2, +3 &
                                            , +4, +5, +6 &
                                            , +4, +5, +6 &
                                            , +4, +5, +6 &
                                            , +7, +8, +9 &
                                            ], shape = [6, 3], order = [2, 1]), RKG)
#else
#error      "Unrecognized interface."
#endif
            Weight_ref = [2, 3, 1]
            ArrayCompact_ref = ArrayVerbose([1, 3, 6],:)

            if (dim == 2_IK) then
                ArrayVerbose = transpose(ArrayVerbose)
                ArrayCompact_ref = transpose(ArrayCompact_ref)
            end if

            csize_ref = size(ArrayCompact_ref, dim)

#if         getCompact_ENABLED
            ArrayCompact = getCompact(ArrayVerbose, dim)
            csize = size(ArrayCompact, dim)
#elif       setCompact_ENABLED
            ArrayCompact = ArrayVerbose
            allocate(Weight(size(ArrayCompact, dim)))
            call setCompact(ArrayCompact, Weight, csize, dim)
#endif

            assertion = assertion .and. csize == csize_ref
            call report()
            call test%assert(assertion, SK_"The compact size of a non-empty verbose array must be properly set.", int(__LINE__, IK))

            assertion = .true._LK
            do i = 1, csize
                if (dim == 1_IK) then
                    assertion = assertion .and. all(ArrayCompact(csize,:) IS_EQUAL ArrayCompact_ref(csize,:))
                else
                    assertion = assertion .and. all(ArrayCompact(:,csize) IS_EQUAL ArrayCompact_ref(:,csize))
                end if
                call report()
                call test%assert(assertion, SK_"A non-empty verbose array must be properly condensed.", int(__LINE__, IK))
#if             setCompact_ENABLED
                assertion = assertion .and. Weight(i) == Weight_ref(i)
                call report()
                call test%assert(assertion, SK_"The Weight of a non-empty verbose array must be properly set.", int(__LINE__, IK))
#endif
            end do

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
#if         setCompact_ENABLED
            if (allocated(Weight)) deallocate(Weight)
#endif
            if (allocated(Weight_ref)) deallocate(Weight_ref)
            if (allocated(ArrayCompact)) deallocate(ArrayCompact)
            if (allocated(ArrayVerbose)) deallocate(ArrayVerbose)
            if (allocated(ArrayCompact_ref)) deallocate(ArrayCompact_ref)
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "ArrayVerbose       ", ArrayVerbose
                write(test%disp%unit,"(*(g0,:,', '))") "ArrayCompact       ", ArrayCompact
                write(test%disp%unit,"(*(g0,:,', '))") "ArrayCompact_ref   ", ArrayCompact_ref
#if             setCompact_ENABLED
                write(test%disp%unit,"(*(g0,:,', '))") "Weight_ref         ", Weight_ref
                write(test%disp%unit,"(*(g0,:,', '))") "Weight             ", Weight
#endif
                write(test%disp%unit,"(*(g0,:,', '))") "csize_ref          ", csize_ref
                write(test%disp%unit,"(*(g0,:,', '))") "csize              ", csize
                write(test%disp%unit,"(*(g0,:,', '))") "dim                ", dim
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef GET_SLICE
#undef IS_EQUAL
#undef GET_SIZE
#undef ALL