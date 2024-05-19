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
!>  This module contains implementations of the procedures of [test_pm_arrayInit](@ref test_pm_arrayInit).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the `logical` operators.
#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif
        ! Define the indexing rules.
#if     D0_ENABLED
#define ALL
#elif   D1_ENABLED
#define SET_ADIM(X) X(:)
#elif   D2_ENABLED
#define SET_ADIM(X) X(:,:)
#elif   D3_ENABLED
#define SET_ADIM(X) X(:,:,:)
#else
#error  "Unrecognized interface."
#endif
        ! Set the dimension of `Core`.
#if     Arr_ENABLED
#define SET_CDIM(X) SET_ADIM(X)
#elif   Sca_ENABLED
#define SET_CDIM(X) X
#else
#error  "Unrecognized interface."
#endif
        integer :: itest
        ! Declare objects.
#if     D0_ENABLED
        character(1,SKG)                :: halo
        character(:,SKG), allocatable   :: Array, array_ref, Core
#else
#if     SK_ENABLED
        character(2,SKG)                :: halo
        character(2,SKG), allocatable   :: &
#elif   IK_ENABLED
        integer(IKG)                    :: halo
        integer(IKG)    , allocatable   :: &
#elif   LK_ENABLED
        logical(LKG)                    :: halo
        logical(LKG)    , allocatable   :: &
#elif   CK_ENABLED
        complex(CKG)                    :: halo
        complex(CKG)    , allocatable   :: &
#elif   RK_ENABLED
        real(RKG)                       :: halo
        real(RKG)       , allocatable   :: &
#else
#error "Unrecognized interface."
#endif
        SET_ADIM(Array), SET_ADIM(array_ref), SET_CDIM(Core)
#endif
        ! Define indexing objects.
#if     D0_ENABLED || D1_ENABLED
        integer(IK) :: Asize, Coffset, Csize
#else
        integer(IK) :: Asize(rank(Array)), Coffset(rank(Array)), Csize(rank(Array))
#endif
        assertion = .true._LK
        do itest = 1, 100

            if (allocated(Core)) deallocate(Core)
            if (allocated(Array)) deallocate(Array)
            if (allocated(array_ref)) deallocate(array_ref)

            call setUnifRand(halo)
            call setUnifRand(Asize, 0_IK, int(50. / max(1,rank(Array))))
            call setUnifRand(Csize, 0_IK, Asize)
            ! Allocate `Array`.
#if         D0_ENABLED
            allocate(character(Asize,SKG) :: array_ref, Array)
            array_ref(:) = repeat(halo, Asize)
#elif       D1_ENABLED
            allocate(array_ref(Asize), Array(Asize), source = halo)
#elif       D2_ENABLED
            allocate(array_ref(Asize(1), Asize(2)), Array(Asize(1), Asize(2)), source = halo)
#elif       D3_ENABLED
            allocate(array_ref(Asize(1), Asize(2), Asize(3)), Array(Asize(1), Asize(2), Asize(3)), source = halo)
#else
#error      "Unrecognized interface."
#endif
            ! Allocate `Core`.
#if         Arr_ENABLED && D0_ENABLED
            allocate(character(Csize,SKG) :: Core)
#elif       Arr_ENABLED && D1_ENABLED
            allocate(Core(Csize))
#elif       Arr_ENABLED && D2_ENABLED
            allocate(Core(Csize(1), Csize(2)))
#elif       Arr_ENABLED && D3_ENABLED
            allocate(Core(Csize(1), Csize(2), Csize(3)))
#else
            allocate(Core, source = halo)
#endif
            call setUnifRand(Core)
            call setUnifRand(Coffset, 0_IK, Asize - Csize)
#if         getCoreHalo_ENABLED && Arr_ENABLED
            Array = getCoreHalo(Asize, Core, halo, Coffset)
#elif       setCoreHalo_ENABLED && Arr_ENABLED
            call setCoreHalo(Array, Core, halo, Coffset)
#elif       getCoreHalo_ENABLED && Sca_ENABLED
            Array = getCoreHalo(Asize, Core, halo, Coffset, Csize)
#elif       setCoreHalo_ENABLED && Sca_ENABLED
            call setCoreHalo(Array, Core, halo, Coffset, Csize)
#else
#error      "Unrecognized interface."
#endif
            call setCoreHalo_ref()
            call report()

            if (allocated(Core)) deallocate(Core)
            if (allocated(Array)) deallocate(Array)
            if (allocated(array_ref)) deallocate(array_ref)

            call setUnifRand(halo)
            call setUnifRand(Asize, 0_IK, int(50. / max(1,rank(Array))))
            call setUnifRand(Csize, 0_IK, Asize)
            ! Allocate `Array`.
#if         D0_ENABLED
            allocate(character(Asize,SKG) :: array_ref, Array)
            array_ref(:) = repeat(halo, Asize)
#elif       D1_ENABLED
            allocate(array_ref(Asize), Array(Asize), source = halo)
#elif       D2_ENABLED
            allocate(array_ref(Asize(1), Asize(2)), Array(Asize(1), Asize(2)), source = halo)
#elif       D3_ENABLED
            allocate(array_ref(Asize(1), Asize(2), Asize(3)), Array(Asize(1), Asize(2), Asize(3)), source = halo)
#else
#error      "Unrecognized interface."
#endif
            ! Allocate `Core`.
#if         Arr_ENABLED && D0_ENABLED
            allocate(character(Csize,SKG) :: Core)
#elif       Arr_ENABLED && D1_ENABLED
            allocate(Core(Csize))
#elif       Arr_ENABLED && D2_ENABLED
            allocate(Core(Csize(1), Csize(2)))
#elif       Arr_ENABLED && D3_ENABLED
            allocate(Core(Csize(1), Csize(2), Csize(3)))
#else
            allocate(Core, source = halo)
#endif
            call setUnifRand(Core)
            call setUnifRand(Coffset, 0_IK, Asize - Csize)
#if         getCoreHalo_ENABLED && Arr_ENABLED
            Array = getCoreHalo(Asize, Core, halo, Coffset)
#elif       setCoreHalo_ENABLED && Arr_ENABLED
            call setCoreHalo(Array, Core, halo, Coffset)
#elif       getCoreHalo_ENABLED && Sca_ENABLED
            Array = getCoreHalo(Asize, Core, halo, Coffset, Csize)
#elif       setCoreHalo_ENABLED && Sca_ENABLED
            call setCoreHalo(Array, Core, halo, Coffset, Csize)
#else
#error      "Unrecognized interface."
#endif
            call setCoreHalo_ref()
            call report()

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine setCoreHalo_ref()
#if         D0_ENABLED
            if (len(Core) == 1) then
                array_ref(coffset + 1 : coffset + csize) = repeat(Core, csize)
            else
                array_ref(coffset + 1 : coffset + csize) = Core
            end if
#elif       D1_ENABLED
            array_ref(coffset + 1_IK : coffset + csize) = Core
#elif       D2_ENABLED
            array_ref(Coffset(1) + 1_IK : Coffset(1) + Csize(1), Coffset(2) + 1_IK : Coffset(2) + Csize(2)) = Core
#elif       D3_ENABLED
            array_ref(Coffset(1) + 1_IK : Coffset(1) + Csize(1), Coffset(2) + 1_IK : Coffset(2) + Csize(2), Coffset(3) + 1_IK : Coffset(3) + Csize(3)) = Core
#else
#error      "Unrecognized interface."
#endif
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report()
            assertion = assertion .and. logical(ALL(Array IS_EQUAL array_ref), LK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call test%disp%skip()
                call test%disp%show("Coffset")
                call test%disp%show( Coffset )
                call test%disp%show("Csize")
                call test%disp%show( Csize )
                call test%disp%show("halo")
                call test%disp%show( halo )
                call test%disp%show("Core")
                call test%disp%show( Core )
                call test%disp%show("Array")
                call test%disp%show( Array )
                call test%disp%show("array_ref")
                call test%disp%show( array_ref )
                call test%disp%skip()
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"The output `array` must be constructed correctly.", int(__LINE__, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef SET_ADIM
#undef SET_CDIM
#undef IS_EQUAL
#undef ALL
