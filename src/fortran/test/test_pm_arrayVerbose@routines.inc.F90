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
!>  This module contains implementations of the tests of the procedures under the generic interfaces of [pm_arrayVerbose](@ref pm_arrayVerbose).
!>
!>  \todo
!>  \phigh The tests in this file still benefit from expansion and improvement.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Define the equivalence operator.
#if     LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif
        ! Define the slicing rules.
#if     SK_ENABLED && D0_ENABLED
#define GET_SIZE(Array, DIM) len(Array, kind = IK)
#define GET_SLICE(i) i:i
#elif   D1_ENABLED
#define GET_SIZE(Array, DIM) size(Array, DIM, IK)
#define GET_SLICE(i) i
#elif   !D2_ENABLED
#error  "Unrecognized interface."
#endif
        integer(IK) :: i, iw, dim, counter, csize, csize_ref
        integer(IK), allocatable :: weight(:)
#if     SK_ENABLED && D0_ENABLED
        character(:,SKG), allocatable :: arrayVerbose, arrayCompact
#elif   SK_ENABLED && D1_ENABLED
        character(2,SKG), allocatable :: arrayVerbose(:), arrayCompact(:)
#elif   IK_ENABLED && D1_ENABLED
        integer(IKG), allocatable :: arrayVerbose(:), arrayCompact(:)
#elif   LK_ENABLED && D1_ENABLED
        logical(LKG), allocatable :: arrayVerbose(:), arrayCompact(:)
#elif   CK_ENABLED && D1_ENABLED
        complex(CKG), allocatable :: arrayVerbose(:), arrayCompact(:)
#elif   RK_ENABLED && D1_ENABLED
        real(RKG), allocatable :: arrayVerbose(:), arrayCompact(:)
#elif   SK_ENABLED && D2_ENABLED
        character(2,SKG), allocatable :: arrayVerbose(:,:), arrayCompact(:,:)
#elif   IK_ENABLED && D2_ENABLED
        integer(IKG), allocatable :: arrayVerbose(:,:), arrayCompact(:,:)
#elif   LK_ENABLED && D2_ENABLED
        logical(LKG), allocatable :: arrayVerbose(:,:), arrayCompact(:,:)
#elif   CK_ENABLED && D2_ENABLED
        complex(CKG), allocatable :: arrayVerbose(:,:), arrayCompact(:,:)
#elif   RK_ENABLED && D2_ENABLED
        real(RKG), allocatable :: arrayVerbose(:,:), arrayCompact(:,:)
#else
#error  "Unrecognized interface."
#endif
        assertion = .true._LK

        !%%%%%%%%%%%%%%%%%%%%%%%
#if     D0_ENABLED || D1_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%

        dim = 1

        block

            call reset()
#if         SK_ENABLED && D0_ENABLED
            allocate(character(0,SKG) :: arrayCompact)
#else
            allocate(arrayCompact(0))
#endif
            allocate(weight(GET_SIZE(arrayCompact, dim)))
            call setUnifRand(weight, -5_IK, 10_IK)
            csize_ref = sum(weight, mask = weight > 0_IK)

#if         getVerbose_ENABLED
            arrayVerbose = getVerbose(arrayCompact, weight, csize_ref)
            csize = GET_SIZE(arrayVerbose, dim)
#elif       setVerbose_ENABLED
#endif

            assertion = assertion .and. csize == csize_ref
            call report()
            call test%assert(assertion, SK_"The verbose size of an empty compact array must be properly set.", int(__LINE__, IK))

            assertion = assertion .and. all(shape(arrayVerbose) == csize_ref)
            call report()
            call test%assert(assertion, SK_"The shape of the output verbose array must be properly set.", int(__LINE__, IK))

            assertion = .true._LK
            counter = 0_IK
            do i = 1, size(weight)
                do iw = 1, weight(i)
                    counter = counter + 1_IK
                    assertion = assertion .and. arrayVerbose(GET_SLICE(counter)) IS_EQUAL arrayCompact(GET_SLICE(i))
                    call report()
                    call test%assert(assertion, SK_"An empty compact array must be properly expanded.", int(__LINE__, IK))
                end do
            end do

        end block

        block

            call reset()
#if         SK_ENABLED && D0_ENABLED
            allocate(character(10,SKG) :: arrayCompact)
#else
            allocate(arrayCompact(10))
#endif
            call setUnifRand(arrayCompact)
            allocate(weight(GET_SIZE(arrayCompact, dim)))
            call setUnifRand(weight, -5_IK, 10_IK)
            csize_ref = sum(weight, mask = weight > 0_IK)

#if         getVerbose_ENABLED
            arrayVerbose = getVerbose(arrayCompact, weight, csize_ref)
            csize = GET_SIZE(arrayVerbose, dim)
#elif       setVerbose_ENABLED
#endif

            assertion = assertion .and. csize == csize_ref
            call report()
            call test%assert(assertion, SK_"The verbose size of a non-empty compact array must be properly set.", int(__LINE__, IK))

            assertion = assertion .and. all(shape(arrayVerbose) == csize_ref)
            call report()
            call test%assert(assertion, SK_"The shape of the output verbose array must be properly set.", int(__LINE__, IK))

            assertion = .true._LK
            counter = 0_IK
            do i = 1, size(weight)
                do iw = 1, weight(i)
                    counter = counter + 1_IK
                    assertion = assertion .and. arrayVerbose(GET_SLICE(counter)) IS_EQUAL arrayCompact(GET_SLICE(i))
                    call report()
                    call test%assert(assertion, SK_"A non-empty compact array must be properly expanded.", int(__LINE__, IK))
                end do
            end do

        end block

        !%%%%%%%%%
#elif   D2_ENABLED
        !%%%%%%%%%

        do dim = 1, 2

            block

                call reset()
                allocate(arrayCompact(0,0))
                allocate(weight(size(arrayCompact, dim)))
                call setUnifRand(weight, -5_IK, 10_IK)
                csize_ref = sum(weight, mask = weight > 0_IK)

#if             getVerbose_ENABLED
                arrayVerbose = getVerbose(arrayCompact, weight, csize_ref, dim)
                csize = size(arrayVerbose, dim)
#elif           setVerbose_ENABLED
#endif

                assertion = assertion .and. csize == csize_ref
                call report()
                call test%assert(assertion, SK_"The verbose size of an empty compact array must be properly set.", int(__LINE__, IK))

                if (dim == 1_IK) then
                    assertion = assertion .and. all(shape(arrayVerbose) == [csize_ref, size(arrayCompact, 2)])
                else
                    assertion = assertion .and. all(shape(arrayVerbose) == [size(arrayCompact, 1), csize_ref])
                end if
                call report()
                call test%assert(assertion, SK_"The shape of the output verbose array must be properly set.", int(__LINE__, IK))

                assertion = .true._LK
                counter = 0_IK
                do i = 1, size(weight)
                    do iw = 1, weight(i)
                        counter = counter + 1_IK
                        if (dim == 1_IK) then
                            assertion = assertion .and. all(arrayVerbose(counter,:) IS_EQUAL arrayCompact(i,:))
                        else
                            assertion = assertion .and. all(arrayVerbose(:,counter) IS_EQUAL arrayCompact(:,i))
                        end if
                        call report()
                        call test%assert(assertion, SK_"An empty compact array must be properly expanded.", int(__LINE__, IK))
                    end do
                end do

            end block

            block

                call reset()
                allocate(arrayCompact(getUnifRand(0,10), getUnifRand(0,10)))
                call setUnifRand(arrayCompact)
                allocate(weight(size(arrayCompact, dim)))
                call setUnifRand(weight, -5_IK, 10_IK)
                csize_ref = sum(weight, mask = weight > 0_IK)

#if             getVerbose_ENABLED
                arrayVerbose = getVerbose(arrayCompact, weight, csize_ref, dim)
                csize = size(arrayVerbose, dim)
#elif           setVerbose_ENABLED
#endif

                assertion = assertion .and. csize == csize_ref
                call report()
                call test%assert(assertion, SK_"The verbose size of a non-empty compact array must be properly set.", int(__LINE__, IK))

                if (dim == 1_IK) then
                    assertion = assertion .and. all(shape(arrayVerbose) == [csize_ref, size(arrayCompact, 2)])
                else
                    assertion = assertion .and. all(shape(arrayVerbose) == [size(arrayCompact, 1), csize_ref])
                end if
                call report()
                call test%assert(assertion, SK_"The shape of the output verbose array must be properly set.", int(__LINE__, IK))

                assertion = .true._LK
                counter = 0_IK
                do i = 1, size(weight)
                    do iw = 1, weight(i)
                        counter = counter + 1_IK
                        if (dim == 1_IK) then
                            assertion = assertion .and. all(arrayVerbose(counter,:) IS_EQUAL arrayCompact(i,:))
                        else
                            assertion = assertion .and. all(arrayVerbose(:,counter) IS_EQUAL arrayCompact(:,i))
                        end if
                        call report()
                        call test%assert(assertion, SK_"A non-empty compact array must be properly expanded.", int(__LINE__, IK))
                    end do
                end do

            end block

        end do

#else
        !%%%%%%%%%%%%%%%%%%%%%%%%
#error  "Unrecognized interface."
        !%%%%%%%%%%%%%%%%%%%%%%%%
#endif

    contains

        subroutine reset()
            if (allocated(weight)) deallocate(weight)
            if (allocated(arrayVerbose)) deallocate(arrayVerbose)
            if (allocated(arrayCompact)) deallocate(arrayCompact)
        end subroutine

        subroutine report()
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "arrayCompact       ", arrayCompact
                write(test%disp%unit,"(*(g0,:,', '))") "arrayVerbose       ", arrayVerbose
                write(test%disp%unit,"(*(g0,:,', '))") "arrayCompact       ", arrayCompact
                write(test%disp%unit,"(*(g0,:,', '))") "weight             ", weight
                write(test%disp%unit,"(*(g0,:,', '))") "shape(arrayVerbose)", shape(arrayVerbose)
                write(test%disp%unit,"(*(g0,:,', '))") "shape(arrayCompact)", shape(arrayCompact)
                write(test%disp%unit,"(*(g0,:,', '))") "size(weight)       ", size(weight)
                write(test%disp%unit,"(*(g0,:,', '))") "csize_ref          ", csize_ref
                write(test%disp%unit,"(*(g0,:,', '))") "csize              ", csize
                write(test%disp%unit,"(*(g0,:,', '))") "dim                ", dim
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
        end subroutine

#undef  GET_SLICE
#undef  IS_EQUAL
#undef  GET_SIZE
#undef  ALL
