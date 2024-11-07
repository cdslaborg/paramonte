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
!>  This module contains implementations of the tests of the procedures under the generic interfaces of [pm_arrayComplement](@ref pm_arrayComplement).
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if     getComplementRange_D1_IK_ENABLED
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        character(*, SK), parameter :: PROCEDURE_NAME = "@getComplementRange()"

        integer(IKG)                :: start, stop, step
        integer(IKG), allocatable   :: setA(:), setB(:), Complement(:), Complement_ref(:)

        assertion = .true.

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setA = [integer(IKG) ::]
        Complement_ref = [integer(IKG) ::]
        start = 2_IKG
        stop = 1_IKG
        step = 1_IKG

        Complement = getComplementRange(setA, start, stop, step)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        call setSorted(setA)
        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setA = int([-1,2,0,4,5], IKG)
        Complement_ref = [integer(IKG) ::]
        start = 5_IKG
        stop = 2_IKG
        step = 1_IKG

        Complement = getComplementRange(setA, start, stop, step)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        call setSorted(setA)
        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        start = -2_IKG
        stop = 5_IKG
        step = 1_IKG
        setA = [integer(IKG) ::]
        Complement_ref = getRange(start, stop, step)

        Complement = getComplementRange(setA, start, stop, step)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        call setSorted(setA)
        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setA = int([-1,2,0,4,5], IKG)
        Complement_ref = int([-2,1,3], IKG)
        start = -2_IKG
        stop = 5_IKG
        step = 1_IKG

        Complement = getComplementRange(setA, start, stop, step)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        call setSorted(setA)
        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setA = int([-1,2,0,4,5], IKG)
        Complement_ref = int([-2,1,3,6], IKG)
        start = -2_IKG
        stop = 6_IKG
        step = 1_IKG

        Complement = getComplementRange(setA, start, stop, step)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        call setSorted(setA)
        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setA = int([-1,2,0,4,5], IKG)
        Complement_ref = int([-2,6], IKG)
        start = -2_IKG
        stop = 6_IKG
        step = 2_IKG

        Complement = getComplementRange(setA, start, stop, step)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        call setSorted(setA)
        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setA = int([-1,2,0,4,5], IKG)
        Complement_ref = int([7, 1, -2], IKG)
        start = 7_IKG
        stop = -2_IKG
        step = -3_IKG

        Complement = getComplementRange(setA, start, stop, step)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        call setSorted(setA)
        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setA = int([-1,2,0,4,6], IKG)
        Complement_ref = int([1, 3, 5], IKG)
        start = 0_IKG
        stop = 5_IKG
        step = 1_IKG

        Complement = getComplementRange(setA, start, stop, step)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        call setSorted(setA)
        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setA = int([-1,2,0,4,6], IKG); setA = setA(size(setA):1:-1)
        Complement_ref = int([1, 3, 5], IKG)
        start = 0_IKG
        stop = 5_IKG
        step = 1_IKG

        Complement = getComplementRange(setA, start, stop, step)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        call setSorted(setA)
        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        setA = int([-1,2,0,4,6], IKG); setA = setA(size(setA):1:-1)
        Complement_ref = [integer(IKG)::]
        start = 0_IKG
        stop = 6_IKG
        step = 2_IKG

        Complement = getComplementRange(setA, start, stop, step)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        call setSorted(setA)
        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
        call report(int(__LINE__, IK))

        Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
        call report(int(__LINE__, IK))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        block
            use pm_distUnif, only: getUnifRand
            integer(IK) :: i
            do i = 1, 255_IK

                start = getUnifRand(-53_IKG, 53_IKG)
                stop = getUnifRand(-53_IKG, 53_IKG)
                do
                    step = getUnifRand(-53_IKG, 53_IKG)
                    if (step /= 0_IKG) exit
                end do
                setA = getUnifRand(-53_IKG, 53_IKG, getUnifRand(0_IK, 126_IK))
                setB = getRange(start, stop, step)
                Complement_ref = getComplement(setA, setB)

                Complement = getComplementRange(setA, start, stop, step)
                call report(int(__LINE__, IK))

                Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .false._LK)
                call report(int(__LINE__, IK))

                Complement = getComplementRange(setA, start, stop, step, sorted = .false._LK, unique = .true._LK)
                call report(int(__LINE__, IK))

                call setSorted(setA)
                if (step < 0_IKG) setA = setA(size(setA):1:-1)

                Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .false._LK)
                call report(int(__LINE__, IK))

                setA = getUnique(setA)
                Complement = getComplementRange(setA, start, stop, step, sorted = .true._LK, unique = .true._LK)
                call report(int(__LINE__, IK))

            end do
        end block

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report(line)
            integer(IK), intent(in) :: line
            integer(IK) :: i
            assertion = assertion .and. size(Complement, kind = IK) == size(Complement_ref, kind = IK)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")
                write(test%disp%unit,"(*(g0,:,', '))") "size(Complement_ref)   ", size(Complement_ref, kind = IK)
                write(test%disp%unit,"(*(g0,:,', '))") "size(Complement    )   ", size(Complement    , kind = IK)
                write(test%disp%unit,"(*(g0,:,', '))") "Complement_ref         ", Complement_ref
                write(test%disp%unit,"(*(g0,:,', '))") "Complement             ", Complement
                if (allocated(setB)) then
                write(test%disp%unit,"(*(g0,:,', '))") "setB                   ", setB
                write(test%disp%unit,"(*(g0,:,', '))") "size(setB)             ", size(setB)
                end if
                write(test%disp%unit,"(*(g0,:,', '))") "size(setA)             ", size(setA)
                write(test%disp%unit,"(*(g0,:,', '))") "setA                   ", setA
                write(test%disp%unit,"(*(g0,:,', '))") "start                  ", start
                write(test%disp%unit,"(*(g0,:,', '))") "stop                   ", stop
                write(test%disp%unit,"(*(g0,:,', '))") "step                   ", step 
                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, PROCEDURE_NAME//SK_": The size of the complement of setA with respect to specified range must be computed correctly.", line)

            do i = 1, size(Complement)
                assertion = assertion .and. Complement(i) == Complement_ref(i)
                if (test%traceable .and. .not. assertion) then
                    ! LCOV_EXCL_START
                    write(test%disp%unit,"(*(g0,:,', '))")
                    write(test%disp%unit,"(*(g0,:,', '))") "start          ", start
                    write(test%disp%unit,"(*(g0,:,', '))") "stop           ", stop
                    write(test%disp%unit,"(*(g0,:,', '))") "step           ", step
                    write(test%disp%unit,"(*(g0,:,', '))") "setA           ", setA
                    write(test%disp%unit,"(*(g0,:,', '))") "Complement_ref ", Complement_ref
                    write(test%disp%unit,"(*(g0,:,', '))") "Complement     ", Complement
                    write(test%disp%unit,"(*(g0,:,', '))")
                    ! LCOV_EXCL_STOP
                end if
                call test%assert(assertion, PROCEDURE_NAME//SK_": The complement of setA with respect to specified range must be computed correctly.", line)
            end do

        end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#else
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
#define IS_EQUAL .eqv.
#else
#define IS_EQUAL ==
#endif

#if     getComplement_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@getComplement()"
#elif   setComplement_ENABLED
        character(*, SK), parameter :: PROCEDURE_NAME = "@setComplement()"
#endif

#if     getComplement_D0_SK_ENABLED || setComplement_D0_SK_ENABLED
#define ALL
        character(:,SKG), allocatable :: setA, setB, Complement, Complement_ref
#elif   getComplement_D1_SK_ENABLED || setComplement_D1_SK_ENABLED
        character(2,SKG), dimension(:), allocatable :: setA, setB, Complement, Complement_ref
#elif   getComplement_D1_IK_ENABLED || setComplement_D1_IK_ENABLED
        integer(IKG)    , dimension(:), allocatable :: setA, setB, Complement, Complement_ref
#elif   getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
        logical(LKG)    , dimension(:), allocatable :: setA, setB, Complement, Complement_ref
#elif   getComplement_D1_CK_ENABLED || setComplement_D1_CK_ENABLED
        complex(CKG)    , dimension(:), allocatable :: setA, setB, Complement, Complement_ref
#elif   getComplement_D1_RK_ENABLED || setComplement_D1_RK_ENABLED
        real(RKG)       , dimension(:), allocatable :: setA, setB, Complement, Complement_ref
#else
#error  "Unrecognized interface."
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        assertion = .true._LK
        call runTestsWith()
        call runTestsWith(iseq = iseq)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine runTestsWith (iseq)

            logical(LK) , external   , optional :: iseq

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getComplement_D0_SK_ENABLED || setComplement_D0_SK_ENABLED
            setA = ""
            setB = ""
            Complement_ref = ""
#elif       getComplement_D1_SK_ENABLED
            allocate(character(2,SKG) :: setA(0), Complement_ref(0), setB(0))
#elif       setComplement_D1_SK_ENABLED
            allocate(setA(0), Complement_ref(0), setB(0))
#elif       getComplement_D1_IK_ENABLED || setComplement_D1_IK_ENABLED
            allocate(setA(0), Complement_ref(0), setB(0))
#elif       getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
            allocate(setA(0), Complement_ref(0), setB(0))
#elif       getComplement_D1_CK_ENABLED || setComplement_D1_CK_ENABLED
            allocate(setA(0), Complement_ref(0), setB(0))
#elif       getComplement_D1_RK_ENABLED || setComplement_D1_RK_ENABLED
            allocate(setA(0), Complement_ref(0), setB(0))
#endif

            call report(iseq)
            call report(iseq, sorted = .true._LK, unique = .true._LK)
            call report(iseq, sorted = .true._LK, unique = .false._LK)
            call report(iseq, sorted = .false._LK, unique = .true._LK)
            call report(iseq, sorted = .false._LK, unique = .false._LK)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getComplement_D0_SK_ENABLED || setComplement_D0_SK_ENABLED
            setA = "ABCD"
            setB = "ABCD"
            Complement_ref = ""
#elif       getComplement_D1_SK_ENABLED
            setA = ["AA", "BB", "CC", "DD"]
            setB = ["AA", "BB", "CC", "DD"]
            allocate(Complement_ref(0))
#elif       getComplement_D1_IK_ENABLED || setComplement_D1_IK_ENABLED
            setA = int([1, 2, 3, 4], IKG)
            setB = int([1, 2, 3, 4], IKG)
            allocate(Complement_ref(0))
#elif       getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
            setA = logical([.true., .true., .false., .false.], LKG)
            setB = logical([.true., .true., .false., .false.], LKG)
            allocate(Complement_ref(0))
#elif       getComplement_D1_CK_ENABLED || setComplement_D1_CK_ENABLED
            setA = cmplx([1, 2, 3, 4], kind = CKG)
            setB = cmplx([1, 2, 3, 4], kind = CKG)
            allocate(Complement_ref(0))
#elif       getComplement_D1_RK_ENABLED || setComplement_D1_RK_ENABLED
            setA = real([1, 2, 3, 4], RKG)
            setB = real([1, 2, 3, 4], RKG)
            allocate(Complement_ref(0))
#endif

            call report(iseq)
#if         !(getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED)
            call report(iseq, sorted = .true._LK, unique = .true._LK)
            call report(iseq, sorted = .false._LK, unique = .true._LK)
#endif
            call report(iseq, sorted = .true._LK, unique = .false._LK)
            call report(iseq, sorted = .false._LK, unique = .false._LK)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getComplement_D0_SK_ENABLED || setComplement_D0_SK_ENABLED
            setA = "ABCD"
            setB = ""
#elif       getComplement_D1_SK_ENABLED
            setA = ["AA", "BB", "CC", "DD"]
            allocate(setB(0))
#elif       getComplement_D1_IK_ENABLED || setComplement_D1_IK_ENABLED
            allocate(setB(0))
            setA = int([1, 2, 3, 4], IKG)
#elif       getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
            allocate(setB(0))
            setA = logical([.true., .true., .false., .false.], LKG)
#elif       getComplement_D1_CK_ENABLED || setComplement_D1_CK_ENABLED
            allocate(setB(0))
            setA = cmplx([1, 2, 3, 4], kind = CKG)
#elif       getComplement_D1_RK_ENABLED || setComplement_D1_RK_ENABLED
            allocate(setB(0))
            setA = real([1, 2, 3, 4], RKG)
#endif
            Complement_ref = setB

            call report(iseq)
#if         !(getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED)
            call report(iseq, sorted = .true._LK, unique = .true._LK)
            call report(iseq, sorted = .false._LK, unique = .true._LK)
#endif
            call report(iseq, sorted = .true._LK, unique = .false._LK)
            call report(iseq, sorted = .false._LK, unique = .false._LK)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getComplement_D0_SK_ENABLED || setComplement_D0_SK_ENABLED
            setA = ""
            setB = "ABCD"
#elif       getComplement_D1_SK_ENABLED
            allocate(setA(0))
            setB = ["AA", "BB", "CC", "DD"]
#elif       getComplement_D1_IK_ENABLED || setComplement_D1_IK_ENABLED
            allocate(setA(0))
            setB = int([1, 2, 3, 4], IKG)
#elif       getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
            allocate(setA(0))
            setB = logical([.true., .true., .false., .false.], LKG)
#elif       getComplement_D1_CK_ENABLED || setComplement_D1_CK_ENABLED
            allocate(setA(0))
            setB = cmplx([1, 2, 3, 4], kind = CKG)
#elif       getComplement_D1_RK_ENABLED || setComplement_D1_RK_ENABLED
            allocate(setA(0))
            setB = real([1, 2, 3, 4], RKG)
#endif
            Complement_ref = setB

            call report(iseq)
#if         !(getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED)
            call report(iseq, sorted = .true._LK, unique = .true._LK)
            call report(iseq, sorted = .false._LK, unique = .true._LK)
#endif
            call report(iseq, sorted = .true._LK, unique = .false._LK)
            call report(iseq, sorted = .false._LK, unique = .false._LK)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getComplement_D0_SK_ENABLED || setComplement_D0_SK_ENABLED
            setA = "ABCD"
            setB = "ABCDE"
            Complement_ref = setB(len(setB):len(setB))
#elif       getComplement_D1_SK_ENABLED
            setA = ["AA", "BB", "CC", "DD"]
            setB = ["AA", "BB", "CC", "DD", "EE"]
            Complement_ref = setB(size(setB):size(setB))
#elif       getComplement_D1_IK_ENABLED || setComplement_D1_IK_ENABLED
            setA = int([1, 2, 3, 4], IKG)
            setB = int([1, 2, 3, 4, 5], IKG)
            Complement_ref = setB(size(setB):size(setB))
#elif       getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
            setA = logical([.true., .true., .false., .false.], LKG)
            setB = logical([.true., .true., .false., .false., .false.], LKG)
            Complement_ref = [logical(LKG)::]
#elif       getComplement_D1_CK_ENABLED || setComplement_D1_CK_ENABLED
            setA = cmplx([1, 2, 3, 4], kind = CKG)
            setB = cmplx([1, 2, 3, 4, 5], kind = CKG)
            Complement_ref = setB(size(setB):size(setB))
#elif       getComplement_D1_RK_ENABLED || setComplement_D1_RK_ENABLED
            setA = real([1, 2, 3, 4], RKG)
            setB = real([1, 2, 3, 4, 5], RKG)
            Complement_ref = setB(size(setB):size(setB))
#endif

            call report(iseq)
#if         !(getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED)
            call report(iseq, sorted = .true._LK, unique = .true._LK)
            call report(iseq, sorted = .false._LK, unique = .true._LK)
#endif
            call report(iseq, sorted = .true._LK, unique = .false._LK)
            call report(iseq, sorted = .false._LK, unique = .false._LK)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getComplement_D0_SK_ENABLED || setComplement_D0_SK_ENABLED
            setA = "ABCD"
            setB = "!ABCD"
            Complement_ref = setB(1:1)
#elif       getComplement_D1_SK_ENABLED
            setA = ["AA", "BB", "CC", "DD"]
            setB = ["!!", "AA", "BB", "CC", "DD"]
            Complement_ref = setB(1:1)
#elif       getComplement_D1_IK_ENABLED || setComplement_D1_IK_ENABLED
            setA = int([1, 2, 3, 4], IKG)
            setB = int([0, 1, 2, 3, 4], IKG)
            Complement_ref = setB(1:1)
#elif       getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
            setA = logical([.true.], LKG)
            setB = logical([.true., .false.], LKG)
            Complement_ref = logical([.false.], LKG)
#elif       getComplement_D1_CK_ENABLED || setComplement_D1_CK_ENABLED
            setA = cmplx([1, 2, 3, 4], kind = CKG)
            setB = cmplx([0, 1, 2, 3, 4], kind = CKG)
            Complement_ref = setB(1:1)
#elif       getComplement_D1_RK_ENABLED || setComplement_D1_RK_ENABLED
            setA = real([1, 2, 3, 4], RKG)
            setB = real([0, 1, 2, 3, 4], RKG)
            Complement_ref = setB(1:1)
#endif

            call report(iseq)
            call report(iseq, sorted = .true._LK, unique = .true._LK)
            call report(iseq, sorted = .false._LK, unique = .true._LK)
            call report(iseq, sorted = .true._LK, unique = .false._LK)
            call report(iseq, sorted = .false._LK, unique = .false._LK)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call reset()

#if         getComplement_D0_SK_ENABLED || setComplement_D0_SK_ENABLED
            setA = "ABCD"
            setB = "!A?BCDEE"
            Complement_ref = "!?EE"
#elif       getComplement_D1_SK_ENABLED
            setA = ["AA", "BB", "CC", "DD"]
            setB = ["!!", "AA", "??", "BB", "CC", "DD", "EE", "EE"]
            Complement_ref = ["!!", "??", "EE", "EE"]
#elif       getComplement_D1_IK_ENABLED || setComplement_D1_IK_ENABLED
            setA = int([1, 2, 3, 4], IKG)
            setB = int([0, 2, 1, 3, 4, 5, 6, 6], IKG)
            Complement_ref = int([0, 5, 6, 6], IKG)
#elif       getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
            setA = logical([.true.], LKG)
            setB = logical([.true., .true., .false., .false., .false.], LKG)
            Complement_ref = logical([.false., .false., .false.], LKG)
#elif       getComplement_D1_CK_ENABLED || setComplement_D1_CK_ENABLED
            setA = cmplx([1, 2, 3, 4], kind = CKG)
            setB = cmplx([0, 2, 1, 3, 4, 5, 6, 6], kind = CKG)
            Complement_ref = cmplx([0, 5, 6, 6], kind = CKG)
#elif       getComplement_D1_RK_ENABLED || setComplement_D1_RK_ENABLED
            setA = real([1, 2, 3, 4], RKG)
            setB = real([0, 2, 1, 3, 4, 5, 6, 6], RKG)
            Complement_ref = real([0, 5, 6, 6], kind = RKG)
#endif

            call report(iseq)
            call report(iseq, sorted = .false._LK, unique = .false._LK)
            call report(iseq, sorted = .false._LK, unique = .true._LK)
            call setSorted(setA)
            call setSorted(setB)
            call report(iseq, sorted = .true._LK, unique = .false._LK)
            setA = getUnique(setA)
            setB = getUnique(setB)
#if         getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
            Complement_ref = Complement_ref(1:1)
#else
            Complement_ref = Complement_ref(1:3)
#endif
            call report(iseq, sorted = .true._LK, unique = .true._LK)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function iseq(elementA, elementB) result(equivalent)
            logical(LK) :: equivalent
#if         getComplement_D0_SK_ENABLED || setComplement_D0_SK_ENABLED
            character(*,SKG), intent(in) :: elementA, elementB
#elif       getComplement_D1_SK_ENABLED || setComplement_D1_SK_ENABLED
            character(*,SKG), intent(in) :: elementA, elementB
#elif       getComplement_D1_IK_ENABLED || setComplement_D1_IK_ENABLED
            integer(IKG)    , intent(in) :: elementA, elementB
#elif       getComplement_D1_LK_ENABLED || setComplement_D1_LK_ENABLED
            logical(LKG)    , intent(in) :: elementA, elementB
#elif       getComplement_D1_CK_ENABLED || setComplement_D1_CK_ENABLED
            complex(CKG)    , intent(in) :: elementA, elementB
#elif       getComplement_D1_RK_ENABLED || setComplement_D1_RK_ENABLED
            real(RKG)       , intent(in) :: elementA, elementB
#endif
            equivalent = elementA IS_EQUAL elementB
        end function

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine report   ( iseq & ! LCOV_EXCL_LINE
                            , sorted & ! LCOV_EXCL_LINE
                            , unique & ! LCOV_EXCL_LINE
                            )

            logical(LK) , external  , optional  :: iseq
            logical(LK) , intent(in), optional  :: sorted
            logical(LK) , intent(in), optional  :: unique

            if (present(sorted) .neqv. present(unique)) error stop "The condition `present(sorted) .neqv. present(unique)` must hold."

            if (present(iseq)) then
                if (present(sorted)) then
#if                 getComplement_ENABLED
                    Complement = getComplement(setA, setB, sorted, unique, iseq)
#elif               setComplement_ENABLED
#error              "Unrecognized interface."
#endif
                else
#if                 getComplement_ENABLED
                    Complement = getComplement(setA, setB, iseq)
#elif               setComplement_ENABLED
#endif
                end if
            else
                if (present(sorted)) then
#if                 getComplement_ENABLED
                    Complement = getComplement(setA, setB, sorted, unique)
#elif               setComplement_ENABLED
#error              "Unrecognized interface."
#endif
                else
#if                 getComplement_ENABLED
                    Complement = getComplement(setA, setB)
#elif               setComplement_ENABLED
#endif
                end if
            end if

            ! Report test results if needed.

            assertion = assertion .and. ALL(Complement IS_EQUAL Complement_ref)
            if (test%traceable .and. .not. assertion) then

                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0,:,', '))")

                write(test%disp%unit,"(*(g0,:,', '))") "setB               ", Complement
                write(test%disp%unit,"(*(g0,:,', '))") "setA               ", Complement
                write(test%disp%unit,"(*(g0,:,', '))") "Complement         ", Complement_ref
                write(test%disp%unit,"(*(g0,:,', '))") "present(unique)    ", present(unique)
                write(test%disp%unit,"(*(g0,:,', '))") "present(sorted)    ", present(sorted)
                write(test%disp%unit,"(*(g0,:,', '))") "present(iseq)      ", present(iseq)

                if (present(sorted)) then
                write(test%disp%unit,"(*(g0,:,', '))") "sorted             ", sorted
                end if

                if (present(unique)) then
                write(test%disp%unit,"(*(g0,:,', '))") "unique             ", unique
                end if

                write(test%disp%unit,"(*(g0,:,', '))")
                ! LCOV_EXCL_STOP

            end if

            call test%assert(assertion, PROCEDURE_NAME//SK_": The complement of setA with respect to setB must be properly set.", int(__LINE__, IK))

        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reset()
            if (allocated(setA)) deallocate(setA)
            if (allocated(setB)) deallocate(setB)
            if (allocated(Complement)) deallocate(Complement)
            if (allocated(Complement_ref)) deallocate(Complement_ref)
        end subroutine reset

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#endif

#undef  IS_EQUAL
#undef  ALL
