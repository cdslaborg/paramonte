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

!>  \brief This module contains tests of the module [pm_err](@ref pm_err).
!>  \author Amir Shahmoradi

module test_pm_err

    use pm_err
    use pm_io, only: TAB
    use pm_test, only: test_type, LK
    use pm_strASCII, only: LF
    implicit none

    private
    public :: setTest
    type(test_type) :: test

    character(*, SK), parameter :: mc_prefix =  TAB//TAB//"ParaMonte"
    character(*, SK), parameter :: mc_msg    =  "What does a fish know about the water in which it swims all its life?    Albert Einstein\n" // &
                                            "Everything should be made as simple as possible, but not simpler.    Albert Einstein\n" // &
                                            "The absence of evidence is not evidence for absence.    Carl Sagan\n" // &
                                            "If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton\n" // &
                                            "I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)
        call test%run(test_note_1, SK_"test_note_1")
        call test%run(test_note_2, SK_"test_note_2")
        call test%run(test_warn_1, SK_"test_warn_1")
        call test%run(test_warn_2, SK_"test_warn_2")
        call test%run(test_abort_1, SK_"test_abort_1")
        call test%run(test_abort_2, SK_"test_abort_2")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_note_1() result(assertion)

        use pm_container, only: css_pdt
        use pm_kind, only: IK
        use pm_val2str, only: getStr

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical(LK) :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 5_IK
        type(css_pdt) , allocatable :: outputList_ref(:)
        type(css_pdt) , allocatable :: outputList(:)

        assertion = .true._LK

        if (allocated(outputList)) deallocate(outputList); allocate(outputList(NLINE))
        if (allocated(outputList_ref)) deallocate(outputList_ref); allocate(outputList_ref(NLINE))

        outputList_ref(1)%val = "        ParaMonte - NOTE: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        outputList_ref(2)%val = "        ParaMonte - NOTE: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        outputList_ref(3)%val = "        ParaMonte - NOTE: The absence of evidence is not evidence for absence.    Carl Sagan"
        outputList_ref(4)%val = "        ParaMonte - NOTE: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        outputList_ref(5)%val = "        ParaMonte - NOTE: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = test%dir%out//"/test_pm_err@test_note_1."//getStr(test%image%id)//".out", status = "replace")

        call note   ( msg = mc_msg &
                    , prefix  = mc_prefix &
                    , newline  = "\n" &
                    , outputUnit = fileUnit &
                    , marginTop = 0_IK &
                    , marginBot = 0_IK &
                    )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(outputList(i)%val)) deallocate(outputList(i)%val)
            allocate(character(132) :: outputList(i)%val)

            read(fileUnit,"(A132)", iostat = iostat) outputList(i)%val
            assertion = iostat == 0_IK
            call test%assert(assertion, desc = "`iostat == 0_IK`.")

            outputList(i)%val = trim(outputList(i)%val)

            assertionCurrent = outputList(i)%val == outputList_ref(i)%val
            assertion = assertion .and. assertionCurrent

            if (test%traceable .and. .not. assertionCurrent) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0))")
                write(test%disp%unit,"(*(g0))") "outputList_ref(",getStr(i),")%val = ", outputList_ref(i)%val
                write(test%disp%unit,"(*(g0))") "outputList    (",getStr(i),")%val = ", outputList(i)%val
                write(test%disp%unit,"(*(g0))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, desc = "`outputList(i)%val == outputList_ref(i)%val`.")

        end do

    end function test_note_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_note_2() result(assertion)

        use pm_container, only: css_pdt
        use pm_kind, only: IK
        use pm_arrayReplace, only: getReplaced
        use pm_val2str, only: getStr

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical(LK) :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 7_IK
        type(css_pdt) , allocatable :: outputList_ref(:)
        type(css_pdt) , allocatable :: outputList(:)

        assertion = .true._LK

        if (allocated(outputList)) deallocate(outputList); allocate(outputList(NLINE))
        if (allocated(outputList_ref)) deallocate(outputList_ref); allocate(outputList_ref(NLINE))

        outputList_ref(1)%val = ""
        outputList_ref(2)%val = " - NOTE: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        outputList_ref(3)%val = " - NOTE: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        outputList_ref(4)%val = " - NOTE: The absence of evidence is not evidence for absence.    Carl Sagan"
        outputList_ref(5)%val = " - NOTE: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        outputList_ref(6)%val = " - NOTE: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"
        outputList_ref(7)%val = ""

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = test%dir%out//"/test_pm_err@test_note_2."//getStr(test%image%id)//".out", status = "replace")

        call note   ( msg = getReplaced(mc_msg, "\n", LF) &
                    , outputUnit = fileUnit &
                    , newline  = LF &
                    )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(outputList(i)%val)) deallocate(outputList(i)%val)
            allocate(character(132) :: outputList(i)%val)

            read(fileUnit,"(A132)", iostat = iostat) outputList(i)%val
            assertion = iostat == 0_IK
            call test%assert(assertion, desc = "`iostat == 0_IK`.")

            outputList(i)%val = trim(outputList(i)%val)

            assertionCurrent = outputList(i)%val == outputList_ref(i)%val
            assertion = assertion .and. assertionCurrent

            if (test%traceable .and. .not. assertionCurrent) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0))")
                write(test%disp%unit,"(*(g0))") "outputList_ref(",getStr(i),")%val = ", outputList_ref(i)%val
                write(test%disp%unit,"(*(g0))") "outputList    (",getStr(i),")%val = ", outputList(i)%val
                write(test%disp%unit,"(*(g0))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, desc = "`outputList(i)%val == outputList_ref(i)%val`.")

        end do

    end function test_note_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_warn_1() result(assertion)

        use pm_container, only: css_pdt
        use pm_kind, only: IK
        use pm_val2str, only: getStr

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical(LK) :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 5_IK
        type(css_pdt) , allocatable :: outputList_ref(:)
        type(css_pdt) , allocatable :: outputList(:)

        assertion = .true._LK

        if (allocated(outputList)) deallocate(outputList); allocate(outputList(NLINE))
        if (allocated(outputList_ref)) deallocate(outputList_ref); allocate(outputList_ref(NLINE))

        outputList_ref(1)%val = "        ParaMonte - WARNING: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        outputList_ref(2)%val = "        ParaMonte - WARNING: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        outputList_ref(3)%val = "        ParaMonte - WARNING: The absence of evidence is not evidence for absence.    Carl Sagan"
        outputList_ref(4)%val = "        ParaMonte - WARNING: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        outputList_ref(5)%val = "        ParaMonte - WARNING: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = test%dir%out//"/test_pm_err@test_warn_1."//getStr(test%image%id)//".out", status = "replace")

        call setWarned  ( msg = mc_msg &
                        , prefix  = mc_prefix &
                        , newline  = "\n" &
                        , outputUnit = fileUnit &
                        , marginTop = 0_IK &
                        , marginBot = 0_IK &
                        )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(outputList(i)%val)) deallocate(outputList(i)%val)
            allocate(character(132) :: outputList(i)%val)

            read(fileUnit,"(A132)", iostat = iostat) outputList(i)%val
            assertion = iostat == 0_IK
            call test%assert(assertion, desc = "`iostat == 0_IK`.")

            outputList(i)%val = trim(outputList(i)%val)

            assertionCurrent = outputList(i)%val == outputList_ref(i)%val
            assertion = assertion .and. assertionCurrent

            if (test%traceable .and. .not. assertionCurrent) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0))")
                write(test%disp%unit,"(*(g0))") "outputList_ref(",getStr(i),")%val = ", outputList_ref(i)%val
                write(test%disp%unit,"(*(g0))") "outputList    (",getStr(i),")%val = ", outputList(i)%val
                write(test%disp%unit,"(*(g0))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, desc = "`outputList(i)%val == outputList_ref(i)%val`.")

        end do

    end function test_warn_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_warn_2() result(assertion)

        use pm_container, only: css_pdt
        use pm_kind, only: IK
        use pm_arrayReplace, only: getReplaced
        use pm_val2str, only: getStr

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical(LK) :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 7_IK
        type(css_pdt) , allocatable :: outputList_ref(:)
        type(css_pdt) , allocatable :: outputList(:)

        assertion = .true._LK

        if (allocated(outputList)) deallocate(outputList); allocate(outputList(NLINE))
        if (allocated(outputList_ref)) deallocate(outputList_ref); allocate(outputList_ref(NLINE))

        outputList_ref(1)%val = ""
        outputList_ref(2)%val = " - WARNING: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        outputList_ref(3)%val = " - WARNING: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        outputList_ref(4)%val = " - WARNING: The absence of evidence is not evidence for absence.    Carl Sagan"
        outputList_ref(5)%val = " - WARNING: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        outputList_ref(6)%val = " - WARNING: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"
        outputList_ref(7)%val = ""

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = test%dir%out//"/test_pm_err@test_warn_2."//getStr(test%image%id)//".out", status = "getReplaced")

        call setWarned  ( msg = getReplaced(mc_msg, "\n", LF) &
                        , outputUnit = fileUnit &
                        , newline  = LF &
                        )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(outputList(i)%val)) deallocate(outputList(i)%val)
            allocate(character(132) :: outputList(i)%val)

            read(fileUnit,"(A132)", iostat = iostat) outputList(i)%val
            assertion = iostat == 0_IK
            call test%assert(assertion, desc = "`iostat == 0_IK`.")

            outputList(i)%val = trim(outputList(i)%val)

            assertionCurrent = outputList(i)%val == outputList_ref(i)%val
            assertion = assertion .and. assertionCurrent

            if (test%traceable .and. .not. assertionCurrent) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0))")
                write(test%disp%unit,"(*(g0))") "outputList_ref(",getStr(i),")%val = ", outputList_ref(i)%val
                write(test%disp%unit,"(*(g0))") "outputList    (",getStr(i),")%val = ", outputList(i)%val
                write(test%disp%unit,"(*(g0))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, desc = "`outputList(i)%val == outputList_ref(i)%val`.")

        end do

    end function test_warn_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_abort_1() result(assertion)

        use pm_container, only: css_pdt
        use pm_kind, only: IK
        use pm_val2str, only: getStr

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical(LK) :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 6_IK
        type(css_pdt) , allocatable :: outputList_ref(:)
        type(css_pdt) , allocatable :: outputList(:)
        type(err_type) :: Err

        assertion = .true._LK
        mv_isTestingMode = .true._LK

        if (allocated(outputList)) deallocate(outputList); allocate(outputList(NLINE))
        if (allocated(outputList_ref)) deallocate(outputList_ref); allocate(outputList_ref(NLINE))

        outputList_ref(1)%val = ""
        outputList_ref(2)%val = "        ParaMonte - FATAL: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        outputList_ref(3)%val = "        ParaMonte - FATAL: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        outputList_ref(4)%val = "        ParaMonte - FATAL: The absence of evidence is not evidence for absence.    Carl Sagan"
        outputList_ref(5)%val = "        ParaMonte - FATAL: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        outputList_ref(6)%val = "        ParaMonte - FATAL: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = test%dir%out//"/test_pm_err@test_abort_1."//getStr(test%image%id)//".out", status = "replace")

        err%msg = mc_msg
        call setAborted ( Err = Err & ! LCOV_EXCL_LINE
                        , newline  = "\n" & ! LCOV_EXCL_LINE
                        , prefix  = mc_prefix & ! LCOV_EXCL_LINE
                        , outputUnit = fileUnit & ! LCOV_EXCL_LINE
                        , renabled = .true._LK & ! LCOV_EXCL_LINE
                        )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(outputList(i)%val)) deallocate(outputList(i)%val)
            allocate(character(132) :: outputList(i)%val)

            read(fileUnit,"(A132)", iostat = iostat) outputList(i)%val
            assertion = iostat == 0_IK
            call test%assert(assertion, desc = "`iostat == 0_IK`.")

            outputList(i)%val = trim(outputList(i)%val)

            assertionCurrent = outputList(i)%val == outputList_ref(i)%val
            assertion = assertion .and. assertionCurrent

            if (test%traceable .and. .not. assertionCurrent) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0))")
                write(test%disp%unit,"(*(g0))") "outputList_ref(",getStr(i),")%val = ", outputList_ref(i)%val
                write(test%disp%unit,"(*(g0))") "outputList    (",getStr(i),")%val = ", outputList(i)%val
                write(test%disp%unit,"(*(g0))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, desc = "`outputList(i)%val == outputList_ref(i)%val`.")

        end do

    end function test_abort_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the effects of an input non-null error code `err%stat`.
    !> Test the effects of missing arguments `prefix`, `returnEnabled`, and `newline`.
    function test_abort_2() result(assertion)

        use pm_container, only: css_pdt
        use pm_kind, only: IK
        use pm_val2str, only: getStr
        use pm_arrayReplace, only: getReplaced

        implicit none

        integer(IK) :: fileUnit, i, iostat
        logical(LK) :: assertion, assertionCurrent
        integer(IK), parameter :: NLINE = 7_IK
        type(css_pdt) , allocatable :: outputList_ref(:)
        type(css_pdt) , allocatable :: outputList(:)
        type(err_type) :: Err

        assertion = .true._LK
        mv_isTestingMode = .true._LK

        if (allocated(outputList)) deallocate(outputList); allocate(outputList(NLINE))
        if (allocated(outputList_ref)) deallocate(outputList_ref); allocate(outputList_ref(NLINE))

        outputList_ref(1)%val = ""
        outputList_ref(2)%val = " - FATAL: What does a fish know about the water in which it swims all its life?    Albert Einstein"
        outputList_ref(3)%val = " - FATAL: Everything should be made as simple as possible, but not simpler.    Albert Einstein"
        outputList_ref(4)%val = " - FATAL: The absence of evidence is not evidence for absence.    Carl Sagan"
        outputList_ref(5)%val = " - FATAL: If I have seen further, it is by standing on the shoulders of giants.    Isaac Newton"
        outputList_ref(6)%val = " - FATAL: I don't pretend to understand the universe - it's much bigger than I am.    Thomas Carlyle"
        outputList_ref(7)%val = " - FATAL: Error Code: 123."

        open(newunit = fileUnit, status = "scratch")
        !open(newunit = fileUnit, file = test%dir%out//"/test_pm_err@test_abort_2."//getStr(test%image%id)//".out", status = "replace")

        err%msg = getReplaced(mc_msg, "\n", LF)
        err%stat = 123_IK
        call setAborted ( Err = Err &
                        , outputUnit = fileUnit &
                        )

        rewind(fileUnit)

        do i = 1, NLINE

            if(allocated(outputList(i)%val)) deallocate(outputList(i)%val)
            allocate(character(132) :: outputList(i)%val)

            read(fileUnit,"(A132)", iostat = iostat) outputList(i)%val
            assertion = iostat == 0_IK
            call test%assert(assertion, desc = "`iostat == 0_IK`.")

            outputList(i)%val = trim(outputList(i)%val)

            assertionCurrent = outputList(i)%val == outputList_ref(i)%val
            assertion = assertion .and. assertionCurrent

            if (test%traceable .and. .not. assertionCurrent) then
                ! LCOV_EXCL_START
                write(test%disp%unit,"(*(g0))")
                write(test%disp%unit,"(*(g0))") "outputList_ref(",getStr(i),")%val = '", outputList_ref(i)%val, "'"
                write(test%disp%unit,"(*(g0))") "outputList    (",getStr(i),")%val = '", outputList(i)%val, "'"
                write(test%disp%unit,"(*(g0))")
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, desc = "`outputList(i)%val == outputList_ref(i)%val`.")

        end do

    end function test_abort_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_err ! LCOV_EXCL_LINE