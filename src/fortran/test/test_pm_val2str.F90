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

!>  \brief This module contains tests of the module [pm_val2Str](@ref pm_val2str).
!>  \author Amir Shahmoradi

module test_pm_val2str

    use pm_val2str
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface

    module function test_getStr_SK_1 () result(assertion); logical(LK) :: assertion; end function
    module function test_getStr_LK_1 () result(assertion); logical(LK) :: assertion; end function
#if RK3_ENABLED
    module function test_getStr_RK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK2_ENABLED
    module function test_getStr_RK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if RK1_ENABLED
    module function test_getStr_RK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK3_ENABLED
    module function test_getStr_CK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK2_ENABLED
    module function test_getStr_CK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if CK1_ENABLED
    module function test_getStr_CK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if IK4_ENABLED
    module function test_getStr_IK4_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if IK3_ENABLED
    module function test_getStr_IK3_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if IK2_ENABLED
    module function test_getStr_IK2_1 () result(assertion); logical(LK) :: assertion; end function
#endif
#if IK1_ENABLED
    module function test_getStr_IK1_1 () result(assertion); logical(LK) :: assertion; end function
#endif

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%run(test_getStr_SK_1  , SK_"test_getStr_SK_1")
        call test%run(test_getStr_LK_1  , SK_"test_getStr_LK_1")

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     RK3_ENABLED
        call test%run(test_getStr_RK3_1, SK_"test_getStr_RK3_1")
#endif
#if     RK2_ENABLED
        call test%run(test_getStr_RK2_1, SK_"test_getStr_RK2_1")
#endif
#if     RK1_ENABLED
        call test%run(test_getStr_RK1_1, SK_"test_getStr_RK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK3_ENABLED
        call test%run(test_getStr_CK3_1, SK_"test_getStr_CK3_1")
#endif
#if     CK2_ENABLED
        call test%run(test_getStr_CK2_1, SK_"test_getStr_CK2_1")
#endif
#if     CK1_ENABLED
        call test%run(test_getStr_CK1_1, SK_"test_getStr_CK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK4_ENABLED
        call test%run(test_getStr_IK4_1, SK_"test_getStr_IK4_1")
#endif
#if     IK3_ENABLED
        call test%run(test_getStr_IK3_1, SK_"test_getStr_IK3_1")
#endif
#if     IK2_ENABLED
        call test%run(test_getStr_IK2_1, SK_"test_getStr_IK2_1")
#endif
#if     IK1_ENABLED
        call test%run(test_getStr_IK1_1  , SK_"test_getStr_IK1_1")
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    function test_getStr_IK3_1() result(assertion)
!        use pm_kind, only: IK => IK3
!        implicit none
!        logical                     :: assertion
!        character(*), parameter     :: this_ref = "12345"
!        integer(IK) , parameter     :: this_int = 12345_IK
!        character(:), allocatable   :: this
!        this = getStr(this_int)
!        assertion = this == this_ref .and. len(this) == len(this_ref)
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = '", this_ref, "'"
!            write(test%disp%unit,"(*(g0))") "this_str  = " , this_int
!            write(test%disp%unit,"(*(g0))") "this      = '", this, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_getStr_IK3_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_IK3_2() result(assertion)
!        use pm_kind, only: IK => IK3
!        implicit none
!        logical                     :: assertion
!        character(*), parameter     :: this_ref = "0000012345"
!        integer(IK) , parameter     :: this_int = 12345_IK
!        character(:), allocatable   :: this
!        this = getStr(this_int,"(1I10.10)")
!        assertion = this == this_ref .and. len(this) == len(this_ref)
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = '", this_ref, "'"
!            write(test%disp%unit,"(*(g0))") "this_str  = " , this_int
!            write(test%disp%unit,"(*(g0))") "this      = '", this, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_getStr_IK3_2
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_IK3_3() result(assertion)
!        use pm_kind, only: IK => IK3
!        implicit none
!        logical                     :: assertion
!        character(*), parameter     :: this_ref = "12345     "
!        integer(IK) , parameter     :: this_int = 12345_IK
!        character(:), allocatable   :: this
!        this = getStr(this_int,len=10_IK)
!        assertion = this == this_ref .and. len(this) == len(this_ref)
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = '", this_ref, "'"
!            write(test%disp%unit,"(*(g0))") "this_str  = " , this_int
!            write(test%disp%unit,"(*(g0))") "this      = '", this, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_getStr_IK3_3
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_IK4_1() result(assertion)
!        use pm_kind, only: IK4
!        logical                         :: assertion
!        character(*, SK), parameter     :: this_ref = "98765432112345"
!        integer(IK4)   , parameter     :: this_int = 98765432112345_IK4
!        character(:, SK), allocatable   :: this
!        this = getStr(this_int)
!        assertion = this == this_ref .and. len(this) == len(this_ref)
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = '", this_ref, "'"
!            write(test%disp%unit,"(*(g0))") "this_str  = " , this_int
!            write(test%disp%unit,"(*(g0))") "this      = '", this, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_getStr_IK4_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_IK4_2() result(assertion)
!        use pm_kind, only: IK4
!        implicit none
!        logical                         :: assertion
!        character(*, SK), parameter     :: this_ref = "00000098765432112345"
!        integer(IK4)   , parameter     :: this_int = 98765432112345_IK4
!        character(:, SK), allocatable   :: this
!        this = getStr(this_int,"(1I20.20)")
!        assertion = this == this_ref .and. len(this) == len(this_ref)
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = '", this_ref, "'"
!            write(test%disp%unit,"(*(g0))") "this_str  = " , this_int
!            write(test%disp%unit,"(*(g0))") "this      = '", this, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_getStr_IK4_2
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_IK4_3() result(assertion)
!        use pm_kind, only: IK4
!        logical                         :: assertion
!        character(*, SK), parameter     :: this_ref = "98765432112345      "
!        integer(IK4)   , parameter     :: this_int = 98765432112345_IK4
!        character(:, SK), allocatable   :: this
!        this = getStr(this_int,len=20_IK)
!        assertion = this == this_ref .and. len(this) == len(this_ref)
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = '", this_ref, "'"
!            write(test%disp%unit,"(*(g0))") "this_str  = " , this_int
!            write(test%disp%unit,"(*(g0))") "this      = '", this, "'"
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_IK4_3
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK1_1() result(assertion)
!        use pm_kind, only: RK => RK1
!        use pm_str2val, only: real32
!        logical                     :: assertion
!        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
!        real(RK)                    :: this
!        this = real32( getStr(this_ref) )
!        assertion = this == this_ref
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK1_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK1_2() result(assertion)
!        use pm_kind, only: RK => RK1
!        use pm_str2val, only: real32
!        logical                     :: assertion
!        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
!        real(RK)                    :: this
!        this = real32( getStr(this_ref,"(g0)") )
!        assertion = this == this_ref
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK1_2
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK1_3() result(assertion)
!        use pm_kind, only: RK => RK1
!        use pm_str2val, only: real32
!        implicit none
!        logical                     :: assertion
!        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
!        real(RK)                    :: this
!        this = real32( getStr(this_ref,"(g0)",20) )
!        assertion = this == this_ref
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_getStr_RK1_3
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !> The input `len` can be larger than the length of the constructed string.
!    function test_getStr_RK1_4() result(assertion)
!        use pm_kind, only: RK => RK1, IK => IK3
!        use pm_str2val, only: real32
!        logical                     :: assertion
!        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
!        real(RK)                    :: this
!        this = real32( getStr(this_ref,"(g0)",263_IK) )
!        assertion = this == this_ref
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK1_4
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK2_1() result(assertion)
!        use pm_kind, only: RK => RK2
!        use pm_str2val, only: real64
!        logical                     :: assertion
!        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
!        real(RK)                    :: this
!        this = real64( getStr(this_ref) )
!        assertion = this == this_ref
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK2_2() result(assertion)
!        use pm_kind, only: RK => RK2
!        use pm_str2val, only: real64
!        logical                     :: assertion
!        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
!        real(RK)                    :: this
!        this = real64( getStr(this_ref,"(g0)") )
!        assertion = this == this_ref
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_2
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK2_3() result(assertion)
!        use pm_kind, only: RK => RK2, IK => IK3
!        use pm_str2val, only: real64
!        logical                     :: assertion
!        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
!        real(RK)                    :: this
!        this = real64( getStr(this_ref,"(g0)",30_IK) )
!        assertion = this == this_ref
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_3
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !> The input `len` can be larger than the length of the constructed string.
!    function test_getStr_RK2_4() result(assertion)
!        use pm_kind, only: RK => RK2, IK => IK3
!        use pm_str2val, only: real64
!        logical                     :: assertion
!        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
!        real(RK)                    :: this
!        this = real64( getStr(this_ref,"(g0)",263_IK) )
!        assertion = this == this_ref
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_4
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK2_D1_1() result(assertion)
!        use pm_kind, only: RK => RK2
!        logical                     :: assertion
!        integer(IK) , parameter     :: nthis = 2
!        real(RK)    , parameter     :: this_ref(nthis) = [1.23456798e-30_RK, 2.32456798e+30_RK]
!        real(RK)                    :: this(nthis)
!        character(:), allocatable   :: string
!        assertion = .true.
!        string = getStr(this_ref)
!        read(string,*) this
!        assertion = all(this == this_ref)
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_D1_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK2_D1_2() result(assertion)
!        use pm_kind, only: RK => RK2
!        logical                     :: assertion
!        integer(IK) , parameter     :: nthis = 2
!        real(RK)    , parameter     :: this_ref(nthis) = [1.23456798e-30_RK, 2.32456798e+30_RK]
!        real(RK)                    :: this(nthis)
!        character(:), allocatable   :: string
!        assertion = .true.
!        string = getStr(this_ref,"(*(g0,:,' '))")
!        read(string,*) this
!        assertion = all(this == this_ref)
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_D1_2
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK2_D1_3() result(assertion)
!        use pm_kind, only: RK => RK2, IK => IK3
!        logical                     :: assertion
!        integer(IK) , parameter     :: nthis = 2
!        real(RK)    , parameter     :: this_ref(nthis) = [1.23456798e-30_RK, 2.32456798e+30_RK]
!        real(RK)                    :: this(nthis)
!        character(:), allocatable   :: string
!        assertion = .true.
!        string = getStr(this_ref,"(*(g0,:,' '))",63_IK)
!        read(string,*) this
!        assertion = all(this == this_ref)
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_getStr_RK2_D1_3
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !> The input `len` can be larger than the length of the constructed string.
!    function test_getStr_RK2_D1_4() result(assertion)
!        use pm_kind, only: RK => RK2, IK => IK3
!        logical                     :: assertion
!        integer(IK) , parameter     :: nthis = 2_IK
!        real(RK)    , parameter     :: this_ref(nthis) = [1.23456798e-30_RK, 2.32456798e+30_RK]
!        real(RK)                    :: this(nthis)
!        character(:), allocatable   :: string
!        assertion = .true.
!        string = getStr(this_ref,"(*(g0,:,' '))",263_IK)
!        read(string,*) this
!        assertion = all(this == this_ref)
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_D1_4
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK2_D2_1() result(assertion)
!        use pm_kind, only: RK => RK2
!        logical                     :: assertion
!        integer(IK) , parameter     :: nthis = 2
!        real(RK)    , parameter     :: this_ref(nthis,nthis) = reshape([1.23456798e-30_RK, 2.32456798e+30_RK, 1.23456798e-30_RK, 2.32456798e+30_RK], shape=shape(this_ref))
!        real(RK)                    :: this(nthis,nthis)
!        character(:), allocatable   :: string
!        assertion = .true.
!        string = getStr(this_ref)
!        read(string,*) this
!        assertion = all(this == this_ref)
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_D2_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK2_D2_2() result(assertion)
!        use pm_kind, only: RK => RK2
!        logical                     :: assertion
!        integer(IK) , parameter     :: nthis = 2
!        real(RK)    , parameter     :: this_ref(nthis,nthis) = reshape([1.23456798e-30_RK, 2.32456798e+30_RK, 1.23456798e-30_RK, 2.32456798e+30_RK], shape=shape(this_ref))
!        real(RK)                    :: this(nthis,nthis)
!        character(:), allocatable   :: string
!        assertion = .true.
!        string = getStr(this_ref,"(*(g0,:,' '))")
!        read(string,*) this
!        assertion = all(this == this_ref)
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_D2_2
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_RK2_D2_3() result(assertion)
!        use pm_kind, only: RK => RK2, IK => IK3
!        logical                     :: assertion
!        integer(IK) , parameter     :: nthis = 2
!        real(RK)    , parameter     :: this_ref(nthis,nthis) = reshape([1.23456798e-30_RK, 2.32456798e+30_RK, 1.23456798e-30_RK, 2.32456798e+30_RK], shape=shape(this_ref))
!        real(RK)                    :: this(nthis,nthis)
!        character(:), allocatable   :: string
!        assertion = .true.
!        string = getStr(this_ref,"(*(g0,:,' '))",128_IK)
!        read(string,*) this
!        assertion = all(this == this_ref)
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_D2_3
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !> The input `len` can be larger than the length of the constructed string.
!    function test_getStr_RK2_D2_4() result(assertion)
!        use pm_kind, only: RK => RK2, IK => IK3
!        logical                     :: assertion
!        integer(IK) , parameter     :: nthis = 2
!        real(RK)    , parameter     :: this_ref(nthis,nthis) = reshape([1.23456798e-30_RK, 2.32456798e+30_RK, 1.23456798e-30_RK, 2.32456798e+30_RK], shape=shape(this_ref))
!        real(RK)                    :: this(nthis,nthis)
!        character(:), allocatable   :: string
!        assertion = .true.
!        string = getStr(this_ref,"(*(g0,:,' '))",256_IK)
!        read(string,*) this
!        assertion = all(this == this_ref)
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "this_ref  = ", this_ref
!            write(test%disp%unit,"(*(g0))") "this      = ", this
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_RK2_D2_4
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_LK_1() result(assertion)
!        implicit none
!        logical                     :: assertion
!        character(:), allocatable   :: strout
!        logical                     :: logic
!        logic = .true.
!        strout = getStr(logic)
!        assertion = strout == "T"
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "logic     = ", logic
!            write(test%disp%unit,"(*(g0))") "strout    = ", strout
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_LK_2() result(assertion)
!        implicit none
!        logical                     :: assertion
!        character(:), allocatable   :: strout
!        logical                     :: logic
!        logic = .false.
!        strout = getStr(logic)
!        assertion = strout == "F"
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "logic     = ", logic
!            write(test%disp%unit,"(*(g0))") "strout    = ", strout
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_IntStr_type_1() result(assertion)
!        use pm_kind, only: IK => IK3
!        logical                     :: assertion
!        type(intStr_type)           :: IntStr
!        IntStr%str = getStr(123_IK)
!        assertion = IntStr%str == "123"
!    end function test_IntStr_type_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_IntStr_type_2() result(assertion)
!        use pm_kind, only: IK => IK4
!        logical                     :: assertion
!        type(intStr_type)           :: IntStr
!        IntStr%str = getStr(123_IK)
!        assertion = IntStr%str == "123"
!    end function test_IntStr_type_2
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_1() result(assertion)
!        use pm_kind, only: RK => RK1
!        implicit none
!        logical                     :: assertion
!        type(RealStr_type)          :: RealStr
!        RealStr%str = getStr(123._RK,"(F10.4)",15)
!        assertion = RealStr%str == "  123.0000     "
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "Number=123_RK1"
!            write(test%disp%unit,"(*(g0))") "RealStr%str = '", RealStr%str, "'"
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_getStr_2() result(assertion)
!        use pm_kind, only: RK => RK2
!        implicit none
!        logical                     :: assertion
!        type(RealStr_type)          :: RealStr
!        RealStr%str = getStr(123._RK,"(F10.4)",15)
!        assertion = RealStr%str == "  123.0000     "
!        if (test%traceable .and. .not. assertion) then
!            ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "Number=123_real64"
!            write(test%disp%unit,"(*(g0))") "RealStr%str = '", RealStr%str, "'"
!            write(test%disp%unit,"(*(g0))")
!            ! LCOV_EXCL_STOP
!        end if
!    end function test_getStr_2
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_padVal2Str_1() result(assertion)
!        use pm_kind, only: IK
!        implicit none
!        logical                     :: assertion
!        integer(IK) , parameter     :: paddedLen = 30_IK
!        character(*), parameter     :: string_nonPadded = "ParaMonte"
!        character(*), parameter     :: stringPadded_ref = "ParaMonte....................."
!        character(*), parameter     :: symbol = "."
!        character(:), allocatable   :: charPadded
!        charPadded = padChar(string_nonPadded, symbol, paddedLen)
!        assertion = charPadded == stringPadded_ref .and. len(charPadded) == len(stringPadded_ref)
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "string_nonPadded  : '", string_nonPadded, "'"
!            write(test%disp%unit,"(*(g0))") "stringPadded_ref  : '", stringPadded_ref, "'"
!            write(test%disp%unit,"(*(g0))") "charPadded      : '", charPadded, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_padVal2Str_1
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_padVal2Str_2() result(assertion)
!        use pm_kind, only: IK
!        implicit none
!        logical                     :: assertion
!        integer(IK) , parameter     :: paddedLen = 9_IK
!        character(*), parameter     :: string_nonPadded = "ParaMonte"
!        character(*), parameter     :: stringPadded_ref = "ParaMonte"
!        character(*), parameter     :: symbol = "."
!        character(:), allocatable   :: charPadded
!        charPadded = padChar(string_nonPadded, symbol, paddedLen)
!        assertion = charPadded == stringPadded_ref .and. len(charPadded) == len(stringPadded_ref)
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "string_nonPadded  : '", string_nonPadded, "'"
!            write(test%disp%unit,"(*(g0))") "stringPadded_ref  : '", stringPadded_ref, "'"
!            write(test%disp%unit,"(*(g0))") "charPadded      : '", charPadded, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_padVal2Str_2
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !>  \brief
!    !> When `len(string) > paddedLen`, the full string must be returned without any padding.
!    function test_padVal2Str_3() result(assertion)
!        use pm_kind, only: IK
!        implicit none
!        logical                     :: assertion
!        integer(IK) , parameter     :: paddedLen = 5_IK
!        character(*), parameter     :: string_nonPadded = "ParaMonte"
!        character(*), parameter     :: stringPadded_ref = "ParaM"
!        character(*), parameter     :: symbol = "."
!        character(:), allocatable   :: charPadded
!        charPadded = padChar(string_nonPadded, symbol, paddedLen)
!        assertion = charPadded == stringPadded_ref .and. len(charPadded) == paddedLen .and. charPadded == stringPadded_ref
!        if (test%traceable .and. .not. assertion) then
!        ! LCOV_EXCL_START
!            write(test%disp%unit,"(*(g0))")
!            write(test%disp%unit,"(*(g0))") "string_nonPadded  : '", string_nonPadded, "'"
!            write(test%disp%unit,"(*(g0))") "stringPadded_ref  : '", stringPadded_ref, "'"
!            write(test%disp%unit,"(*(g0))") "charPadded      : '", charPadded, "'"
!            write(test%disp%unit,"(*(g0))")
!        end if
!        ! LCOV_EXCL_STOP
!    end function test_padVal2Str_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_val2str ! LCOV_EXCL_LINE