program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayPad, only: getPadded

    implicit none

    character(:, SK), allocatable   :: string_SK    
    character(2, SK), allocatable   :: Array_SK(:)  ! Can be any processor-supported kind.
    logical(LK)     , allocatable   :: Array_LK(:)  ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: Array_IK(:)  ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: Array_CK(:)  ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: Array_RK(:)  ! Can be any processor-supported kind.

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Pad an array to the left and right with symbols.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()

    allocate(Array_SK(3:8), source = [SK_"AA", SK_"BB", SK_"CC", SK_"DD", SK_"EE", SK_"FF"])
    allocate(Array_IK(3:8), source = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK])
    allocate(Array_RK(3:8), source = [1._RK, 2._RK, 3._RK, 4._RK, 5._RK, 6._RK])
    allocate(Array_CK(3:8), source = [(1._CK, -1._CK), (2._CK, -2._CK), (3._CK, -3._CK), (4._CK, -4._CK), (5._CK, -5._CK), (6._CK, -6._CK)])
    allocate(Array_LK(3:8), source = [.true._LK, .true._LK, .true._LK, .true._LK, .true._LK, .true._LK])

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Pad a character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = SK_"In God We Trust"
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("len(string_SK)")
    call disp%show( len(string_SK) )
    call disp%show("string_SK = getPadded(string_SK, lpsize = 5_IK, rpsize = 5_IK, lpfill = SK_'-', rpfill = SK_'+')")
                    string_SK = getPadded(string_SK, lpsize = 5_IK, rpsize = 5_IK, lpfill = SK_'-', rpfill = SK_'+')
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("len(string_SK)")
    call disp%show( len(string_SK) )
    call disp%skip()

    string_SK = SK_"In Science We Invest"
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("len(string_SK)")
    call disp%show( len(string_SK) )
    call disp%show("string_SK = getPadded(string_SK, lpsize = 4_IK, rpsize = 4_IK, lpfill = SK_' ', rpfill = SK_' ', lmsize = 4_IK, rmsize = 8_IK, lmfill = SK_'*', rmfill = SK_'*')")
                    string_SK = getPadded(string_SK, lpsize = 4_IK, rpsize = 4_IK, lpfill = SK_' ', rpfill = SK_' ', lmsize = 4_IK, rmsize = 8_IK, lmfill = SK_'*', rmfill = SK_'*')
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("len(string_SK)")
    call disp%show( len(string_SK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Pad a character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("lbound(Array_SK), ubound(Array_SK)")
    call disp%show([lbound(Array_SK), ubound(Array_SK)])
    call disp%show("Array_SK = getPadded(Array_SK, lpsize = 5_IK, rpsize = 5_IK, lpfill = SK_'  ', rpfill = SK_'++')")
                    Array_SK = getPadded(Array_SK, lpsize = 5_IK, rpsize = 5_IK, lpfill = SK_'  ', rpfill = SK_'++')
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("lbound(Array_SK), ubound(Array_SK)")
    call disp%show([lbound(Array_SK), ubound(Array_SK)])

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Pad a logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_LK")
    call disp%show( Array_LK )
    call disp%show("lbound(Array_LK), ubound(Array_LK)")
    call disp%show([lbound(Array_LK), ubound(Array_LK)])
    call disp%show("Array_LK = getPadded(Array_LK, lpsize = 5_IK, rpsize = 5_IK, lpfill = .false._LK, rpfill = .false._LK)")
                    Array_LK = getPadded(Array_LK, lpsize = 5_IK, rpsize = 5_IK, lpfill = .false._LK, rpfill = .false._LK)
    call disp%show("Array_LK")
    call disp%show( Array_LK )
    call disp%show("lbound(Array_LK), ubound(Array_LK)")
    call disp%show([lbound(Array_LK), ubound(Array_LK)])

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Pad an integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("lbound(Array_IK), ubound(Array_IK)")
    call disp%show([lbound(Array_IK), ubound(Array_IK)])
    call disp%show("Array_IK = getPadded(Array_IK, lpsize = 5_IK, rpsize = 5_IK, lpfill = -0_IK, rpfill = +0_IK)")
                    Array_IK = getPadded(Array_IK, lpsize = 5_IK, rpsize = 5_IK, lpfill = -0_IK, rpfill = +0_IK)
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("lbound(Array_IK), ubound(Array_IK)")
    call disp%show([lbound(Array_IK), ubound(Array_IK)])

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Pad a complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("lbound(Array_CK), ubound(Array_CK)")
    call disp%show([lbound(Array_CK), ubound(Array_CK)])
    call disp%show("Array_CK = getPadded(Array_CK, lpsize = 5_IK, rpsize = 5_IK, lpfill = (-999._CK,-999._CK), rpfill = (+999._CK,+999._CK))")
                    Array_CK = getPadded(Array_CK, lpsize = 5_IK, rpsize = 5_IK, lpfill = (-999._CK,-999._CK), rpfill = (+999._CK,+999._CK))
    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("lbound(Array_CK), ubound(Array_CK)")
    call disp%show([lbound(Array_CK), ubound(Array_CK)])

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%show("! Pad a real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("lbound(Array_RK), ubound(Array_RK)")
    call disp%show([lbound(Array_RK), ubound(Array_RK)])
    call disp%show("Array_RK = getPadded(Array_RK, lpsize = 5_IK, rpsize = 5_IK, lpfill = -0._RK, rpfill = +0._RK)")
                    Array_RK = getPadded(Array_RK, lpsize = 5_IK, rpsize = 5_IK, lpfill = -0._RK, rpfill = +0._RK)
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("lbound(Array_RK), ubound(Array_RK)")
    call disp%show([lbound(Array_RK), ubound(Array_RK)])

contains

    subroutine reset()
        if (allocated(Array_SK)) deallocate(Array_SK)
        if (allocated(Array_IK)) deallocate(Array_IK)
        if (allocated(Array_RK)) deallocate(Array_RK)
        if (allocated(Array_CK)) deallocate(Array_CK)
        if (allocated(Array_LK)) deallocate(Array_LK)
    end subroutine

end program example