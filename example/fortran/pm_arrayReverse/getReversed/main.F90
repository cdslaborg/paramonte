program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayReverse, only: getReversed

    implicit none

    character(:, SK), allocatable   :: string1_SK   , string2_SK 
    character(9, SK), allocatable   :: Array1_SK(:) , Array2_SK(:) ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: Array1_IK(:) , Array2_IK(:) ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: Array1_CK(:) , Array2_CK(:) ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: Array1_RK(:) , Array2_RK(:) ! Can be any processor-supported kind.
    logical(LK)     , allocatable   :: Array1_LK(:) , Array2_LK(:)

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Reverse an array via a function call.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string1_SK = "ABCDEF"
    Array1_SK = ["AA", "BB", "CC", "DD", "EE", "FF"]
    Array1_IK = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK]
    Array1_RK = [1._RK, 2._RK, 3._RK, 4._RK, 5._RK, 6._RK]
    Array1_CK = [(1._CK, -1._CK), (2._CK, -2._CK), (3._CK, -3._CK), (4._CK, -4._CK), (5._CK, -5._CK), (6._CK, -6._CK)]
    Array1_LK = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK, .true._LK]

    string2_SK = string1_SK
    Array2_SK = Array1_SK
    Array2_IK = Array1_IK
    Array2_RK = Array1_RK
    Array2_CK = Array1_CK
    Array2_LK = Array1_LK

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Reverse character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string2_SK")
    call disp%show( string2_SK, deliml = SK_"""" )
    call disp%show("getReversed(string2_SK)")
    call disp%show( getReversed(string2_SK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Reverse character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_SK")
    call disp%show( Array2_SK, deliml = SK_"""" )
    call disp%show("getReversed(Array2_SK)")
    Array2_SK = getReversed(Array2_SK) ! bypass gfortran 10.3 bug.
    call disp%show( Array2_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Reverse logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_LK")
    call disp%show( Array2_LK )
    call disp%show("getReversed(Array2_LK)")
    call disp%show( getReversed(Array2_LK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Reverse integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_IK")
    call disp%show( Array2_IK )
    call disp%show("getReversed(Array2_IK)")
    call disp%show( getReversed(Array2_IK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Reverse complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_CK")
    call disp%show( Array2_CK )
    call disp%show("getReversed(Array2_CK)")
    call disp%show( getReversed(Array2_CK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Reverse real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_RK")
    call disp%show( Array2_RK )
    call disp%show("getReversed(Array2_RK)")
    call disp%show( getReversed(Array2_RK) )

end program example