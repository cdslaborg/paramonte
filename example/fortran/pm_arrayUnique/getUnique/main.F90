program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayUnique, only: getUnique

    implicit none

    character(:, SK), allocatable   :: string_SK    
    character(9, SK), allocatable   :: Array_SK(:) ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: Array_IK(:) ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: Array_CK(:) ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: Array_RK(:) ! Can be any processor-supported kind.

    type(display_type)              :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find all unique elements in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ParaMonte is a Monte Carlo Library."
    Array_SK = ["ParaMonte", "PARAMONTE", "paramonte", "ParaMonte", "ParaMonte", "Paramonte"]
    Array_IK = [1_IK, 1_IK, 4_IK, 4_IK]
    Array_RK = [1._RK, 2._RK, 4._RK, 4._RK]
    Array_CK = [(1._CK, -1._CK), (1._CK, -2._CK), (1._CK, -1._CK), (4._CK, -4._CK)]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("getUnique(string_SK)")
    call disp%show( getUnique(string_SK), deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("getUnique(Array_SK)")
    Array_SK = getUnique(Array_SK) ! This is a bug in gfortran as of version 11.
    call disp%show( Array_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("getUnique(Array_IK)")
    call disp%show( getUnique(Array_IK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("getUnique(Array_CK)")
    call disp%show( getUnique(Array_CK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("getUnique(Array_RK)")
    call disp%show( getUnique(Array_RK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find all unique elements according to the user-specified equivalence criterion.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_RK = [1.01_RK, 1.04_RK, 0.98_RK, 1.0_RK, 1.02_RK, 2._RK]

    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("getUnique(Array_RK, iseq = iseq_RK)")
    call disp%show( getUnique(Array_RK, iseq = iseq_RK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique case-insensitive instances within the character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ABBAbbA"

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("getUnique(string_SK, iseq = iseq_SK)")
    call disp%show( getUnique(string_SK, iseq = iseq_SK), deliml = SK_"""" )

contains

    pure function iseq_RK(element1, element2) result(iseq)
        real(RK)    , intent(in)    :: element1, element2
        logical(LK)                 :: iseq
        iseq = abs(abs(element1) - abs(element2)) < 0.05_RK
    end function

    pure function iseq_SK(element1, element2) result(iseq)
        use pm_strASCII, only: getStrLower
        character(*, SK)    , intent(in)    :: element1, element2
        logical(LK)                     :: iseq
        iseq = getStrLower(element1) == getStrLower(element2)
    end function

end program example