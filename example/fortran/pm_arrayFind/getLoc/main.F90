program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayFind, only: getLoc

    implicit none

    integer(IK)     , allocatable   :: instance(:)  , Loc(:)                ! Must be of default kind IK

    character(:, SK), allocatable   :: string_SK    , stringPattern_SK      
    character(9, SK), allocatable   :: array_SK(:)  , arrayPattern_SK(:)    ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: array_IK(:)  , arrayPattern_IK(:)    ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: array_CK(:)  , arrayPattern_CK(:)    ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: array_RK(:)  , arrayPattern_RK(:)    ! Can be any processor-supported kind.
    logical(LK)     , allocatable   :: array_LK(:)  , arrayPattern_LK(:)

    type(display_type)              :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find all instances of pattern in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ParaMonte is a Machine Learning Library "
    array_SK = ["ParaMonte", "XXXXXXXXX", "is       ", "XXXXXXXXX", "a        ", "XXXXXXXXX", "Monte    ", "XXXXXXXXX", "Carlo    ", "XXXXXXXXX", "Library. ", "XXXXXXXXX"]
    array_IK = [1_IK, 0_IK, 2_IK, 0_IK, 3_IK, 0_IK, 4_IK]
    array_RK = [1._RK, 0._RK, 2._RK, 0._RK, 3._RK, 0._RK, 4._RK]
    array_CK = [(1._CK, -1._CK), (0._CK, -0._CK), (2._CK, -2._CK), (0._CK, -0._CK), (3._CK, -3._CK), (0._CK, -0._CK), (4._CK, -4._CK)]
    array_LK = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK, .true._LK, .false._LK]

    stringPattern_SK = " "
    arrayPattern_SK = ["XXXXXXXXX"]
    arrayPattern_IK = [0_IK]
    arrayPattern_RK = [0._RK]
    arrayPattern_CK = [(0._CK, -0._CK)]
    arrayPattern_LK = [.true._LK]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("Loc = getLoc(string_SK, stringPattern_SK)")
                    Loc = getLoc(string_SK, stringPattern_SK)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("arrayPattern_SK")
    call disp%show( arrayPattern_SK, deliml = SK_"""" )
    call disp%show("Loc = getLoc(array_SK, arrayPattern_SK)")
                    Loc = getLoc(array_SK, arrayPattern_SK)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_LK")
    call disp%show( array_LK )
    call disp%show("arrayPattern_LK")
    call disp%show( arrayPattern_LK )
    call disp%show("Loc = getLoc(array_LK, arrayPattern_LK)")
                    Loc = getLoc(array_LK, arrayPattern_LK)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_IK")
    call disp%show( array_IK )
    call disp%show("arrayPattern_IK")
    call disp%show( arrayPattern_IK )
    call disp%show("Loc = getLoc(array_IK, arrayPattern_IK)")
                    Loc = getLoc(array_IK, arrayPattern_IK)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_CK")
    call disp%show( array_CK )
    call disp%show("arrayPattern_CK")
    call disp%show( arrayPattern_CK )
    call disp%show("Loc = getLoc(array_CK, arrayPattern_CK)")
                    Loc = getLoc(array_CK, arrayPattern_CK)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_RK")
    call disp%show( array_RK )
    call disp%show("arrayPattern_RK")
    call disp%show( arrayPattern_RK )
    call disp%show("Loc = getLoc(array_RK, arrayPattern_RK)")
                    Loc = getLoc(array_RK, arrayPattern_RK)
    call disp%show("Loc")
    call disp%show( Loc )



    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find only particular instances of pattern in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "A_A_A_A_A_A_A_A_A"
    array_SK = ["A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A"]
    array_IK = [0_IK, 1_IK, 0_IK, 2_IK, 3_IK, 0_IK, 4_IK, 5_IK, 0_IK, 0_IK]
    array_RK = [0._RK, 1._RK, 0._RK, 2._RK, 3._RK, 0._RK, 4._RK, 5._RK, 0._RK, 0._RK]
    array_CK = [(0._CK, -0._CK), (1._CK, -1._CK), (0._CK, -0._CK), (2._CK, -2._CK), (3._CK, -3._CK), (0._CK, -0._CK), (4._CK, -4._CK), (5._CK, -5._CK), (0._CK, -0._CK), (0._CK, -0._CK)]
    array_LK = [.false._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .false._LK]

    stringPattern_SK = "_"
    arrayPattern_SK = ["_"]
    arrayPattern_IK = [0_IK]
    arrayPattern_RK = [0._RK]
    arrayPattern_CK = [(0._CK, -0._CK)]
    arrayPattern_LK = [.false._LK]

    instance = [-3, 2, -4] ! remove at the second occurrence from the beginning and the third occurrence from the end. Duplicate indices are ignored.

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(string_SK, stringPattern_SK, instance = instance)")
                    Loc = getLoc(string_SK, stringPattern_SK, instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find vector `pattern` from character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("arrayPattern_SK")
    call disp%show( arrayPattern_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_SK, arrayPattern_SK, instance = instance)")
                    Loc = getLoc(array_SK, arrayPattern_SK, instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find character array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("arrayPattern_SK(1)")
    call disp%show( arrayPattern_SK(1), deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_SK, arrayPattern_SK(1), instance = instance)")
                    Loc = getLoc(array_SK, arrayPattern_SK(1), instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find logical array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_LK")
    call disp%show( array_LK )
    call disp%show("arrayPattern_LK")
    call disp%show( arrayPattern_LK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_LK, arrayPattern_LK, instance = instance)")
                    Loc = getLoc(array_LK, arrayPattern_LK, instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find logical array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_LK")
    call disp%show( array_LK )
    call disp%show("arrayPattern_LK(1)")
    call disp%show( arrayPattern_LK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_LK, arrayPattern_LK(1), instance = instance)")
                    Loc = getLoc(array_LK, arrayPattern_LK(1), instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find integer array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_IK")
    call disp%show( array_IK )
    call disp%show("arrayPattern_IK")
    call disp%show( arrayPattern_IK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_IK, arrayPattern_IK, instance = instance)")
                    Loc = getLoc(array_IK, arrayPattern_IK, instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find integer array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_IK")
    call disp%show( array_IK )
    call disp%show("arrayPattern_IK(1)")
    call disp%show( arrayPattern_IK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_IK, arrayPattern_IK(1), instance = instance)")
                    Loc = getLoc(array_IK, arrayPattern_IK(1), instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find complex array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_CK")
    call disp%show( array_CK )
    call disp%show("arrayPattern_CK")
    call disp%show( arrayPattern_CK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_CK, arrayPattern_CK, instance = instance)")
                    Loc = getLoc(array_CK, arrayPattern_CK, instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find complex array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_CK")
    call disp%show( array_CK )
    call disp%show("arrayPattern_CK(1)")
    call disp%show( arrayPattern_CK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_CK, arrayPattern_CK(1), instance = instance)")
                    Loc = getLoc(array_CK, arrayPattern_CK(1), instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find real array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_RK")
    call disp%show( array_RK )
    call disp%show("arrayPattern_RK")
    call disp%show( arrayPattern_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_RK, arrayPattern_RK, instance = instance)")
                    Loc = getLoc(array_RK, arrayPattern_RK, instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find real array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_RK")
    call disp%show( array_RK )
    call disp%show("arrayPattern_RK(1)")
    call disp%show( arrayPattern_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_RK, arrayPattern_RK(1), instance = instance)")
                    Loc = getLoc(array_RK, arrayPattern_RK(1), instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Adjust the blindness to exclusively find non-overlapping instances of `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    string_SK = "AAAAAAAA"
    stringPattern_SK = "AAA"

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("Loc = getLoc(string_SK, stringPattern_SK, blindness = 1_IK) ! The default blindness episode after each detection is `blindness = 1_IK`")
                    Loc = getLoc(string_SK, stringPattern_SK, blindness = 1_IK)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("Loc = getLoc(string_SK, stringPattern_SK, blindness = size(stringPattern_SK, kind = IK)) ! Find only non-overlapping patterns.")
                        Loc = getLoc(string_SK, stringPattern_SK, blindness = len(stringPattern_SK, kind = IK))
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("Loc = getLoc(string_SK, stringPattern_SK, blindness = 2_IK) ! Find instances with jumps of size 2.")
                    Loc = getLoc(string_SK, stringPattern_SK, blindness = 2_IK)
    call disp%show("Loc")
    call disp%show( Loc )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find specific instances with a user-defined equivalence test.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find case-insensitive instances of vector `pattern` within the character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ABBAbbA"
    stringPattern_SK = "bb"

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("Loc = getLoc(string_SK, stringPattern_SK, iseq = iseq_SK)")
                    Loc = getLoc(string_SK, stringPattern_SK, iseq = iseq_SK)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find specific instances of vector `pattern` within the real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    array_RK = [0._RK, 1.01_RK, 1.04_RK, 0.98_RK, 1.0_RK, 1.02_RK]
    arrayPattern_RK = [-1._RK, 1._RK]
    instance = [-2_IK]

    call disp%show("array_RK")
    call disp%show( array_RK )
    call disp%show("arrayPattern_RK")
    call disp%show( arrayPattern_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_RK, arrayPattern_RK, iseq = iseq_vec_RK, instance = instance)")
                    Loc = getLoc(array_RK, arrayPattern_RK, iseq = iseq_vec_RK, instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find specific instances of scalar `pattern` within the real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("array_RK")
    call disp%show( array_RK )
    call disp%show("arrayPattern_RK(1)")
    call disp%show( arrayPattern_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("Loc = getLoc(array_RK, arrayPattern_RK(1), iseq = iseq_RK, instance = instance)")
                    Loc = getLoc(array_RK, arrayPattern_RK(1), iseq = iseq_RK, instance = instance)
    call disp%show("Loc")
    call disp%show( Loc )

contains

    pure function iseq_SK(ArraySegment, pattern) result(equivalent)
        use pm_strASCII, only: getStrLower
        character(*, SK), intent(in)    :: pattern, ArraySegment
        logical(LK)                     :: equivalent
        equivalent = pattern == getStrLower(ArraySegment)
    end function

    function iseq_RK(arraysegment, pattern) result(equivalent)
        real(RK)        , intent(in)    :: pattern, arraySegment
        logical(LK)                     :: equivalent
        equivalent = abs(abs(pattern) - abs(arraySegment)) < 0.05_RK
    end function

    function iseq_vec_RK(ArraySegment, pattern, lenPattern) result(equivalent)
        integer(IK)     , intent(in)    :: lenPattern
        real(RK)        , intent(in)    :: pattern(lenPattern), ArraySegment(lenPattern)
        logical(LK)                     :: equivalent
        equivalent = all(abs(abs(pattern) - abs(ArraySegment)) < 0.05_RK)
    end function

end program example