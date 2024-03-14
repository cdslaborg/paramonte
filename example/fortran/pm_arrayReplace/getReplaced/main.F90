program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayReplace, only: getReplaced

    implicit none

    integer(IK)     , allocatable   :: instance(:) ! Must be of default kind IK
    character(:, SK), allocatable   :: string_SK    , stringPattern_SK  , stringReplacement_SK      
    character(9, SK), allocatable   :: Array_SK(:)  , ArrayPattern_SK(:), ArrayReplacement_SK(:)    ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: Array_IK(:)  , ArrayPattern_IK(:), ArrayReplacement_IK(:)    ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: Array_CK(:)  , ArrayPattern_CK(:), ArrayReplacement_CK(:)    ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: Array_RK(:)  , ArrayPattern_RK(:), ArrayReplacement_RK(:)    ! Can be any processor-supported kind.
    logical(LK)     , allocatable   :: Array_LK(:)  , ArrayPattern_LK(:), ArrayReplacement_LK(:)

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace all instances of pattern in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ParaMonte is a Monte Carlo Library."
    Array_SK = ["ParaMonte", "is       ", "a        ", "Monte    ", "Carlo    ", "Library. "]
    Array_IK = [1_IK, 2_IK, 3_IK, 4_IK]
    Array_RK = [1._RK, 2._RK, 3._RK, 4._RK]
    Array_CK = [(1._CK, -1._CK), (2._CK, -2._CK), (3._CK, -3._CK), (4._CK, -4._CK)]
    Array_LK = [.false._LK, .true._LK, .true._LK, .false._LK]

    stringPattern_SK = "a Monte Carlo"
    ArrayPattern_SK = ["a        ", "Monte    ", "Carlo    "]
    ArrayPattern_IK = [3_IK]
    ArrayPattern_RK = [3._RK]
    ArrayPattern_CK = [(3._CK, -3._CK)]
    ArrayPattern_LK = [.true._LK, .true._LK]

    stringReplacement_SK = "now a Machine Learning"
    ArrayReplacement_SK = ["now      ", "a        ", "Machine  ", "Learning "]
    ArrayReplacement_IK = [-1_IK, 0_IK, 1_IK]
    ArrayReplacement_RK = [-1._RK, 0._RK, 1._RK]
    ArrayReplacement_CK = [(-1._CK, +1._CK), (0._CK, 0._CK), (+1._CK, -1._CK)]
    ArrayReplacement_LK = [.false._LK]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("stringReplacement_SK")
    call disp%show( stringReplacement_SK, deliml = SK_"""" )
    call disp%show("getReplaced(string_SK, stringPattern_SK, stringReplacement_SK)")
    call disp%show( getReplaced(string_SK, stringPattern_SK, stringReplacement_SK), deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("ArrayPattern_SK")
    call disp%show( ArrayPattern_SK, deliml = SK_"""" )
    call disp%show("ArrayReplacement_SK")
    call disp%show( ArrayReplacement_SK, deliml = SK_"""" )
    call disp%show("getReplaced(Array_SK, ArrayPattern_SK, ArrayReplacement_SK)")
    call disp%show("getReplaced(Array_SK, ArrayPattern_SK, ArrayReplacement_SK)") ! gfortran 11 bug: does not allow allocatable deferred length array because the interface is assumed-length allocatable.

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_LK")
    call disp%show( Array_LK )
    call disp%show("ArrayPattern_LK")
    call disp%show( ArrayPattern_LK )
    call disp%show("ArrayReplacement_LK")
    call disp%show( ArrayReplacement_LK )
    call disp%show("getReplaced(Array_LK, ArrayPattern_LK, ArrayReplacement_LK)")
    call disp%show( getReplaced(Array_LK, ArrayPattern_LK, ArrayReplacement_LK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("ArrayPattern_IK")
    call disp%show( ArrayPattern_IK )
    call disp%show("ArrayReplacement_IK")
    call disp%show( ArrayReplacement_IK )
    call disp%show("getReplaced(Array_IK, ArrayPattern_IK, ArrayReplacement_IK)")
    call disp%show( getReplaced(Array_IK, ArrayPattern_IK, ArrayReplacement_IK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("ArrayPattern_CK")
    call disp%show( ArrayPattern_CK )
    call disp%show("ArrayReplacement_CK")
    call disp%show( ArrayReplacement_CK )
    call disp%show("getReplaced(Array_CK, ArrayPattern_CK, ArrayReplacement_CK)")
    call disp%show( getReplaced(Array_CK, ArrayPattern_CK, ArrayReplacement_CK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("ArrayPattern_RK")
    call disp%show( ArrayPattern_RK )
    call disp%show("ArrayReplacement_RK")
    call disp%show( ArrayReplacement_RK )
    call disp%show("getReplaced(Array_RK, ArrayPattern_RK, ArrayReplacement_RK)")
    call disp%show( getReplaced(Array_RK, ArrayPattern_RK, ArrayReplacement_RK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace only particular instances of pattern in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "1ParaMonte-2ParaMonte-3ParaMonte-4ParaMonte-5ParaMonte"
   !Array_SK = ["ParaMonte", "ParaMonte", "ParaMonte", "ParaMonte", "ParaMonte"]
    Array_IK = [1_IK, 1_IK, 1_IK, 1_IK, 1_IK]
    Array_RK = [1._RK, 1._RK, 1._RK, 1._RK, 1._RK]
    Array_CK = [(+1._CK, -1._CK), (+1._CK, -1._CK), (+1._CK, -1._CK), (+1._CK, -1._CK), (+1._CK, -1._CK)]

    stringPattern_SK = "ParaMonte"
   !ArrayPattern_SK = "ParaMonte"
    ArrayPattern_IK = [1_IK]
    ArrayPattern_RK = [1._RK]
    ArrayPattern_CK = [(+1._CK, -1._CK)]

    stringReplacement_SK = ""
   !ArrayReplacement_SK = ["XXXXXXXXX"]
    ArrayReplacement_IK = [-2_IK,-3_IK]
    ArrayReplacement_RK = [-2._RK,-3._RK]
    ArrayReplacement_CK = [(-2._CK, +2._CK), (-3._CK, +3._CK)]

    ! Replace only the second occurrence from the beginning as well as the first and the second occurrences from the end.

    instance = [-1_IK, -2_IK, 2_IK]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace specific instances within the character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("stringReplacement_SK")
    call disp%show( stringReplacement_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("getReplaced(string_SK, stringPattern_SK, stringReplacement_SK, instance)")
    call disp%show( getReplaced(string_SK, stringPattern_SK, stringReplacement_SK, instance), deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace specific instances within the character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("ArrayPattern_SK")
    call disp%show( ArrayPattern_SK, deliml = SK_"""" )
    call disp%show("ArrayReplacement_SK")
    call disp%show( ArrayReplacement_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("getReplaced(Array_SK, ArrayPattern_SK, ArrayReplacement_SK, instance)")
    Array_SK = getReplaced(Array_SK, ArrayPattern_SK, ArrayReplacement_SK, instance) ! gfortran 11 bug: does not allow passing function result directly to `disp%show()`.
    call disp%show( Array_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace specific instances within the integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("ArrayPattern_IK")
    call disp%show( ArrayPattern_IK )
    call disp%show("ArrayReplacement_IK")
    call disp%show( ArrayReplacement_IK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("getReplaced(Array_IK, ArrayPattern_IK, ArrayReplacement_IK, instance)")
    call disp%show( getReplaced(Array_IK, ArrayPattern_IK, ArrayReplacement_IK, instance) )

    call disp%skip()
    call disp%show("Array_IK = [integer(IK) :: 5, 8, 4, 7, 5, 5, 5]")
                    Array_IK = [integer(IK) :: 5, 8, 4, 7, 5, 5, 5]
    call disp%show("ArrayPattern_IK = [5_IK]")
                    ArrayPattern_IK = [5_IK]
    call disp%show("ArrayReplacement_IK")
                    ArrayReplacement_IK = [integer(IK) ::]
    call disp%show("instance = [integer(IK) :: -2, 3, 2, -13, -11, 6, 4, -7, 8]")
                    instance = [integer(IK) :: -2, 3, 2, -13, -11, 6, 4, -7, 8]
    call disp%show("getReplaced(Array_IK, ArrayPattern_IK(1), ArrayReplacement_IK, instance, sorted = .false._LK, unique = .false._LK)")
    call disp%show( getReplaced(Array_IK, ArrayPattern_IK(1), ArrayReplacement_IK, instance, sorted = .false._LK, unique = .false._LK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace specific instances within the logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_LK = [logical(LK) :: .false., .true.,  .true.,  .false., .false., .true.,  .false., .false., .false., .true.,  .false., .false., .false., .true.,  .false.]")
                    Array_LK = [logical(LK) :: .false., .true.,  .true.,  .false., .false., .true.,  .false., .false., .false., .true.,  .false., .false., .false., .true.,  .false.]
    call disp%show("ArrayPattern_LK = [logical(LK) :: .true.]")
                    ArrayPattern_LK = [logical(LK) :: .true.]
    call disp%show("ArrayReplacement_LK")
                    ArrayReplacement_LK = [logical(LK) ::]
    call disp%show("instance = [integer(IK) :: 3]")
                    instance = [integer(IK) :: 3]
    call disp%show("getReplaced(Array_LK, ArrayPattern_LK, ArrayReplacement_LK, instance, sorted = .false._LK, unique = .false._LK)")
    call disp%show( getReplaced(Array_LK, ArrayPattern_LK, ArrayReplacement_LK, instance, sorted = .false._LK, unique = .false._LK) )

    call disp%skip()
    call disp%show("Array_LK = [logical(LK) :: .false.]")
                    Array_LK = [logical(LK) :: .false.]
    call disp%show("ArrayPattern_LK = [logical(LK) :: .false.]")
                    ArrayPattern_LK = [logical(LK) :: .false.]
    call disp%show("ArrayReplacement_LK = [logical(LK) ::]")
                    ArrayReplacement_LK = [logical(LK) ::]
    call disp%show("instance = [integer(IK) :: 0]")
                    instance = [integer(IK) :: 0]
    call disp%show("getReplaced(Array_LK, ArrayPattern_LK(1), ArrayReplacement_LK, instance, sorted = .false._LK, unique = .false._LK)")
    call disp%show( getReplaced(Array_LK, ArrayPattern_LK(1), ArrayReplacement_LK, instance, sorted = .false._LK, unique = .false._LK) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace specific instances within the complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("ArrayPattern_CK")
    call disp%show( ArrayPattern_CK )
    call disp%show("ArrayReplacement_CK")
    call disp%show( ArrayReplacement_CK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("getReplaced(Array_CK, ArrayPattern_CK, ArrayReplacement_CK, instance)")
    call disp%show( getReplaced(Array_CK, ArrayPattern_CK, ArrayReplacement_CK, instance) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace specific instances within the real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("ArrayPattern_RK")
    call disp%show( ArrayPattern_RK )
    call disp%show("ArrayReplacement_RK")
    call disp%show( ArrayReplacement_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("getReplaced(Array_RK, ArrayPattern_RK, ArrayReplacement_RK, instance)")
    call disp%show( getReplaced(Array_RK, ArrayPattern_RK, ArrayReplacement_RK, instance) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace specific instances with a user-defined equality test.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace specific instances within the real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_RK = [1.01_RK, 1.04_RK, 0.98_RK, 1.0_RK, 1.02_RK]
    ArrayPattern_RK = [-1._RK, 1._RK]
    ArrayReplacement_RK = [0._RK,0._RK]
    instance = [-2]

    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("ArrayPattern_RK")
    call disp%show( ArrayPattern_RK )
    call disp%show("ArrayReplacement_RK")
    call disp%show( ArrayReplacement_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("getReplaced(Array_RK, ArrayPattern_RK, ArrayReplacement_RK, iseq = iseq_RK, instance = instance)")
    call disp%show( getReplaced(Array_RK, ArrayPattern_RK, ArrayReplacement_RK, iseq = iseq_RK, instance = instance) )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Replace case-insensitive instances within the character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ABBAbbA"
    stringPattern_SK = "bb"
    stringReplacement_SK = "x"

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("stringReplacement_SK")
    call disp%show( stringReplacement_SK, deliml = SK_"""" )
    call disp%show("getReplaced(string_SK, stringPattern_SK, stringReplacement_SK, iseq = iseq_SK)")
    call disp%show( getReplaced(string_SK, stringPattern_SK, stringReplacement_SK, iseq = iseq_SK), deliml = SK_"""" )

contains

    function iseq_RK(Segment, pattern, lenPattren) result(equivalent)
        integer(IK)     , intent(in)    :: lenPattren
        real(RK)        , intent(in)    :: Segment(lenPattren)
        real(RK)        , intent(in)    :: pattern(lenPattren)
        logical(LK)                     :: equivalent
        equivalent = all(abs(abs(pattern) - abs(Segment)) < 0.05_RK)
    end function

    function iseq_SK(segment, pattern) result(equivalent)
        use pm_strASCII, only: getStrLower
        character(*, SK), intent(in)    :: segment, pattern
        logical(LK)                     :: equivalent
        equivalent = getStrLower(segment) == pattern
    end function

end program example