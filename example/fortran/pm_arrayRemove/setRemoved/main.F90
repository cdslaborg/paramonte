program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayRemove, only: setRemoved

    implicit none

    integer(IK) , allocatable   :: instance(:) ! Must be of default kind IK

    character(:, SK), allocatable   :: stringRef_SK     , string_SK     , stringPattern_SK      
    character(9, SK), allocatable   :: ArrayRef_SK(:)   , Array_SK(:)   , ArrayPattern_SK(:)    ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: ArrayRef_IK(:)   , Array_IK(:)   , ArrayPattern_IK(:)    ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: ArrayRef_CK(:)   , Array_CK(:)   , ArrayPattern_CK(:)    ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: ArrayRef_RK(:)   , Array_RK(:)   , ArrayPattern_RK(:)    ! Can be any processor-supported kind.
    logical(LK)     , allocatable   :: ArrayRef_LK(:)   , Array_LK(:)   , ArrayPattern_LK(:)

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    stringRef_SK = "ParaMonte is a Machine Learning Library "
    ArrayRef_SK = ["ParaMonte", "XXXXXXXXX", "is       ", "XXXXXXXXX", "a        ", "XXXXXXXXX", "Monte    ", "XXXXXXXXX", "Carlo    ", "XXXXXXXXX", "Library. ", "XXXXXXXXX"]
    ArrayRef_IK = [1_IK, 0_IK, 2_IK, 0_IK, 3_IK, 0_IK, 4_IK]
    ArrayRef_RK = [1._RK, 0._RK, 2._RK, 0._RK, 3._RK, 0._RK, 4._RK]
    ArrayRef_CK = [(1._CK, -1._CK), (0._CK, -0._CK), (2._CK, -2._CK), (0._CK, -0._CK), (3._CK, -3._CK), (0._CK, -0._CK), (4._CK, -4._CK)]
    ArrayRef_LK = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK, .true._LK, .false._LK]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove all instances of pattern in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = stringRef_SK
    Array_SK = ArrayRef_SK
    Array_IK = ArrayRef_IK
    Array_RK = ArrayRef_RK
    Array_CK = ArrayRef_CK
    Array_LK = ArrayRef_LK

    stringPattern_SK = " "
    ArrayPattern_SK = ["XXXXXXXXX"]
    ArrayPattern_IK = [0_IK]
    ArrayPattern_RK = [0._RK]
    ArrayPattern_CK = [(0._CK, -0._CK)]
    ArrayPattern_LK = [.true._LK]

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("call setRemoved(string_SK, stringPattern_SK)")
                    call setRemoved(string_SK, stringPattern_SK)
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("ArrayPattern_SK")
    call disp%show( ArrayPattern_SK, deliml = SK_"""" )
    call disp%show("call setRemoved(Array_SK, ArrayPattern_SK)")
                    call setRemoved(Array_SK, ArrayPattern_SK)
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_LK")
    call disp%show( Array_LK )
    call disp%show("ArrayPattern_LK")
    call disp%show( ArrayPattern_LK )
    call disp%show("call setRemoved(Array_LK, ArrayPattern_LK)")
                    call setRemoved(Array_LK, ArrayPattern_LK)
    call disp%show("Array_LK")
    call disp%show( Array_LK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("ArrayPattern_IK")
    call disp%show( ArrayPattern_IK )
    call disp%show("call setRemoved(Array_IK, ArrayPattern_IK)")
                    call setRemoved(Array_IK, ArrayPattern_IK)
    call disp%show("Array_IK")
    call disp%show( Array_IK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("ArrayPattern_CK")
    call disp%show( ArrayPattern_CK )
    call disp%show("call setRemoved(Array_CK, ArrayPattern_CK)")
                    call setRemoved(Array_CK, ArrayPattern_CK)
    call disp%show("Array_CK")
    call disp%show( Array_CK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("ArrayPattern_RK")
    call disp%show( ArrayPattern_RK )
    call disp%show("call setRemoved(Array_RK, ArrayPattern_RK)")
                    call setRemoved(Array_RK, ArrayPattern_RK)
    call disp%show("Array_RK")
    call disp%show( Array_RK )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove only particular instances of pattern in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    stringRef_SK = "A_A_A_A_A_A_A_A_A"
    ArrayRef_SK = ["A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A", "_", "A"]
    ArrayRef_IK = [0_IK, 1_IK, 0_IK, 2_IK, 3_IK, 0_IK, 4_IK, 5_IK, 0_IK, 0_IK]
    ArrayRef_RK = [0._RK, 1._RK, 0._RK, 2._RK, 3._RK, 0._RK, 4._RK, 5._RK, 0._RK, 0._RK]
    ArrayRef_CK = [(0._CK, -0._CK), (1._CK, -1._CK), (0._CK, -0._CK), (2._CK, -2._CK), (3._CK, -3._CK), (0._CK, -0._CK), (4._CK, -4._CK), (5._CK, -5._CK), (0._CK, -0._CK), (0._CK, -0._CK)]
    ArrayRef_LK = [.false._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .true._LK, .true._LK, .false._LK, .false._LK]

    stringPattern_SK = "_"
    ArrayPattern_SK = ["_"]
    ArrayPattern_IK = [0_IK]
    ArrayPattern_RK = [0._RK]
    ArrayPattern_CK = [(0._CK, -0._CK)]
    ArrayPattern_LK = [.false._LK]

    instance = [-3, 2, -4] ! remove at the second occurrence from the beginning and the third occurrence from the end. Duplicate indices are ignored.

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = stringRef_SK
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(string_SK, stringPattern_SK, instance = instance)")
                    call setRemoved(string_SK, stringPattern_SK, instance = instance)
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove vector `pattern` from character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_SK = ArrayRef_SK
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("ArrayPattern_SK")
    call disp%show( ArrayPattern_SK, deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_SK, ArrayPattern_SK, instance = instance)")
                    call setRemoved(Array_SK, ArrayPattern_SK, instance = instance)
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove character array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_SK = ArrayRef_SK
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("ArrayPattern_SK(1)")
    call disp%show( ArrayPattern_SK(1), deliml = SK_"""" )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_SK, ArrayPattern_SK(1), instance = instance)")
                    call setRemoved(Array_SK, ArrayPattern_SK(1), instance = instance)
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove logical array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_LK = ArrayRef_LK
    call disp%show("Array_LK")
    call disp%show( Array_LK )
    call disp%show("ArrayPattern_LK")
    call disp%show( ArrayPattern_LK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_LK, ArrayPattern_LK, instance = instance)")
                    call setRemoved(Array_LK, ArrayPattern_LK, instance = instance)
    call disp%show("Array_LK")
    call disp%show( Array_LK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove logical array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_LK = ArrayRef_LK
    call disp%show("Array_LK")
    call disp%show( Array_LK )
    call disp%show("ArrayPattern_LK(1)")
    call disp%show( ArrayPattern_LK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_LK, ArrayPattern_LK(1), instance = instance)")
                    call setRemoved(Array_LK, ArrayPattern_LK(1), instance = instance)
    call disp%show("Array_LK")
    call disp%show( Array_LK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove integer array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_IK = ArrayRef_IK
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("ArrayPattern_IK")
    call disp%show( ArrayPattern_IK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_IK, ArrayPattern_IK, instance = instance)")
                    call setRemoved(Array_IK, ArrayPattern_IK, instance = instance)
    call disp%show("Array_IK")
    call disp%show( Array_IK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove integer array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_IK = ArrayRef_IK
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("ArrayPattern_IK(1)")
    call disp%show( ArrayPattern_IK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_IK, ArrayPattern_IK(1), instance = instance)")
                    call setRemoved(Array_IK, ArrayPattern_IK(1), instance = instance)
    call disp%show("Array_IK")
    call disp%show( Array_IK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove complex array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_CK = ArrayRef_CK
    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("ArrayPattern_CK")
    call disp%show( ArrayPattern_CK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_CK, ArrayPattern_CK, instance = instance)")
                    call setRemoved(Array_CK, ArrayPattern_CK, instance = instance)
    call disp%show("Array_CK")
    call disp%show( Array_CK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove complex array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_CK = ArrayRef_CK
    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("ArrayPattern_CK(1)")
    call disp%show( ArrayPattern_CK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_CK, ArrayPattern_CK(1), instance = instance)")
                    call setRemoved(Array_CK, ArrayPattern_CK(1), instance = instance)
    call disp%show("Array_CK")
    call disp%show( Array_CK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove real array with vector `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_RK = ArrayRef_RK
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("ArrayPattern_RK")
    call disp%show( ArrayPattern_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_RK, ArrayPattern_RK, instance = instance)")
                    call setRemoved(Array_RK, ArrayPattern_RK, instance = instance)
    call disp%show("Array_RK")
    call disp%show( Array_RK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove real array with scalar `pattern`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_RK = ArrayRef_RK
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("ArrayPattern_RK(1)")
    call disp%show( ArrayPattern_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_RK, ArrayPattern_RK(1), instance = instance)")
                    call setRemoved(Array_RK, ArrayPattern_RK(1), instance = instance)
    call disp%show("Array_RK")
    call disp%show( Array_RK )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove specific instances with a user-defined equivalence test.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove case-insensitive instances of vector `pattern` within the character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ABBAbbA"
    stringPattern_SK = "bb"

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringPattern_SK")
    call disp%show( stringPattern_SK, deliml = SK_"""" )
    call disp%show("call setRemoved(string_SK, stringPattern_SK, iseq = iseq_SK) ! case-insensitive string removal.")
                    call setRemoved(string_SK, stringPattern_SK, iseq = iseq_SK)
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove specific instances of vector `pattern` within the real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    ArrayRef_RK = [0._RK, 1.01_RK, 1.04_RK, 0.98_RK, 1.0_RK, 1.02_RK]
    ArrayPattern_RK = [-1._RK, 1._RK]
    instance = [-2_IK]

    Array_RK = ArrayRef_RK
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("ArrayPattern_RK")
    call disp%show( ArrayPattern_RK )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_RK, ArrayPattern_RK, iseq = iseq_vec_RK, instance = instance)")
                    call setRemoved(Array_RK, ArrayPattern_RK, iseq = iseq_vec_RK, instance = instance)
    call disp%show("Array_RK")
    call disp%show( Array_RK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remove specific instances of scalar `pattern` within the real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_RK = ArrayRef_RK
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("ArrayPattern_RK(1)")
    call disp%show( ArrayPattern_RK(1) )
    call disp%show("instance")
    call disp%show( instance )
    call disp%show("call setRemoved(Array_RK, ArrayPattern_RK(1), iseq = iseq_RK, instance = instance)")
                    call setRemoved(Array_RK, ArrayPattern_RK(1), iseq = iseq_RK, instance = instance)
    call disp%show("Array_RK")
    call disp%show( Array_RK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Call to setRemoved() preserves the lower bound of the input `Array`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    deallocate(Array_IK)
    ArrayPattern_IK = [0_IK]
    allocate(Array_IK(-7:-1), source = [1_IK, 0_IK, 2_IK, 0_IK, 3_IK, 0_IK, 4_IK])
    call disp%show("[lbound(Array_IK,1), ubound(Array_IK,1)]")
    call disp%show( [lbound(Array_IK,1), ubound(Array_IK,1)] )
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("ArrayPattern_IK(1)")
    call disp%show( ArrayPattern_IK(1) )
    call disp%show("call setRemoved(Array_IK, ArrayPattern_IK(1))")
                    call setRemoved(Array_IK, ArrayPattern_IK(1))
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("[lbound(Array_IK,1), ubound(Array_IK,1)]")
    call disp%show( [lbound(Array_IK,1), ubound(Array_IK,1)] )

contains

    pure function iseq_SK(ArraySegment, pattern) result(equivalent)
        use pm_strASCII, only: getStrLower
        character(*, SK), intent(in)    :: pattern, ArraySegment
        logical(LK)                     :: equivalent
        equivalent = pattern == getStrLower(ArraySegment)
    end function

    pure function iseq_RK(arraysegment, pattern) result(equivalent)
        real(RK)    , intent(in)    :: pattern, arraySegment
        logical(LK)                 :: equivalent
        equivalent = abs(abs(pattern) - abs(arraySegment)) < 0.05_RK
    end function

    pure function iseq_vec_RK(ArraySegment, pattern, lenPattern) result(equivalent)
        integer(IK) , intent(in)    :: lenPattern
        real(RK)    , intent(in)    :: pattern(lenPattern), ArraySegment(lenPattern)
        logical(LK)                 :: equivalent
        equivalent = all(abs(abs(pattern) - abs(ArraySegment)) < 0.05_RK)
    end function

end program example