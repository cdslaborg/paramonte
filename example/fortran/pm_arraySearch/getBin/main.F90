program example

    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arraySearch, only: getBin

    implicit none

    character(:, SK), allocatable   :: string_SK, strval_SK 
    character(2, SK), allocatable   :: Array_SK(:)          ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: Array_IK(:)          ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: Array_RK(:)          ! Can be any processor-supported kind.

    character(2,SK) :: value_SK  ! Can be any processor-supported kind.
    integer(IK)     :: value_IK  ! Can be any processor-supported kind.
    real(RK)        :: value_RK  ! Can be any processor-supported kind.

    integer(IK)     :: bin ! Must be of default kind IK

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the index of the largest element in the ascending-ordered `array` that is smaller than or equal to `value` via Binary Search.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "abcdefghi"
    Array_SK = ["aa", "bb", "dd", "ee", "ff", "gg", "hh", "ii", "jj"]
    Array_RK = [0._RK, 1._RK, 2._RK, 4._RK, 5._RK, 6._RK]
    Array_IK = [0_IK, 1_IK, 2_IK, 4_IK, 5_IK, 6_IK]

    strval_SK = "c"
    value_SK = "cc"
    value_IK = 3_IK
    value_RK = 3.5_RK

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the index of character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("strval_SK")
    call disp%show( strval_SK, deliml = SK_"""" )
    call disp%show("bin = getBin(string_SK, strval_SK) ! index of the largest value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(string_SK, strval_SK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    strval_SK = "cz"
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("strval_SK")
    call disp%show( strval_SK, deliml = SK_"""" )
    call disp%show("bin = getBin(string_SK, strval_SK) ! index of the largest value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(string_SK, strval_SK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    strval_SK = "ca"
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("strval_SK")
    call disp%show( strval_SK, deliml = SK_"""" )
    call disp%show("bin = getBin(string_SK, strval_SK) ! index of the largest value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(string_SK, strval_SK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the index of character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("value_SK")
    call disp%show( value_SK, deliml = SK_"""" )
    call disp%show("bin = getBin(Array_SK, value_SK) ! index of the largest value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(Array_SK, value_SK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the index of integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("value_IK")
    call disp%show( value_IK )
    call disp%show("bin = getBin(Array_IK, value_IK) ! index of the largest value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(Array_IK, value_IK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the index of real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("value_RK")
    call disp%show( value_RK )
    call disp%show("bin = getBin(Array_RK, value_RK) ! index of the largest value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(Array_RK, value_RK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! The index of a `value` that is out of range.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    value_IK = Array_IK(1) - 1_IK
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("value_IK")
    call disp%show( value_IK )
    call disp%show("bin = getBin(Array_IK, value_IK) ! index of the largest value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(Array_IK, value_IK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    value_IK = Array_IK(size(Array_IK)) + 1_IK
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("value_IK")
    call disp%show( value_IK )
    call disp%show("bin = getBin(Array_IK, value_IK) ! index of the largest value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(Array_IK, value_IK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the index of the smallest element in the descending-ordered `array` that is")
    call disp%show("! larger than or equal to `value` via Binary Search with a custom comparison function.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the index of a `value` in a descending-ordered input `Array`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    value_IK = 2_IK
    Array_IK = Array_IK(size(Array_IK):1:-1)
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("value_IK")
    call disp%show( value_IK )
    call disp%show("bin = getBin(Array_IK, value_IK, isLess_IK) ! index of the smallest value in `array` that is larger than or equal to the input `value`.")
                    bin = getBin(Array_IK, value_IK, isLess_IK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find the index of a `value` in an input `array` whose absolute value is ascending-ordered.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_RK = [( Array_RK(bin)*(-1)**bin, bin = 1, size(Array_RK, kind = IK) )]
    call disp%skip()
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("value_RK")
    call disp%show( value_RK )
    call disp%show("bin = getBin(Array_RK, value_RK, isLess_RK) ! index of the largest absolute value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(Array_RK, value_RK, isLess_RK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Be mindful of duplicate values in the input `Array`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    value_IK = 0_IK
    Array_IK = int([0, 0, 3, 3, 3], kind = IK)
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("value_IK")
    call disp%show( value_IK )
    call disp%show("bin = getBin(Array_IK, value_IK) ! index of the largest value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(Array_IK, value_IK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

    value_IK = 3_IK
    Array_IK = int([0, 0, 3, 3, 3], kind = IK)
    call disp%skip()
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("value_IK")
    call disp%show( value_IK )
    call disp%show("bin = getBin(Array_IK, value_IK) ! index of the largest value in `array` that is smaller than or equal to the input `value`.")
                    bin = getBin(Array_IK, value_IK)
    call disp%show("bin")
    call disp%show( bin )
    call disp%skip()

contains

    !>  \brief
    !> Custom-comparison function to find the index of a `value` in an array sorted in descending-order.
    pure function isLess_IK(value, segment) result(less)
        use pm_kind, only: LK
        integer(IK) , intent(in)    :: value, segment
        logical(LK)                 :: less
        less = value > segment
    end function

    !>  \brief
    !> Custom-comparison function to find the index of a `value` in an array whose absolute value is sorted in ascending-order.
    pure function isLess_RK(value, segment) result(less)
        use pm_kind, only: LK
        real(RK)    , intent(in)    :: value, segment
        logical(LK)                 :: less
        less = value < abs(segment)
    end function

end program example