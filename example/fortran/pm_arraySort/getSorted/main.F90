program example

    use pm_kind, only: SK, IK, LK
    use pm_arraySort, only: isort
    use pm_arraySort, only: getSorted
    use pm_distUnif, only: getUnifRand
    use pm_arrayRemap, only: getRemapped
    use pm_arrayResize, only: setResized
    use pm_io, only: display_type

    implicit none

    integer(IK), allocatable :: index(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Sort a string of characters of arbitrary kind in arbitrary orders.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(:), allocatable :: array
        call disp%skip()
        call disp%show("!Sort array in ascending order with case-sensitivity.")
        call disp%show("array = 'ParaMonte'")
                        array = 'ParaMonte'
        call disp%show("index = getSorted(array, isort)")
                        index = getSorted(array, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("getRemapped(array, index)")
        call disp%show( getRemapped(array, index) , deliml = SK_"""" )
        call disp%show("array = getSorted(array)")
                        array = getSorted(array)
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%skip()

        call disp%skip()
        call disp%show("!Sort array in ascending order without case-sensitivity via a custom-designed input comparison function.")
        array = "ParaMonte"
        call disp%skip()
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%show("index = getSorted(array, isSortedChar, isort)")
                        index = getSorted(array, isSortedChar, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("getRemapped(array, index)")
        call disp%show( getRemapped(array, index) , deliml = SK_"""" )
        call disp%show("array = getSorted(array, isSortedChar)")
                        array = getSorted(array, isSortedChar)
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Sort array of strings of the same length of arbitrary kind in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        character(10), allocatable :: array(:)
        array = [ "ParaMonte " &
                , "V.2       " &
                , "is        " &
                , "a         " &
                , "Parallel  " &
                , "Monte     " &
                , "Carlo     " &
                , "and       " &
                , "Machine   " &
                , "Learning  " &
                , "Library.  " &
                ]
        call disp%skip()
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%show("index = getSorted(array, isort)")
                        index = getSorted(array, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("array(index)")
        call disp%show( array(index) , deliml = SK_"""" )
        call disp%show("array = getSorted(array)")
                        array = getSorted(array)
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Sort array of integer values of arbitrary kind in arbitrary orders.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        integer, allocatable :: array(:)
        call disp%skip()
        call disp%show("!Sort integer values in ascending order.")
        call disp%show("array = getUnifRand(-9, 9, s1 = getUnifRand(5_IK, 10_IK))")
                        array = getUnifRand(-9, 9, s1 = getUnifRand(5_IK, 10_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("index = getSorted(array, isort)")
                        index = getSorted(array, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("array(index)")
        call disp%show( array(index) )
        call disp%show("array = getSorted(array)")
                        array = getSorted(array)
        call disp%show("array")
        call disp%show( array )
        call disp%skip()

        call disp%skip()
        call disp%show("!Sort integer values in descending order via a custom-designed input comparison function.")
        call disp%skip()
        call disp%show("array")
        call disp%show( array )
        call disp%show("index = getSorted(array, isSortedInteger, isort)")
                        index = getSorted(array, isSortedInteger, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("array(index)")
        call disp%show( array(index) )
        call disp%show("array = getSorted(array, isSortedInteger)")
                        array = getSorted(array, isSortedInteger)
        call disp%show("array")
        call disp%show( array )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Sort array of logical values of arbitrary kind in arbitrary orders.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        logical, allocatable :: array(:)
        call disp%skip()
        call disp%show("!Sort logical values in ascending order.")
        call disp%show("array = getUnifRand(.false., .true., s1 = getUnifRand(5_IK, 10_IK))")
                        array = getUnifRand(.false., .true., s1 = getUnifRand(5_IK, 10_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("index = getSorted(array, isort)")
                        index = getSorted(array, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("array(index)")
        call disp%show( array(index) )
        call disp%show("array = getSorted(array)")
                        array = getSorted(array)
        call disp%show("array")
        call disp%show( array )
        call disp%skip()

        call disp%skip()
        call disp%show("!Sort logical values in descending order via a custom-designed input comparison function.")
        call disp%skip()
        call disp%show("array")
        call disp%show( array )
        call disp%show("index = getSorted(array, isSortedLogical, isort)")
                        index = getSorted(array, isSortedLogical, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("array(index)")
        call disp%show( array(index) )
        call disp%show("array = getSorted(array, isSortedLogical)")
                        array = getSorted(array, isSortedLogical)
        call disp%show("array")
        call disp%show( array )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Sort array of complex values of arbitrary kind in arbitrary orders.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        complex, allocatable :: array(:)
        call disp%skip()
        call disp%show("!Sort complex values in ascending order.")
        call disp%show("array = getUnifRand((-9., -9), (9., 9.), s1 = getUnifRand(3_IK, 5_IK))")
                        array = getUnifRand((-9., -9), (9., 9.), s1 = getUnifRand(3_IK, 5_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("index = getSorted(array, isort)")
                        index = getSorted(array, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("array(index)")
        call disp%show( array(index) )
        call disp%show("array = getSorted(array)")
                        array = getSorted(array)
        call disp%show("array")
        call disp%show( array )
        call disp%skip()

        call disp%skip()
        call disp%show("!Sort complex values in ascending order of their modulus via a custom-designed input comparison function.")
        call disp%skip()
        call disp%show("array")
        call disp%show( array )
        call disp%show("index = getSorted(array, isSortedComplex, isort)")
                        index = getSorted(array, isSortedComplex, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("array(index)")
        call disp%show( array(index) )
        call disp%show("array = getSorted(array, isSortedComplex, isort)")
                        array = getSorted(array, isSortedComplex, isort)
        call disp%show("array")
        call disp%show( array )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Sort array of real values of arbitrary kind in arbitrary orders.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real, allocatable :: array(:)
        call disp%skip()
        call disp%show("!Sort real values in ascending order.")
        call disp%show("array = getUnifRand(-9., 9., s1 = getUnifRand(4_IK, 8_IK))")
                        array = getUnifRand(-9., 9., s1 = getUnifRand(4_IK, 8_IK))
        call disp%show("array")
        call disp%show( array )
        call disp%show("index = getSorted(array, isort)")
                        index = getSorted(array, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("array(index)")
        call disp%show( array(index) )
        call disp%show("array = getSorted(array)")
                        array = getSorted(array)
        call disp%show("array")
        call disp%show( array )
        call disp%skip()

        call disp%skip()
        call disp%show("!Sort real values in ascending order of their modulus via a custom-designed input comparison function.")
        call disp%skip()
        call disp%show("array")
        call disp%show( array )
        call disp%show("index = getSorted(array, isSortedReal, isort)")
                        index = getSorted(array, isSortedReal, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("array(index)")
        call disp%show( array(index) )
        call disp%show("array = getSorted(array, isSortedReal, isort)")
                        array = getSorted(array, isSortedReal, isort)
        call disp%show("array")
        call disp%show( array )
        call disp%skip()
    end block

#if !__GFORTRAN__ && 0
    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Sort array of string containers in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_container, only: css_pdt ! \bug Intel ifort 2021.6 cannot handle PDT aliases such as `strc`.
        type(css_pdt), allocatable :: array(:)
        array = [ css_pdt("ParaMonte") &
                , css_pdt("V.2") &
                , css_pdt("is") &
                , css_pdt("a") &
                , css_pdt("Parallel") &
                , css_pdt("Monte") &
                , css_pdt("Carlo") &
                , css_pdt("and") &
                , css_pdt("Machine") &
                , css_pdt("Learning") &
                , css_pdt("Library.") &
                ]
        call disp%skip()
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%show("index = getSorted(array, isort)")
                        index = getSorted(array, isort)
        call disp%show("index")
        call disp%show( index )
        call disp%show("getRemapped(array, index)")
        call disp%show( getRemapped(array, index) , deliml = SK_"""" )
        call disp%show("array = getSorted(array)")
                        array = getSorted(array)
        call disp%show("array")
        call disp%show( array , deliml = SK_"""" )
        call disp%skip()
    end block
#endif

contains

    function isSortedInteger(a,b) result(sorted)
        use pm_kind, only: LK
        integer         , intent(in)    :: a, b
        logical(LK)                     :: sorted
        sorted = a > b
    end function

    function isSortedLogical(a,b) result(sorted)
        use pm_logicalCompare, only: operator(>)
        use pm_kind, only: LK
        logical         , intent(in)    :: a, b
        logical(LK)                     :: sorted
        sorted = a > b
    end function

    function isSortedComplex(a,b) result(sorted)
        use pm_kind, only: LK
        complex         , intent(in)    :: a, b
        logical(LK)                     :: sorted
        sorted = abs(a) < abs(b)
    end function

    function isSortedReal(a,b) result(sorted)
        use pm_kind, only: LK
        real            , intent(in)    :: a, b
        logical(LK)                     :: sorted
        sorted = abs(a) < abs(b)
    end function

    function isSortedChar(a,b) result(sorted)
        use pm_strASCII, only: getStrLower
        use pm_kind, only: LK
        character(1)    , intent(in)    :: a, b
        logical(LK)                     :: sorted
        sorted = getStrLower(a) < getStrLower(b)
    end function

end program example