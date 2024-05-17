program example

    use pm_io, only: display_type
    use pm_kind, only: SK, IK, LK, CK, RK ! all other processor kinds are also supported.
    use pm_distUnif, only: setUnifRand
    use pm_arrayRank, only: getRankModified

    implicit none

    integer(IK) , parameter :: NP = 5_IK

    ! Vector of indices of the sorted array.

    integer(IK)                         :: i
    integer(IK)         , allocatable   :: rank(:)

    character(:, SK)    , allocatable   :: string_SK
    character(9, SK)    , allocatable   :: vector_SK(:)
    real(RK)                            :: vector_RK(NP)
    integer(IK)                         :: vector_IK(NP)
    logical(LK)                         :: vector_LK(NP)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Define the unsorted arrays.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !call setUnifRand(vector_SK, SK_"aaaaaaaaa", SK_"zzzzzzzzz")
    call setUnifRand(vector_RK, -0.5_RK, +0.5_RK)
    call setUnifRand(vector_IK, 1_IK, 2_IK**NP)
    call setUnifRand(vector_LK)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!rank the characters of a string in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = SK_"ParaMonte"
    call allocateRank(len(string_SK))

    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("rank = getRankModified(string_SK)")
                    rank = getRankModified(string_SK)
    call disp%show("rank")
    call disp%show( rank )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!rank array of strings of the same length in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    vector_SK = [ "ParaMonte" &
                , "V.2      " &
                , "is       " &
                , "a        " &
                , "Parallel " &
                , "Monte    " &
                , "Carlo    " &
                , "and      " &
                , "a        " &
                , "Machine  " &
                , "Learning " &
                , "Library. " &
                ]

    call allocateRank(size(vector_SK))
    call disp%skip()
    call disp%show("vector_SK")
    call disp%show( vector_SK, deliml = SK_"""" )
    call disp%show("rank = getRankModified(vector_SK)")
                    rank = getRankModified(vector_SK)
    call disp%show("rank")
    call disp%show( rank )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!rank vector of integers of arbitrary kinds in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call allocateRank(size(vector_IK))
    call disp%skip()
    call disp%show("vector_IK")
    call disp%show( vector_IK )
    call disp%show("rank = getRankModified(vector_IK)")
                    rank = getRankModified(vector_IK)
    call disp%show("rank")
    call disp%show( rank )
    call disp%skip()

    call allocateRank(size(vector_IK))
    call disp%skip()
    call disp%show("vector_IK = [1_IK, 2_IK, 3_IK, 2_IK, 1_IK]")
                    vector_IK = [1_IK, 2_IK, 3_IK, 2_IK, 1_IK]
    call disp%show("rank = getRankModified(vector_IK)")
                    rank = getRankModified(vector_IK)
    call disp%show("rank")
    call disp%show( rank )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!rank vector of logicals of arbitrary kinds in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call allocateRank(size(vector_LK))
    call disp%skip()
    call disp%show("vector_LK")
    call disp%show( vector_LK )
    call disp%show("rank = getRankModified(vector_LK)")
                    rank = getRankModified(vector_LK)
    call disp%show("rank")
    call disp%show( rank )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!rank arrays of reals of arbitrary kinds in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call allocateRank(size(vector_RK))
    call disp%skip()
    call disp%show("vector_RK")
    call disp%show( vector_RK )
    call disp%show("rank = getRankModified(vector_RK)")
                    rank = getRankModified(vector_RK)
    call disp%show("rank")
    call disp%show( rank )
    call disp%skip()

    block
        real, allocatable :: vector_RK(:)
        vector_RK = [1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 5.0, 5.0]
        call allocateRank(size(vector_RK))
        call disp%skip()
        call disp%show("vector_RK")
        call disp%show( vector_RK )
        call disp%show("rank = getRankModified(vector_RK)")
                        rank = getRankModified(vector_RK)
        call disp%show("rank")
        call disp%show( rank )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!rank according to an input user-defined comparison function.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    vector_IK = int([(i, i = 1_IK, NP)], kind = IK)
    call allocateRank(size(vector_IK))
    call disp%skip()
    call disp%show("!rank in DESCENDING (decreasing) order via an input custom-designed `isSorted()` function.")
    call disp%skip()
    call disp%show("vector_IK")
    call disp%show( vector_IK )
    call disp%show("rank = getRankModified(vector_IK, isSorted_IK)")
                    rank = getRankModified(vector_IK, isSorted_IK)
    call disp%show("rank")
    call disp%show( rank )
    call disp%skip()

    call random_number(vector_RK); vector_RK  = vector_RK  - 0.5_RK
    call allocateRank(size(vector_RK))
    call disp%skip()
    call disp%show("!rank in ascending order solely based on the magnitude of numbers using a custom comparison function.")
    call disp%skip()
    call disp%show("vector_RK")
    call disp%show( vector_RK )
    call disp%show("rank = getRankModified(vector_RK, isSorted_RK)")
                    rank = getRankModified(vector_RK, isSorted_RK)
    call disp%show("rank")
    call disp%show( rank )
    call disp%skip()

    string_SK = "ParaMonte"
    call allocateRank(len(string_SK))
    call disp%skip()
    call disp%show("!rank string in ascending order without case-sensitivity via a custom-designed input comparison function.")
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("rank = getRankModified(string_SK, isSortedChar1)")
                    rank = getRankModified(string_SK, isSortedChar1)
    call disp%show("rank")
    call disp%show( rank )
    call disp%skip()

    string_SK = "2211"
    call allocateRank(len(string_SK))
    call disp%skip()
    call disp%show("!rank string in ascending order without case-sensitivity via a custom-designed input comparison function.")
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("rank = getRankModified(string_SK, isSortedChar2)")
                    rank = getRankModified(string_SK, isSortedChar2)
    call disp%show("rank")
    call disp%show( rank )
    call disp%skip()

    ! gfortran 12 / ifort still cannot digest PDTs.
#if PDT_ENABLED
    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!rank array of strings of varying length in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_container, only: css_pdt
        type(css_pdt) , allocatable   :: cssvec(:)
        cssvec =[ css_pdt("ParaMonte") &
                , css_pdt("V.2") &
                , css_pdt("is") &
                , css_pdt("a") &
                , css_pdt("Parallel") &
                , css_pdt("Monte") &
                , css_pdt("Carlo") &
                , css_pdt("and") &
                , css_pdt("a") &
                , css_pdt("Machine") &
                , css_pdt("Learning") &
                , css_pdt("Library.") &
                ]
        call allocateRank(size(cssvec))
        call disp%skip()
        call disp%show("cssvec")
        call disp%show( cssvec, deliml = SK_"""" )
        call disp%show("rank = getRankModified(cssvec)")
                        rank = getRankModified(cssvec)
        call disp%show("rank")
        call disp%show( rank )
        call disp%skip()
    end block
#endif

contains

    function isSorted_IK(a,b) result(isSorted)
        use pm_kind, only: LK, IK
        integer(IK)    , intent(in)    :: a, b
        logical(LK)                     :: isSorted
        isSorted = a > b
    end function

    function isSorted_RK(a,b) result(isSorted)
        use pm_kind, only: LK, IK
        integer(IK)    , intent(in)    :: a, b
        logical(LK)                     :: isSorted
        isSorted = abs(a) < abs(b)
    end function

    function isSortedChar1(a,b) result(isSorted)
        use pm_strASCII, only: getStrLower
        use pm_kind, only: LK, SK
        character(1, SK), intent(in)    :: a, b
        logical(LK)                     :: isSorted
        isSorted = getStrLower(a) < getStrLower(b)
    end function

    function isSortedChar2(a,b) result(isSorted)
        use pm_kind, only: LK, SK
        character(1, SK), intent(in)    :: a, b
        logical(LK)                     :: isSorted
        isSorted = a < b
    end function

    subroutine allocateRank(size)
        integer, intent(in) :: size
        if (allocated(rank)) deallocate(rank)
        allocate(rank(size))
    end subroutine

end program example