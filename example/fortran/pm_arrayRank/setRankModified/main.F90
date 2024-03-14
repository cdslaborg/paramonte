program example

    use pm_io, only: display_type
    use pm_kind, only: SK, IK, LK, CK, RK ! all other processor kinds are also supported.
    use pm_distUnif, only: setUnifRand
    use pm_arrayRank, only: setRankModified
    use pm_container, only: css_pdt

    implicit none

    integer(IK) , parameter :: NP = 5_IK

    ! Vector of indices of the sorted array.

    integer(IK)                         :: i
    integer(IK)         , allocatable   :: Rank(:)
    
    character(:, SK)    , allocatable   :: string_SK
    character(9, SK)    , allocatable   :: Vector_SK(:)
    real(RK)                            :: Vector_RK(NP)
    integer(IK)                         :: Vector_IK(NP)
    logical(LK)                         :: Vector_LK(NP)
    type(css_pdt) , allocatable   :: cssvec(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Define the unsorted arrays.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !call setUnifRand(Vector_SK, SK_"aaaaaaaaa", SK_"zzzzzzzzz")
    call setUnifRand(Vector_RK, -0.5_RK, +0.5_RK)
    call setUnifRand(Vector_IK, 1_IK, 2_IK**NP)
    call setUnifRand(Vector_LK)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Rank the characters of a string in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = SK_"ParaMonte"
    call allocateRank(len(string_SK))

    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("call setRankModified(Rank, string_SK)")
                    call setRankModified(Rank, string_SK)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Rank array of strings of the same length in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Vector_SK = [ "ParaMonte" &
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

    call allocateRank(size(Vector_SK))
    call disp%skip()
    call disp%show("Vector_SK")
    call disp%show( Vector_SK, deliml = SK_"""" )
    call disp%show("call setRankModified(Rank, Vector_SK)")
                    call setRankModified(Rank, Vector_SK)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Rank vector of integers of arbitrary kinds in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call allocateRank(size(Vector_IK))
    call disp%skip()
    call disp%show("Vector_IK")
    call disp%show( Vector_IK )
    call disp%show("call setRankModified(Rank, Vector_IK)")
                    call setRankModified(Rank, Vector_IK)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()

    call allocateRank(size(Vector_IK))
    call disp%skip()
    call disp%show("Vector_IK = [1_IK, 2_IK, 3_IK, 2_IK, 1_IK]")
                    Vector_IK = [1_IK, 2_IK, 3_IK, 2_IK, 1_IK]
    call disp%show("call setRankModified(Rank, Vector_IK)")
                    call setRankModified(Rank, Vector_IK)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Rank vector of logicals of arbitrary kinds in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call allocateRank(size(Vector_LK))
    call disp%skip()
    call disp%show("Vector_LK")
    call disp%show( Vector_LK )
    call disp%show("call setRankModified(Rank, Vector_LK)")
                    call setRankModified(Rank, Vector_LK)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Rank arrays of reals of arbitrary kinds in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call allocateRank(size(Vector_RK))
    call disp%skip()
    call disp%show("Vector_RK")
    call disp%show( Vector_RK )
    call disp%show("call setRankModified(Rank, Vector_RK)")
                    call setRankModified(Rank, Vector_RK)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()

    block
        real, allocatable :: Vector_RK(:)
        Vector_RK = [1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 5.0, 5.0]
        call allocateRank(size(Vector_RK))
        call disp%skip()
        call disp%show("Vector_RK")
        call disp%show( Vector_RK )
        call disp%show("call setRankModified(Rank, Vector_RK)")
                        call setRankModified(Rank, Vector_RK)
        call disp%show("Rank")
        call disp%show( Rank )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Rank according to an input user-defined comparison function.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Vector_IK = int([(i, i = 1_IK, NP)], kind = IK)
    call allocateRank(size(Vector_IK))
    call disp%skip()
    call disp%show("!Rank in DESCENDING (decreasing) order via an input custom-designed `isSorted()` function.")
    call disp%skip()
    call disp%show("Vector_IK")
    call disp%show( Vector_IK )
    call disp%show("call setRankModified(Rank, Vector_IK, isSorted_IK)")
                    call setRankModified(Rank, Vector_IK, isSorted_IK)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()

    call random_number(Vector_RK); Vector_RK  = Vector_RK  - 0.5_RK
    call allocateRank(size(Vector_RK))
    call disp%skip()
    call disp%show("!Rank in ascending order solely based on the magnitude of numbers using a custom comparison function.")
    call disp%skip()
    call disp%show("Vector_RK")
    call disp%show( Vector_RK )
    call disp%show("call setRankModified(Rank, Vector_RK, isSorted_RK)")
                    call setRankModified(Rank, Vector_RK, isSorted_RK)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()

    string_SK = "ParaMonte"
    call allocateRank(len(string_SK))
    call disp%skip()
    call disp%show("!Rank string in ascending order without case-sensitivity via a custom-designed input comparison function.")
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("call setRankModified(Rank, string_SK, isSortedChar1)")
                    call setRankModified(Rank, string_SK, isSortedChar1)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()

    string_SK = "2211"
    call allocateRank(len(string_SK))
    call disp%skip()
    call disp%show("!Rank string in ascending order without case-sensitivity via a custom-designed input comparison function.")
    call disp%skip()
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("call setRankModified(Rank, string_SK, isSortedChar2)")
                    call setRankModified(Rank, string_SK, isSortedChar2)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()

    ! gfortran 12 still cannot digest PDTs.
#if !__GFORTRAN__
    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Rank array of strings of varying length in ascending order.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    cssvec = [ css_pdt("ParaMonte") &
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
    call disp%show("call setRankModified(Rank, cssvec)")
                    call setRankModified(Rank, cssvec)
    call disp%show("Rank")
    call disp%show( Rank )
    call disp%skip()
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
        if (allocated(Rank)) deallocate(Rank)
        allocate(Rank(size))
    end subroutine

end program example