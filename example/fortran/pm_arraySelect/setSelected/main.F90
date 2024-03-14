program example

    use pm_io, only: display_type
    use pm_kind, only: IK, RK32, RK64, RK128
    use pm_kind, only: IK8, IK16, IK32, IK64
    use pm_kind, only: LK, SK
    use pm_arraySelect, only: setSelected
    use pm_container, only: strc => css_pdt

    implicit none

    integer(IK) , parameter :: NP = 6_IK

    ! 1-dimensional array of real values.

    real(RK32 )                     :: array_RK32 (NP), selection_RK32
    real(RK64 )                     :: array_RK64 (NP), selection_RK64
    real(RK128)                     :: array_RK128(NP), selection_RK128

    ! 1-dimensional array of integer values.

    integer(IK8 )                   :: array_IK8  (NP), selection_IK8
    integer(IK16)                   :: array_IK16 (NP), selection_IK16
    integer(IK32)                   :: array_IK32 (NP), selection_IK32
    integer(IK64)                   :: array_IK64 (NP), selection_IK64

    ! Character array

    character(:, SK), allocatable   :: string_SK
    character(1, SK)                :: selectionString_SK

    character(10,SK), allocatable   :: array_SK(:)
    character(10,SK)                :: selectionarray_SK

    type(strc)      , allocatable   :: array_PSSK(:)
    type(strc)                      :: selection_PSSK

    integer(IK)                     :: i, rank

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Define the unselected arrays.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call random_number(array_RK32 ); array_RK32  = array_RK32  - 0.5_RK32
    call random_number(array_RK64 ); array_RK64  = array_RK64  - 0.5_RK64
    call random_number(array_RK128); array_RK128 = array_RK128 - 0.5_RK128

    array_IK8  = int(array_RK128 * huge(0_IK8 ), kind = IK8 )
    array_IK16 = int(array_RK128 * huge(0_IK16), kind = IK16)
    array_IK32 = int(array_RK128 * huge(0_IK32), kind = IK32)
    array_IK64 = int(array_RK128 * huge(0_IK64), kind = IK64)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Select arrays in ascending order.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Select the `rank`th smallest element in the input `integer` Array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    rank = 3_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IK8")
    call disp%show( array_IK8 )
    call disp%show("call setSelected(selection_IK8, array_IK8, rank)")
                    call setSelected(selection_IK8, array_IK8, rank)
    call disp%show("selection_IK8")
    call disp%show( selection_IK8 )
    call disp%show("array_IK8")
    call disp%show( array_IK8 )
    call disp%skip()

    rank = 4_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IK16")
    call disp%show( array_IK16 )
    call disp%show("call setSelected(selection_IK16, array_IK16, rank)")
                    call setSelected(selection_IK16, array_IK16, rank)
    call disp%show("selection_IK16")
    call disp%show( selection_IK16 )
    call disp%show("array_IK16")
    call disp%show( array_IK16 )
    call disp%skip()

    rank = 1_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IK32")
    call disp%show( array_IK32 )
    call disp%show("call setSelected(selection_IK32, array_IK32, rank)")
                    call setSelected(selection_IK32, array_IK32, rank)
    call disp%show("selection_IK32")
    call disp%show( selection_IK32 )
    call disp%show("array_IK32")
    call disp%show( array_IK32 )
    call disp%skip()

    rank = 5_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IK64")
    call disp%show( array_IK64 )
    call disp%show("call setSelected(selection_IK64, array_IK64, rank)")
                    call setSelected(selection_IK64, array_IK64, rank)
    call disp%show("selection_IK64")
    call disp%show( selection_IK64 )
    call disp%show("array_IK64")
    call disp%show( array_IK64 )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Select the `rank`th smallest element in the input `real` Array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    rank = 3_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_RK32")
    call disp%show( array_RK32 )
    call disp%show("call setSelected(selection_RK32, array_RK32, rank)")
                    call setSelected(selection_RK32, array_RK32, rank)
    call disp%show("selection_RK32")
    call disp%show( selection_RK32 )
    call disp%show("array_RK32")
    call disp%show( array_RK32 )
    call disp%skip()

    rank = 2_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_RK64")
    call disp%show( array_RK64 )
    call disp%show("call setSelected(selection_RK64, array_RK64, rank)")
                    call setSelected(selection_RK64, array_RK64, rank)
    call disp%show("selection_RK64")
    call disp%show( selection_RK64 )
    call disp%show("array_RK64")
    call disp%show( array_RK64 )
    call disp%skip()

    rank = size(array_RK128, kind = IK)
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_RK128")
    call disp%show( array_RK128 )
    call disp%show("call setSelected(selection_RK128, array_RK128, rank)")
                    call setSelected(selection_RK128, array_RK128, rank)
    call disp%show("selection_RK128")
    call disp%show( selection_RK128 )
    call disp%show("array_RK128")
    call disp%show( array_RK128 )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Select the `rank`th smallest element in the input Fortran `string` (character) array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    array_SK =  [ "ParaMonte " &
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

    rank = 3_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%show("call setSelected(Selectionarray_SK, array_SK, rank)")
                    call setSelected(Selectionarray_SK, array_SK, rank)
    call disp%show("Selectionarray_SK")
    call disp%show( Selectionarray_SK, deliml = SK_"""" )
    call disp%show("array_SK")
    call disp%show( array_SK, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Select the `rank`th smallest element in the input Fortran `string` (character) scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    rank = 3_IK
    string_SK = "ParaMonte"
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("call setSelected(selectionString_SK, string_SK, rank)")
                    call setSelected(selectionString_SK, string_SK, rank)
    call disp%show("selectionString_SK")
    call disp%show( selectionString_SK, deliml = SK_"""" )
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Select the `rank`th smallest element in the input array of containers of Fortran `string` (character) scalars.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

#if !__GFORTRAN__ && 0
    array_PSSK = [ strc("ParaMonte") &
                , strc("V.2") &
                , strc("is") &
                , strc("a") &
                , strc("Parallel") &
                , strc("Monte") &
                , strc("Carlo") &
                , strc("and") &
                , strc("Machine") &
                , strc("Learning") &
                , strc("Library.") &
                ]

    rank = 3_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_PSSK")
    call disp%show( array_PSSK, deliml = SK_"""" )
    call disp%show("call setSelected(selection_PSSK, array_PSSK, rank)")
                    call setSelected(selection_PSSK, array_PSSK, rank)
    call disp%show("selection_PSSK")
    call disp%show( selection_PSSK, deliml = SK_"""" )
    call disp%show("array_PSSK")
    call disp%show( array_PSSK, deliml = SK_"""" )
    call disp%skip()
#endif

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Select according to an input user-defined comparison function.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    array_IK8 = int([((-2)**i, i = 1, NP)], kind = IK8)
    rank = 3_IK
    call disp%skip()
    call disp%show("!Select the `rank`th largest element via an input custom-designed `isSorted()` function.")
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IK8")
    call disp%show( array_IK8 )
    call disp%show("call setSelected(selection_IK8, array_IK8, rank, isSorted_IK8)")
                    call setSelected(selection_IK8, array_IK8, rank, isSorted_IK8)
    call disp%show("selection_IK8")
    call disp%show( selection_IK8 )
    call disp%show("array_IK8")
    call disp%show( array_IK8 )
    call disp%skip()

    call random_number(array_RK32); array_RK32  = array_RK32  - 0.5_RK32
    call disp%skip()
    call disp%show("!Select the `rank`th smallest element solely based on the magnitude of numbers using a custom comparison function.")
    call disp%skip()
    rank = 3_IK
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_RK32")
    call disp%show( array_RK32 )
    call disp%show("call setSelected(selection_RK32, array_RK32, rank, isSorted_RK32)")
                    call setSelected(selection_RK32, array_RK32, rank, isSorted_RK32)
    call disp%show("selection_RK32")
    call disp%show( selection_RK32 )
    call disp%show("array_RK32")
    call disp%show( array_RK32 )
    call disp%skip()

    string_SK = "ParaMonte"
    call disp%skip()
    call disp%show("!Select the `rank`th smallest element with case-sensitivity (default behavior).")
    call disp%skip()
    rank = 3_IK
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("call setSelected(selectionString_SK, string_SK, rank)")
                    call setSelected(selectionString_SK, string_SK, rank)
    call disp%show("selectionString_SK")
    call disp%show( selectionString_SK, deliml = SK_"""" )
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%skip()

    string_SK = "ParaMonte"
    call disp%skip()
    call disp%show("!Select the `rank`th smallest element WITHOUT case-sensitivity via a custom-designed input comparison function.")
    call disp%skip()
    rank = 3_IK
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("call setSelected(selectionString_SK, string_SK, rank, isSorted_SK)")
                    call setSelected(selectionString_SK, string_SK, rank, isSorted_SK)
    call disp%show("selectionString_SK")
    call disp%show( selectionString_SK, deliml = SK_"""" )
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Expedite the selection process for a partially sorted array via optional arguments `lb` or `ub`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    array_IK8 = [-3_IK8, -1_IK8, -2_IK8, 1_IK8, 2_IK8, 3_IK8] ! all elements after index `3` are sorted (`ub = 3_IK`).
    rank = 3_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IK8")
    call disp%show( array_IK8 )
    call disp%show("call setSelected(selection_IK8, array_IK8, rank, ub = 3_IK) ! all elements after index `3` are sorted (`ub = 3_IK`).")
                    call setSelected(selection_IK8, array_IK8, rank, ub = 3_IK)
    call disp%show("selection_IK8")
    call disp%show( selection_IK8 )
    call disp%show("array_IK8")
    call disp%show( array_IK8 )
    call disp%skip()

    array_IK8 = [-3_IK8, -1_IK8, 3_IK8, -2_IK8, 2_IK8, 1_IK8]
    call disp%skip()
    call disp%show("!Select the second, fourth, and fifth smallest elements of the array sequentially in order.")
    call disp%skip()
    block
        integer(IK) :: i, RankList(0:3) = [0_IK, 2_IK, 4_IK, 5_IK]
        call disp%show("RankList")
        call disp%show( RankList )
        do i = 1, size(RankList(1:))
            call disp%skip()
            call disp%show("RankList(i)")
            call disp%show( RankList(i) )
            call disp%show("array_IK8")
            call disp%show( array_IK8 )
            call disp%show("call setSelected(selection_IK8, array_IK8, rank = RankList(i), lb = RankList(i-1)+1) ! all elements after `lb` need sorting.")
                            call setSelected(selection_IK8, array_IK8, rank = RankList(i), lb = RankList(i-1)+1)
            call disp%show("selection_IK8")
            call disp%show( selection_IK8 )
            call disp%show("array_IK8")
            call disp%show( array_IK8 )
        end do
    end block
    call disp%skip()

contains

    function isSorted_IK8(a,b) result(isSorted)
        use pm_kind, only: LK, IK8
        integer(IK8)    , intent(in)    :: a, b
        logical(LK)                     :: isSorted
        isSorted = a > b
    end function

    function isSorted_RK32(a,b) result(isSorted)
        use pm_kind, only: LK, IK8
        integer(IK8)    , intent(in)    :: a, b
        logical(LK)                     :: isSorted
        isSorted = abs(a) < abs(b)
    end function

    function isSorted_SK(a,b) result(isSorted)
        use pm_strASCII, only: getStrLower
        use pm_kind, only: LK, SK
        character(1, SK), intent(in)    :: a, b
        logical(LK)                     :: isSorted
        isSorted = getStrLower(a) < getStrLower(b)
    end function

end program example