program example

    use pm_io, only: display_type
    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKL, RKD, RKH
    use pm_kind, only: IKL, IKS, IKD, IKH
    use pm_arraySelect, only: setSelected
    use pm_container, only: strc => css_pdt

    implicit none

    integer(IK) , parameter :: NP = 6_IK

    ! 1-dimensional array of real values.

    real(RKL) :: array_RKL (NP), selection_RKL
    real(RKD) :: array_RKD (NP), selection_RKD
    real(RKH) :: array_RKH(NP), selection_RKH

    ! 1-dimensional array of integer values.

    integer(IKL) :: array_IKL  (NP), selection_IKL
    integer(IKS) :: array_IKS (NP), selection_IKS
    integer(IKD) :: array_IKD (NP), selection_IKD
    integer(IKH) :: array_IKH (NP), selection_IKH

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

    call random_number(array_RKL); array_RKL = array_RKL - 0.5_RKL
    call random_number(array_RKD); array_RKD = array_RKD - 0.5_RKD
    call random_number(array_RKH); array_RKH = array_RKH - 0.5_RKH

    array_IKL = int(array_RKH * huge(0_IKL), kind = IKL)
    array_IKS = int(array_RKH * huge(0_IKS), kind = IKS)
    array_IKD = int(array_RKH * huge(0_IKD), kind = IKD)
    array_IKH = int(array_RKH * huge(0_IKH), kind = IKH)

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
    call disp%show("array_IKL")
    call disp%show( array_IKL )
    call disp%show("call setSelected(selection_IKL, array_IKL, rank)")
                    call setSelected(selection_IKL, array_IKL, rank)
    call disp%show("selection_IKL")
    call disp%show( selection_IKL )
    call disp%show("array_IKL")
    call disp%show( array_IKL )
    call disp%skip()

    rank = 4_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IKS")
    call disp%show( array_IKS )
    call disp%show("call setSelected(selection_IKS, array_IKS, rank)")
                    call setSelected(selection_IKS, array_IKS, rank)
    call disp%show("selection_IKS")
    call disp%show( selection_IKS )
    call disp%show("array_IKS")
    call disp%show( array_IKS )
    call disp%skip()

    rank = 1_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IKD")
    call disp%show( array_IKD )
    call disp%show("call setSelected(selection_IKD, array_IKD, rank)")
                    call setSelected(selection_IKD, array_IKD, rank)
    call disp%show("selection_IKD")
    call disp%show( selection_IKD )
    call disp%show("array_IKD")
    call disp%show( array_IKD )
    call disp%skip()

    rank = 5_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IKH")
    call disp%show( array_IKH )
    call disp%show("call setSelected(selection_IKH, array_IKH, rank)")
                    call setSelected(selection_IKH, array_IKH, rank)
    call disp%show("selection_IKH")
    call disp%show( selection_IKH )
    call disp%show("array_IKH")
    call disp%show( array_IKH )
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
    call disp%show("array_RKL")
    call disp%show( array_RKL )
    call disp%show("call setSelected(selection_RKL, array_RKL, rank)")
                    call setSelected(selection_RKL, array_RKL, rank)
    call disp%show("selection_RKL")
    call disp%show( selection_RKL )
    call disp%show("array_RKL")
    call disp%show( array_RKL )
    call disp%skip()

    rank = 2_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_RKD")
    call disp%show( array_RKD )
    call disp%show("call setSelected(selection_RKD, array_RKD, rank)")
                    call setSelected(selection_RKD, array_RKD, rank)
    call disp%show("selection_RKD")
    call disp%show( selection_RKD )
    call disp%show("array_RKD")
    call disp%show( array_RKD )
    call disp%skip()

    rank = size(array_RKH, kind = IK)
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_RKH")
    call disp%show( array_RKH )
    call disp%show("call setSelected(selection_RKH, array_RKH, rank)")
                    call setSelected(selection_RKH, array_RKH, rank)
    call disp%show("selection_RKH")
    call disp%show( selection_RKH )
    call disp%show("array_RKH")
    call disp%show( array_RKH )
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

    array_IKL = int([((-2)**i, i = 1, NP)], kind = IKL)
    rank = 3_IK
    call disp%skip()
    call disp%show("!Select the `rank`th largest element via an input custom-designed `isSorted()` function.")
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IKL")
    call disp%show( array_IKL )
    call disp%show("call setSelected(selection_IKL, array_IKL, rank, isSorted_IKL)")
                    call setSelected(selection_IKL, array_IKL, rank, isSorted_IKL)
    call disp%show("selection_IKL")
    call disp%show( selection_IKL )
    call disp%show("array_IKL")
    call disp%show( array_IKL )
    call disp%skip()

    call random_number(array_RKL); array_RKL  = array_RKL  - 0.5_RKL
    call disp%skip()
    call disp%show("!Select the `rank`th smallest element solely based on the magnitude of numbers using a custom comparison function.")
    call disp%skip()
    rank = 3_IK
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_RKL")
    call disp%show( array_RKL )
    call disp%show("call setSelected(selection_RKL, array_RKL, rank, isSorted_RKL)")
                    call setSelected(selection_RKL, array_RKL, rank, isSorted_RKL)
    call disp%show("selection_RKL")
    call disp%show( selection_RKL )
    call disp%show("array_RKL")
    call disp%show( array_RKL )
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

    array_IKL = [-3_IKL, -1_IKL, -2_IKL, 1_IKL, 2_IKL, 3_IKL] ! all elements after index `3` are sorted (`ub = 3_IK`).
    rank = 3_IK
    call disp%skip()
    call disp%show("rank")
    call disp%show( rank )
    call disp%show("array_IKL")
    call disp%show( array_IKL )
    call disp%show("call setSelected(selection_IKL, array_IKL, rank, ub = 3_IK) ! all elements after index `3` are sorted (`ub = 3_IK`).")
                    call setSelected(selection_IKL, array_IKL, rank, ub = 3_IK)
    call disp%show("selection_IKL")
    call disp%show( selection_IKL )
    call disp%show("array_IKL")
    call disp%show( array_IKL )
    call disp%skip()

    array_IKL = [-3_IKL, -1_IKL, 3_IKL, -2_IKL, 2_IKL, 1_IKL]
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
            call disp%show("array_IKL")
            call disp%show( array_IKL )
            call disp%show("call setSelected(selection_IKL, array_IKL, rank = RankList(i), lb = RankList(i-1)+1) ! all elements after `lb` need sorting.")
                            call setSelected(selection_IKL, array_IKL, rank = RankList(i), lb = RankList(i-1)+1)
            call disp%show("selection_IKL")
            call disp%show( selection_IKL )
            call disp%show("array_IKL")
            call disp%show( array_IKL )
        end do
    end block
    call disp%skip()

contains

    function isSorted_IKL(a,b) result(isSorted)
        use pm_kind, only: LK, IKL
        integer(IKL)    , intent(in)    :: a, b
        logical(LK)                     :: isSorted
        isSorted = a > b
    end function

    function isSorted_RKL(a,b) result(isSorted)
        use pm_kind, only: LK, IKL
        integer(IKL)    , intent(in)    :: a, b
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