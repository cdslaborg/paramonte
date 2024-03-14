program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayInsert, only: setInserted

    implicit none

    integer(IK) :: i
    integer(IK) , allocatable   :: index(:) ! Must be of default kind IK

    character(:, SK), allocatable   :: string_SK    , stringNew_SK  , stringInsertion_SK    
    character(9, SK), allocatable   :: Array_SK(:)  , ArrayNew_SK(:), ArrayInsertion_SK(:)  ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: Array_IK(:)  , ArrayNew_IK(:), ArrayInsertion_IK(:)  ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: Array_CK(:)  , ArrayNew_CK(:), ArrayInsertion_CK(:)  ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: Array_RK(:)  , ArrayNew_RK(:), ArrayInsertion_RK(:)  ! Can be any processor-supported kind.
    logical(LK)     , allocatable   :: Array_LK(:)  , ArrayNew_LK(:), ArrayInsertion_LK(:)

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert an array-like `insertion` at the specified locations in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "OOOOOOOOO"
    Array_SK = ["O", "O", "O", "O", "O", "O", "O", "O", "O"]
    Array_IK = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK, 7_IK, 8_IK, 9_IK]
    Array_RK = [1._RK, 2._RK, 3._RK, 4._RK, 5._RK, 6._RK, 7._RK, 8._RK, 9._RK]
    Array_CK = [(1._CK, -1._CK), (2._CK, -2._CK), (3._CK, -3._CK), (4._CK, -4._CK), (5._CK, -5._CK), (6._CK, -6._CK), (7._CK, -7._CK), (8._CK, -8._CK), (9._CK, -9._CK)]
    Array_LK = [.false._LK, .false._LK, .false._LK, .false._LK, .false._LK, .false._LK, .false._LK, .false._LK, .false._LK]

    stringInsertion_SK = "+%"
    ArrayInsertion_SK = ["++", "%%"]
    ArrayInsertion_IK = [0_IK, -1_IK]
    ArrayInsertion_RK = [0._RK, -1._RK]
    ArrayInsertion_CK = [(0._CK, -0._CK), (-1._CK, +1._CK)]
    ArrayInsertion_LK = [.true._LK, .true._LK]

    index = int([1, 1, 1, 2, 4, 4, -1, size(Array_SK) + 1], kind = IK) ! Duplicate indices result in multiple instances being inserted in the same location.

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert string `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(character(len(string_SK)+size(index)*len(stringInsertion_SK)) :: stringNew_SK)
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringInsertion_SK")
    call disp%show( stringInsertion_SK, deliml = SK_"""" )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(stringNew_SK, string_SK, stringInsertion_SK, index = index)")
                    call setInserted(stringNew_SK, string_SK, stringInsertion_SK, index = index)
    call disp%show("stringNew_SK")
    call disp%show( stringNew_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert character `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(ArrayNew_SK(size(Array_SK) + size(index) * size(ArrayInsertion_SK)))
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("ArrayInsertion_SK")
    call disp%show( ArrayInsertion_SK, deliml = SK_"""" )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(ArrayNew_SK, Array_SK, ArrayInsertion_SK, index = index)")
                    call setInserted(ArrayNew_SK, Array_SK, ArrayInsertion_SK, index = index)
    call disp%show("ArrayNew_SK")
    call disp%show( ArrayNew_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert logical `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(ArrayNew_LK(size(Array_LK) + size(index) * size(ArrayInsertion_LK)))
    call disp%show("Array_LK")
    call disp%show( Array_LK )
    call disp%show("ArrayInsertion_LK")
    call disp%show( ArrayInsertion_LK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(ArrayNew_LK, Array_LK, ArrayInsertion_LK, index = index)")
                    call setInserted(ArrayNew_LK, Array_LK, ArrayInsertion_LK, index = index)
    call disp%show("ArrayNew_LK")
    call disp%show( ArrayNew_LK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert integer `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(ArrayNew_IK(size(Array_IK) + size(index) * size(ArrayInsertion_IK)))
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("ArrayInsertion_IK")
    call disp%show( ArrayInsertion_IK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(ArrayNew_IK, Array_IK, ArrayInsertion_IK, index = index)")
                    call setInserted(ArrayNew_IK, Array_IK, ArrayInsertion_IK, index = index)
    call disp%show("ArrayNew_IK")
    call disp%show( ArrayNew_IK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert complex `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(ArrayNew_CK(size(Array_CK) + size(index) * size(ArrayInsertion_CK)))
    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("ArrayInsertion_CK")
    call disp%show( ArrayInsertion_CK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(ArrayNew_CK, Array_CK, ArrayInsertion_CK, index = index)")
                    call setInserted(ArrayNew_CK, Array_CK, ArrayInsertion_CK, index = index)
    call disp%show("ArrayNew_CK")
    call disp%show( ArrayNew_CK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert real `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(ArrayNew_RK(size(Array_RK) + size(index) * size(ArrayInsertion_RK)))
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("ArrayInsertion_RK")
    call disp%show( ArrayInsertion_RK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(ArrayNew_RK, Array_RK, ArrayInsertion_RK, index = index)")
                    call setInserted(ArrayNew_RK, Array_RK, ArrayInsertion_RK, index = index)
    call disp%show("ArrayNew_RK")
    call disp%show( ArrayNew_RK )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert an scalar `insertion` at the specified locations in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call reset()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert string `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(character(len(string_SK)+size(index)) :: stringNew_SK)
    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("stringInsertion_SK(1:1)")
    call disp%show( stringInsertion_SK(1:1), deliml = SK_"""" )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(stringNew_SK, string_SK, stringInsertion_SK(1:1), index = index)")
                    call setInserted(stringNew_SK, string_SK, stringInsertion_SK(1:1), index = index)
    call disp%show("stringNew_SK")
    call disp%show( stringNew_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert character `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(ArrayNew_SK(size(Array_SK) + size(index)))
    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("ArrayInsertion_SK(1)")
    call disp%show( ArrayInsertion_SK(1), deliml = SK_"""" )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(ArrayNew_SK, Array_SK, ArrayInsertion_SK(1), index = index)")
                    call setInserted(ArrayNew_SK, Array_SK, ArrayInsertion_SK(1), index = index)
    call disp%show("ArrayNew_SK")
    call disp%show( ArrayNew_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert logical `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(ArrayNew_LK(size(Array_LK) + size(index)))
    call disp%show("Array_LK")
    call disp%show( Array_LK )
    call disp%show("ArrayInsertion_LK")
    call disp%show( ArrayInsertion_LK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(ArrayNew_LK, Array_LK, ArrayInsertion_LK(1), index = index)")
                    call setInserted(ArrayNew_LK, Array_LK, ArrayInsertion_LK(1), index = index)
    call disp%show("ArrayNew_LK")
    call disp%show( ArrayNew_LK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert integer `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(ArrayNew_IK(size(Array_IK) + size(index)))
    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("ArrayInsertion_IK")
    call disp%show( ArrayInsertion_IK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(ArrayNew_IK, Array_IK, ArrayInsertion_IK(1), index = index)")
                    call setInserted(ArrayNew_IK, Array_IK, ArrayInsertion_IK(1), index = index)
    call disp%show("ArrayNew_IK")
    call disp%show( ArrayNew_IK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert complex `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(ArrayNew_CK(size(Array_CK) + size(index)))
    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("ArrayInsertion_CK")
    call disp%show( ArrayInsertion_CK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(ArrayNew_CK, Array_CK, ArrayInsertion_CK(1), index = index)")
                    call setInserted(ArrayNew_CK, Array_CK, ArrayInsertion_CK(1), index = index)
    call disp%show("ArrayNew_CK")
    call disp%show( ArrayNew_CK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Insert real `insertion`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(ArrayNew_RK(size(Array_RK) + size(index)))
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("ArrayInsertion_RK")
    call disp%show( ArrayInsertion_RK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setInserted(ArrayNew_RK, Array_RK, ArrayInsertion_RK(1), index = index)")
                    call setInserted(ArrayNew_RK, Array_RK, ArrayInsertion_RK(1), index = index)
    call disp%show("ArrayNew_RK")
    call disp%show( ArrayNew_RK )

contains

    subroutine reset()
        if (allocated(stringNew_SK)) deallocate(stringNew_SK)
        if (allocated(ArrayNew_SK)) deallocate(ArrayNew_SK)
        if (allocated(ArrayNew_IK)) deallocate(ArrayNew_IK)
        if (allocated(ArrayNew_RK)) deallocate(ArrayNew_RK)
        if (allocated(ArrayNew_CK)) deallocate(ArrayNew_CK)
        if (allocated(ArrayNew_LK)) deallocate(ArrayNew_LK)
    end subroutine

end program example