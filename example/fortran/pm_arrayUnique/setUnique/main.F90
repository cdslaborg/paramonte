program example

    use pm_kind, only: LK
    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayUnique, only: setUnique
    use pm_container, only: Vector => cvi_type

    implicit none

    type(Vector)    , allocatable   :: index(:)                         ! The jagged-array of locations of occurrences of the unique elements.
    integer(IK)     , allocatable   :: Count(:)                         ! Must be of default kind IK
    character(:, SK), allocatable   :: string_SK    , stringUnique_SK   
    character(9, SK), allocatable   :: Array_SK(:)  , UniqueArray_SK(:) ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: Array_IK(:)  , UniqueArray_IK(:) ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: Array_CK(:)  , UniqueArray_CK(:) ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: Array_RK(:)  , UniqueArray_RK(:) ! Can be any processor-supported kind.

    type(display_type)              :: disp

    disp = display_type(file = "main.out.F90")

    string_SK = "ParaMonte is a Monte Carlo Library."
    Array_SK = ["ParaMonte", "PARAMONTE", "paramonte", "ParaMonte", "ParaMonte", "Paramonte"]
    Array_IK = [1_IK, 1_IK, 3_IK, 4_IK, 4_IK, 4_IK]
    Array_RK = [1._RK, 1._RK, 3._RK, 4._RK, 4._RK, 4._RK]
    Array_CK = [(1._CK, -1._CK), (1._CK, -1._CK), (3._CK, -3._CK), (4._CK, -4._CK), (4._CK, -4._CK), (4._CK, -4._CK)]


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find all unique elements along with their frequency and locations of occurrences in array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in character scalar along with their frequency and locations of occurrences.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("call setUnique(string_SK, stringUnique_SK, Count, index)")
                    call setUnique(string_SK, stringUnique_SK, Count, index)
    call disp%show("stringUnique_SK")
    call disp%show( stringUnique_SK, deliml = SK_"""" )
    call disp%show("Count")
    call disp%show( Count )
    call disp%show("index")
    call disp%show( index )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in character array along with their frequency and locations of occurrences.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_SK")
    call disp%show( Array_SK, deliml = SK_"""" )
    call disp%show("call setUnique(Array_SK, UniqueArray_SK, Count, index)")
                    call setUnique(Array_SK, UniqueArray_SK, Count, index)
    call disp%show("UniqueArray_SK")
    call disp%show( UniqueArray_SK, deliml = SK_"""" )
    call disp%show("Count")
    call disp%show( Count )
    call disp%show("index")
    call disp%show( index )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in integer array along with their frequency and locations of occurrences.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_IK")
    call disp%show( Array_IK )
    call disp%show("call setUnique(Array_IK, UniqueArray_IK, Count, index)")
                    call setUnique(Array_IK, UniqueArray_IK, Count, index)
    call disp%show("UniqueArray_IK")
    call disp%show( UniqueArray_IK )
    call disp%show("Count")
    call disp%show( Count )
    call disp%show("index")
    call disp%show( index )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in complex array along with their frequency and locations of occurrences.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array_CK")
    call disp%show( Array_CK )
    call disp%show("call setUnique(Array_CK, UniqueArray_CK, Count, index)")
                    call setUnique(Array_CK, UniqueArray_CK, Count, index)
    call disp%show("UniqueArray_CK")
    call disp%show( UniqueArray_CK )
    call disp%show("Count")
    call disp%show( Count )
    call disp%show("index")
    call disp%show( index )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in real array along with their frequency and locations of occurrences.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_RK = Array_RK
    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("call setUnique(Array_RK, UniqueArray_RK, Count, index)")
                    call setUnique(Array_RK, UniqueArray_RK, Count, index)
    call disp%show("UniqueArray_RK")
    call disp%show( UniqueArray_RK )
    call disp%show("Count")
    call disp%show( Count )
    call disp%show("index")
    call disp%show( index )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Find all unique elements according to the user-specified equivalence criterion.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique elements in real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    Array_RK = [1.01_RK, 1.04_RK, 0.98_RK, 1.0_RK, 1.02_RK, 2._RK]

    call disp%show("Array_RK")
    call disp%show( Array_RK )
    call disp%show("call setUnique(Array_RK, UniqueArray_RK, Count, iseq = iseq_RK)")
                    call setUnique(Array_RK, UniqueArray_RK, Count, iseq = iseq_RK)
    call disp%show("UniqueArray_RK")
    call disp%show( UniqueArray_RK )
    call disp%show("Count")
    call disp%show( Count )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Unique case-insensitive instances within the character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    string_SK = "ABBAbbA"

    call disp%show("string_SK")
    call disp%show( string_SK, deliml = SK_"""" )
    call disp%show("call setUnique(string_SK, stringUnique_SK, Count, iseq = iseq_SK)")
                    call setUnique(string_SK, stringUnique_SK, Count, iseq = iseq_SK)
    call disp%show("stringUnique_SK")
    call disp%show( stringUnique_SK, deliml = SK_"""" )
    call disp%show("Count")
    call disp%show( Count )

contains

    pure function iseq_RK(element1, element2) result(iseq)
        real(RK)    , intent(in)    :: element1, element2
        logical(LK)                 :: iseq
        iseq = abs(abs(element1) - abs(element2)) < 0.05_RK
    end function

    pure function iseq_SK(element1, element2) result(iseq)
        use pm_strASCII, only: getStrLower
        character(*, SK)    , intent(in)    :: element1, element2
        logical(LK)                     :: iseq
        iseq = getStrLower(element1) == getStrLower(element2)
    end function

end program example