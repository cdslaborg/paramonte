program example

    use pm_kind, only: SK ! default string kind
    use pm_kind, only: LK
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayRemap, only: setRemapped, reverse

    implicit none

    integer(IK)     , allocatable   :: index(:) ! Must be of default kind IK
    character(:, SK), allocatable   :: string1_SK   , string2_SK    
    character(9, SK), allocatable   :: Array1_SK(:) , Array2_SK(:)  ! Can be any processor-supported kind.
    integer(IK)     , allocatable   :: Array1_IK(:) , Array2_IK(:)  ! Can be any processor-supported kind.
    complex(CK)     , allocatable   :: Array1_CK(:) , Array2_CK(:)  ! Can be any processor-supported kind.
    real(RK)        , allocatable   :: Array1_RK(:) , Array2_RK(:)  ! Can be any processor-supported kind.
    logical(LK)     , allocatable   :: Array1_LK(:) , Array2_LK(:)

    type(display_type)              :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap an ALLOCATABLE array IN-PLACE.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    index = [2_IK, 1_IK, 4_IK, 3_IK, 3_IK, 2_IK]

    string1_SK = "ABCDEF"
    Array1_SK = ["AA", "BB", "CC", "DD", "EE", "FF"]
    Array1_IK = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK]
    Array1_RK = [1._RK, 2._RK, 3._RK, 4._RK, 5._RK, 6._RK]
    Array1_CK = [(1._CK, -1._CK), (2._CK, -2._CK), (3._CK, -3._CK), (4._CK, -4._CK), (5._CK, -5._CK), (6._CK, -6._CK)]
    Array1_LK = [.false._LK, .true._LK, .false._LK, .true._LK, .false._LK, .true._LK]

    string2_SK = string1_SK
    Array2_SK = Array1_SK
    Array2_IK = Array1_IK
    Array2_RK = Array1_RK
    Array2_CK = Array1_CK
    Array2_LK = Array1_LK

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string2_SK")
    call disp%show( string2_SK, deliml = SK_"""" )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(string2_SK, index)")
                    call setRemapped(string2_SK, index)
    call disp%show("string2_SK")
    call disp%show( string2_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_SK")
    call disp%show( Array2_SK, deliml = SK_"""" )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array2_SK, index)")
                    call setRemapped(Array2_SK, index)
    call disp%show("Array2_SK")
    call disp%show( Array2_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_LK")
    call disp%show( Array2_LK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array2_LK, index)")
                    call setRemapped(Array2_LK, index)
    call disp%show("Array2_LK")
    call disp%show( Array2_LK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_IK")
    call disp%show( Array2_IK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array2_IK, index)")
                    call setRemapped(Array2_IK, index)
    call disp%show("Array2_IK")
    call disp%show( Array2_IK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_CK")
    call disp%show( Array2_CK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array2_CK, index)")
                    call setRemapped(Array2_CK, index)
    call disp%show("Array2_CK")
    call disp%show( Array2_CK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_RK")
    call disp%show( Array2_RK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array2_RK, index)")
                    call setRemapped(Array2_RK, index)
    call disp%show("Array2_RK")
    call disp%show( Array2_RK )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap an ALLOCATABLE array IN-PLACE in backward direction.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    index = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK]
    string2_SK = string1_SK
    Array2_SK = Array1_SK
    Array2_IK = Array1_IK
    Array2_RK = Array1_RK
    Array2_CK = Array1_CK
    Array2_LK = Array1_LK

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string2_SK")
    call disp%show( string2_SK, deliml = SK_"""" )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(string2_SK, index, reverse)")
                    call setRemapped(string2_SK, index, reverse)
    call disp%show("string2_SK")
    call disp%show( string2_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_SK")
    call disp%show( Array2_SK, deliml = SK_"""" )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array2_SK, index, reverse)")
                    call setRemapped(Array2_SK, index, reverse)
    call disp%show("Array2_SK")
    call disp%show( Array2_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_LK")
    call disp%show( Array2_LK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array2_LK, index, reverse)")
                    call setRemapped(Array2_LK, index, reverse)
    call disp%show("Array2_LK")
    call disp%show( Array2_LK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_IK")
    call disp%show( Array2_IK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array2_IK, index, reverse)")
                    call setRemapped(Array2_IK, index, reverse)
    call disp%show("Array2_IK")
    call disp%show( Array2_IK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_CK")
    call disp%show( Array2_CK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array2_CK, index, reverse)")
                    call setRemapped(Array2_CK, index, reverse)
    call disp%show("Array2_CK")
    call disp%show( Array2_CK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array2_RK")
    call disp%show( Array2_RK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array2_RK, index, reverse)")
                    call setRemapped(Array2_RK, index, reverse)
    call disp%show("Array2_RK")
    call disp%show( Array2_RK )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap input array and return it in a new array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap character scalar.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("string1_SK")
    call disp%show( string1_SK, deliml = SK_"""" )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(string1_SK, index, reverse, arrayNew = string2_SK)")
                    call setRemapped(string1_SK, index, reverse, arrayNew = string2_SK)
    call disp%show("string1_SK")
    call disp%show( string1_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap character array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array1_SK")
    call disp%show( Array1_SK, deliml = SK_"""" )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array1_SK, index, reverse, arrayNew = Array2_SK)")
                    call setRemapped(Array1_SK, index, reverse, arrayNew = Array2_SK)
    call disp%show("Array1_SK")
    call disp%show( Array2_SK, deliml = SK_"""" )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap logical array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array1_LK")
    call disp%show( Array1_LK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array1_LK, index, reverse, arrayNew = Array2_LK)")
                    call setRemapped(Array1_LK, index, reverse, arrayNew = Array2_LK)
    call disp%show("Array2_LK")
    call disp%show( Array2_LK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap integer array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array1_IK")
    call disp%show( Array1_IK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array1_IK, index, reverse, arrayNew = Array2_IK)")
                    call setRemapped(Array1_IK, index, reverse, arrayNew = Array2_IK)
    call disp%show("Array2_IK")
    call disp%show( Array2_IK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap complex array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array1_CK")
    call disp%show( Array1_CK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array1_CK, index, reverse, arrayNew = Array2_CK)")
                    call setRemapped(Array1_CK, index, reverse, arrayNew = Array2_CK)
    call disp%show("Array2_CK")
    call disp%show( Array2_CK )

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remap real array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("Array1_RK")
    call disp%show( Array1_RK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array1_RK, index, reverse, arrayNew = Array2_RK)")
                    call setRemapped(Array1_RK, index, reverse, arrayNew = Array2_RK)
    call disp%show("Array2_RK")
    call disp%show( Array2_RK )


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Remapping preserves the lower bound of the input allocatable array.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    index = [2_IK, 1_IK, 4_IK, 3_IK, 3_IK, 2_IK] - 4
    deallocate(Array1_IK); 
    allocate(Array1_IK(-5 : -5 + size(index) - 1))
    Array1_IK(:) = [1_IK, 2_IK, 3_IK, 4_IK, 5_IK, 6_IK]

    call disp%show("[lbound(Array1_IK,1), ubound(Array1_IK,1)]")
    call disp%show( [lbound(Array1_IK,1), ubound(Array1_IK,1)] )
    call disp%show("Array1_IK")
    call disp%show( Array1_IK )
    call disp%show("index")
    call disp%show( index )
    call disp%show("call setRemapped(Array1_IK, index)")
                    call setRemapped(Array1_IK, index)
    call disp%show("Array1_IK")
    call disp%show( Array1_IK )
    call disp%show("[lbound(Array1_IK,1), ubound(Array1_IK,1)]")
    call disp%show( [lbound(Array1_IK,1), ubound(Array1_IK,1)] )

end program example