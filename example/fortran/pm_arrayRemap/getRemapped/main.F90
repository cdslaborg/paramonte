program example

    use pm_kind, only: SK ! All kinds are supported.
    use pm_kind, only: LK ! All kinds are supported.
    use pm_kind, only: IK ! All kinds are supported.
    use pm_kind, only: CK ! All kinds are supported.
    use pm_kind, only: RK ! All kinds are supported.
    use pm_io, only: display_type
    use pm_arrayRemap, only: getRemapped, reverse

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
    call disp%show("string2_SK = getRemapped(string2_SK, index)")
                    string2_SK = getRemapped(string2_SK, index)
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
    call disp%show("Array2_SK = getRemapped(Array2_SK, index)")
                    Array2_SK = getRemapped(Array2_SK, index)
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
    call disp%show("Array2_LK = getRemapped(Array2_LK, index)")
                    Array2_LK = getRemapped(Array2_LK, index)
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
    call disp%show("Array2_IK = getRemapped(Array2_IK, index)")
                    Array2_IK = getRemapped(Array2_IK, index)
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
    call disp%show("Array2_CK = getRemapped(Array2_CK, index)")
                    Array2_CK = getRemapped(Array2_CK, index)
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
    call disp%show("Array2_RK = getRemapped(Array2_RK, index)")
                    Array2_RK = getRemapped(Array2_RK, index)
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
    call disp%show("string2_SK = getRemapped(string2_SK, index, reverse)")
                    string2_SK = getRemapped(string2_SK, index, reverse)
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
    call disp%show("Array2_SK = getRemapped(Array2_SK, index, reverse)")
                    Array2_SK = getRemapped(Array2_SK, index, reverse)
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
    call disp%show("Array2_LK = getRemapped(Array2_LK, index, reverse)")
                    Array2_LK = getRemapped(Array2_LK, index, reverse)
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
    call disp%show("Array2_IK = getRemapped(Array2_IK, index, reverse)")
                    Array2_IK = getRemapped(Array2_IK, index, reverse)
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
    call disp%show("Array2_CK = getRemapped(Array2_CK, index, reverse)")
                    Array2_CK = getRemapped(Array2_CK, index, reverse)
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
    call disp%show("Array2_RK = getRemapped(Array2_RK, index, reverse)")
                    Array2_RK = getRemapped(Array2_RK, index, reverse)
    call disp%show("Array2_RK")
    call disp%show( Array2_RK )

end program example