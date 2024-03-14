program example

    use pm_kind, only: SK, IK, LK, CK, RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_arrayCompact, only: getCompact

    implicit none

    integer(IK) , parameter :: ND = 2_IK, NP = 8_IK
    type(display_type)      :: disp

    disp = display_type(file = "main.out.F90")

    ! Condense 1D array.

    block

        character(:, SK), allocatable   :: string
        character(2, SK), allocatable   :: Array_SK(:)
        integer(IK)     , allocatable   :: Array_IK(:)
        logical(LK)     , allocatable   :: Array_LK(:)
        complex(CK)     , allocatable   :: Array_CK(:)
        real(RK)        , allocatable   :: Array_RK(:)

        string = "AABCCCDC"
        Array_SK = [ "AA", "AA", "BB", "CC", "CC", "CC", "DD", "CC" ]
        Array_LK = [ .false., .false., .true., .false., .false., .false., .true., .false. ]
        Array_IK = [ 1, 1, 2, 3, 3, 3, 4, 3 ]
        Array_CK = [ 1, 1, 2, 3, 3, 3, 4, 3 ]
        Array_RK = [ 1, 1, 2, 3, 3, 3, 4, 3 ]

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Condense a 1D array.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%skip()
        call disp%show("string")
        call disp%show( string , deliml = SK_"""" )
        call disp%show("getCompact(string)")
        call disp%show( getCompact(string) , deliml = SK_"""" )

        call disp%skip()
        call disp%show("Array_SK")
        call disp%show( Array_SK , deliml = SK_"""" )
        call disp%show("getCompact(Array_SK)")
        call disp%show( getCompact(Array_SK) , deliml = SK_"""" )

        call disp%skip()
        call disp%show("Array_IK")
        call disp%show( Array_IK )
        call disp%show("getCompact(Array_IK)")
        call disp%show( getCompact(Array_IK) )

        call disp%skip()
        call disp%show("Array_LK")
        call disp%show( Array_LK )
        call disp%show("getCompact(Array_LK)")
        call disp%show( getCompact(Array_LK) )

        call disp%skip()
        call disp%show("Array_CK")
        call disp%show( Array_CK )
        call disp%show("getCompact(Array_CK)")
        call disp%show( getCompact(Array_CK) )

        call disp%skip()
        call disp%show("Array_RK")
        call disp%show( Array_RK )
        call disp%show("getCompact(Array_RK)")
        call disp%show( getCompact(Array_RK) )

    end block

    ! Flatten 2D array.

    block

        character(2, SK), allocatable   :: Array_SK(:,:)
        integer(IK)     , allocatable   :: Array_IK(:,:)
        logical(LK)     , allocatable   :: Array_LK(:,:)
        complex(CK)     , allocatable   :: Array_CK(:,:)
        real(RK)        , allocatable   :: Array_RK(:,:)

        Array_SK = transpose(reshape([ "AA", "AA", "BB", "CC", "CC", "CC", "DD", "CC", "EE", "EE", "FF", "GG", "GG", "GG", "HH", "GG" ], shape = [NP, ND]))
        Array_LK = transpose(reshape([ .false., .false., .true., .false., .false., .false., .true., .false., .true., .true., .false., .true., .true., .true., .false., .true. ], shape = [NP, ND]))
        Array_IK = transpose(reshape([ 1, 1, 2, 3, 3, 3, 4, 3, 5, 5, 6, 7, 7, 7, 8, 7 ], shape = [NP, ND]))
        Array_CK = transpose(reshape([ 1, 1, 2, 3, 3, 3, 4, 3, 5, 5, 6, 7, 7, 7, 8, 7 ], shape = [NP, ND]))
        Array_RK = transpose(reshape([ 1, 1, 2, 3, 3, 3, 4, 3, 5, 5, 6, 7, 7, 7, 8, 7 ], shape = [NP, ND]))

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Condense a 2D array along the desired axis `dim`.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%skip()
        call disp%show("Array_SK")
        call disp%show( Array_SK , deliml = SK_"""" )
        call disp%show("getCompact(Array_SK, dim = 2_IK)")
        call disp%show( getCompact(Array_SK, dim = 2_IK) , deliml = SK_"""" )

        call disp%skip()
        call disp%show("Array_IK")
        call disp%show( Array_IK )
        call disp%show("getCompact(Array_IK, dim = 2_IK)")
        call disp%show( getCompact(Array_IK, dim = 2_IK) )

        call disp%skip()
        call disp%show("Array_LK")
        call disp%show( Array_LK )
        call disp%show("getCompact(Array_LK, dim = 2_IK)")
        call disp%show( getCompact(Array_LK, dim = 2_IK) )

        call disp%skip()
        call disp%show("Array_CK")
        call disp%show( Array_CK )
        call disp%show("getCompact(Array_CK, dim = 2_IK)")
        call disp%show( getCompact(Array_CK, dim = 2_IK) )

        call disp%skip()
        call disp%show("Array_RK")
        call disp%show( Array_RK )
        call disp%show("getCompact(Array_RK, dim = 2_IK)")
        call disp%show( getCompact(Array_RK, dim = 2_IK) )

    end block

end program example