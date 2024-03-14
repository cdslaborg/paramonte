program example

    use pm_kind, only: SK, IK, LK, CK, RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_arrayCompact, only: setCompact
    use pm_arrayVerbose, only: getVerbose

    implicit none

    integer(IK) , parameter :: ND = 2_IK, NP = 8_IK
    integer(IK)             :: Weight(NP)
    integer(IK)             :: newsize
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
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(string, Weight, newsize)")
                        call setCompact(string, Weight, newsize)
        call disp%show("string(1:newsize)")
        call disp%show( string(1:newsize) , deliml = SK_"""" )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(string(1:newsize), Weight(1:newsize), sum(Weight(1:newsize)))")
        call disp%show( getVerbose(string(1:newsize), Weight(1:newsize), sum(Weight(1:newsize))) , deliml = SK_"""" )

        call disp%skip()
        call disp%show("Array_SK")
        call disp%show( Array_SK , deliml = SK_"""" )
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(Array_SK, Weight, newsize)")
                        call setCompact(Array_SK, Weight, newsize)
        call disp%show("Array_SK(1:newsize)")
        call disp%show( Array_SK(1:newsize) , deliml = SK_"""" )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(Array_SK(1:newsize), Weight(1:newsize), sum(Weight(1:newsize)))")
        call disp%show( getVerbose(Array_SK(1:newsize), Weight(1:newsize), sum(Weight(1:newsize))) , deliml = SK_"""" )

        call disp%skip()
        call disp%show("Array_IK")
        call disp%show( Array_IK )
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(Array_IK, Weight, newsize)")
                        call setCompact(Array_IK, Weight, newsize)
        call disp%show("Array_IK(1:newsize)")
        call disp%show( Array_IK(1:newsize) )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(Array_IK(1:newsize), Weight(1:newsize), sum(Weight(1:newsize)))")
        call disp%show( getVerbose(Array_IK(1:newsize), Weight(1:newsize), sum(Weight(1:newsize))) )

        call disp%skip()
        call disp%show("Array_LK")
        call disp%show( Array_LK )
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(Array_LK, Weight, newsize)")
                        call setCompact(Array_LK, Weight, newsize)
        call disp%show("Array_LK(1:newsize)")
        call disp%show( Array_LK(1:newsize) )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(Array_LK(1:newsize), Weight(1:newsize), sum(Weight(1:newsize)))")
        call disp%show( getVerbose(Array_LK(1:newsize), Weight(1:newsize), sum(Weight(1:newsize))) )

        call disp%skip()
        call disp%show("Array_CK")
        call disp%show( Array_CK )
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(Array_CK, Weight, newsize)")
                        call setCompact(Array_CK, Weight, newsize)
        call disp%show("Array_CK(1:newsize)")
        call disp%show( Array_CK(1:newsize) )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(Array_CK(1:newsize), Weight(1:newsize), sum(Weight(1:newsize)))")
        call disp%show( getVerbose(Array_CK(1:newsize), Weight(1:newsize), sum(Weight(1:newsize))) )

        call disp%skip()
        call disp%show("Array_RK")
        call disp%show( Array_RK )
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(Array_RK, Weight, newsize)")
                        call setCompact(Array_RK, Weight, newsize)
        call disp%show("Array_RK(1:newsize)")
        call disp%show( Array_RK(1:newsize) )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(Array_RK(1:newsize), Weight(1:newsize), sum(Weight(1:newsize)))")
        call disp%show( getVerbose(Array_RK(1:newsize), Weight(1:newsize), sum(Weight(1:newsize))) )

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
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(Array_SK, Weight, newsize, dim = 2_IK)")
                        call setCompact(Array_SK, Weight, newsize, dim = 2_IK)
        call disp%show("Array_SK(:,1:newsize)")
        call disp%show( Array_SK(:,1:newsize) , deliml = SK_"""" )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(Array_SK(:,1:newsize), Weight(1:newsize), sum(Weight(1:newsize)), dim = 2_IK)")
        call disp%show( getVerbose(Array_SK(:,1:newsize), Weight(1:newsize), sum(Weight(1:newsize)), dim = 2_IK) , deliml = SK_"""" )

        call disp%skip()
        call disp%show("Array_IK")
        call disp%show( Array_IK )
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(Array_IK, Weight, newsize, dim = 2_IK)")
                        call setCompact(Array_IK, Weight, newsize, dim = 2_IK)
        call disp%show("Array_IK(:,1:newsize)")
        call disp%show( Array_IK(:,1:newsize) )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(Array_IK(:,1:newsize), Weight(1:newsize), sum(Weight(1:newsize)), dim = 2_IK)")
        call disp%show( getVerbose(Array_IK(:,1:newsize), Weight(1:newsize), sum(Weight(1:newsize)), dim = 2_IK) )

        call disp%skip()
        call disp%show("Array_LK")
        call disp%show( Array_LK )
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(Array_LK, Weight, newsize, dim = 2_IK)")
                        call setCompact(Array_LK, Weight, newsize, dim = 2_IK)
        call disp%show("Array_LK(:,1:newsize)")
        call disp%show( Array_LK(:,1:newsize) )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(Array_LK(:,1:newsize), Weight(1:newsize), sum(Weight(1:newsize)), dim = 2_IK)")
        call disp%show( getVerbose(Array_LK(:,1:newsize), Weight(1:newsize), sum(Weight(1:newsize)), dim = 2_IK) )

        call disp%skip()
        call disp%show("Array_CK")
        call disp%show( Array_CK )
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(Array_CK, Weight, newsize, dim = 2_IK)")
                        call setCompact(Array_CK, Weight, newsize, dim = 2_IK)
        call disp%show("Array_CK(:,1:newsize)")
        call disp%show( Array_CK(:,1:newsize) )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(Array_CK(:,1:newsize), Weight(1:newsize), sum(Weight(1:newsize)), dim = 2_IK)")
        call disp%show( getVerbose(Array_CK(:,1:newsize), Weight(1:newsize), sum(Weight(1:newsize)), dim = 2_IK) )

        call disp%skip()
        call disp%show("Array_RK")
        call disp%show( Array_RK )
        call disp%show("size(Weight)")
        call disp%show( size(Weight) )
        call disp%show("call setCompact(Array_RK, Weight, newsize, dim = 2_IK)")
                        call setCompact(Array_RK, Weight, newsize, dim = 2_IK)
        call disp%show("Array_RK(:,1:newsize)")
        call disp%show( Array_RK(:,1:newsize) )
        call disp%show("Weight(1:newsize)")
        call disp%show( Weight(1:newsize) )
        call disp%show("getVerbose(Array_RK(:,1:newsize), Weight(1:newsize), sum(Weight(1:newsize)), dim = 2_IK)")
        call disp%show( getVerbose(Array_RK(:,1:newsize), Weight(1:newsize), sum(Weight(1:newsize)), dim = 2_IK) )

    end block

end program example