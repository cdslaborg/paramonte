program example

    use pm_kind, only: SK, IK, LK, CK, RK ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_arrayVerbose, only: getVerbose

    implicit none

    integer(IK) , parameter :: ND = 2_IK, NP = 4_IK
    integer(IK)             :: Weight(NP) = [-1_IK, 2_IK, 0_IK, 1_IK]
    type(display_type)      :: disp

    disp = display_type(file = "main.out.F90")

    ! Flatten 1D array.

    block

        character(:, SK), allocatable   :: WeightedString
        character(2, SK), allocatable   :: WeightedArray_SK(:)
        integer(IK)     , allocatable   :: WeightedArray_IK(:)
        logical(LK)     , allocatable   :: WeightedArray_LK(:)
        complex(CK)     , allocatable   :: WeightedArray_CK(:)
        real(RK)        , allocatable   :: WeightedArray_RK(:)

        WeightedString = "ABCD"
        WeightedArray_SK = [ "AA", "BB", "CC", "DD" ]
        WeightedArray_IK = [ 1, 2, 3, 4 ]
        WeightedArray_LK = [ .false., .true., .false., .true. ]
        WeightedArray_CK = [ 1, 2, 3, 4 ]
        WeightedArray_RK = [ 1, 2, 3, 4 ]

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Flatten weighted 1D array.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedString")
        call disp%show( WeightedString , deliml = SK_"""" )
        call disp%show("getVerbose(WeightedString, Weight, sum(Weight, mask = Weight > 0_IK))")
        call disp%show( getVerbose(WeightedString, Weight, sum(Weight, mask = Weight > 0_IK)) , deliml = SK_"""" )

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedArray_SK")
        call disp%show( WeightedArray_SK , deliml = SK_"""" )
        call disp%show("getVerbose(WeightedArray_SK, Weight, sum(Weight, mask = Weight > 0_IK))")
        call disp%show( getVerbose(WeightedArray_SK, Weight, sum(Weight, mask = Weight > 0_IK)) , deliml = SK_"""" )

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedArray_IK")
        call disp%show( WeightedArray_IK )
        call disp%show("getVerbose(WeightedArray_IK, Weight, sum(Weight, mask = Weight > 0_IK))")
        call disp%show( getVerbose(WeightedArray_IK, Weight, sum(Weight, mask = Weight > 0_IK)) )

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedArray_LK")
        call disp%show( WeightedArray_LK )
        call disp%show("getVerbose(WeightedArray_LK, Weight, sum(Weight, mask = Weight > 0_IK))")
        call disp%show( getVerbose(WeightedArray_LK, Weight, sum(Weight, mask = Weight > 0_IK)) )

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedArray_CK")
        call disp%show( WeightedArray_CK )
        call disp%show("getVerbose(WeightedArray_CK, Weight, sum(Weight, mask = Weight > 0_IK))")
        call disp%show( getVerbose(WeightedArray_CK, Weight, sum(Weight, mask = Weight > 0_IK)) )

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedArray_RK")
        call disp%show( WeightedArray_RK )
        call disp%show("getVerbose(WeightedArray_RK, Weight, sum(Weight, mask = Weight > 0_IK))")
        call disp%show( getVerbose(WeightedArray_RK, Weight, sum(Weight, mask = Weight > 0_IK)) )

    end block

    ! Flatten 2D array.

    block

        character(2, SK), allocatable   :: WeightedArray_SK(:,:)
        integer(IK)     , allocatable   :: WeightedArray_IK(:,:)
        logical(LK)     , allocatable   :: WeightedArray_LK(:,:)
        complex(CK)     , allocatable   :: WeightedArray_CK(:,:)
        real(RK)        , allocatable   :: WeightedArray_RK(:,:)

        WeightedArray_SK = transpose(reshape([ "AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH" ], shape = [NP, ND]))
        WeightedArray_LK = transpose(reshape([ .false., .true., .false., .true., .false., .true., .false., .true. ], shape = [NP, ND]))
        WeightedArray_IK = transpose(reshape([ 1, 2, 3, 4, 5, 6, 7, 8 ], shape = [NP, ND]))
        WeightedArray_CK = transpose(reshape([ 1, 2, 3, 4, 5, 6, 7, 8 ], shape = [NP, ND]))
        WeightedArray_RK = transpose(reshape([ 1, 2, 3, 4, 5, 6, 7, 8 ], shape = [NP, ND]))

        call disp%skip()
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Flatten weighted 2D array along the desired axis `dim`.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip()

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedArray_SK")
        call disp%show( WeightedArray_SK , deliml = SK_"""" )
        call disp%show("getVerbose(WeightedArray_SK, Weight, sum(Weight, mask = Weight > 0_IK), dim = 2_IK)")
        call disp%show( getVerbose(WeightedArray_SK, Weight, sum(Weight, mask = Weight > 0_IK), dim = 2_IK) , deliml = SK_"""" )

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedArray_IK")
        call disp%show( WeightedArray_IK )
        call disp%show("getVerbose(WeightedArray_IK, Weight, sum(Weight, mask = Weight > 0_IK), dim = 2_IK)")
        call disp%show( getVerbose(WeightedArray_IK, Weight, sum(Weight, mask = Weight > 0_IK), dim = 2_IK) )

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedArray_LK")
        call disp%show( WeightedArray_LK )
        call disp%show("getVerbose(WeightedArray_LK, Weight, sum(Weight, mask = Weight > 0_IK), dim = 2_IK)")
        call disp%show( getVerbose(WeightedArray_LK, Weight, sum(Weight, mask = Weight > 0_IK), dim = 2_IK) )

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedArray_CK")
        call disp%show( WeightedArray_CK )
        call disp%show("getVerbose(WeightedArray_CK, Weight, sum(Weight, mask = Weight > 0_IK), dim = 2_IK)")
        call disp%show( getVerbose(WeightedArray_CK, Weight, sum(Weight, mask = Weight > 0_IK), dim = 2_IK) )

        call disp%skip()
        call disp%show("Weight")
        call disp%show( Weight )
        call disp%show("WeightedArray_RK")
        call disp%show( WeightedArray_RK )
        call disp%show("getVerbose(WeightedArray_RK, Weight, sum(Weight, mask = Weight > 0_IK), dim = 2_IK)")
        call disp%show( getVerbose(WeightedArray_RK, Weight, sum(Weight, mask = Weight > 0_IK), dim = 2_IK) )

    end block

end program example