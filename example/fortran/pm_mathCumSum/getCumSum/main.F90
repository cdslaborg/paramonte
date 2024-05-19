program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_mathCumSum, only: getCumSum, forward, backward, nothing, reverse

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: IKG => IK
        integer(IKG), allocatable :: array(:)
        call disp%skip
        call disp%show("array = [1, 2, 3, 4]")
                        array = [1, 2, 3, 4]
        call disp%show("getCumSum(array)")
        call disp%show( getCumSum(array) )
        call disp%show("getCumSum(array, action = reverse)")
        call disp%show( getCumSum(array, action = reverse) )
        call disp%show("getCumSum(array, direction = backward)")
        call disp%show( getCumSum(array, direction = backward) )
        call disp%show("getCumSum(array, direction = backward, action = reverse)")
        call disp%show( getCumSum(array, direction = backward, action = reverse) )
        call disp%skip
    end block

    block
        use pm_kind, only: CKG => CKS
        complex(CKG), allocatable :: array(:)
        call disp%skip
        call disp%show("array = cmplx([real(CKG) :: 1, 2, 3, 4], -[real(CKG) :: 1, 2, 3, 4], CKG)")
                        array = cmplx([real(CKG) :: 1, 2, 3, 4], -[real(CKG) :: 1, 2, 3, 4], CKG)
        call disp%show("getCumSum(array)")
        call disp%show( getCumSum(array) )
        call disp%show("getCumSum(array, action = reverse)")
        call disp%show( getCumSum(array, action = reverse) )
        call disp%show("getCumSum(array, direction = backward)")
        call disp%show( getCumSum(array, direction = backward) )
        call disp%show("getCumSum(array, direction = backward, action = reverse)")
        call disp%show( getCumSum(array, direction = backward, action = reverse) )
        call disp%skip
    end block

    block
        use pm_kind, only: RKG => RKS
        real(RKG), allocatable :: array(:)
        call disp%skip
        call disp%show("array = [1, 2, 3, 4]")
                        array = [1, 2, 3, 4]
        call disp%show("getCumSum(array)")
        call disp%show( getCumSum(array) )
        call disp%show("getCumSum(array, action = reverse)")
        call disp%show( getCumSum(array, action = reverse) )
        call disp%show("getCumSum(array, direction = backward)")
        call disp%show( getCumSum(array, direction = backward) )
        call disp%show("getCumSum(array, direction = backward, action = reverse)")
        call disp%show( getCumSum(array, direction = backward, action = reverse) )
        call disp%skip
    end block

    block
        use pm_kind, only: RKG => RKH
        real(RKG), allocatable :: array(:)
        call disp%skip
        call disp%show("array = [1, 2, 3, 4]")
                        array = [1, 2, 3, 4]
        call disp%show("getCumSum(array)")
        call disp%show( getCumSum(array) )
        call disp%show("getCumSum(array, action = reverse)")
        call disp%show( getCumSum(array, action = reverse) )
        call disp%show("getCumSum(array, direction = backward)")
        call disp%show( getCumSum(array, direction = backward) )
        call disp%show("getCumSum(array, direction = backward, action = reverse)")
        call disp%show( getCumSum(array, direction = backward, action = reverse) )
        call disp%skip
    end block

end program example