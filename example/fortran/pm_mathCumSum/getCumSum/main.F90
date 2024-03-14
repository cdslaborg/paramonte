program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_mathCumSum, only: getCumSum, forward, backward, nothing, reverse

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: IKC => IK
        integer(IKC), allocatable :: array(:)
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
        use pm_kind, only: CKC => CK32
        complex(CKC), allocatable :: array(:)
        call disp%skip
        call disp%show("array = cmplx([real(CKC) :: 1, 2, 3, 4], -[real(CKC) :: 1, 2, 3, 4], CKC)")
                        array = cmplx([real(CKC) :: 1, 2, 3, 4], -[real(CKC) :: 1, 2, 3, 4], CKC)
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
        use pm_kind, only: RKC => RK32
        real(RKC), allocatable :: array(:)
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
        use pm_kind, only: RKC => RKH
        real(RKC), allocatable :: array(:)
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