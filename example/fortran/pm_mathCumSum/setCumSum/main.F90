program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_arrayResize, only: setResized
    use pm_mathCumSum, only: setCumSum, forward, backward, nothing, reverse

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: IKG => IK
        integer(IKG), allocatable :: array(:), cumsum(:)
        call disp%skip
        call disp%show("array = [1, 2, 3, 4]")
                        array = [1, 2, 3, 4]
        call disp%show("call setResized(cumsum, size(array, 1, IK))")
                        call setResized(cumsum, size(array, 1, IK))
        call disp%show("call setCumSum(cumsum, array)")
                        call setCumSum(cumsum, array)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, forward, nothing)")
                        call setCumSum(cumsum, array, forward, nothing)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, forward, reverse)")
                        call setCumSum(cumsum, array, forward, reverse)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, backward, nothing)")
                        call setCumSum(cumsum, array, backward, nothing)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, backward, reverse)")
                        call setCumSum(cumsum, array, backward, reverse)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%skip
    end block

    block
        use pm_kind, only: IKG => IK
        integer(IKG), allocatable :: array(:), cumsum(:), reference(:)
        call disp%skip
        call disp%show("reference = [1, 2, 3, 4]")
                        reference = [1, 2, 3, 4]
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array)")
                        call setCumSum(array)
        call disp%show("array")
        call disp%show( array )
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, forward, nothing)")
                        call setCumSum(array, forward, nothing)
        call disp%show("array")
        call disp%show( array )
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, forward, reverse)")
                        call setCumSum(array, forward, reverse)
        call disp%show("array")        
        call disp%show( array )        
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, backward, nothing)")
                        call setCumSum(array, backward, nothing)
        call disp%show("array")        
        call disp%show( array )        
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, backward, reverse)")
                        call setCumSum(array, backward, reverse)
        call disp%show("array")
        call disp%show( array )
        call disp%skip
    end block

    block
        use pm_kind, only: CKG => CKS
        complex(CKG), allocatable :: array(:), cumsum(:)
        call disp%skip
        call disp%show("array = cmplx([real(CKG) :: 1, 2, 3, 4], -[real(CKG) :: 1, 2, 3, 4], CKG)")
                        array = cmplx([real(CKG) :: 1, 2, 3, 4], -[real(CKG) :: 1, 2, 3, 4], CKG)
        call disp%show("call setResized(cumsum, size(array, 1, IK))")
                        call setResized(cumsum, size(array, 1, IK))
        call disp%show("call setCumSum(cumsum, array)")
                        call setCumSum(cumsum, array)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, forward, nothing)")
                        call setCumSum(cumsum, array, forward, nothing)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, forward, reverse)")
                        call setCumSum(cumsum, array, forward, reverse)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, backward, nothing)")
                        call setCumSum(cumsum, array, backward, nothing)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, backward, reverse)")
                        call setCumSum(cumsum, array, backward, reverse)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%skip
    end block

    block
        use pm_kind, only: CKG => CKS
        complex(CKG), allocatable :: array(:), cumsum(:), reference(:)
        call disp%skip
        call disp%show("reference = cmplx([real(CKG) :: 1, 2, 3, 4], -[real(CKG) :: 1, 2, 3, 4], CKG)")
                        reference = cmplx([real(CKG) :: 1, 2, 3, 4], -[real(CKG) :: 1, 2, 3, 4], CKG)
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array)")
                        call setCumSum(array)
        call disp%show("array")
        call disp%show( array )
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, forward, nothing)")
                        call setCumSum(array, forward, nothing)
        call disp%show("array")
        call disp%show( array )
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, forward, reverse)")
                        call setCumSum(array, forward, reverse)
        call disp%show("array")        
        call disp%show( array )        
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, backward, nothing)")
                        call setCumSum(array, backward, nothing)
        call disp%show("array")        
        call disp%show( array )        
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, backward, reverse)")
                        call setCumSum(array, backward, reverse)
        call disp%show("array")
        call disp%show( array )
        call disp%skip
    end block

    block
        use pm_kind, only: RKG => RKH
        real(RKG), allocatable :: array(:), cumsum(:)
        call disp%skip
        call disp%show("array = [real(RKG) :: 1, 2, 3, 4]")
                        array = [real(RKG) :: 1, 2, 3, 4]
        call disp%show("call setResized(cumsum, size(array, 1, IK))")
                        call setResized(cumsum, size(array, 1, IK))
        call disp%show("call setCumSum(cumsum, array)")
                        call setCumSum(cumsum, array)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, forward, nothing)")
                        call setCumSum(cumsum, array, forward, nothing)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, forward, reverse)")
                        call setCumSum(cumsum, array, forward, reverse)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, backward, nothing)")
                        call setCumSum(cumsum, array, backward, nothing)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%show("call setCumSum(cumsum, array, backward, reverse)")
                        call setCumSum(cumsum, array, backward, reverse)
        call disp%show("cumsum")
        call disp%show( cumsum )
        call disp%skip
    end block

    block
        use pm_kind, only: RKG => RKH
        real(RKG), allocatable :: array(:), cumsum(:), reference(:)
        call disp%skip
        call disp%show("reference = [1, 2, 3, 4]")
                        reference = [1, 2, 3, 4]
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array)")
                        call setCumSum(array)
        call disp%show("array")
        call disp%show( array )
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, forward, nothing)")
                        call setCumSum(array, forward, nothing)
        call disp%show("array")
        call disp%show( array )
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, forward, reverse)")
                        call setCumSum(array, forward, reverse)
        call disp%show("array")        
        call disp%show( array )        
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, backward, nothing)")
                        call setCumSum(array, backward, nothing)
        call disp%show("array")        
        call disp%show( array )        
        call disp%show("array = reference")
                        array = reference
        call disp%show("call setCumSum(array, backward, reverse)")
                        call setCumSum(array, backward, reverse)
        call disp%show("array")
        call disp%show( array )
        call disp%skip
    end block

end program example