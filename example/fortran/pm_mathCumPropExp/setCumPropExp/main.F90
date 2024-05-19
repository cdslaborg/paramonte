program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_mathCumPropExp, only: setCumPropExp
    use pm_mathCumPropExp, only: forward, backward, reverse, nothing, sequence, selection
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: setUnifRand

    implicit none


    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: RKG => RKH ! all real kinds are supported.
        real, allocatable :: array(:), cumPropExp(:)
        call disp%skip
        call disp%show("array = log([1., 2., 3., 4.])")
                        array = log([1., 2., 3., 4.])
        call disp%show("array")
        call disp%show( array )
        call disp%show("call setResized(cumPropExp, size(array, 1, IK))")
                        call setResized(cumPropExp, size(array, 1, IK))
        call disp%show("call setCumPropExp(cumPropExp, array, maxval(array), sequence)")
                        call setCumPropExp(cumPropExp, array, maxval(array), sequence)
        call disp%show("cumPropExp")
        call disp%show( cumPropExp )
        call disp%show("call setCumPropExp(cumPropExp, array, maxval(array), selection)")
                        call setCumPropExp(cumPropExp, array, maxval(array), selection)
        call disp%show("cumPropExp")
        call disp%show( cumPropExp )
        call disp%skip
        call disp%show("array = array(size(array):1:-1)")
                        array = array(size(array):1:-1)
        call disp%show("array")
        call disp%show( array )
        call disp%show("call setCumPropExp(cumPropExp, array, maxval(array), sequence, backward, nothing)")
                        call setCumPropExp(cumPropExp, array, maxval(array), sequence, backward, nothing)
        call disp%show("cumPropExp")
        call disp%show( cumPropExp )
        call disp%skip
        call disp%show("call setCumPropExp(cumPropExp, array, maxval(array), sequence, backward, reverse)")
                        call setCumPropExp(cumPropExp, array, maxval(array), sequence, backward, reverse)
        call disp%show("cumPropExp")
        call disp%show( cumPropExp )
        call disp%skip
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute and return the cumulative sum in-place, within the input `Array`.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKG => RKH ! all real kinds are supported.
        real, allocatable :: array(:)
        call disp%skip
        call disp%show("array = log([1., 2., 3., 4.])")
                        array = log([1., 2., 3., 4.])
        call disp%show("array")
        call disp%show( array )
        call disp%show("call setCumPropExp(array, maxval(array), sequence)")
                        call setCumPropExp(array, maxval(array), sequence)
        call disp%show("array")
        call disp%show( array )
        call disp%show("call setCumPropExp(array, maxval(array), selection)")
                        call setCumPropExp(array, maxval(array), selection)
        call disp%show("array")
        call disp%show( array )
        call disp%skip
        call disp%show("array = array(size(array):1:-1)")
                        array = array(size(array):1:-1)
        call disp%show("array")
        call disp%show( array )
        call disp%show("call setCumPropExp(array, maxval(array), sequence, backward, nothing)")
                        call setCumPropExp(array, maxval(array), sequence, backward, nothing)
        call disp%show("array")
        call disp%show( array )
        call disp%skip
        call disp%show("call setCumPropExp(array, maxval(array), sequence, backward, reverse)")
                        call setCumPropExp(array, maxval(array), sequence, backward, reverse)
        call disp%show("array")
        call disp%show( array )
        call disp%skip
    end block

end program example