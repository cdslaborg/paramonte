program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_mathCumPropExp, only: getCumPropExp
    use pm_mathCumPropExp, only: forward, backward, reverse, nothing, sequence, selection
    use pm_distUnif, only: setUnifRand

    implicit none


    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: RKC => RKH ! all real kinds are supported.
        real, allocatable :: array(:), cumPropExp(:)
        call disp%skip
        call disp%show("array = log([1., 2., 3., 4.])")
                        array = log([1., 2., 3., 4.])
        call disp%show("array")
        call disp%show( array )
        call disp%show("cumPropExp = getCumPropExp(array, maxval(array))")
                        cumPropExp = getCumPropExp(array, maxval(array))
        call disp%show("cumPropExp")
        call disp%show( cumPropExp )
        call disp%show("cumPropExp = getCumPropExp(array, maxval(array), sequence)")
                        cumPropExp = getCumPropExp(array, maxval(array), sequence)
        call disp%show("cumPropExp")
        call disp%show( cumPropExp )
        call disp%show("cumPropExp = getCumPropExp(array, maxval(array), selection)")
                        cumPropExp = getCumPropExp(array, maxval(array), selection)
        call disp%show("cumPropExp")
        call disp%show( cumPropExp )
        call disp%skip
        call disp%show("array = array(size(array):1:-1)")
                        array = array(size(array):1:-1)
        call disp%show("array")
        call disp%show( array )
        call disp%show("cumPropExp = getCumPropExp(array, maxval(array), direction = backward)")
                        cumPropExp = getCumPropExp(array, maxval(array), direction = backward)
        call disp%show("cumPropExp")
        call disp%show( cumPropExp )
        call disp%skip
        call disp%show("cumPropExp = getCumPropExp(array, maxval(array), direction = backward, action = reverse)")
                        cumPropExp = getCumPropExp(array, maxval(array), direction = backward, action = reverse)
        call disp%show("cumPropExp")
        call disp%show( cumPropExp )
        call disp%show("cumPropExp = getCumPropExp(array, maxval(array), action = reverse)")
                        cumPropExp = getCumPropExp(array, maxval(array), action = reverse)
        call disp%show("cumPropExp")
        call disp%show( cumPropExp )
        call disp%skip
    end block

end program example