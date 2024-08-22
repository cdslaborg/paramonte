program example

    use pm_kind, only: RKG => RKS ! any real kind is supported.
    use pm_kind, only: IKG => IKS ! any integer kind is supported.
    use pm_kind, only: SK, IK, LK, RKD
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type
    use pm_mathRound, only: pnint
    use pm_err, only: setAsserted

    implicit none

    integer(IK) :: i, ival, irep
    real(RKG), allocatable :: val(:)
    integer(IKG), allocatable :: whole(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    do i = 0, 1
        call disp%skip()
        call disp%show("i")
        call disp%show( i )
        call disp%show("val = (-1)**i * [5., 5.2, 5.5, 5.8, 6.]")
                        val = (-1)**i * [5., 5.2, 5.5, 5.8, 6.]
        do ival = 1, size(val)
            call disp%skip()
            call disp%show("ival")
            call disp%show( ival )
            call disp%show("val(ival)")
            call disp%show( val(ival) )
            call disp%show("whole = pnint([(val(ival), irep = 1, 20)])")
                            whole = pnint([(val(ival), irep = 1, 20)])
            call disp%show("whole")
            call disp%show( whole )
            call disp%show("call setAsserted(all(abs(val(ival) - whole) < 1))")
                            call setAsserted(all(abs(val(ival) - whole) < 1))
            call disp%skip()
        end do
    end do

end program example