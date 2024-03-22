program example

    use pm_kind, only: IKS, IKD
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_mathFactoring, only: getFactoring
    use pm_distLogUnif, only: getLogUnifRand

    implicit none

    integer(IK) :: i, n
    integer(IK), allocatable :: Factoring(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the factoring of a scalar or array of integers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    do i = 1, 10
        call disp%skip()
        call disp%show("n = getLogUnifRand(2, huge(1))")
                        n = getLogUnifRand(2, huge(1))
        call disp%show("n")
        call disp%show( n )
        call disp%show("Factoring = getFactoring(n)")
                        Factoring = getFactoring(n)
        call disp%show("Factoring")
        call disp%show( Factoring )
        call disp%show("n - product(Factoring) ! By definition, it must be zero.")
        call disp%show( n - product(Factoring) )
        call disp%skip()
    end do

end program example