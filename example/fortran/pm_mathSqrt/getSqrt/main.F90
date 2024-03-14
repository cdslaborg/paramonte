program example

    use pm_kind, only: IKC => IK32 ! any integer kind is supported.
    use pm_kind, only: SK, IK, LK, RKD
    use pm_io, only: display_type
    use pm_distLogUnif, only: getLogUnifRand
    use pm_mathSqrt, only: getSqrt, linear, binary
    use pm_err, only: setAsserted

    implicit none

    integer(IK) :: i
    integer(IKC) :: posint, intSqrt

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    do i = 1, 10
        call disp%skip()
        call disp%show("posint = getLogUnifRand(1, huge(1))")
                        posint = getLogUnifRand(1, huge(1))
        call disp%show("posint")
        call disp%show( posint )
        call disp%show("intSqrt = getSqrt(posint, binary)")
                        intSqrt = getSqrt(posint, binary)
        call disp%show("[intSqrt, floor(sqrt(real(posint, RKD)))]")
        call disp%show( [intSqrt, floor(sqrt(real(posint, RKD)))] )
        call disp%show("call setAsserted(intSqrt == floor(sqrt(real(posint, RKD))))")
                        call setAsserted(intSqrt == floor(sqrt(real(posint, RKD))))
        call disp%skip()
    end do

    do posint = huge(0_IKC), huge(0_IKC) - 10_IKC, -1_IKC
        call disp%skip()
        call disp%show("posint")
        call disp%show( posint )
        call disp%show("intSqrt = getSqrt(posint, binary)")
                        intSqrt = getSqrt(posint, binary)
        call disp%show("[intSqrt, floor(sqrt(real(posint, RKD)))]")
        call disp%show( [intSqrt, floor(sqrt(real(posint, RKD)))] )
        call disp%show("call setAsserted(intSqrt == floor(sqrt(real(posint, RKD))))")
                        call setAsserted(intSqrt == floor(sqrt(real(posint, RKD))))
        call disp%skip()
    end do

    do posint = 0_IKC, 10_IKC
        call disp%skip()
        call disp%show("posint")
        call disp%show( posint )
        call disp%show("intSqrt = getSqrt(posint, binary)")
                        intSqrt = getSqrt(posint, binary)
        call disp%show("[intSqrt, getSqrt(posint, linear), floor(sqrt(real(posint, RKD)))]")
        call disp%show( [intSqrt, getSqrt(posint, linear), floor(sqrt(real(posint, RKD)))] )
        call disp%show("call setAsserted(intSqrt == floor(sqrt(real(posint, RKD))))")
                        call setAsserted(intSqrt == floor(sqrt(real(posint, RKD))))
        call disp%skip()
    end do

end program example