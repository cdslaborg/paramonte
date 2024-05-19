program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_mathSubAdd, only: operator(.subadd.)
    use pm_mathErf, only: setErfInv
    use pm_arrayResize, only: setResized
    use pm_io, only: display_type

    implicit none

    real(RKH), allocatable :: erfinv(:), x(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("x = [real(RKH) :: 0., 0.5]")
                    x = [real(RKH) :: 0., 0.5]
    call disp%show("x")
    call disp%show( x )
    call disp%show("call setResized(erfinv, size(x, 1, IK))")
                    call setResized(erfinv, size(x, 1, IK))
    call disp%show("call setErfInv(erfinv, x, abserr = epsilon(x))")
                    call setErfInv(erfinv, x, abserr = epsilon(x))
    call disp%show("erfinv")
    call disp%show( erfinv )
    call disp%show("erf(erfinv)")
    call disp%show( erf(erfinv) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = .subadd. 1. - .subadd. epsilon(0.)")
                    x = .subadd. 1. - .subadd. epsilon(0.)
    call disp%show("x")
    call disp%show( x )
    call disp%show("call setResized(erfinv, size(x, 1, IK))")
                    call setResized(erfinv, size(x, 1, IK))
    call disp%show("call setErfInv(erfinv, x, abserr = epsilon(x))")
                    call setErfInv(erfinv, x, abserr = epsilon(x))
    call disp%show("erfinv")
    call disp%show( erfinv )
    call disp%show("erf(erfinv)")
    call disp%show( erf(erfinv) )
    call disp%skip()

    call disp%skip()
    call disp%show("x = [-.99, -.75, -0.5, -0.1, 0., .1, .5, .75, .99]")
                    x = [-.99, -.75, -0.5, -0.1, 0., .1, .5, .75, .99]
    call disp%show("x")
    call disp%show( x )
    call disp%show("call setResized(erfinv, size(x, 1, IK))")
                    call setResized(erfinv, size(x, 1, IK))
    call disp%show("call setErfInv(erfinv, x, abserr = epsilon(x))")
                    call setErfInv(erfinv, x, abserr = epsilon(x))
    call disp%show("erfinv")
    call disp%show( erfinv )
    call disp%show("erf(erfinv)")
    call disp%show( erf(erfinv) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Incomplete Beta function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_kind, only: RKB, RKG => RKH
        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 2000
        integer :: fileUnit, i
        real(RKG) :: erfval(NP), erfinv(NP)

        call setLinSpace(erfval, -1._RKG + sqrt(epsilon(0._RKG)), 1._RKG - sqrt(epsilon(0._RKG)))

        open(newunit = fileUnit, file = "setErfInv.RK.txt")
        call setErfInv(erfinv, erfval, abserr = 100 * epsilon(0._RKG))
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))") real(erfval(i), RKG), real(erfinv(i), RKG)
        end do
        close(fileUnit)

        block
            use pm_kind, only: RKG => RKS
            open(newunit = fileUnit, file = "setErfInv.RKS.abserr.txt")
            call setErfInv(erfinv, erfval, abserr = real(100 * epsilon(0._RKG), RKH))
            do i = 1, NP
                write(fileUnit, "(*(g0,:,' '))") real(erfval(i), RKG), real(erfval(i), RKG) - erf(real(erfinv(i), RKB))
            end do
            close(fileUnit)
        end block

        block
            use pm_kind, only: RKG => RKD
            open(newunit = fileUnit, file = "setErfInv.RKD.abserr.txt")
            call setErfInv(erfinv, erfval, abserr = real(100 * epsilon(0._RKG), RKH))
            do i = 1, NP
                write(fileUnit, "(*(g0,:,' '))") real(erfval(i), RKG), real(erfval(i), RKG) - erf(real(erfinv(i), RKB))
            end do
            close(fileUnit)
        end block

        block
            use pm_kind, only: RKG => RKH
            open(newunit = fileUnit, file = "setErfInv.RKH.abserr.txt")
            call setErfInv(erfinv, erfval, abserr = real(100 * epsilon(0._RKG), RKH))
            do i = 1, NP
                write(fileUnit, "(*(g0,:,' '))") real(erfval(i), RKG), real(erfval(i), RKG) - erf(real(erfinv(i), RKB))
            end do
            close(fileUnit)
        end block

    end block

end program example