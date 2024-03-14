program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_mathErf, only: getErfInv
    use pm_io, only: display_type

    implicit none

    real(RKH), allocatable :: erfinv(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("erfinv = [real(RKH) :: getErfInv(0.), getErfInv(0.d0)]")
                    erfinv = [real(RKH) :: getErfInv(0.), getErfInv(0.d0)]
    call disp%show("erfinv")
    call disp%show( erfinv )
    call disp%show("erf(erfinv)")
    call disp%show( erf(erfinv) )
    call disp%skip()

    call disp%skip()
    call disp%show("erfinv = [getErfInv(-1. + epsilon(0.))]")
                    erfinv = [getErfInv(-1. + epsilon(0.))]
    call disp%show("erfinv")
    call disp%show( erfinv )
    call disp%show("erf(erfinv)")
    call disp%show( erf(erfinv) )
    call disp%skip()

    call disp%skip()
    call disp%show("erfinv = [getErfInv(+1. - epsilon(0.))]")
                    erfinv = [getErfInv(+1. - epsilon(0.))]
    call disp%show("erfinv")
    call disp%show( erfinv )
    call disp%show("erf(erfinv)")
    call disp%show( erf(erfinv) )
    call disp%skip()

    call disp%skip()
    call disp%show("erfinv = getErfInv([real :: -.99, -.75, -0.5, -0.1, 0., .1, .5, .75, .99])")
                    erfinv = getErfInv([real :: -.99, -.75, -0.5, -0.1, 0., .1, .5, .75, .99])
    call disp%show("erfinv")
    call disp%show( erfinv )
    call disp%show("erf(erfinv)")
    call disp%show( erf(erfinv) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array of the regularized Incomplete Beta function for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_kind, only: RKB, RKC => RKH, RKL => RK32, RK32, RK64
        use pm_arraySpace, only: setLinSpace
        integer(IK) , parameter :: NP = 2000
        integer :: fileUnit, i
        real(RKC) :: erfval(NP)
        call setLinSpace(erfval, -1._RKC + epsilon(0._RKL), 1._RKC - epsilon(0._RKL))
        open(newunit = fileUnit, file = "getErfInv.RK.txt")
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))") real(erfval(i), RK64), getErfInv(real(erfval(i), RK64))
        end do
        close(fileUnit)
        open(newunit = fileUnit, file = "getErfInv.RK32.abserr.txt")
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))") real(erfval(i), RK32), real(erfval(i), RK32) - erf(real(getErfInv(real(erfval(i), RK32)), RKB))
        end do
        close(fileUnit)
        open(newunit = fileUnit, file = "getErfInv.RK64.abserr.txt")
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))") real(erfval(i), RK64), real(erfval(i), RK64) - erf(real(getErfInv(real(erfval(i), RK64)), RKB))
        end do
        close(fileUnit)
        open(newunit = fileUnit, file = "getErfInv.RK128.abserr.txt")
        do i = 1, NP
            write(fileUnit, "(*(g0,:,' '))") real(erfval(i), RKH), real(erfval(i), RKH) - erf(real(getErfInv(real(erfval(i), RKH)), RKB))
        end do
        close(fileUnit)
    end block

end program example