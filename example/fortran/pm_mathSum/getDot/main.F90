program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_distUnif, only: setUnifRand
    use pm_arrayRank, only: getRankDense
    use pm_arraySpace, only: setLinSpace
    use pm_mathSum, only: getDot, iteration, recursion, kahanbabu, fablocked, nablocked
    use pm_io, only: display_type

    implicit none

    real(RKH) :: truth
    real(RKH), allocatable :: dotres(:), relerr(:)
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKG => RKS
        integer(IK), parameter :: lenx = 10**7
        real(TKG) :: x(lenx), y(lenx), lb, ub
        call disp%skip()
        call disp%show("lenx")
        call disp%show( lenx )
        call disp%show("lb = real(lenx, TKG); ub = 1._TKG")
                        lb = real(lenx, TKG); ub = 1._TKG
        call disp%show("call setLinSpace(x, lb, ub) ! call setUnifRand(x)")
                        call setLinSpace(x, lb, ub) ! call setUnifRand(x)
        call disp%show("y = 1 !/ x")
                        y = 1 !/ x
        call disp%show("truth = (real(ub, TKG) + real(lb, RKH)) * size(x, 1, IK) / 2 ! reference high-precision value for comparison")
                        truth = (real(ub, TKG) + real(lb, RKH)) * size(x, 1, IK) / 2 ! reference high-precision value for comparison
        call disp%show("truth")
        call disp%show( truth )
        call disp%show("dotres = [getDot(x, y), getDot(x, y, iteration), getDot(x, y, recursion), getDot(x, y, kahanbabu), getDot(x, y, fablocked), getDot(x, y, nablocked)]")
                        dotres = [getDot(x, y), getDot(x, y, iteration), getDot(x, y, recursion), getDot(x, y, kahanbabu), getDot(x, y, fablocked), getDot(x, y, nablocked)]
        call disp%show("dotres")
        call disp%show( dotres )
        call disp%show("relerr = abs(truth - dotres) / truth")
                        relerr = abs(truth - dotres) / truth
        call disp%show("relerr")
        call disp%show( relerr )
        call disp%show("getRankDense(relerr)")
        call disp%show( getRankDense(relerr) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKD
        integer(IK), parameter :: lenx = 10**8
        real(TKG) :: x(lenx), y(lenx), lb, ub
        call disp%skip()
        call disp%show("lenx")
        call disp%show( lenx )
        call disp%show("lb = real(lenx, TKG); ub = 1._TKG")
                        lb = real(lenx, TKG); ub = 1._TKG
        call disp%show("call setLinSpace(x, lb, ub) ! call setUnifRand(x)")
                        call setLinSpace(x, lb, ub) ! call setUnifRand(x)
        call disp%show("y = 1 !/ x")
                        y = 1 !/ x
        call disp%show("truth = (real(ub, TKG) + real(lb, RKH)) * size(x, 1, IK) / 2 ! reference high-precision value for comparison")
                        truth = (real(ub, TKG) + real(lb, RKH)) * size(x, 1, IK) / 2 ! reference high-precision value for comparison
        call disp%show("truth")
        call disp%show( truth )
        call disp%show("dotres = [getDot(x, y), getDot(x, y, iteration), getDot(x, y, recursion), getDot(x, y, kahanbabu), getDot(x, y, fablocked), getDot(x, y, nablocked)]")
                        dotres = [getDot(x, y), getDot(x, y, iteration), getDot(x, y, recursion), getDot(x, y, kahanbabu), getDot(x, y, fablocked), getDot(x, y, nablocked)]
        call disp%show("dotres")
        call disp%show( dotres )
        call disp%show("relerr = abs(truth - dotres) / truth")
                        relerr = abs(truth - dotres) / truth
        call disp%show("relerr")
        call disp%show( relerr )
        call disp%show("getRankDense(relerr)")
        call disp%show( getRankDense(relerr) )
        call disp%skip()
    end block

    block
        call disp%skip()
        call disp%show("[getDot([real ::], [real ::]), getDot([real :: 1], [real :: 1]), getDot([real :: 1, 1], [real :: 1, 1]), getDot([real :: 1, 1, 1], [real :: 1, 1, 1])]")
        call disp%show( [getDot([real ::], [real ::]), getDot([real :: 1], [real :: 1]), getDot([real :: 1, 1], [real :: 1, 1]), getDot([real :: 1, 1, 1], [real :: 1, 1, 1])] )
        call disp%skip()
    end block

end program example