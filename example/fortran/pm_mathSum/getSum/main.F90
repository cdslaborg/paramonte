program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_distUnif, only: setUnifRand
    use pm_arrayRank, only: getRankDense
    use pm_arraySpace, only: setLinSpace
    use pm_mathSum, only: getSum, iteration, recursion, kahanbabu, fablocked, nablocked
    use pm_io, only: display_type

    implicit none

    real(RKH) :: truth
    real(RKH), allocatable :: sumres(:), relerr(:)
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKG => RKS
        integer(IK), parameter :: lenx = 10**7
        real(TKG) :: x(lenx), lb, ub
        call disp%skip()
        call disp%show("lenx")
        call disp%show( lenx )
        call disp%show("lb = real(lenx, TKG); ub = 1._TKG")
                        lb = real(lenx, TKG); ub = 1._TKG
        call disp%show("call setLinSpace(x, lb, ub) ! call setUnifRand(x)")
                        call setLinSpace(x, lb, ub) ! call setUnifRand(x)
        call disp%show("truth = (real(ub, TKG) + real(lb, RKH)) * size(x, 1, IK) / 2 ! reference high-precision value for comparison")
                        truth = (real(ub, TKG) + real(lb, RKH)) * size(x, 1, IK) / 2 ! reference high-precision value for comparison
        call disp%show("truth")
        call disp%show( truth )
        call disp%show("sumres = [getSum(x), getSum(x, iteration), getSum(x, recursion), getSum(x, kahanbabu), getSum(x, fablocked), getSum(x, nablocked)]")
                        sumres = [getSum(x), getSum(x, iteration), getSum(x, recursion), getSum(x, kahanbabu), getSum(x, fablocked), getSum(x, nablocked)]
        call disp%show("sumres")
        call disp%show( sumres )
        call disp%show("relerr = abs(truth - sumres) / truth")
                        relerr = abs(truth - sumres) / truth
        call disp%show("relerr")
        call disp%show( relerr )
        call disp%show("getRankDense(relerr)")
        call disp%show( getRankDense(relerr) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => RKD
        integer(IK), parameter :: lenx = 10**8
        real(TKG) :: x(lenx), lb, ub
        call disp%skip()
        call disp%show("lenx")
        call disp%show( lenx )
        call disp%show("lb = real(lenx, TKG); ub = 1._TKG")
                        lb = real(lenx, TKG); ub = 1._TKG
        call disp%show("call setLinSpace(x, lb, ub) ! call setUnifRand(x)")
                        call setLinSpace(x, lb, ub) ! call setUnifRand(x)
        call disp%show("truth = (real(ub, TKG) + real(lb, RKH)) * size(x, 1, IK) / 2 ! reference high-precision value for comparison")
                        truth = (real(ub, TKG) + real(lb, RKH)) * size(x, 1, IK) / 2 ! reference high-precision value for comparison
        call disp%show("truth")
        call disp%show( truth )
        call disp%show("sumres = [getSum(x), getSum(x, iteration), getSum(x, recursion), getSum(x, kahanbabu), getSum(x, fablocked), getSum(x, nablocked)]")
                        sumres = [getSum(x), getSum(x, iteration), getSum(x, recursion), getSum(x, kahanbabu), getSum(x, fablocked), getSum(x, nablocked)]
        call disp%show("sumres")
        call disp%show( sumres )
        call disp%show("relerr = abs(truth - sumres) / truth")
                        relerr = abs(truth - sumres) / truth
        call disp%show("relerr")
        call disp%show( relerr )
        call disp%show("getRankDense(relerr)")
        call disp%show( getRankDense(relerr) )
        call disp%skip()
    end block

    block
        call disp%skip()
        call disp%show("[getSum([real ::]), getSum([real :: 1]), getSum([real :: 1, 1]), getSum([real :: 1, 1, 1])]")
        call disp%show( [getSum([real ::]), getSum([real :: 1]), getSum([real :: 1, 1]), getSum([real :: 1, 1, 1])] )
        call disp%skip()
    end block

end program example