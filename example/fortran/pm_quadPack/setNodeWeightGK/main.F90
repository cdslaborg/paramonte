program example

    use pm_kind, only: SK, IK, RK => RKH ! testing with highest real precision available. all other real kinds are also supported.
    use pm_quadpack, only: nodeK15, weightK15, weightG7
    use pm_quadpack, only: nodeK21, weightK21, weightG10
    use pm_quadpack, only: nodeK31, weightK31, weightG15
    use pm_quadpack, only: nodeK41, weightK41, weightG20
    use pm_quadpack, only: nodeK51, weightK51, weightG25
    use pm_quadpack, only: nodeK61, weightK61, weightG30
    use pm_quadpack, only: setNodeWeightGK
    use pm_io, only: display_type

    implicit none

    integer(IK) , parameter     :: MAX_NPG = 30_IK  !   MAXimum Number of Points in the Gauss quadrature rule considered in this example.
    real(RK)                    :: nodeK(MAX_NPG+1), weightK(MAX_NPG+1), weightG(MAX_NPG+1)
    real(RK)    , allocatable   :: RefNodeK(:), RefWeightK(:), RefWeightG(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 7-15 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(RefWeightG, source = weightG7)
    allocate(RefWeightK, source = weightK15)
    allocate(RefNodeK, source = nodeK15)
    call setNodeWeight(npg = 7)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 10-21 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(RefWeightG, source = weightG10)
    allocate(RefWeightK, source = weightK21)
    allocate(RefNodeK, source = nodeK21)
    call setNodeWeight(npg = 10)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 15-31 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(RefWeightG, source = weightG15)
    allocate(RefWeightK, source = weightK31)
    allocate(RefNodeK, source = nodeK31)
    call setNodeWeight(npg = 15)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 20-41 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(RefWeightG, source = weightG20)
    allocate(RefWeightK, source = weightK41)
    allocate(RefNodeK, source = nodeK41)
    call setNodeWeight(npg = 20)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 25-51 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(RefWeightG, source = weightG25)
    allocate(RefWeightK, source = weightK51)
    allocate(RefNodeK, source = nodeK51)
    call setNodeWeight(npg = 25)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 30-61 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    allocate(RefWeightG, source = weightG30)
    allocate(RefWeightK, source = weightK61)
    allocate(RefNodeK, source = nodeK61)
    call setNodeWeight(npg = 30)

contains

    subroutine setNodeWeight(npg)
        use pm_val2str, only: getStr
        integer(IK) , intent(in)    :: npg  !   number of points in the Gauss quadrature rule.
        integer(IK)                 :: npk  !   number of extra points by the Kronrod extension rule ( = # points in the Gauss rule + 1).
        integer(IK)                 :: i,j
        call disp%skip()
        call disp%show("npk = "//getStr(npg)//" + 1")
                        npk = npg + 1_IK
        call disp%show("call setNodeWeightGK(nodeK(1:npk), weightK(1:npk), weightG(1:npk/2))")
                        call setNodeWeightGK(nodeK(1:npk), weightK(1:npk), weightG(1:npk/2))
        call disp%show("[nodeK, weightK, weightG]")
        do i = 1, npk
            if (mod(i, 2_IK) == 0_IK) then
                call disp%show( [nodeK(i), weightK(i), weightG(i/2)] )
            else
                call disp%show( [nodeK(i), weightK(i)] )
            end if
        end do
        call disp%skip()
        call disp%show("[nodeK - exact, weightK - exact, weightG - exact]")
        do i = 1, npk
            j = npk - i + 1
            if (mod(i, 2_IK) == 0_IK) then
                call disp%show( [nodeK(i) - RefNodeK(j), weightK(i) - RefWeightK(j), weightG(i/2) - RefWeightG(size(RefWeightG) - i/2 + 1)] )
            else
                call disp%show( [nodeK(i) - RefNodeK(j), weightK(i) - RefWeightK(j)] )
            end if
        end do
        call disp%skip()
        deallocate(RefNodeK, RefWeightG, RefWeightK)
    end subroutine

end program example