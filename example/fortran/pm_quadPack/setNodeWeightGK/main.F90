program example

    use pm_kind, only: SK, IK, RKG => RKH ! testing with highest real precision available. all other real kinds are also supported.
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
    real(RKG)                   :: nodeK(MAX_NPG+1), weightK(MAX_NPG+1), weightG(MAX_NPG+1)
    real(RKG)   , allocatable   :: refNodeK(:), refWeightK(:), refWeightG(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 7-15 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    refWeightG = real(weightG7, RKG)
    refWeightK = real(weightK15, RKG)
    refNodeK = real(nodeK15, RKG)
    call setNodeWeight(npg = 7)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 10-21 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    refWeightG = real(weightG10, RKG)
    refWeightK = real(weightK21, RKG)
    refNodeK = real(nodeK21, RKG)
    call setNodeWeight(npg = 10)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 15-31 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    refWeightG = real(weightG15, RKG)
    refWeightK = real(weightK31, RKG)
    refNodeK = real(nodeK31, RKG)
    call setNodeWeight(npg = 15)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 20-41 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    refWeightG = real(weightG20, RKG)
    refWeightK = real(weightK41, RKG)
    refNodeK = real(nodeK41, RKG)
    call setNodeWeight(npg = 20)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 25-51 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    refWeightG = real(weightG25, RKG)
    refWeightK = real(weightK51, RKG)
    refNodeK = real(nodeK51, RKG)
    call setNodeWeight(npg = 25)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the 30-61 Gauss-Kronrod nodes and weights.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    refWeightG = real(weightG30, RKG)
    refWeightK = real(weightK61, RKG)
    refNodeK = real(nodeK61, RKG)
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
                call disp%show( [nodeK(i) - refNodeK(j), weightK(i) - refWeightK(j), weightG(i/2) - refWeightG(size(refWeightG) - i/2 + 1)] )
            else
                call disp%show( [nodeK(i) - refNodeK(j), weightK(i) - refWeightK(j)] )
            end if
        end do
        call disp%skip()
        deallocate(refNodeK, refWeightG, refWeightK)
    end subroutine

end program example