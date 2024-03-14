program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_matrixPack, only: isMatPack, rfpack
    use pm_err, only: setAsserted
    use pm_val2str, only: getStr

    implicit none

    integer(IK) :: i, shape(2)
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    do i = 1, 30
        call disp%skip
        call disp%show("shape = getUnifRand(1_IK, 7_IK, 2_IK)")
                        shape = getUnifRand(1_IK, 7_IK, 2_IK)
        call disp%show("shape")
        call disp%show( shape )
        call disp%show("isMatPack(rfpack, shape)")
        call disp%show( isMatPack(rfpack, shape) )
    end do

end program example