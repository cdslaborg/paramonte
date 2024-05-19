program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distKolm, only: setKolmRand
    use pm_distUnif, only: getUnifRand
    use pm_arrayResize, only: setResized

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: rand(:)
        call disp%skip()
        call disp%show("call setResized(rand, 3_IK)")
                        call setResized(rand, 3_IK)
        call disp%show("call setKolmRand(rand(1), getUnifRand(0._TKG, 1._TKG))")
                        call setKolmRand(rand(1), getUnifRand(0._TKG, 1._TKG))
        call disp%show("rand(1)")
        call disp%show( rand(1) )
        call disp%show("call setKolmRand(rand, getUnifRand(0._TKG, 1._TKG, 3_IK))")
                        call setKolmRand(rand, getUnifRand(0._TKG, 1._TKG, 3_IK))
        call disp%show("rand")
        call disp%show( rand )
        call disp%skip()
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_io, only: getErrTableWrite
        use pm_kind, only: TKG => RKH
        real(TKG) :: rand(1500)
        call setKolmRand(rand, getUnifRand(0._TKG, 1._TKG, size(rand, 1, IK)))
        if (0 /= getErrTableWrite(SK_"setKolmRand.RK.txt", rand, header = SK_"rand")) error stop "table output failed."
    end block

end program example