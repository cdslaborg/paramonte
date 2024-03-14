program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distKolm, only: getKolmRand
    use pm_distUnif, only: getUnifRand

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: rand(:)
        call disp%skip()
        call disp%show("rand = [getKolmRand(getUnifRand(0._TKC, 1._TKC))]")
                        rand = [getKolmRand(getUnifRand(0._TKC, 1._TKC))]
        call disp%show("rand")
        call disp%show( rand )
        call disp%show("rand = [getKolmRand(getUnifRand(0._TKC, 1._TKC, 3_IK))]")
                        rand = [getKolmRand(getUnifRand(0._TKC, 1._TKC, 3_IK))]
        call disp%show("rand")
        call disp%show( rand )
        call disp%skip()
    end block

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_io, only: getErrTableWrite
        use pm_kind, only: TKC => RKH
        real(TKC) :: rand(1500)
        rand = getKolmRand(getUnifRand(0._TKC, 1._TKC, size(rand, 1, IK)))
        if (0 /= getErrTableWrite(SK_"getKolmRand.RK.txt", rand, header = SK_"rand")) error stop "table output failed."
    end block

end program example