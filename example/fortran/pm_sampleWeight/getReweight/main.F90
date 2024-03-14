program example

    use pm_kind, only: SK, IK, LK
    use pm_sampleWeight, only: getReweight
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    integer(IK) :: lenwei, itry, ntry = 10
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: IKC => IK ! all processor kinds are supported.
        integer(IKC), allocatable :: weight(:)
        integer(IKC) :: skip
        do itry = 1, ntry
            call disp%show("lenwei = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                            lenwei = getUnifRand(0, 9); skip = getUnifRand(1, 4)
            call disp%show("lenwei")
            call disp%show( lenwei )
            call disp%show("skip")
            call disp%show( skip )
            call disp%show("weight = getUnifRand(-1, 9, lenwei)")
                            weight = getUnifRand(-1, 9, lenwei)
            call disp%show("weight")
            call disp%show( weight )
            call disp%show("weight = getReweight(weight, skip)")
                            weight = getReweight(weight, skip)
            call disp%show("weight")
            call disp%show( weight )
            call disp%skip()
        end do
    end block

    block
        use pm_kind, only: RKC => RKS ! all processor kinds are supported.
        real(RKC), allocatable :: weight(:)
        real(RKC) :: skip
        do itry = 1, ntry
            call disp%show("lenwei = getUnifRand(0, 9); skip = getUnifRand(.1_RKC, 4._RKC)")
                            lenwei = getUnifRand(0, 9); skip = getUnifRand(.1_RKC, 4._RKC)
            call disp%show("lenwei")
            call disp%show( lenwei )
            call disp%show("skip")
            call disp%show( skip )
            call disp%show("weight = getUnifRand(-1._RKC, 9._RKC, lenwei)")
                            weight = getUnifRand(-1._RKC, 9._RKC, lenwei)
            call disp%show("weight")
            call disp%show( weight )
            call disp%show("weight = getReweight(weight, skip)")
                            weight = getReweight(weight, skip)
            call disp%show("weight")
            call disp%show( weight )
            call disp%skip()
        end do
    end block

end program example