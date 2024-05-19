program example

    use pm_kind, only: SK, IK, LK
    use pm_sampleWeight, only: setReweight
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    integer(IK) :: lenwei, itry, ntry = 10
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: IKG => IK ! all processor kinds are supported.
        integer(IKG), allocatable :: weight(:)
        integer(IKG) :: skip
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
            call disp%show("call setReweight(weight, skip)")
                            call setReweight(weight, skip)
            call disp%show("weight")
            call disp%show( weight )
            call disp%skip()
        end do
    end block

    block
        use pm_kind, only: RKG => RKS ! all processor kinds are supported.
        real(RKG), allocatable :: weight(:)
        real(RKG) :: skip
        do itry = 1, ntry
            call disp%show("lenwei = getUnifRand(0, 9); skip = getUnifRand(.1_RKG, 4._RKG)")
                            lenwei = getUnifRand(0, 9); skip = getUnifRand(.1_RKG, 4._RKG)
            call disp%show("lenwei")
            call disp%show( lenwei )
            call disp%show("skip")
            call disp%show( skip )
            call disp%show("weight = getUnifRand(-1._RKG, 9._RKG, lenwei)")
                            weight = getUnifRand(-1._RKG, 9._RKG, lenwei)
            call disp%show("weight")
            call disp%show( weight )
            call disp%show("call setReweight(weight, skip)")
                            call setReweight(weight, skip)
            call disp%show("weight")
            call disp%show( weight )
            call disp%skip()
        end do
    end block

end program example