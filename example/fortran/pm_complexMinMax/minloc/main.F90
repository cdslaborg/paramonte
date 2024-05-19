program example

    use pm_kind, only: SK, IK
    use pm_complexMinMax, only: minloc
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    integer(IK) :: ndim, nsam
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: CKG => CKS
        complex(CKG), allocatable :: array(:,:)
        call disp%skip()
        call disp%show("ndim = 3; nsam = 5")
                        ndim = 3; nsam = 5
        call disp%show("array = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), ndim, nsam)")
                        array = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), ndim, nsam)
        call disp%show("array")
        call disp%show( array )
        call disp%show("minloc(array)")
        call disp%show( minloc(array) )
        call disp%show("minloc(array, dim = 1_IK)")
        call disp%show( minloc(array, dim = 1_IK) )
        call disp%show("minloc(array, dim = 2_IK)")
        call disp%show( minloc(array, dim = 2_IK) )
        call disp%show("minloc([complex(CKG)::])")
        call disp%show( minloc([complex(CKG)::]) )
        call disp%show("minloc([complex(CKG)::], dim = 1_IK)")
        call disp%show( minloc([complex(CKG)::], dim = 1_IK) )
        call disp%show("minloc([complex(CKG)::], dim = 2_IK)")
        call disp%show( minloc([complex(CKG)::], dim = 2_IK) )
        call disp%skip()
    end block

    block
        use pm_kind, only: CKG => CKD
        complex(CKG), allocatable :: array(:,:)
        call disp%skip()
        call disp%show("ndim = 3; nsam = 5")
                        ndim = 3; nsam = 5
        call disp%show("array = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), ndim, nsam)")
                        array = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), ndim, nsam)
        call disp%show("array")
        call disp%show( array )
        call disp%show("minloc(array)")
        call disp%show( minloc(array) )
        call disp%show("minloc(array, dim = 1_IK)")
        call disp%show( minloc(array, dim = 1_IK) )
        call disp%show("minloc(array, dim = 2_IK)")
        call disp%show( minloc(array, dim = 2_IK) )
        call disp%show("minloc([complex(CKG)::])")
        call disp%show( minloc([complex(CKG)::]) )
        call disp%show("minloc([complex(CKG)::], dim = 1_IK)")
        call disp%show( minloc([complex(CKG)::], dim = 1_IK) )
        call disp%show("minloc([complex(CKG)::], dim = 2_IK)")
        call disp%show( minloc([complex(CKG)::], dim = 2_IK) )
        call disp%skip()
    end block

    block
        use pm_kind, only: CKG => CKH
        complex(CKG), allocatable :: array(:,:)
        call disp%skip()
        call disp%show("ndim = 3; nsam = 5")
                        ndim = 3; nsam = 5
        call disp%show("array = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), ndim, nsam)")
                        array = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), ndim, nsam)
        call disp%show("array")
        call disp%show( array )
        call disp%show("minloc(array)")
        call disp%show( minloc(array) )
        call disp%show("minloc(array, dim = 1_IK)")
        call disp%show( minloc(array, dim = 1_IK) )
        call disp%show("minloc(array, dim = 2_IK)")
        call disp%show( minloc(array, dim = 2_IK) )
        call disp%show("minloc([complex(CKG)::])")
        call disp%show( minloc([complex(CKG)::]) )
        call disp%show("minloc([complex(CKG)::], dim = 1_IK)")
        call disp%show( minloc([complex(CKG)::], dim = 1_IK) )
        call disp%show("minloc([complex(CKG)::], dim = 2_IK)")
        call disp%show( minloc([complex(CKG)::], dim = 2_IK) )
        call disp%skip()
    end block

end program example