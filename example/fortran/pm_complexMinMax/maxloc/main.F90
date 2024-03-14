program example

    use pm_kind, only: SK, IK
    use pm_complexMinMax, only: maxloc
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    integer(IK) :: ndim, nsam
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: CKC => CKS
        complex(CKC), allocatable :: array(:,:)
        call disp%skip()
        call disp%show("ndim = 3; nsam = 5")
                        ndim = 3; nsam = 5
        call disp%show("array = getUnifRand((-9._CKC, -9._CKC), (+9._CKC, +9._CKC), ndim, nsam)")
                        array = getUnifRand((-9._CKC, -9._CKC), (+9._CKC, +9._CKC), ndim, nsam)
        call disp%show("array")
        call disp%show( array )
        call disp%show("maxloc(array)")
        call disp%show( maxloc(array) )
        call disp%show("maxloc(array, dim = 1_IK)")
        call disp%show( maxloc(array, dim = 1_IK) )
        call disp%show("maxloc(array, dim = 2_IK)")
        call disp%show( maxloc(array, dim = 2_IK) )
        call disp%show("maxloc([complex(CKC)::])")
        call disp%show( maxloc([complex(CKC)::]) )
        call disp%show("maxloc([complex(CKC)::], dim = 1_IK)")
        call disp%show( maxloc([complex(CKC)::], dim = 1_IK) )
        call disp%show("maxloc([complex(CKC)::], dim = 2_IK)")
        call disp%show( maxloc([complex(CKC)::], dim = 2_IK) )
        call disp%skip()
    end block

    block
        use pm_kind, only: CKC => CKD
        complex(CKC), allocatable :: array(:,:)
        call disp%skip()
        call disp%show("ndim = 3; nsam = 5")
                        ndim = 3; nsam = 5
        call disp%show("array = getUnifRand((-9._CKC, -9._CKC), (+9._CKC, +9._CKC), ndim, nsam)")
                        array = getUnifRand((-9._CKC, -9._CKC), (+9._CKC, +9._CKC), ndim, nsam)
        call disp%show("array")
        call disp%show( array )
        call disp%show("maxloc(array)")
        call disp%show( maxloc(array) )
        call disp%show("maxloc(array, dim = 1_IK)")
        call disp%show( maxloc(array, dim = 1_IK) )
        call disp%show("maxloc(array, dim = 2_IK)")
        call disp%show( maxloc(array, dim = 2_IK) )
        call disp%show("maxloc([complex(CKC)::])")
        call disp%show( maxloc([complex(CKC)::]) )
        call disp%show("maxloc([complex(CKC)::], dim = 1_IK)")
        call disp%show( maxloc([complex(CKC)::], dim = 1_IK) )
        call disp%show("maxloc([complex(CKC)::], dim = 2_IK)")
        call disp%show( maxloc([complex(CKC)::], dim = 2_IK) )
        call disp%skip()
    end block

    block
        use pm_kind, only: CKC => CKH
        complex(CKC), allocatable :: array(:,:)
        call disp%skip()
        call disp%show("ndim = 3; nsam = 5")
                        ndim = 3; nsam = 5
        call disp%show("array = getUnifRand((-9._CKC, -9._CKC), (+9._CKC, +9._CKC), ndim, nsam)")
                        array = getUnifRand((-9._CKC, -9._CKC), (+9._CKC, +9._CKC), ndim, nsam)
        call disp%show("array")
        call disp%show( array )
        call disp%show("maxloc(array)")
        call disp%show( maxloc(array) )
        call disp%show("maxloc(array, dim = 1_IK)")
        call disp%show( maxloc(array, dim = 1_IK) )
        call disp%show("maxloc(array, dim = 2_IK)")
        call disp%show( maxloc(array, dim = 2_IK) )
        call disp%show("maxloc([complex(CKC)::])")
        call disp%show( maxloc([complex(CKC)::]) )
        call disp%show("maxloc([complex(CKC)::], dim = 1_IK)")
        call disp%show( maxloc([complex(CKC)::], dim = 1_IK) )
        call disp%show("maxloc([complex(CKC)::], dim = 2_IK)")
        call disp%show( maxloc([complex(CKC)::], dim = 2_IK) )
        call disp%skip()
    end block

end program example