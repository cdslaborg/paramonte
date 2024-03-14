program example

    use pm_arrayVerbose, only: getVerbose
    use pm_arrayRefine, only: getRefined
    use pm_arrayChoice, only: getChoice
    use pm_distUnif, only: getUnifRand
    use pm_distBern, only: isHead
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    integer(IK) :: nsam, itry, ntry = 10
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("! Refine 1D array.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        block
            use pm_kind, only: TKC => SK ! all kinds are supported.
            character(:,TKC), allocatable :: array
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("nsam")
                call disp%show( nsam )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = getUnifRand(repeat('A', nsam), repeat('Z', nsam))")
                                array = getUnifRand(repeat('A', nsam), repeat('Z', nsam))
                call disp%show("array")
                call disp%show( array , deliml = TKC_"""" )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0))")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0)) , deliml = TKC_"""" )
                call disp%show("array = getRefined(array, weight, skip)")
                                array = getRefined(array, weight, skip)
                call disp%show("array")
                call disp%show( array , deliml = TKC_"""" )
                call disp%skip()
            end do
        end block

        block
            use pm_kind, only: TKC => SK ! all kinds are supported.
            character(2,TKC), allocatable :: array(:)
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("nsam")
                call disp%show( nsam )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = getUnifRand('AA', 'ZZ', nsam)")
                                array = getUnifRand('AA', 'ZZ', nsam)
                call disp%show("array")
                call disp%show( array , deliml = TKC_"""" )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0))")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0)) , deliml = TKC_"""" )
                call disp%show("array = getRefined(array, weight, skip)")
                                array = getRefined(array, weight, skip)
                call disp%show("array")
                call disp%show( array , deliml = TKC_"""" )
                call disp%skip()
            end do
        end block

        block
            use pm_kind, only: TKC => IK ! all kinds are supported.
            integer(TKC), allocatable :: array(:)
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("nsam")
                call disp%show( nsam )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = getUnifRand(0, 9, nsam)")
                                array = getUnifRand(0, 9, nsam)
                call disp%show("array")
                call disp%show( array )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0))")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0)) )
                call disp%show("array = getRefined(array, weight, skip)")
                                array = getRefined(array, weight, skip)
                call disp%show("array")
                call disp%show( array )
                call disp%skip()
            end do
        end block

        block
            use pm_kind, only: TKC => LK ! all kinds are supported.
            logical(TKC), allocatable :: array(:)
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("nsam")
                call disp%show( nsam )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = getUnifRand(.false., .true., nsam)")
                                array = getUnifRand(.false., .true., nsam)
                call disp%show("array")
                call disp%show( array )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0))")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0)) )
                call disp%show("array = getRefined(array, weight, skip)")
                                array = getRefined(array, weight, skip)
                call disp%show("array")
                call disp%show( array )
                call disp%skip()
            end do
        end block

        block
            use pm_kind, only: TKC => CKS ! all kinds are supported.
            complex(TKC), allocatable :: array(:)
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("nsam")
                call disp%show( nsam )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = cmplx(getUnifRand(0, 9, nsam), getUnifRand(0, 9, nsam), TKC)")
                                array = cmplx(getUnifRand(0, 9, nsam), getUnifRand(0, 9, nsam), TKC)
                call disp%show("array")
                call disp%show( array )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0))")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0)) )
                call disp%show("array = getRefined(array, weight, skip)")
                                array = getRefined(array, weight, skip)
                call disp%show("array")
                call disp%show( array )
                call disp%skip()
            end do
        end block

        block
            use pm_kind, only: TKC => RKS ! all kinds are supported.
            real(TKC), allocatable :: array(:)
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("nsam")
                call disp%show( nsam )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = getUnifRand(0, 9, nsam)")
                                array = getUnifRand(0, 9, nsam)
                call disp%show("array")
                call disp%show( array )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0))")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0)) )
                call disp%show("array = getRefined(array, weight, skip)")
                                array = getRefined(array, weight, skip)
                call disp%show("array")
                call disp%show( array )
                call disp%skip()
            end do
        end block

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%show("! Refine 2D array.")
    call disp%show("!%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        integer(IK) :: dim, ndim

        block
            use pm_kind, only: TKC => SK ! all kinds are supported.
            character(2,TKC), allocatable :: array(:,:), arref(:,:)
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("dim = 2; ndim = merge(0, getUnifRand(1, 3), isHead(.1)); nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                dim = 2; ndim = merge(0, getUnifRand(1, 3), isHead(.1)); nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("[dim, ndim, nsam]")
                call disp%show( [dim, ndim, nsam] )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = getUnifRand('AA', 'ZZ', ndim, nsam)")
                                array = getUnifRand('AA', 'ZZ', ndim, nsam)
                call disp%show("array")
                call disp%show( array , deliml = TKC_"""" )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0), dim)")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0), dim) , deliml = TKC_"""" )
                call disp%show("arref = getRefined(array, dim, weight, skip)")
                                arref = getRefined(array, dim, weight, skip)
                call disp%show("arref")
                call disp%show( arref , deliml = TKC_"""" )
                call disp%show("arref = getRefined(transpose(array), 3_IK - dim, weight, skip)")
                                arref = getRefined(transpose(array), 3_IK - dim, weight, skip)
                call disp%show("arref")
                call disp%show( arref , deliml = TKC_"""" )
                call disp%skip()
            end do
        end block

        block
            use pm_kind, only: TKC => IK ! all kinds are supported.
            integer(TKC), allocatable :: array(:,:), arref(:,:)
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("dim = 2; ndim = merge(0, getUnifRand(1, 3), isHead(.1)); nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                dim = 2; ndim = merge(0, getUnifRand(1, 3), isHead(.1)); nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("[dim, ndim, nsam]")
                call disp%show( [dim, ndim, nsam] )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = getUnifRand(0, 9, ndim, nsam)")
                                array = getUnifRand(0, 9, ndim, nsam)
                call disp%show("array")
                call disp%show( array )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0), dim)")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0), dim) )
                call disp%show("arref = getRefined(array, dim, weight, skip)")
                                arref = getRefined(array, dim, weight, skip)
                call disp%show("arref")
                call disp%show( arref )
                call disp%show("arref = getRefined(transpose(array), 3_IK - dim, weight, skip)")
                                arref = getRefined(transpose(array), 3_IK - dim, weight, skip)
                call disp%show("arref")
                call disp%show( arref )
                call disp%skip()
            end do
        end block

        block
            use pm_kind, only: TKC => LK ! all kinds are supported.
            logical(TKC), allocatable :: array(:,:), arref(:,:)
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("dim = 2; ndim = merge(0, getUnifRand(1, 3), isHead(.1)); nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                dim = 2; ndim = merge(0, getUnifRand(1, 3), isHead(.1)); nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("[dim, ndim, nsam]")
                call disp%show( [dim, ndim, nsam] )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = getUnifRand(.false., .true., ndim, nsam)")
                                array = getUnifRand(.false., .true., ndim, nsam)
                call disp%show("array")
                call disp%show( array )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0), dim)")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0), dim) )
                call disp%show("arref = getRefined(array, dim, weight, skip)")
                                arref = getRefined(array, dim, weight, skip)
                call disp%show("arref")
                call disp%show( arref )
                call disp%show("arref = getRefined(transpose(array), 3_IK - dim, weight, skip)")
                                arref = getRefined(transpose(array), 3_IK - dim, weight, skip)
                call disp%show("arref")
                call disp%show( arref )
                call disp%skip()
            end do
        end block

        block
            use pm_kind, only: TKC => CKS ! all kinds are supported.
            complex(TKC), allocatable :: array(:,:), arref(:,:)
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("dim = 2; ndim = merge(0, getUnifRand(1, 3), isHead(.1)); nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                dim = 2; ndim = merge(0, getUnifRand(1, 3), isHead(.1)); nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("[dim, ndim, nsam]")
                call disp%show( [dim, ndim, nsam] )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = cmplx(getUnifRand(0, 9, ndim, nsam), getUnifRand(0, 9, ndim, nsam), TKC)")
                                array = cmplx(getUnifRand(0, 9, ndim, nsam), getUnifRand(0, 9, ndim, nsam), TKC)
                call disp%show("array")
                call disp%show( array )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0), dim)")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0), dim) )
                call disp%show("arref = getRefined(array, dim, weight, skip)")
                                arref = getRefined(array, dim, weight, skip)
                call disp%show("arref")
                call disp%show( arref )
                call disp%show("arref = getRefined(transpose(array), 3_IK - dim, weight, skip)")
                                arref = getRefined(transpose(array), 3_IK - dim, weight, skip)
                call disp%show("arref")
                call disp%show( arref )
                call disp%skip()
            end do
        end block

        block
            use pm_kind, only: TKC => RKS ! all kinds are supported.
            real(TKC), allocatable :: array(:,:), arref(:,:)
            integer(IK), allocatable :: weight(:)
            integer(IK) :: skip
            do itry = 1, ntry
                call disp%show("dim = 2; ndim = merge(0, getUnifRand(1, 3), isHead(.1)); nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)")
                                dim = 2; ndim = merge(0, getUnifRand(1, 3), isHead(.1)); nsam = getUnifRand(0, 9); skip = getUnifRand(1, 4)
                call disp%show("[dim, ndim, nsam]")
                call disp%show( [dim, ndim, nsam] )
                call disp%show("skip")
                call disp%show( skip )
                call disp%show("weight = getUnifRand(-1, 9, nsam)")
                                weight = getUnifRand(-1, 9, nsam)
                call disp%show("weight")
                call disp%show( weight )
                call disp%show("array = getUnifRand(0, 9, ndim, nsam)")
                                array = getUnifRand(0, 9, ndim, nsam)
                call disp%show("array")
                call disp%show( array )
                call disp%show("getVerbose(array, weight, sum(weight, mask = weight > 0), dim)")
                call disp%show( getVerbose(array, weight, sum(weight, mask = weight > 0), dim) )
                call disp%show("arref = getRefined(array, dim, weight, skip)")
                                arref = getRefined(array, dim, weight, skip)
                call disp%show("arref")
                call disp%show( arref )
                call disp%show("arref = getRefined(transpose(array), 3_IK - dim, weight, skip)")
                                arref = getRefined(transpose(array), 3_IK - dim, weight, skip)
                call disp%show("arref")
                call disp%show( arref )
                call disp%skip()
            end do
        end block

    end block

end program example