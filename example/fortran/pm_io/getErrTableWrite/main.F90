program example

    use iso_fortran_env, only: output_unit, input_unit, error_unit
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: getErrTableWrite, trans
    use pm_io, only: getContentsFrom
    use pm_distUnif, only: getUnifRand

    implicit none

    character(:, SK), allocatable :: file
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: SKG => SK
        character(2,SKG), allocatable :: table(:,:)
        call disp%skip
        call disp%show("table = getUnifRand('aa', 'zz', 3_IK, 6_IK)")
                        table = getUnifRand('aa', 'zz', 3_IK, 6_IK)
        call disp%show("table")
        call disp%show( table , deliml = SK_"""" )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, deliml = SK_'''')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, deliml = SK_'''')) error stop 'table write failed.'
        call disp%show("getContentsFrom(file, del = .true._LK)")
        call disp%show( getContentsFrom(file, del = .true._LK) )
        call disp%skip
    end block

    block
        use pm_kind, only: IKG => IK
        integer(IKG), allocatable :: table(:,:)
        call disp%skip
        call disp%show("table = getUnifRand(-1, 1, 4_IK, 2_IK)")
                        table = getUnifRand(-1, 1, 4_IK, 2_IK)
        call disp%show("table")
        call disp%show( table )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, header = SK_'col1 col2', sep = SK_' ')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, header = SK_'col1 col2', sep = SK_' ')) error stop 'table write failed.'
        call disp%show("getContentsFrom(file, del = .true._LK)")
        call disp%show( getContentsFrom(file, del = .true._LK) )
        call disp%show("if (0 /= getErrTableWrite(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table write failed.'
        call disp%show("getContentsFrom(file, del = .true._LK)")
        call disp%show( getContentsFrom(file, del = .true._LK) )
        call disp%skip
    end block

    block
        use pm_kind, only: LKG => LK
        logical(LKG), allocatable :: table(:,:)
        call disp%skip
        call disp%show("table = getUnifRand(.false., .true., 4_IK, 2_IK)")
                        table = getUnifRand(.false., .true., 4_IK, 2_IK)
        call disp%show("table")
        call disp%show( table )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, header = SK_'col1 col2', sep = SK_' ')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, header = SK_'col1 col2', sep = SK_' ')) error stop 'table write failed.'
        call disp%show("getContentsFrom(file, del = .true._LK)")
        call disp%show( getContentsFrom(file, del = .true._LK) )
        call disp%show("if (0 /= getErrTableWrite(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table write failed.'
        call disp%show("getContentsFrom(file, del = .true._LK)")
        call disp%show( getContentsFrom(file, del = .true._LK) )
        call disp%skip
    end block

    block
        use pm_kind, only: CKG => CKS
        complex(CKG), allocatable :: table(:,:)
        call disp%skip
        call disp%show("table = getUnifRand((-1., -1.), (+1., +1.), 4_IK, 2_IK)")
                        table = getUnifRand((-1., -1.), (+1., +1.), 4_IK, 2_IK)
        call disp%show("table")
        call disp%show( table )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, header = SK_'col1 col2', sep = SK_' ')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, header = SK_'col1 col2', sep = SK_' ')) error stop 'table write failed.'
        call disp%show("getContentsFrom(file, del = .true._LK)")
        call disp%show( getContentsFrom(file, del = .true._LK) )
        call disp%show("if (0 /= getErrTableWrite(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table write failed.'
        call disp%show("getContentsFrom(file, del = .true._LK)")
        call disp%show( getContentsFrom(file, del = .true._LK) )
        call disp%skip
    end block

    block
        use pm_kind, only: RKG => RKS
        real(RKG), allocatable :: table(:,:)
        call disp%skip
        call disp%show("table = getUnifRand(-1, 1, 4_IK, 2_IK)")
                        table = getUnifRand(-1, 1, 4_IK, 2_IK)
        call disp%show("table")
        call disp%show( table )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, header = SK_'col1 col2', sep = SK_' ')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, header = SK_'col1 col2', sep = SK_' ')) error stop 'table write failed.'
        call disp%show("getContentsFrom(file, del = .true._LK)")
        call disp%show( getContentsFrom(file, del = .true._LK) )
        call disp%show("if (0 /= getErrTableWrite(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table write failed.'
        call disp%show("getContentsFrom(file, del = .true._LK)")
        call disp%show( getContentsFrom(file, del = .true._LK) )
        call disp%skip
    end block

end program example