program example

    use iso_fortran_env, only: output_unit, input_unit, error_unit
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_io, only: LEN_IOMSG, trans
    use pm_io, only: getErrTableWrite
    use pm_io, only: getErrTableRead
    use pm_io, only: getContentsFrom
    use pm_distUnif, only: getUnifRand

    implicit none

    character(LEN_IOMSG, SK) :: iomsg
    character(:, SK), allocatable :: file, header
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: SKC => SK
        character(2,SKC), allocatable :: table(:,:)
        call disp%skip
        call disp%show("table = getUnifRand('aa', 'zz', 3_IK, 6_IK)")
                        table = getUnifRand('aa', 'zz', 3_IK, 6_IK)
        call disp%show("table")
        call disp%show( table , deliml = SK_"""" )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, deliml = SKC_'''')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, deliml = SKC_'''')) error stop 'table write failed.'
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table)) error stop 'table write failed.'")
                        if (0 /= getErrTableRead(file, table)) error stop 'table write failed.'
        call disp%show("table")
        call disp%show( table , deliml = SK_"""" )
        call disp%skip
    end block

    block
        use pm_kind, only: IKC => IK
        integer(IKC), allocatable :: table(:,:)
        call disp%skip
        call disp%show("table = getUnifRand(-1, 1, 4_IK, 2_IK)")
                        table = getUnifRand(-1, 1, 4_IK, 2_IK)
        call disp%show("table")
        call disp%show( table )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, header = 'col1 col2', sep = SK_' ')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, header = 'col1 col2', sep = SK_' ')) error stop 'table write failed.'
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, header = header)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, header = header)) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, header = header, sep = SK_' ')) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, header = header, sep = SK_' ')) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%skip
        call disp%show("if (0 /= getErrTableWrite(file, table, trans, sep = '|', roff = 2_IK)) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, trans, sep = '|', roff = 2_IK)) error stop 'table write failed.'
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, sep = SK_'|', roff = 2_IK, iomsg = iomsg)) error stop trim(iomsg)")
                        if (0 /= getErrTableRead(file, table, sep = SK_'|', roff = 2_IK, iomsg = iomsg)) error stop trim(iomsg)
        call disp%show("table")
        call disp%show( table )
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%skip
    end block

    block
        use pm_kind, only: LKC => LK
        logical(LKC), allocatable :: table(:,:)
        call disp%skip
        call disp%show("table = getUnifRand(.false., .true., 4_IK, 2_IK)")
                        table = getUnifRand(.false., .true., 4_IK, 2_IK)
        call disp%show("table")
        call disp%show( table )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, header = 'col1 col2', sep = SK_' ')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, header = 'col1 col2', sep = SK_' ')) error stop 'table write failed.'
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, header = header)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, header = header)) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, header = header, sep = SK_' ')) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, header = header, sep = SK_' ')) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%skip
        call disp%show("if (0 /= getErrTableWrite(file, table, trans, sep = '|', roff = 2_IK)) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, trans, sep = '|', roff = 2_IK)) error stop 'table write failed.'
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'
        call disp%show("table")
        call disp%show( table )
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%skip
    end block

    block
        use pm_kind, only: CKC => CK32
        complex(CKC), allocatable :: table(:,:)
        call disp%skip
        call disp%show("table = getUnifRand((-1., -1.), (+1., +1.), 4_IK, 2_IK)")
                        table = getUnifRand((-1., -1.), (+1., +1.), 4_IK, 2_IK)
        call disp%show("table")
        call disp%show( table )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, header = 'col1 col2', sep = SK_' ')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, header = 'col1 col2', sep = SK_' ')) error stop 'table write failed.'
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, header = header)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, header = header)) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, header = header, sep = SK_' ')) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, header = header, sep = SK_' ')) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%skip
        call disp%show("if (0 /= getErrTableWrite(file, table, trans, sep = '|', roff = 2_IK)) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, trans, sep = '|', roff = 2_IK)) error stop 'table write failed.'
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, sep = SK_'|', roff = 2_IK, iomsg = iomsg)) error stop trim(iomsg)")
                        if (0 /= getErrTableRead(file, table, sep = SK_'|', roff = 2_IK, iomsg = iomsg)) error stop trim(iomsg)
        call disp%show("table")
        call disp%show( table )
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%skip
    end block

    block
        use pm_kind, only: RKC => RKS
        real(RKC), allocatable :: table(:,:)
        call disp%skip
        call disp%show("table = getUnifRand(-1, 1, 4_IK, 2_IK)")
                        table = getUnifRand(-1, 1, 4_IK, 2_IK)
        call disp%show("table")
        call disp%show( table )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, header = 'col1 col2', sep = SK_' ')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, header = 'col1 col2', sep = SK_' ')) error stop 'table write failed.'
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, header = header)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, header = header)) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, header = header, sep = SK_' ')) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, header = header, sep = SK_' ')) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%skip
        call disp%show("if (0 /= getErrTableWrite(file, table, trans, sep = '|', roff = 2_IK)) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, trans, sep = '|', roff = 2_IK)) error stop 'table write failed.'
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'
        call disp%show("table")
        call disp%show( table )
        call disp%show("deallocate(table)")
                        deallocate(table)
        call disp%show("if (0 /= getErrTableRead(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'")
                        if (0 /= getErrTableRead(file, table, trans, sep = SK_'|', roff = 2_IK)) error stop 'table read failed.'
        call disp%show("header")
        call disp%show( header , deliml = SK_"""" )
        call disp%show("table")
        call disp%show( table )
        call disp%skip
    end block

end program example