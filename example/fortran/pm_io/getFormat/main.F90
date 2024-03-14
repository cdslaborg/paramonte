program example

    use pm_kind, only: SK, IK, LK, RKS, RKH
    use pm_io, only: display_type
    use pm_io, only: getFormat

    implicit none

    character(:, SK), allocatable :: format

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("format = getFormat()")
                    format = getFormat()
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("format = getFormat([''])")
                    format = getFormat([''])
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("format = getFormat([1])")
                    format = getFormat([1])
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("format = getFormat([.false.])")
                    format = getFormat([.false.])
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("format = getFormat([cmplx(0, kind = RKS)])")
                    format = getFormat([cmplx(0, kind = RKS)])
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("write(disp%unit, format) cmplx(3.1415926535897932384626433832795028841971_RKH, kind = RKH)")
                    write(disp%unit, format) cmplx(3.1415926535897932384626433832795028841971_RKH, kind = RKH)
    call disp%show("format = getFormat([cmplx(0, kind = RKH)])")
                    format = getFormat([cmplx(0, kind = RKH)])
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("write(disp%unit, format) cmplx(3.1415926535897932384626433832795028841971_RKH)")
                    write(disp%unit, format) cmplx(3.1415926535897932384626433832795028841971_RKH)
    call disp%show("format = getFormat(deliml = '', subsep = SK_'', delimr = 'i', subcount = 2_IK, ndigit = 1_IK, signed = .true._LK) ! math-style complex value.")
                    format = getFormat(deliml = '', subsep = SK_'', delimr = 'i', subcount = 2_IK, ndigit = 1_IK, signed = .true._LK)
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("write(disp%unit, format) cmplx(3.1415926535897932384626433832795028841971_RKH)")
                    write(disp%unit, format) cmplx(3.1415926535897932384626433832795028841971_RKH)
    call disp%show("format = getFormat([real(0, kind = RKS)])")
                    format = getFormat([real(0, kind = RKS)])
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("write(disp%unit, format) 'pi', 3.1415926535897932384626433832795028841971_RKH")
                    write(disp%unit, format) 'pi', 3.1415926535897932384626433832795028841971_RKH
    call disp%show("format = getFormat([real(0, kind = RKH)])")
                    format = getFormat([real(0, kind = RKH)])
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("write(disp%unit, format) 'pi', 3.1415926535897932384626433832795028841971_RKH")
                    write(disp%unit, format) 'pi', 3.1415926535897932384626433832795028841971_RKH
    call disp%skip

    call disp%skip
    call disp%show("format = getFormat(prefix = SK_'ParaMonte: ', sep = SK_' = ', ndigit = 8_IK, signed = .true._LK)")
                    format = getFormat(prefix = SK_'ParaMonte: ', sep = SK_' = ', ndigit = 8_IK, signed = .true._LK)
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("write(disp%unit, format) 'huge(0._RKH)', huge(0._RKH)")
    call disp%skip

    call disp%skip
    call disp%show("format = getFormat(mold = [0._RKH], prefix = SK_'ParaMonte: ', sep = SK_' = ')")
                    format = getFormat(mold = [0._RKH], prefix = SK_'ParaMonte: ', sep = SK_' = ')
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("write(disp%unit, format) 'huge(0._RKH)', huge(0._RKH)")
                    write(disp%unit, format) 'huge(0._RKH)', huge(0._RKH)
    call disp%skip

    call disp%skip
    call disp%show("format = getFormat()")
                    format = getFormat()
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%skip

    call disp%skip
    call disp%show("format = getFormat(mold = [SK_''], prefix = SK_'ParaMonte: ', subcount = 2, deliml = SK_'''', subsep = SK_' ') ! string format with delimited subfields")
                    format = getFormat(mold = [SK_''], prefix = SK_'ParaMonte: ', subcount = 2, deliml = SK_'''', subsep = SK_' ')
    call disp%show("format")
    call disp%show( format , deliml = '''' )
    call disp%show("write(disp%unit, format) 'monte', 'carlo', 'machine', 'learning', 'parallel', 'library'")
                    write(disp%unit, format) 'monte', 'carlo', 'machine', 'learning', 'parallel', 'library'
    call disp%skip

end program example