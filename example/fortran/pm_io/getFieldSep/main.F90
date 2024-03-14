program example

    use iso_fortran_env, only: output_unit, input_unit, error_unit
    use pm_str, only: NLC
    use pm_val2str, only: getStr
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_io, only: getErrTableWrite
    use pm_io, only: getContentsFrom
    use pm_io, only: getFieldSep, csv, fld
    use pm_io, only: setContentsTo

    implicit none

    integer(IK) :: nfield
    character(:, SK), allocatable :: file
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Search for single-character field separators.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

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
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps = SKC_';|, ')")
        call disp%show( getFieldSep(file, seps = SKC_';|, ') , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_';|, ', form = csv)")
        call disp%show( getFieldSep(file, seps = SKC_';|, ', form = csv) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_';|, ', form = csv, nfield = nfield)")
        call disp%show( getFieldSep(file, seps = SKC_';|, ', form = csv, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = repeat(SK_'""a"", '//getStr([1, 2, 3])//NLC, 2_IK))")
                        call setContentsTo(file, contents = repeat(SK_"""a"", "//getStr([1, 2, 3])//NLC, 2_IK))
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps = SKC_', ')")
        call disp%show( getFieldSep(file, seps = SKC_', ') , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_' ,')")
        call disp%show( getFieldSep(file, seps = SKC_' ,') , deliml = SK_"""" )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = repeat(SK_'(1, -1), (2, -2), (3, -3)'//NLC, 2_IK))")
                        call setContentsTo(file, contents = repeat(SK_'(1, -1),   (2, -2), (3, -3)'//NLC, 2_IK))
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps = SKC_', ', nfield = nfield)")
        call disp%show( getFieldSep(file, seps = SKC_', ', nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%show("getFieldSep(file, seps = SKC_' ,', nfield = nfield)")
        call disp%show( getFieldSep(file, seps = SKC_' ,', nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%show("getFieldSep(file, seps = SKC_', ', form = fld, nfield = nfield)")
        call disp%show( getFieldSep(file, seps = SKC_', ', form = fld, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%show("getFieldSep(file, seps = SKC_' ,', form = fld, nfield = nfield)")
        call disp%show( getFieldSep(file, seps = SKC_' ,', form = fld, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = SK_'""a,"",'//getStr([1, 2, 3])//NLC//SK_'""a"",'//getStr([1, 2, 3])//NLC)")
                        call setContentsTo(file, contents = SK_"""a,"","//getStr([1, 2, 3])//NLC//SK_"""a"","//getStr([1, 2, 3])//NLC)
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps = SKC_', ')")
        call disp%show( getFieldSep(file, seps = SKC_', ') , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_', ', form = csv) ! quoted strings are respected in csv format.")
        call disp%show( getFieldSep(file, seps = SKC_', ', form = csv) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_' ,')")
        call disp%show( getFieldSep(file, seps = SKC_' ,') , deliml = SK_"""" )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = SK_'""a""   '//getStr([1, 2, 3])//NLC//SK_'""a"" '//getStr([1, 2, 3])//NLC)")
                        call setContentsTo(file, contents = SK_"""a""   "//getStr([1, 2, 3])//NLC//SK_"""a"" "//getStr([1, 2, 3])//NLC)
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps = SKC_' ')")
        call disp%show( getFieldSep(file, seps = SKC_' ') , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_' ', form = fld) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (fld) format.")
        call disp%show( getFieldSep(file, seps = SKC_' ', form = fld) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_' ', form = fld, nfield = nfield) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (fld) format.")
        call disp%show( getFieldSep(file, seps = SKC_' ', form = fld, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = SK_'""a'//NLC//'"" '//getStr([1, 2, 3])//NLC//SK_'""a"" '//getStr([1, 2, 3])//NLC) ! double-line csv field")
                        call setContentsTo(file, contents = SK_"""a"//NLC//""" "//getStr([1, 2, 3])//NLC//SK_"""a"" "//getStr([1, 2, 3])//NLC)
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps = SKC_',')")
        call disp%show( getFieldSep(file, seps = SKC_',') , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_',', form = csv) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (csv) format.")
        call disp%show( getFieldSep(file, seps = SKC_',', form = csv) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_',', form = csv, nfield = nfield) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (csv) format.")
        call disp%show( getFieldSep(file, seps = SKC_',', form = csv, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%show("getFieldSep(file, seps = SKC_',', form = csv, nfield = nfield) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (csv) format.")
        call disp%show( getFieldSep(file, seps = SKC_',', form = csv, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = SK_'""a'//NLC//NLC//'""   '//getStr([1, 2, 3])//NLC//SK_'""a"" '//getStr([1, 2, 3])//NLC) ! three-line FLD field")
                        call setContentsTo(file, contents = SK_"""a"//NLC//NLC//"""   "//getStr([1, 2, 3])//NLC//SK_"""a"" "//getStr([1, 2, 3])//NLC)
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps = SKC_' ')")
        call disp%show( getFieldSep(file, seps = SKC_' ') , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_' ', form = fld) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (fld) format.")
        call disp%show( getFieldSep(file, seps = SKC_' ', form = fld) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps = SKC_' ', form = fld, nfield = nfield) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (fld) format.")
        call disp%show( getFieldSep(file, seps = SKC_' ', form = fld, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%show("getFieldSep(file, seps = SKC_',', form = fld, nfield = nfield) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (fld) format.")
        call disp%show( getFieldSep(file, seps = SKC_',', form = fld, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%skip
    end block

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Retrieve whole string field separators.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    block
        use pm_io, only: css_type
        use pm_kind, only: SKC => SK
        character(2,SKC), allocatable :: table(:,:)
        type(css_type), allocatable :: seps(:)
        call disp%skip
        call disp%show("table = getUnifRand('aa', 'zz', 3_IK, 6_IK)")
                        table = getUnifRand('aa', 'zz', 3_IK, 6_IK)
        call disp%show("table")
        call disp%show( table , deliml = SK_"""" )
        call disp%show("file = 'temp.temp'")
                        file = 'temp.temp'
        call disp%show("if (0 /= getErrTableWrite(file, table, deliml = SKC_'''')) error stop 'table write failed.'")
                        if (0 /= getErrTableWrite(file, table, deliml = SKC_'''')) error stop 'table write failed.'
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("seps = [css_type(SK_'|||'), css_type(SK_',')]")
                        seps = [css_type(SK_'|||'), css_type(SK_',')]
        call disp%show("seps")
        call disp%show( seps , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps)")
        call disp%show( getFieldSep(file, seps) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps, form = csv)")
        call disp%show( getFieldSep(file, seps, form = csv) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps, form = csv, nfield = nfield)")
        call disp%show( getFieldSep(file, seps, form = csv, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = repeat(SK_'""a"", '//getStr([1, 2, 3])//NLC, 2_IK))")
                        call setContentsTo(file, contents = repeat(SK_"""a"", "//getStr([1, 2, 3])//NLC, 2_IK))
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("seps = css_type([SK_',', SK_' '], trimmed = .true._LK)")
                        seps = css_type([SK_',', SK_' '], trimmed = .true._LK)
        call disp%show("seps")
        call disp%show( seps , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps)")
        call disp%show( getFieldSep(file, seps) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps(2:1:-1))")
        call disp%show( getFieldSep(file, seps(2:1:-1)) , deliml = SK_"""" )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = repeat(SK_'(1, -1), (2, -2), (3, -3)'//NLC, 2_IK))")
                        call setContentsTo(file, contents = repeat(SK_'(1, -1),   (2, -2), (3, -3)'//NLC, 2_IK))
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps, nfield = nfield)")
        call disp%show( getFieldSep(file, seps, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%show("getFieldSep(file, seps(2:1:-1), nfield = nfield)")
        call disp%show( getFieldSep(file, seps(2:1:-1), nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%show("getFieldSep(file, seps, form = fld, nfield = nfield)")
        call disp%show( getFieldSep(file, seps, form = fld, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%show("getFieldSep(file, seps(2:1:-1), form = fld, nfield = nfield)")
        call disp%show( getFieldSep(file, seps(2:1:-1), form = fld, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = SK_'""a,"",'//getStr([1, 2, 3])//NLC//SK_'""a"",'//getStr([1, 2, 3])//NLC)")
                        call setContentsTo(file, contents = SK_"""a,"","//getStr([1, 2, 3])//NLC//SK_"""a"","//getStr([1, 2, 3])//NLC)
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps)")
        call disp%show( getFieldSep(file, seps) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps, form = csv) ! quoted strings are respected in csv format.")
        call disp%show( getFieldSep(file, seps, form = csv) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps(2:1:-1))")
        call disp%show( getFieldSep(file, seps(2:1:-1)) , deliml = SK_"""" )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = SK_'""a""   '//getStr([1, 2, 3])//NLC//SK_'""a"" '//getStr([1, 2, 3])//NLC)")
                        call setContentsTo(file, contents = SK_"""a""   "//getStr([1, 2, 3])//NLC//SK_"""a"" "//getStr([1, 2, 3])//NLC)
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps(2:2))")
        call disp%show( getFieldSep(file, seps(2:2)) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps(2:2), form = fld) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (fld) format.")
        call disp%show( getFieldSep(file, seps(2:2), form = fld) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps(2:2), form = fld, nfield = nfield) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (fld) format.")
        call disp%show( getFieldSep(file, seps(2:2), form = fld, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = SK_'""a'//NLC//'"" '//getStr([1, 2, 3])//NLC//SK_'""a"" '//getStr([1, 2, 3])//NLC) ! double-line csv field")
                        call setContentsTo(file, contents = SK_"""a"//NLC//""" "//getStr([1, 2, 3])//NLC//SK_"""a"" "//getStr([1, 2, 3])//NLC)
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps(1:1))")
        call disp%show( getFieldSep(file, seps(1:1)) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps(1:1), form = csv) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (csv) format.")
        call disp%show( getFieldSep(file, seps(1:1), form = csv) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps(1:1), form = csv, nfield = nfield) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (csv) format.")
        call disp%show( getFieldSep(file, seps(1:1), form = csv, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%show("getFieldSep(file, seps(1:1), form = csv, nfield = nfield) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (csv) format.")
        call disp%show( getFieldSep(file, seps(1:1), form = csv, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%skip
        call disp%show("call setContentsTo(file, contents = SK_'""a'//NLC//NLC//'""   '//getStr([1, 2, 3])//NLC//SK_'""a"" '//getStr([1, 2, 3])//NLC) ! three-line FLD field")
                        call setContentsTo(file, contents = SK_"""a"//NLC//NLC//"""   "//getStr([1, 2, 3])//NLC//SK_"""a"" "//getStr([1, 2, 3])//NLC)
        call disp%show("getContentsFrom(file)")
        call disp%show( getContentsFrom(file) )
        call disp%show("getFieldSep(file, seps(2:2))")
        call disp%show( getFieldSep(file, seps(2:2)) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps(2:2), form = fld) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (fld) format.")
        call disp%show( getFieldSep(file, seps(2:2), form = fld) , deliml = SK_"""" )
        call disp%show("getFieldSep(file, seps(2:2), form = fld, nfield = nfield) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (fld) format.")
        call disp%show( getFieldSep(file, seps(2:2), form = fld, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%show("getFieldSep(file, seps(1:1), form = fld, nfield = nfield) ! multiple adjacent blanks count as a single separator in Fortran-list-directed (fld) format.")
        call disp%show( getFieldSep(file, seps(1:1), form = fld, nfield = nfield) , deliml = SK_"""" )
        call disp%show("nfield")
        call disp%show( nfield )
        call disp%skip
    end block

end program example