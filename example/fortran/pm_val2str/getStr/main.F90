program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: IK8, IK16, IK32, IK64
    use pm_kind, only: CK32, CK64, CK128
    use pm_kind, only: RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_val2str, only: getStr

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getStr('paramonte')")
    call disp%show( getStr('paramonte') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr('  paramonte')")
    call disp%show( getStr('  paramonte') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr('  paramonte  ')")
    call disp%show( getStr('  paramonte  ') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(['paramonte', 'library  '])")
    call disp%show( getStr(['paramonte', 'library  ']) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(['paramonte  ', '  library  '])")
    call disp%show( getStr(['paramonte  ', '  library  ']) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(-1_IK8)")
    call disp%show( getStr(-1_IK8) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(.true.)")
    call disp%show( getStr(.true.) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(.false.)")
    call disp%show( getStr(.false.) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(123_IK16)")
    call disp%show( getStr(123_IK16) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(1234_IK32, signed = .true._LK)")
    call disp%show( getStr(1234_IK32, signed = .true._LK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(1234_IK32, length = 10_IK, signed = .true._LK)")
    call disp%show( getStr(1234_IK32, length = 10_IK, signed = .true._LK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(987654321_IK64, length = 15_IK)")
    call disp%show( getStr(987654321_IK64, length = 15_IK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(987654321_IK64, length = 15_IK, format = SK_'(1I15)')")
    call disp%show( getStr(987654321_IK64, length = 15_IK, format = SK_'(1I15)') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(123._RK32)")
    call disp%show( getStr(123._RK32) )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(123._RK64, format = SK_'(E25.15)') ! use scientific notation with a width of 25.")
    call disp%show( getStr(123._RK64, format = SK_'(E25.15)') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr([9.e1000_RK128, -123456._RK128])")
    call disp%show( getStr([9.e1000_RK128, -123456._RK128]) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr((1._CK128,-1._CK128))")
    call disp%show( getStr((1._CK128,-1._CK128)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr([(1._CK64,-1._CK64), (2._CK64,-2._CK64)], signed = .true._LK)")
    call disp%show( getStr([(1._CK64,-1._CK64), (2._CK64,-2._CK64)], signed = .true._LK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(reshape([(1._CK32,-1._CK32), (2._CK32,-2._CK32), (3._CK32,-3._CK32), (4._CK32,-4._CK32)], shape = [2,2]), signed = .true._LK) ! input object of rank 2")
    call disp%show( getStr(reshape([(1._CK32,-1._CK32), (2._CK32,-2._CK32), (3._CK32,-3._CK32), (4._CK32,-4._CK32)], shape = [2,2]), signed = .true._LK) , deliml = SK_"""" )
    call disp%skip()

end program example