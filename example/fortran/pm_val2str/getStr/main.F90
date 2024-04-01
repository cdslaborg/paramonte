program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: IKL, IKS, IKD, IKH
    use pm_kind, only: CKL, CKD, CKH
    use pm_kind, only: RKL, RKD, RKH
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
    call disp%show("getStr(-1_IKL)")
    call disp%show( getStr(-1_IKL) , deliml = SK_"""" )
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
    call disp%show("getStr(123_IKS)")
    call disp%show( getStr(123_IKS) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(1234_IKD, signed = .true._LK)")
    call disp%show( getStr(1234_IKD, signed = .true._LK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(1234_IKD, length = 10_IK, signed = .true._LK)")
    call disp%show( getStr(1234_IKD, length = 10_IK, signed = .true._LK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(987654321_IKH, length = 15_IK)")
    call disp%show( getStr(987654321_IKH, length = 15_IK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(987654321_IKH, length = 15_IK, format = SK_'(1I15)')")
    call disp%show( getStr(987654321_IKH, length = 15_IK, format = SK_'(1I15)') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(123._RKL)")
    call disp%show( getStr(123._RKL) )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(123._RKD, format = SK_'(E25.15)') ! use scientific notation with a width of 25.")
    call disp%show( getStr(123._RKD, format = SK_'(E25.15)') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr([9.e100_RKH, -123456._RKH])")
    call disp%show( getStr([9.e100_RKH, -123456._RKH]) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr((1._CKH,-1._CKH))")
    call disp%show( getStr((1._CKH,-1._CKH)) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr([(1._CKD,-1._CKD), (2._CKD,-2._CKD)], signed = .true._LK)")
    call disp%show( getStr([(1._CKD,-1._CKD), (2._CKD,-2._CKD)], signed = .true._LK) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getStr(reshape([(1._CKL,-1._CKL), (2._CKL,-2._CKL), (3._CKL,-3._CKL), (4._CKL,-4._CKL)], shape = [2,2]), signed = .true._LK) ! input object of rank 2")
    call disp%show( getStr(reshape([(1._CKL,-1._CKL), (2._CKL,-2._CKL), (3._CKL,-3._CKL), (4._CKL,-4._CKL)], shape = [2,2]), signed = .true._LK) , deliml = SK_"""" )
    call disp%skip()

end program example