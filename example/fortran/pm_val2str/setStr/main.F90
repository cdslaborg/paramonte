program example

    use pm_kind, only: IK, LK
    use pm_kind, only: SK ! All other processor string kinds are also suppored.
    use pm_kind, only: IK8, IK16, IK32, IK64
    use pm_kind, only: CK32, CK64, CK128
    use pm_kind, only: RK32, RK64, RK128
    use pm_io, only: display_type
    use pm_val2str, only: setStr

    implicit none

    integer(IK)         :: lenstr
    character(1023, SK) :: str

    type(display_type)  :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call setStr(str, lenstr, 'paramonte')")
                    call setStr(str, lenstr, 'paramonte') 
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, '  paramonte')")
                    call setStr(str, lenstr, '  paramonte')
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, '  paramonte  ')")
                    call setStr(str, lenstr, '  paramonte  ')
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, ['paramonte', 'library  '])")
                    call setStr(str, lenstr, ['paramonte', 'library  '])
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, ['paramonte  ', '  library  '])")
                    call setStr(str, lenstr, ['paramonte  ', '  library  '])
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, .true.)")
                    call setStr(str, lenstr, .true.)
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, .false.)")
                    call setStr(str, lenstr, .false.)
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, -1_IK8)")
                    call setStr(str, lenstr, -1_IK8)
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, 123_IK16)")
                    call setStr(str, lenstr, 123_IK16)
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, 1234_IK32, signed = .true._LK)")
                    call setStr(str, lenstr, 1234_IK32, signed = .true._LK)
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, 987654321_IK64)")
                    call setStr(str, lenstr, 987654321_IK64)
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, 987654321_IK64, format = '(I15)')")
                    call setStr(str, lenstr, 987654321_IK64, format = '(I15)')
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, 123._RK32)")
                    call setStr(str, lenstr, 123._RK32)
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, 123._RK64, format = '(E25.15)') ! use scientific notation with a width of 25.")
                    call setStr(str, lenstr, 123._RK64, format = '(E25.15)')
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, [9.e1000_RK128, -123456._RK128])")
                    call setStr(str, lenstr, [9.e1000_RK128, -123456._RK128])
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, (1._CK128,-1._CK128))")
                    call setStr(str, lenstr, (1._CK128,-1._CK128))
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, [(1._CK64,-1._CK64), (2._CK64,-2._CK64)], signed = .true._LK)")
                    call setStr(str, lenstr, [(1._CK64,-1._CK64), (2._CK64,-2._CK64)], signed = .true._LK)
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("call setStr(str, lenstr, reshape([(1._CK32,-1._CK32), (2._CK32,-2._CK32), (3._CK32,-3._CK32), (4._CK32,-4._CK32)], shape = [2,2]), signed = .true._LK) ! input object of rank 2")
                    call setStr(str, lenstr, reshape([(1._CK32,-1._CK32), (2._CK32,-2._CK32), (3._CK32,-3._CK32), (4._CK32,-4._CK32)], shape = [2,2]), signed = .true._LK)
    call disp%show("lenstr")
    call disp%show( lenstr )
    call disp%show("str(1:lenstr)")
    call disp%show( str(1:lenstr) , deliml = SK_"""" )
    call disp%skip()

end program example