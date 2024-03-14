program example

    use pm_kind, only: SK, IK ! all string kinds are supported.
    use pm_kind, only: RKH
    use pm_io, only: display_type
    use pm_str, only: getTrimmedTZ
    use pm_val2str, only: getStr

    implicit none

    character(:, SK), allocatable :: str

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the number of digits in a given integer.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getTrimmedTZ('paramonte V. 2.000')")
    call disp%show( getTrimmedTZ('paramonte V. 2.000') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getTrimmedTZ('0000')")
    call disp%show( getTrimmedTZ('0000') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getTrimmedTZ('0001')")
    call disp%show( getTrimmedTZ('0001') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("str = getStr(1._RKH)")
                    str = getStr(1._RKH)
    call disp%show("str")
    call disp%show( str , deliml = SK_"""" )
    call disp%show("getTrimmedTZ(str)")
    call disp%show( getTrimmedTZ(str) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("str = getStr(1._RKH + epsilon(1._RKH))")
                    str = getStr(1._RKH + epsilon(1._RKH))
    call disp%show("str")
    call disp%show( str , deliml = SK_"""" )
    call disp%show("getTrimmedTZ(str)")
    call disp%show( getTrimmedTZ(str) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("str = getStr(1.2_RKH + epsilon(1._RKH))")
                    str = getStr(1.2_RKH + epsilon(1._RKH))
    call disp%show("str")
    call disp%show( str , deliml = SK_"""" )
    call disp%show("getTrimmedTZ(str)")
    call disp%show( getTrimmedTZ(str) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("str = getStr(1.03_RKH + epsilon(1._RKH))")
                    str = getStr(1.03_RKH + epsilon(1._RKH))
    call disp%show("str")
    call disp%show( str , deliml = SK_"""" )
    call disp%show("getTrimmedTZ(str)")
    call disp%show( getTrimmedTZ(str) , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getTrimmedTZ('1.375007620000000000000000')")
    call disp%show( getTrimmedTZ('1.375007620000000000000000') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getTrimmedTZ('.000000000000000000000000')")
    call disp%show( getTrimmedTZ('.000000000000000000000000') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getTrimmedTZ('.000.000...010')")
    call disp%show( getTrimmedTZ('.000.000...010') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getTrimmedTZ('.000 .000.. .010')")
    call disp%show( getTrimmedTZ('.000 .000.. .010') , deliml = SK_"""" )
    call disp%skip()

end program example