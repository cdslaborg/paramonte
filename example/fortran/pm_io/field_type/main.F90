program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: field_type
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_distUnif, only: setUnifRand

    implicit none

    type(display_type) :: disp
    disp = display_type(file = SK_"main.out.F90")

    call disp%skip()
    call disp%show('disp = display_type(file = SK_"main.out.F90", format = field_type(complex = SK_"math"))')
                    disp = display_type(file = SK_"main.out.F90", format = field_type(complex = SK_"math"))
    call disp%show("call disp%show( getUnifRand((-1., -1.), (+1., +1.), 5_IK) )")
                    call disp%show( getUnifRand((-1., -1.), (+1., +1.), 5_IK) )
    call disp%skip()
    call disp%show('disp = display_type(file = SK_"main.out.F90", deliml = field_type(string = SK_"**", real = SK_"("), delimr = field_type(string = SK_"**", real = SK_")"))')
                    disp = display_type(file = SK_"main.out.F90", deliml = field_type(string = SK_"**", real = SK_"("), delimr = field_type(string = SK_"**", real = SK_")"))
    call disp%show("call disp%show('This is a custom (star-) delimited string.')", deliml = SK_"", delimr = SK_"")
                    call disp%show('This is a custom (star-) delimited string.')
    call disp%show("call disp%show( getUnifRand(-1., +1., 5_IK) )", deliml = SK_"", delimr = SK_"")
                    call disp%show( getUnifRand(-1., +1., 5_IK) )
    call disp%skip()

end program example