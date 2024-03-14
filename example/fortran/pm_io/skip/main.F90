program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = SK_"main.out.F90")

    call disp%skip()
    call disp%show('disp = display_type(file = SK_"main.out.F90")')
    call disp%show("call disp%show('before skip')")
                    call disp%show('before skip', deliml = """")
    call disp%show("call disp%skip()")
                    call disp%skip()
    call disp%show("call disp%show('after skip')")
                    call disp%show('after skip', deliml = """")
    call disp%skip()

    call disp%skip()
    call disp%show('disp = display_type(disp%unit)')
                    disp = display_type(disp%unit)
    call disp%show("call disp%show('before skip')")
                    call disp%show('before skip', deliml = """")
    call disp%show("call disp%skip(3)")
                    call disp%skip(3)
    call disp%show("call disp%show('after skip')")
                    call disp%show('after skip', deliml = """")
    call disp%skip()

end program example