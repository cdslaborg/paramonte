program example

    use pm_kind, only: LK, IK, SK
    use pm_io, only: display_type
    use pm_sysPath, only: getPathVerbatim

    implicit none

    logical(LK) :: failed

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getPathVerbatim(SK_'Para/Monte')")
    call disp%show( getPathVerbatim(SK_"Para/Monte") , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathVerbatim(SK_""Para""/'Monte"")")
    call disp%show( getPathVerbatim(SK_"Para""/'Monte") , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathVerbatim(SK_""Para/'\''Monte"")")
    call disp%show( getPathVerbatim(SK_"Para/'\''Monte") , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathVerbatim(SK_""'Para/'\''Monte'"")")
    call disp%show( getPathVerbatim(SK_"'Para/'\''Monte'") , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathVerbatim(SK_'Para\Monte')")
    call disp%show( getPathVerbatim(SK_"Para\Monte") , deliml = SK_"'" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathVerbatim(SK_'""Para\Monte""')")
    call disp%show( getPathVerbatim(SK_"""Para\Monte""") , deliml = SK_"'" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathVerbatim(SK_""Para""\'Monte"")")
    call disp%show( getPathVerbatim(SK_"Para""\'Monte") , deliml = SK_"'" )
    call disp%skip()

    call disp%skip()
    call disp%show("getPathVerbatim(SK_'""')")
    call disp%show( getPathVerbatim(SK_'""') , deliml = SK_"'" )
    call disp%skip()

end program example