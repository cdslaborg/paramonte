program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_paramonte, only: getParaMonteSplash

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getParaMonteSplash()")
    call disp%show( getParaMonteSplash() , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getParaMonteSplash(width = 90_IK, fill = SK_'.', lwfill = SK_'.', rwfill = SK_'.')")
    call disp%show( getParaMonteSplash(width = 90_IK, fill = SK_'.', lwfill = SK_'.', rwfill = SK_'.') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getParaMonteSplash(width = 90_IK, fill = SK_'.', lwsize = 1_IK, twsize = 1_IK, rwsize = 2_IK, bwsize = 3_IK, lwfill = SK_':', twfill = SK_'-', rwfill = SK_'/', bwfill = SK_'%')")
    call disp%show( getParaMonteSplash(width = 90_IK, fill = SK_'.', lwsize = 1_IK, twsize = 1_IK, rwsize = 2_IK, bwsize = 3_IK, lwfill = SK_':', twfill = SK_'-', rwfill = SK_'/', bwfill = SK_'%') , deliml = SK_"""" )
    call disp%skip()

    call disp%skip()
    call disp%show("getParaMonteSplash(width = 30_IK)")
    call disp%show( getParaMonteSplash(width = 30_IK) , deliml = SK_"""" )
    call disp%skip()

end program example