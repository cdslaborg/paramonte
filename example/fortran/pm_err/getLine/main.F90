#ifndef __LINE__
#define __LINE__ 0
#endif
program example

    use pm_kind, only: SK
    use pm_err, only: getLine
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%show("getLine(__LINE__)//SK_': This is the line number in the source file where printing was instructed.'")
    call disp%show( getLine(__LINE__)//SK_': This is the line number in the source file where printing was instructed.', deliml = SK_"""" )
    call disp%skip()

end program example