program example

    use pm_kind, only: SK, IK
    use pm_distUnifPar, only: getUnifParLogPDF
    use pm_io, only: display_type

    implicit none

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("1. / exp(getUnifParLogPDF(log(2.), [1, 2, 3, 4])) ! volume of the support.")
    call disp%show( 1. / exp(getUnifParLogPDF(log(2.), [1, 2, 3, 4])) )
    call disp%skip()

    call disp%skip()
    call disp%show("1. / exp(getUnifParLogPDF(log([1., 2., 3., 4.]))) ! volume of the support.")
    call disp%show( 1. / exp(getUnifParLogPDF(log([1., 2., 3., 4.]))) )
    call disp%skip()

end program example