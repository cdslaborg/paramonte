program example

    use pm_kind, only: SK, IK
    use pm_matrixInit, only: uppLowDia, getMatInit
    use pm_distUnifEll, only: getUnifEllLogPDF
    use pm_io, only: display_type

    implicit none

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("1. / exp(getUnifEllLogPDF(log(1.), [1, 2, 3, 4])) ! volume of the spherical support with unit radius.")
    call disp%show( 1. / exp(getUnifEllLogPDF(log(1.), [1, 2, 3, 4])) )
    call disp%skip()

    call disp%skip()
    call disp%show("1. / exp(getUnifEllLogPDF(log([1., 1., 1., 1.]))) ! volume of the support.")
    call disp%show( 1. / exp(getUnifEllLogPDF(log([1., 1., 1., 1.]))) )
    call disp%skip()

    call disp%skip()
    call disp%show("1. / exp(getUnifEllLogPDF(getMatInit([2, 2], uppLowDia, 0.5, 0.5, 1.))) ! volume of the support.")
    call disp%show( 1. / exp(getUnifEllLogPDF(getMatInit([2, 2], uppLowDia, 0.5, 0.5, 1.))) )
    call disp%skip()

    call disp%skip()
    call disp%show("1. / exp(getUnifEllLogPDF(getMatInit([2, 2], uppLowDia, 0.99, 0.99, 1.))) ! volume of the support.")
    call disp%show( 1. / exp(getUnifEllLogPDF(getMatInit([2, 2], uppLowDia, 0.99, 0.99, 1.))) )
    call disp%skip()

end program example