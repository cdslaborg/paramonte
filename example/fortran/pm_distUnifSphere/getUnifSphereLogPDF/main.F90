program example

    use pm_kind, only: SK, IK
    use pm_matrixInit, only: uppLowDia, getMatInit
    use pm_distUnifSphere, only: getUnifSphereLogPDF
    use pm_io, only: display_type

    implicit none

    type(display_type)      :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("1. / exp(getUnifSphereLogPDF(log(1.), ndim = [1, 2, 3, 4])) ! area of the sphere embedded in ndim dimensions with unit radius.")
    call disp%show( 1. / exp(getUnifSphereLogPDF(log(1.), ndim = [1, 2, 3, 4])) )
    call disp%skip()

end program example