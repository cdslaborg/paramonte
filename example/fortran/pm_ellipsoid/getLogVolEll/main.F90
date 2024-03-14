program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_ellipsoid, only: getLogVolEll

    implicit none

    real, allocatable :: gramian(:,:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("gramian = reshape([1., 0., 0., 1.], [2,2]) ! correlation 0.")
                    gramian = reshape([1., 0., 0., 1.], [2,2])
    call disp%show("exp(getLogVolEll(gramian)) ! correlation 0.")
    call disp%show( exp(getLogVolEll(gramian)) )
    call disp%skip()

    call disp%skip()
    call disp%show("gramian = reshape([1., 0.5, 0.5, 1.], [2,2]) ! correlation 0.5")
                    gramian = reshape([1., 0.5, 0.5, 1.], [2,2])
    call disp%show("exp(getLogVolEll(gramian))")
    call disp%show( exp(getLogVolEll(gramian)) )
    call disp%skip()

    call disp%skip()
    call disp%show("gramian = reshape([1., 0.99, 0.99, 1.], [2,2]) ! correlation 0.99")
                    gramian = reshape([1., 0.99, 0.99, 1.], [2,2])
    call disp%show("exp(getLogVolEll(gramian))")
    call disp%show( exp(getLogVolEll(gramian)) )
    call disp%skip()

end program example