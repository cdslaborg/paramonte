program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_ellipsoid, only: isMemberEll

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("isMemberEll(radius = 1., point = [0., 0.])")
    call disp%show( isMemberEll(radius = 1., point = [0., 0.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMemberEll(radius = 1., point = [1., 0.])")
    call disp%show( isMemberEll(radius = 1., point = [1., 0.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMemberEll(radius = 1., point = [1., 1.])")
    call disp%show( isMemberEll(radius = 1., point = [1., 1.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMemberEll(radius = 1., point = -[1., 1.])")
    call disp%show( isMemberEll(radius = 1., point = -[1., 1.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMemberEll(radius = 1., center = [.5, .5], point = -[1., 1.])")
    call disp%show( isMemberEll(radius = 1., center = [.5, .5], point = -[1., 1.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMemberEll(radius = 1., center = -[.5, .5], point = -[1., 1.])")
    call disp%show( isMemberEll(radius = 1., center = -[.5, .5], point = -[1., 1.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMemberEll(invGram = reshape([1., 0. , 0., 1.], [2,2]), point = [0., 1.]) ! correlation .0")
    call disp%show( isMemberEll(invGram = reshape([1., 0. , 0., 1.], [2,2]), point = [0., 1.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMemberEll(invGram = reshape([5.2632, -4.7368, -4.7368, 5.2632], [2,2]), point = -[1., 1.]) ! correlation .9")
    call disp%show( isMemberEll(invGram = reshape([5.2632, -4.7368, -4.7368, 5.2632], [2,2]), point = -[1., 1.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMemberEll(invGram = reshape([5.2632, -4.7368, -4.7368, 5.2632], [2,2]), center = [.5, .5], point = [1., 0.]) ! correlation .9")
    call disp%show( isMemberEll(invGram = reshape([5.2632, -4.7368, -4.7368, 5.2632], [2,2]), center = [.5, .5], point = [1., 0.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("isMemberEll(invGram = reshape([5.2632, -4.7368, -4.7368, 5.2632], [2,2]), center = -[.5, .5], point = -[1., 1.]) ! correlation .9")
    call disp%show( isMemberEll(invGram = reshape([5.2632, -4.7368, -4.7368, 5.2632], [2,2]), center = -[.5, .5], point = -[1., 1.]) )
    call disp%skip()

end program example