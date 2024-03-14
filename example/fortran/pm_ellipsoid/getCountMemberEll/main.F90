program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_ellipsoid, only: getCountMemberEll

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getCountMemberEll(radius = 1., point = reshape([0., 0., 1., 0., 1., 1., -1., -1.], [2, 4]))")
    call disp%show( getCountMemberEll(radius = 1., point = reshape([0., 0., 1., 0., 1., 1., -1., -1.], [2, 4])) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountMemberEll(radius = 2., center = [-1.5, 0.], point = reshape([0., 0., 1., 0., 1., 1., -1., -1.], [2, 4]))")
    call disp%show( getCountMemberEll(radius = 2., center = [-1.5, 0.], point = reshape([0., 0., 1., 0., 1., 1., -1., -1.], [2, 4])) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountMemberEll(invGram = reshape([1., 0. , 0., 1.], [2,2]), point = reshape([0., 0., 1., 0., 1., 1., -1., -1.], [2, 4])) ! correlation 0.")
    call disp%show( getCountMemberEll(invGram = reshape([1., 0. , 0., 1.], [2,2]), point = reshape([0., 0., 1., 0., 1., 1., -1., -1.], [2, 4])) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountMemberEll(invGram = reshape([5.2632, -4.7368, -4.7368, 5.2632], [2,2]), center = [0., 0.], point = reshape([0., 0., 1., 0., 1., 1., -1., -1.], [2, 4])) ! correlation 0.9")
    call disp%show( getCountMemberEll(invGram = reshape([5.2632, -4.7368, -4.7368, 5.2632], [2,2]), center = [0., 0.], point = reshape([0., 0., 1., 0., 1., 1., -1., -1.], [2, 4])) )
    call disp%skip()

end program example